package edu.mit.csail.cgs.reeder.germ;

import java.io.*;
import java.lang.reflect.Array;
import java.sql.SQLException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqAlignment;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLoader;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.ChromosomeGenerator;
import edu.mit.csail.cgs.projects.readdb.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.*;

public class SproutStorage {
	
	private DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	//private SortedSet<Pair<StrandedPoint,StrandedPoint>> leftPairSet;
	//private Map<Pair<StrandedPoint,StrandedPoint>,Integer> leftPairIndex;
	//private List<Pair<StrandedPoint,StrandedPoint>> leftPairList;
	private Pair<StrandedPoint,StrandedPoint>[] leftPairArray;
	//private SortedSet<Pair<StrandedPoint,StrandedPoint>> rightPairSet;
	//private Map<Pair<StrandedPoint,StrandedPoint>,Integer> rightPairIndex;
	//private List<Pair<StrandedPoint,StrandedPoint>> rightPairList;
	private Pair<StrandedPoint,StrandedPoint>[] rightPairArray;
	
	private Comparator<Pair<StrandedPoint,StrandedPoint>> leftComp = new LeftPairComparator();
	private Comparator<Pair<StrandedPoint,StrandedPoint>> rightComp = new RightPairComparator();


	public static void main(String[] args) throws NotFoundException, SQLException, IOException, ClientException {
		
		Genome g = SproutUtils.parseGenome(args);
		String outfile = Args.parseString(args, "outfile", "");
		List<ChipSeqLocator> locators = Args.parseChipSeq(args,"align");
		String type = Args.parseString(args, "type", "");
		int numgroups = Args.parseInteger(args, "numgroups", 0);
		int num = Args.parseInteger(args, "num", 0);
		SproutStorage storage;
		if (type.equals("")) {
			if (numgroups==0) {
				storage = fromReadDB(g, locators);
			} else {
				storage = fromReadDB(g, locators, num, numgroups);
			}
		} else if (type.equals("self")) {
			storage = fromReadDB(g, locators, true);
		} else {
			storage = fromReadDB(g, locators, false);
		}
		storage.toFile(outfile);
		
		
		/*
		Genome g = SproutUtils.parseGenome(args);
		String infile = Args.parseString(args, "infile", "");
		String outfile = Args.parseString(args, "outfile", "");
		fromFileToFile(g, infile, outfile);
		*/
		
		/*
		String readfile = Args.parseString(args, "readfile", "");
		SproutStorage storage = fromFile(g,readfile,Args.parseInteger(args, "size", 0));
		System.err.println(storage.getLeftN()+"\t"+storage.getRightN());
		*/
		/*
		Genome g = Args.parseGenome(args).cdr();
		String outfile = Args.parseString(args, "outfile", "");
		List<ChipSeqLocator> locators = Args.parseChipSeq(args,"align");
		List<Pair<StrandedPoint,StrandedPoint>> pairlist = fromReadDBNonDup(g, locators);
		PrintStream out = new PrintStream(outfile);
		for (Pair<StrandedPoint,StrandedPoint> pair : pairlist) {
			out.println(pair.car()+"\t"+pair.cdr());
		}
		out.flush();
		out.close();
		*/
	}

	public SproutStorage(Pair<StrandedPoint,StrandedPoint>[] pairs) {
		this.leftPairArray = pairs;
		this.rightPairArray = Arrays.copyOf(pairs, pairs.length);
		System.err.println("storage initialized "+dfm.format(new Date()));
		Arrays.sort(this.leftPairArray,leftComp);
		System.err.println("first sort completed "+dfm.format(new Date()));
		Arrays.sort(this.rightPairArray,rightComp);
		System.err.println("second sort completed "+dfm.format(new Date()));
		System.err.println(getLeftN()+" pairs in storage");
		
		/*
		this.leftPairList = new ArrayList<Pair<StrandedPoint,StrandedPoint>>(pairs);
		//this.rightPairList = new ArrayList<Pair<StrandedPoint,StrandedPoint>>(pairs);
		Collections.sort(this.leftPairList, leftComp);
		Collections.sort(this.rightPairList, rightComp);
		*/
		
		/*
		this.rightPairSet = new TreeSet<Pair<StrandedPoint,StrandedPoint>>(new RightPairComparator());
		this.rightPairSet.addAll(pairs);
		this.leftPairIndex = new HashMap<Pair<StrandedPoint,StrandedPoint>,Integer>();
		this.rightPairIndex = new HashMap<Pair<StrandedPoint,StrandedPoint>,Integer>();
		this.leftPairList = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		this.rightPairList = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		int index = 0;
		for (Pair<StrandedPoint,StrandedPoint> pair : pairs) {
			leftPairList.add(pair);
			leftPairIndex.put(pair, index++);
		}
		index = 0;
		
		for (Pair<StrandedPoint,StrandedPoint> pair : this.rightPairSet) {
			rightPairList.add(pair);
			rightPairIndex.put(pair, index++);
		}
		*/
	}

	public static SproutStorage fromReadDB(Genome g, List<ChipSeqLocator> locators) throws SQLException, NotFoundException, IOException, ClientException {
		SortedSet<Pair<StrandedPoint,StrandedPoint>> pairs = new TreeSet<Pair<StrandedPoint,StrandedPoint>>(new LeftPairComparator());
		Map<Integer,String> revChromIDMap = g.getRevChromIDMap();
		int readdbcount = 0;
		Set<String> aligns = new HashSet<String>();
		Client client = new Client();
		ChipSeqLoader loader = new ChipSeqLoader();
		for (ChipSeqLocator locator : locators) {
			for (ChipSeqAlignment a : loader.loadAlignments(locator,g)) {
				aligns.add(Integer.toString(a.getDBID()));
			}
		}

		ChromosomeGenerator<Genome> cg = new ChromosomeGenerator<Genome>();
		Iterator<Region> chromiter = cg.execute(g);
		while (chromiter.hasNext()) {
			Region chrom = chromiter.next();

			try {
				for (String a : aligns) {
					List<PairedHit> hits = client.getPairedHits(a, g.getChromID(chrom.getChrom()), true, chrom.getStart(), chrom.getEnd(), null, null);
					readdbcount += hits.size();
					for (PairedHit hit : hits) {
						StrandedPoint leftPoint = new StrandedPoint(g, revChromIDMap.get(hit.leftChrom), hit.leftPos, hit.leftStrand ? '+' : '-');
						StrandedPoint rightPoint = new StrandedPoint(g, revChromIDMap.get(hit.rightChrom), hit.rightPos, hit.rightStrand ? '+' : '-');
						if (leftPoint.compareTo(rightPoint)>0) {
							pairs.add(new Pair<StrandedPoint,StrandedPoint>(rightPoint,leftPoint));
						} else {
							pairs.add(new Pair<StrandedPoint,StrandedPoint>(leftPoint,rightPoint));
						}
					}
				}
			} catch (Exception e) {
				System.err.println(e);
			}
		}
		System.err.println("Hits in ReadDB: "+readdbcount);
		System.err.println("Non-duplicate hits: "+pairs.size());
		client.close();
		Pair<StrandedPoint,StrandedPoint>[] tor = (Pair<StrandedPoint,StrandedPoint>[]) Array.newInstance(Pair.class, pairs.size());
		int index = 0;
		for (Pair<StrandedPoint,StrandedPoint> pair : pairs) {
			tor[index++] = pair;
		}
		return new SproutStorage(tor);
	}
	
	public static SproutStorage fromReadDB(Genome g, List<ChipSeqLocator> locators, int num, int outof) throws SQLException, NotFoundException, IOException, ClientException {
		SortedSet<Pair<StrandedPoint,StrandedPoint>> pairs = new TreeSet<Pair<StrandedPoint,StrandedPoint>>(new LeftPairComparator());
		Map<Integer,String> revChromIDMap = g.getRevChromIDMap();
		Map<String,Integer> chromsizemap = g.getChromLengthMap();
		int readdbcount = 0;
		Set<String> aligns = new HashSet<String>();
		Client client = new Client();
		ChipSeqLoader loader = new ChipSeqLoader();
		for (ChipSeqLocator locator : locators) {
			for (ChipSeqAlignment a : loader.loadAlignments(locator,g)) {
				aligns.add(Integer.toString(a.getDBID()));
			}
		}

		List<String> chroms = new ArrayList<String>(chromsizemap.keySet());
		Collections.sort(chroms);
		int pergroup = chroms.size() / outof;
		int start = (num-1)*pergroup;
		int end = num==outof ? chroms.size() : num*pergroup;
		for (int i=start; i<end; i++) {
			Region chrom = SproutUtils.chromRegion(chroms.get(i), g);

			try {
				for (String a : aligns) {
					List<PairedHit> hits = client.getPairedHits(a, g.getChromID(chrom.getChrom()), true, chrom.getStart(), chrom.getEnd(), null, null);
					readdbcount += hits.size();
					for (PairedHit hit : hits) {
						StrandedPoint leftPoint = new StrandedPoint(g, revChromIDMap.get(hit.leftChrom), hit.leftPos, hit.leftStrand ? '+' : '-');
						StrandedPoint rightPoint = new StrandedPoint(g, revChromIDMap.get(hit.rightChrom), hit.rightPos, hit.rightStrand ? '+' : '-');
						if (leftPoint.compareTo(rightPoint)>0) {
							pairs.add(new Pair<StrandedPoint,StrandedPoint>(rightPoint,leftPoint));
						} else {
							pairs.add(new Pair<StrandedPoint,StrandedPoint>(leftPoint,rightPoint));
						}
					}
				}
			} catch (Exception e) {
				System.err.println(e);
			}
			System.err.println(chrom);
		}
		System.err.println("Hits in ReadDB: "+readdbcount);
		System.err.println("Non-duplicate hits: "+pairs.size());
		client.close();
		Pair<StrandedPoint,StrandedPoint>[] tor = (Pair<StrandedPoint,StrandedPoint>[]) Array.newInstance(Pair.class, pairs.size());
		int index = 0;
		for (Pair<StrandedPoint,StrandedPoint> pair : pairs) {
			tor[index++] = pair;
		}
		return new SproutStorage(tor);
	}
	
	public static SproutStorage fromReadDB(Genome g, List<ChipSeqLocator> locators, boolean self) throws SQLException, NotFoundException, IOException, ClientException {
		SortedSet<Pair<StrandedPoint,StrandedPoint>> pairs = new TreeSet<Pair<StrandedPoint,StrandedPoint>>(new LeftPairComparator());
		Map<Integer,String> revChromIDMap = g.getRevChromIDMap();
		int readdbcount = 0;
		Set<String> aligns = new HashSet<String>();
		Client client = new Client();
		ChipSeqLoader loader = new ChipSeqLoader();
		for (ChipSeqLocator locator : locators) {
			for (ChipSeqAlignment a : loader.loadAlignments(locator,g)) {
				aligns.add(Integer.toString(a.getDBID()));
			}
		}

		ChromosomeGenerator<Genome> cg = new ChromosomeGenerator<Genome>();
		Iterator<Region> chromiter = cg.execute(g);
		while (chromiter.hasNext()) {
			Region chrom = chromiter.next();

			try {
				for (String a : aligns) {
					List<PairedHit> hits = client.getPairedHits(a, g.getChromID(chrom.getChrom()), true, chrom.getStart(), chrom.getEnd(), null, null);
					readdbcount += hits.size();
					for (PairedHit hit : hits) {
						StrandedPoint leftPoint = new StrandedPoint(g, revChromIDMap.get(hit.leftChrom), hit.leftPos, hit.leftStrand ? '+' : '-');
						StrandedPoint rightPoint = new StrandedPoint(g, revChromIDMap.get(hit.rightChrom), hit.rightPos, hit.rightStrand ? '+' : '-');
						if (leftPoint.compareTo(rightPoint)>0) {
							StrandedPoint tmp = leftPoint;
							leftPoint = rightPoint;
							rightPoint = tmp;
							if (self && leftPoint.getChrom().equals(rightPoint.getChrom()) && (leftPoint.getStrand()=='-' && rightPoint.getStrand()=='+')) {
								pairs.add(new Pair<StrandedPoint,StrandedPoint>(leftPoint,rightPoint));
							} else if (!self && leftPoint.getChrom().equals(rightPoint.getChrom()) && !(leftPoint.getStrand()=='-' && rightPoint.getStrand()=='+')) {
								pairs.add(new Pair<StrandedPoint,StrandedPoint>(leftPoint,rightPoint));
							}
						} else {
							if (self && leftPoint.getChrom().equals(rightPoint.getChrom()) && (leftPoint.getStrand()=='-' && rightPoint.getStrand()=='+')) {
								pairs.add(new Pair<StrandedPoint,StrandedPoint>(leftPoint,rightPoint));
							} else if (!self && leftPoint.getChrom().equals(rightPoint.getChrom()) && !(leftPoint.getStrand()=='-' && rightPoint.getStrand()=='+')) {
								pairs.add(new Pair<StrandedPoint,StrandedPoint>(leftPoint,rightPoint));
							}
						}
					}
				}
			} catch (Exception e) {
				System.err.println(e);
			}
		}
		System.err.println("Hits in ReadDB: "+readdbcount);
		System.err.println("Non-duplicate hits: "+pairs.size());
		client.close();
		Pair<StrandedPoint,StrandedPoint>[] tor = (Pair<StrandedPoint,StrandedPoint>[]) Array.newInstance(Pair.class, pairs.size());
		int index = 0;
		for (Pair<StrandedPoint,StrandedPoint> pair : pairs) {
			tor[index++] = pair;
		}
		return new SproutStorage(tor);
	}
	
	public static List<Pair<StrandedPoint,StrandedPoint>> fromReadDBNonDup(Genome g, List<ChipSeqLocator> locators) throws SQLException, NotFoundException, IOException, ClientException {
		//SortedSet<Pair<StrandedPoint,StrandedPoint>> pairs = new TreeSet<Pair<StrandedPoint,StrandedPoint>>(new LeftPairComparator());
		List<Pair<StrandedPoint,StrandedPoint>> pairlist = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		Map<Integer,String> revChromIDMap = g.getRevChromIDMap();
		int readdbcount = 0;
		Set<String> aligns = new HashSet<String>();
		Client client = new Client();
		ChipSeqLoader loader = new ChipSeqLoader();
		for (ChipSeqLocator locator : locators) {
			for (ChipSeqAlignment a : loader.loadAlignments(locator,g)) {
				aligns.add(Integer.toString(a.getDBID()));
			}
		}

		ChromosomeGenerator<Genome> cg = new ChromosomeGenerator<Genome>();
		Iterator<Region> chromiter = cg.execute(g);
		while (chromiter.hasNext()) {
			Region chrom = chromiter.next();

			try {
				for (String a : aligns) {
					List<PairedHit> hits = client.getPairedHits(a, g.getChromID(chrom.getChrom()), true, chrom.getStart(), chrom.getEnd(), null, null);
					readdbcount += hits.size();
					for (PairedHit hit : hits) {
						StrandedPoint leftPoint = new StrandedPoint(g, revChromIDMap.get(hit.leftChrom), hit.leftPos, hit.leftStrand ? '+' : '-');
						StrandedPoint rightPoint = new StrandedPoint(g, revChromIDMap.get(hit.rightChrom), hit.rightPos, hit.rightStrand ? '+' : '-');
						if (leftPoint.compareTo(rightPoint)>0) {
							pairlist.add(new Pair<StrandedPoint,StrandedPoint>(rightPoint,leftPoint));
						} else {
							pairlist.add(new Pair<StrandedPoint,StrandedPoint>(leftPoint,rightPoint));
						}
					}
				}
			} catch (Exception e) {
				System.err.println(e);
			}
		}
		System.err.println("Hits in ReadDB: "+readdbcount);
		System.err.println("Non-duplicate hits: "+pairlist.size());
		client.close();
		return pairlist;
	}
	
	public static SproutStorage fromFileEven(Genome g, String file) throws IOException {
		SortedSet<Pair<StrandedPoint,StrandedPoint>> pairs = new TreeSet<Pair<StrandedPoint,StrandedPoint>>(new LeftPairComparator());
		int chrom1, chrom2;
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			StrandedPoint p1 = StrandedPoint.fromString(g, split[0]);
			StrandedPoint p2 = StrandedPoint.fromString(g, split[1]);
			try {
				chrom1 = Integer.valueOf(p1.getChrom());
				chrom2 = Integer.valueOf(p2.getChrom());
				if (chrom1%2==0 && chrom2%2==0) {
					pairs.add(new Pair<StrandedPoint,StrandedPoint>(p1,p2));
				}
			} catch (NumberFormatException e) {
				
			}
		}
		r.close();
		Pair<StrandedPoint,StrandedPoint>[] tor = (Pair<StrandedPoint,StrandedPoint>[]) Array.newInstance(Pair.class, pairs.size());
		int index = 0;
		for (Pair<StrandedPoint,StrandedPoint> pair : pairs) {
			tor[index++] = pair;
		}
		return new SproutStorage(tor);
	}
	
	public static SproutStorage fromFile(Genome g, String file, Region region) throws IOException {
		SortedSet<Pair<StrandedPoint,StrandedPoint>> pairs = new TreeSet<Pair<StrandedPoint,StrandedPoint>>(new LeftPairComparator());
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			StrandedPoint p1 = StrandedPoint.fromString(g, split[0]);
			StrandedPoint p2 = StrandedPoint.fromString(g, split[1]);
			if (p1.compareTo(p1)>0) {
				System.err.println(p1+" > "+p2);
			}
			if (region.contains(p1) && region.contains(p2)) {
				pairs.add(new Pair<StrandedPoint,StrandedPoint>(p1,p2));
			}
		}
		r.close();
		Pair<StrandedPoint,StrandedPoint>[] tor = (Pair<StrandedPoint,StrandedPoint>[]) Array.newInstance(Pair.class, pairs.size());
		int index = 0;
		for (Pair<StrandedPoint,StrandedPoint> pair : pairs) {
			tor[index++] = pair;
		}
		return new SproutStorage(tor);
	}
	
	public static SproutStorage fromFile(Genome g, String file) throws IOException {
		SortedSet<Pair<StrandedPoint,StrandedPoint>> pairs = new TreeSet<Pair<StrandedPoint,StrandedPoint>>(new LeftPairComparator());
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			pairs.add(new Pair<StrandedPoint,StrandedPoint>(StrandedPoint.fromString(g, split[0]),StrandedPoint.fromString(g, split[1])));
		}
		r.close();
		Pair<StrandedPoint,StrandedPoint>[] tor = (Pair<StrandedPoint,StrandedPoint>[]) Array.newInstance(Pair.class, pairs.size());
		int index = 0;
		for (Pair<StrandedPoint,StrandedPoint> pair : pairs) {
			tor[index++] = pair;
		}
		return new SproutStorage(tor);
	}
	
	public static void fromFileToFile(Genome g, String file, String outfile) throws IOException {
		SortedSet<Pair<StrandedPoint,StrandedPoint>> pairs = new TreeSet<Pair<StrandedPoint,StrandedPoint>>(new LeftPairComparator());
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			pairs.add(new Pair<StrandedPoint,StrandedPoint>(StrandedPoint.fromString(g, split[0]),StrandedPoint.fromString(g, split[1])));
		}
		r.close();
		PrintStream out = new PrintStream(outfile);
		for (Pair<StrandedPoint,StrandedPoint> pair : pairs) {
			out.println(pair.car()+"\t"+pair.cdr());
		}
		out.flush();
		out.close();
	}
	
	public static SproutStorage fromFile(Genome g, String file, int size) throws IOException {
		Pair<StrandedPoint,StrandedPoint>[] tor = (Pair<StrandedPoint,StrandedPoint>[]) Array.newInstance(Pair.class, size);
		int index = 0;
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			tor[index++] = new Pair<StrandedPoint,StrandedPoint>(StrandedPoint.fromString(g, split[0]),StrandedPoint.fromString(g, split[1]));
		}
		r.close();
		return new SproutStorage(tor);
	}
	
	public void toFile(String file) throws FileNotFoundException {
		PrintStream out = new PrintStream(file);
		for (int i=0; i<leftPairArray.length; i++) {
			Pair<StrandedPoint,StrandedPoint> pair = leftPairArray[i];
			out.println(pair.car()+"\t"+pair.cdr());
		}
		out.flush();
		out.close();
	}
	
	public int getPileup(int width, Point p) {
		Pair<Integer,Integer> minmaxn = getLeftMinMaxn(new Region(p.getGenome(),p.getChrom(),p.getLocation()-width, p.getLocation()+1));
		int pileup = 0;
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pair = getLeftPair(n);
			if (PairedReadDistribution.isSelfLigation(pair) && pair.cdr().getLocation()>p.getLocation() && pair.car().distance(pair.cdr())<=width) {
				pileup++;
			}
		}
		return pileup;
	}
	
	public Pair<Integer,Integer> getLeftMinMaxn(Region r) {
		StrandedPoint minPoint = new StrandedPoint(r.getGenome(), "", 0, Character.MIN_VALUE);
		StrandedPoint startPoint = new StrandedPoint(r.getGenome(), r.getChrom(), r.getStart(), Character.MIN_VALUE);
		StrandedPoint endPoint = new StrandedPoint(r.getGenome(), r.getChrom(), r.getEnd(), Character.MAX_VALUE);
		/*
		SortedSet<Pair<StrandedPoint,StrandedPoint>> subset = leftPairSet.subSet(new Pair<StrandedPoint,StrandedPoint>(startPoint,minPoint), 
				new Pair<StrandedPoint,StrandedPoint>(endPoint,minPoint));
		if (subset.isEmpty()) {
			return new Pair<Integer,Integer>(-1,-1);
		} else {
			return new Pair<Integer,Integer>(leftPairIndex.get(subset.first()),leftPairIndex.get(subset.last()));
		}
		*/
		int minn = Arrays.binarySearch(leftPairArray, new Pair<StrandedPoint,StrandedPoint>(startPoint,minPoint), leftComp);
		if (minn<0) {
			minn = -minn - 1;
		}
		int maxn = Arrays.binarySearch(leftPairArray, new Pair<StrandedPoint,StrandedPoint>(endPoint,minPoint), leftComp);
		if (maxn<0) {
			maxn = -maxn - 1;
		}
		return new Pair<Integer,Integer>(minn,maxn);
	}
	
	public Pair<Integer,Integer> getRightMinMaxn(Region r) {
		StrandedPoint minPoint = new StrandedPoint(r.getGenome(), "", 0, Character.MIN_VALUE);
		StrandedPoint startPoint = new StrandedPoint(r.getGenome(), r.getChrom(), r.getStart(), Character.MIN_VALUE);
		StrandedPoint endPoint = new StrandedPoint(r.getGenome(), r.getChrom(), r.getEnd(), Character.MAX_VALUE);
		/*
		SortedSet<Pair<StrandedPoint,StrandedPoint>> subset = rightPairSet.subSet(new Pair<StrandedPoint,StrandedPoint>(minPoint,startPoint), 
				new Pair<StrandedPoint,StrandedPoint>(minPoint,endPoint));
		if (subset.isEmpty()) {
			return new Pair<Integer,Integer>(-1,-1);
		} else {
			return new Pair<Integer,Integer>(rightPairIndex.get(subset.first()),rightPairIndex.get(subset.last()));
		}
		*/
		int minn = Arrays.binarySearch(rightPairArray, new Pair<StrandedPoint,StrandedPoint>(minPoint,startPoint), rightComp);
		if (minn<0) {
			minn = -minn - 1;
		}
		int maxn = Arrays.binarySearch(rightPairArray, new Pair<StrandedPoint,StrandedPoint>(minPoint,endPoint), rightComp);
		if (maxn<0) {
			maxn = -maxn - 1;
		}
		return new Pair<Integer,Integer>(minn,maxn);
	}

	public Pair<StrandedPoint,StrandedPoint> getLeftPair(int n) {
		return leftPairArray[n];
	}
	
	public Pair<StrandedPoint,StrandedPoint> getRightPair(int n) {
		return rightPairArray[n];
	}

	public int getLeftN() {
		return leftPairArray.length;
	}
	
	public int getRightN() {
		return leftPairArray.length;
	}

	public static class LeftPairComparator implements Comparator<Pair<StrandedPoint,StrandedPoint>> {

		public int compare(Pair<StrandedPoint, StrandedPoint> arg0,
				Pair<StrandedPoint, StrandedPoint> arg1) {
			int tor = arg0.car().compareToStranded(arg1.car());
			if (tor==0) {
				return arg0.cdr().compareToStranded(arg1.cdr());
			} else {
				return tor;
			}
		}

	}
	
	public static class RightPairComparator implements Comparator<Pair<StrandedPoint,StrandedPoint>> {

		public int compare(Pair<StrandedPoint, StrandedPoint> arg0,
				Pair<StrandedPoint, StrandedPoint> arg1) {
			int tor = arg0.cdr().compareToStranded(arg1.cdr());
			if (tor==0) {
				return arg0.car().compareToStranded(arg1.car());
			} else {
				return tor;
			}
		}

	}

}
