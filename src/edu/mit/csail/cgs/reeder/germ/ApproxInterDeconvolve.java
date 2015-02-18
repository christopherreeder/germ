package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class ApproxInterDeconvolve {

	private static DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	private Genome g;
	private SproutStorage storage;
	private BBFileReader bbfr;
	private List<Double> readDist;
	private int binsize;

	public ApproxInterDeconvolve(String[] args) throws NotFoundException, NumberFormatException, IOException {
		g = SproutUtils.parseGenome(args);
		String margfile = Args.parseString(args, "margfile", "");
		String readdistfile = Args.parseString(args, "readdistfile", "");
		String readfile = Args.parseString(args, "readfile", "");
		int storagesize = Args.parseInteger(args, "storagesize", 0);
		binsize = Args.parseInteger(args, "binsize", 0);

		bbfr = new BBFileReader(margfile);
		storage = SproutStorage.fromFile(g, readfile, storagesize);
		readDist = new ArrayList<Double>();
		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(readdistfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			readDist = new ArrayList<Double>();
			for (int i=0; i<split.length; i++) {
				readDist.add(Double.parseDouble(split[i]));
			}
		}
	}

	/**
	 * Assumes readDist is for + strand reads
	 * @param args
	 * @throws NotFoundException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NotFoundException, IOException {
		String regionfile = Args.parseString(args, "regionfile", "");
		String outfile = Args.parseString(args, "outfile", "");
		PrintStream out = new PrintStream(outfile);

		int inc = 100;
		int index = 0;
		ApproxInterDeconvolve aid = new ApproxInterDeconvolve(args);
		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(regionfile));
		while ((s = r.readLine()) != null) {
			out.println(s);
			split = s.split("\t");
			Region query = Region.fromString(aid.g, split[0]);
			Region region = SproutUtils.chromRegion(query.getChrom(), aid.g);
			double marg = aid.getMarg(query);
			out.println(marg);
			double[] tor = aid.deconvolveStrip(query, region);
			Region reg = null;
			Region maxreg = null;
			double score = 0;
			List<Double> scores = new ArrayList<Double>();
			for (int i=0; i<tor.length; i++) {
				if (tor[i]>0) {
					if (reg==null) {
						reg = new Region(aid.g,region.getChrom(),i*aid.binsize,i*aid.binsize+aid.binsize);
						maxreg = new Region(aid.g,region.getChrom(),i*aid.binsize,i*aid.binsize+aid.binsize);
						score = tor[i];
						scores = new ArrayList<Double>();
						scores.add(tor[i]);
					} else if (i*aid.binsize-reg.getEnd()>aid.binsize) {
						out.print(maxreg+"\t"+score+"\t"+reg+"\t[");
						for (Double d : scores) {
							out.print(d+",");
						}
						out.println("]");
						reg = new Region(aid.g,region.getChrom(),i*aid.binsize,i*aid.binsize+aid.binsize);
						maxreg = new Region(aid.g,region.getChrom(),i*aid.binsize,i*aid.binsize+aid.binsize);
						score = tor[i];
						scores = new ArrayList<Double>();
						scores.add(tor[i]);
					} else {
						reg = new Region(aid.g,reg.getChrom(),reg.getStart(),i*aid.binsize+aid.binsize);
						if (tor[i]>score) {
							maxreg = new Region(aid.g,region.getChrom(),i*aid.binsize,i*aid.binsize+aid.binsize);
							score = tor[i];
						}
						scores.add(tor[i]);
					}
				}
			}
			out.print(maxreg+"\t"+score+"\t"+reg+"\t[");
			for (Double d : scores) {
				out.print(d+",");
			}
			out.println("]");

			index++;
			if (index % inc == 0) {
				System.err.println(index+" regions completed "+dfm.format(new Date()));
			}
		}
		out.flush();
		out.close();

		/*
		Genome g = SproutUtils.parseGenome(args);
		String margfile = Args.parseString(args, "margfile", "");
		String readdistfile = Args.parseString(args, "readdistfile", "");
		String readfile = Args.parseString(args, "readfile", "");
		int storagesize = Args.parseInteger(args, "storagesize", 0);
		String region1 = Args.parseString(args, "region1", "");
		String region2 = Args.parseString(args, "region2", "");
		int binsize = Args.parseInteger(args, "binsize", 0);
		String outfile = Args.parseString(args, "outfile", "");
		PrintStream out = new PrintStream(outfile);

		List<Double> readDist = new ArrayList<Double>();
		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(readdistfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			readDist = new ArrayList<Double>();
			for (int i=0; i<split.length; i++) {
				readDist.add(Double.parseDouble(split[i]));
			}
		}
		Region reg1 = Region.fromString(g, region1);
		Region bigreg1 = reg1.expand(readDist.size()/2, readDist.size()/2);
		Region reg2 = Region.fromString(g, region2);
		Region bigreg2 = reg2.expand(readDist.size()/2, readDist.size()/2);
		double[][] tor = new double[reg1.getWidth()/binsize][reg2.getWidth()/binsize];
		SproutStorage storage = SproutStorage.fromFile(g, readfile, storagesize);
		double[] marg1 = new double[reg1.getWidth()/binsize], marg2 = new double[reg2.getWidth()/binsize];
		BBFileReader bbfr = new BBFileReader(margfile);
		BigWigIterator bwiter = bbfr.getBigWigIterator("chr"+reg1.getChrom(), reg1.getStart(), "chr"+reg1.getChrom(), reg1.getEnd(), false);
		for (int i=0; i<marg1.length; i++) {
			marg1[i] = bwiter.next().getWigValue();
		}
		bwiter = bbfr.getBigWigIterator("chr"+reg2.getChrom(), reg2.getStart(), "chr"+reg2.getChrom(), reg2.getEnd(), false);
		for (int i=0; i<marg2.length; i++) {
			marg2[i] = bwiter.next().getWigValue();
		}

		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg1);
		List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (bigreg2.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}

		for (Pair<StrandedPoint,StrandedPoint> pair : localreads) {
			int x = (pair.car().getLocation()-reg1.getStart())/binsize;
			int y = (pair.cdr().getLocation()-reg2.getStart())/binsize;
			char xstrand = pair.car().getStrand();
			char ystrand = pair.cdr().getStrand();
			for (int i = Math.max(x-readDist.size()/2, 0); i<Math.min(x+readDist.size()/2, tor.length); i++) {
				int distind = x-i+readDist.size()/2;
				if (xstrand=='-') distind = readDist.size() - distind;
				double pru = readDist.get(distind) * marg1[x];
				for (int j = Math.max(y-readDist.size()/2, 0); j<Math.min(y+readDist.size()/2, tor[i].length); j++) {
					distind = y-j+readDist.size()/2;
					if (ystrand=='-') distind = readDist.size() - distind;
					double prv = readDist.get(distind) * marg2[y];
					tor[i][j] += pru * prv;
				}
			}
		}

		for (int i=0; i<tor.length; i++) {
			for (int j=0; j<tor[i].length; j++) {
				if (tor[i][j]>0) {
					out.println(i+"\t"+j+"\t"+tor[i][j]);
				}
			}
		}
		out.flush();
		out.close();
		 */
	}

	public double[][] deconvolve(Region reg1, Region reg2, PrintStream jointinfo) {
		Region bigreg1 = reg1.expand(binsize*readDist.size()/2, binsize*readDist.size()/2);
		Region bigreg2 = reg2.expand(binsize*readDist.size()/2, binsize*readDist.size()/2);

		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg1);
		Set<Pair<StrandedPoint,StrandedPoint>> localreads = new HashSet<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (bigreg2.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}
		minmaxn = storage.getLeftMinMaxn(bigreg2);
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (bigreg1.contains(pairn.cdr())) {
				localreads.add(pairn.swap());
			}
		}
		if (jointinfo!=null) {
			jointinfo.println("read pairs considered: "+localreads.size());
		}
		if (!localreads.isEmpty()) {
			double[] marg1 = new double[reg1.getWidth()/binsize], marg2 = new double[reg2.getWidth()/binsize];
			BigWigIterator bwiter = bbfr.getBigWigIterator("chr"+reg1.getChrom(), reg1.getStart(), "chr"+reg1.getChrom(), reg1.getEnd(), false);
			for (int i=0; i<marg1.length; i++) {
				marg1[i] = bwiter.next().getWigValue();
			}
			bwiter = bbfr.getBigWigIterator("chr"+reg2.getChrom(), reg2.getStart(), "chr"+reg2.getChrom(), reg2.getEnd(), false);
			for (int i=0; i<marg2.length; i++) {
				try {
					if (bwiter.hasNext()) {
						marg2[i] = bwiter.next().getWigValue();
					}
				} catch (Exception e) {
					System.err.println(reg2+"\t"+i);
					e.printStackTrace();
				}
			}
			double[][] tor = new double[reg1.getWidth()/binsize][reg2.getWidth()/binsize];
			for (Pair<StrandedPoint,StrandedPoint> pair : localreads) {
				int x = (pair.car().getLocation()-reg1.getStart())/binsize;
				int y = (pair.cdr().getLocation()-reg2.getStart())/binsize;
				char xstrand = pair.car().getStrand();
				char ystrand = pair.cdr().getStrand();
				for (int i = Math.max(x-readDist.size()/2, 0); i<Math.min(x+readDist.size()/2, tor.length); i++) {
					int distind = x-i+readDist.size()/2;
					if (xstrand=='-') distind = readDist.size() - distind;
					double pru = readDist.get(distind) * marg1[i];
					for (int j = Math.max(y-readDist.size()/2, 0); j<Math.min(y+readDist.size()/2, tor[i].length); j++) {
						distind = y-j+readDist.size()/2;
						if (ystrand=='-') distind = readDist.size() - distind;
						double prv = readDist.get(distind) * marg2[j];
						tor[i][j] += pru * prv;
					}
				}
			}
			return tor;
		} else {
			return null;
		}
	}
	
	public int[][] deconvolveCount(Region reg1, Region reg2, PrintStream jointinfo) {
		Region bigreg1 = reg1.expand(binsize*readDist.size()/2, binsize*readDist.size()/2);
		Region bigreg2 = reg2.expand(binsize*readDist.size()/2, binsize*readDist.size()/2);

		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg1);
		Set<Pair<StrandedPoint,StrandedPoint>> localreads = new HashSet<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (bigreg2.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}
		minmaxn = storage.getLeftMinMaxn(bigreg2);
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (bigreg1.contains(pairn.cdr())) {
				localreads.add(pairn.swap());
			}
		}
		if (jointinfo!=null) {
			jointinfo.println("read pairs considered: "+localreads.size());
		}
		if (!localreads.isEmpty()) {
			double[] marg1 = new double[reg1.getWidth()/binsize], marg2 = new double[reg2.getWidth()/binsize];
			BigWigIterator bwiter = bbfr.getBigWigIterator("chr"+reg1.getChrom(), reg1.getStart(), "chr"+reg1.getChrom(), reg1.getEnd(), false);
			for (int i=0; i<marg1.length; i++) {
				marg1[i] = bwiter.next().getWigValue();
			}
			bwiter = bbfr.getBigWigIterator("chr"+reg2.getChrom(), reg2.getStart(), "chr"+reg2.getChrom(), reg2.getEnd(), false);
			for (int i=0; i<marg2.length; i++) {
				if (bwiter.hasNext()) {
					marg2[i] = bwiter.next().getWigValue();
				}
			}
			int[][] tor = new int[reg1.getWidth()/binsize][reg2.getWidth()/binsize];
			for (Pair<StrandedPoint,StrandedPoint> pair : localreads) {
				int x = (pair.car().getLocation()-reg1.getStart())/binsize;
				int y = (pair.cdr().getLocation()-reg2.getStart())/binsize;
				char xstrand = pair.car().getStrand();
				char ystrand = pair.cdr().getStrand();
				for (int i = Math.max(x-readDist.size()/2, 0); i<Math.min(x+readDist.size()/2, tor.length); i++) {
					int distind = x-i+readDist.size()/2;
					if (xstrand=='-') distind = readDist.size() - distind;
					//double pru = readDist.get(distind) * marg1[i];
					for (int j = Math.max(y-readDist.size()/2, 0); j<Math.min(y+readDist.size()/2, tor[i].length); j++) {
						distind = y-j+readDist.size()/2;
						if (ystrand=='-') distind = readDist.size() - distind;
						//double prv = readDist.get(distind) * marg2[j];
						tor[i][j] ++;
					}
				}
			}
			return tor;
		} else {
			return null;
		}
	}
	
	public int deconvolveNumSamps(Region reg1, Region reg2, PrintStream jointinfo) {
		Region bigreg1 = reg1.expand(binsize*readDist.size()/2, binsize*readDist.size()/2);
		Region bigreg2 = reg2.expand(binsize*readDist.size()/2, binsize*readDist.size()/2);

		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg1);
		Set<Pair<StrandedPoint,StrandedPoint>> localreads = new HashSet<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (bigreg2.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}
		minmaxn = storage.getLeftMinMaxn(bigreg2);
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (bigreg1.contains(pairn.cdr())) {
				localreads.add(pairn.swap());
			}
		}
		if (jointinfo!=null) {
			jointinfo.println("read pairs considered: "+localreads.size());
		}
		return localreads.size();
	}

	public double getMarg(Region reg1) {
		double marg1 = 0;
		BigWigIterator bwiter = bbfr.getBigWigIterator("chr"+reg1.getChrom(), reg1.getStart(), "chr"+reg1.getChrom(), reg1.getEnd(), false);
		for (int i=0; i<reg1.getWidth()/binsize; i++) {
			marg1 += bwiter.next().getWigValue();
		}
		marg1 /= (double)(reg1.getWidth()/binsize);
		return marg1;
	}

	public double[] deconvolveStrip(Region reg1, Region reg2) {
		Region bigreg1 = reg1.expand(binsize*readDist.size()/2, binsize*readDist.size()/2);
		Region bigreg2 = reg2.expand(binsize*readDist.size()/2, binsize*readDist.size()/2);
		double[] tor = new double[reg2.getWidth()/binsize];
		double marg1 = 0;
		double[] marg2 = new double[reg2.getWidth()/binsize];
		BigWigIterator bwiter = bbfr.getBigWigIterator("chr"+reg1.getChrom(), reg1.getStart(), "chr"+reg1.getChrom(), reg1.getEnd(), false);
		for (int i=0; i<reg1.getWidth()/binsize; i++) {
			marg1 += bwiter.next().getWigValue();
		}
		marg1 /= (double)(reg1.getWidth()/binsize);
		bwiter = bbfr.getBigWigIterator("chr"+reg2.getChrom(), reg2.getStart(), "chr"+reg2.getChrom(), reg2.getEnd(), false);
		for (int i=0; i<marg2.length; i++) {
			marg2[i] = bwiter.next().getWigValue();
		}

		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg1);
		Set<Pair<StrandedPoint,StrandedPoint>> localreads = new HashSet<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (bigreg2.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}
		minmaxn = storage.getLeftMinMaxn(bigreg2);
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (bigreg1.contains(pairn.cdr())) {
				localreads.add(pairn.swap());
			}
		}

		int reg1mid = (reg1.getEnd()-reg1.getStart())/(2*binsize);
		for (Pair<StrandedPoint,StrandedPoint> pair : localreads) {
			int x = (pair.car().getLocation()-reg1.getStart())/binsize;
			int y = (pair.cdr().getLocation()-reg2.getStart())/binsize;
			char xstrand = pair.car().getStrand();
			char ystrand = pair.cdr().getStrand();
			int distind = x-reg1mid+readDist.size()/2;
			if (xstrand=='-') distind = readDist.size() - distind;
			double pru = readDist.get(distind) * marg1;
			for (int j = Math.max(y-readDist.size()/2, 0); j<Math.min(y+readDist.size()/2, tor.length); j++) {
				distind = y-j+readDist.size()/2;
				if (ystrand=='-') distind = readDist.size() - distind;
				double prv = readDist.get(distind) * marg2[j];
				tor[j] += pru * prv;
			}
		}
		return tor;
	}

}
