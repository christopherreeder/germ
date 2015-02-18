package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class SproutUtils {

	public static Region chromRegion(String chrom, Genome g) {
		return new Region(g, chrom, 1, g.getChromLength(chrom));
	}
	
	public static List<Pair<Point,Point>> readInters(String file, Region r1, Region r2, Genome g) throws IOException {
		List<Pair<Point,Point>> tor = new ArrayList<Pair<Point,Point>>();
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			Point tmp1 = Point.fromString(g, split[0]);
			Point tmp2 = Point.fromString(g, split[1]);
			if ((r1.contains(tmp1) && r2.contains(tmp2)) || (r1.contains(tmp2) && r2.contains(tmp1))) {
				tor.add(new Pair<Point,Point>(tmp1,tmp2));
			}
		}
		return tor;
	}
	
	public static List<Point> readEvents(String file, Region r1, Genome g) throws NumberFormatException, IOException {
		List<Point> tor = new ArrayList<Point>();
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			Point tmp1 = Point.fromString(g, split[0]);
			if (r1.contains(tmp1)) {
				tor.add(tmp1);
			}
		}
		return tor;
	}
	
	public static List<Pair<Point,Point>> readInters(String file, Genome g) throws IOException {
		List<Pair<Point,Point>> tor = new ArrayList<Pair<Point,Point>>();
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			Point tmp1 = Point.fromString(g, split[0]);
			Point tmp2 = Point.fromString(g, split[1]);
			tor.add(new Pair<Point,Point>(tmp1,tmp2));
		}
		return tor;
	}
	
	public static List<Double> readDoubleList(String file) throws IOException {
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		List<Double> tor = new ArrayList<Double>();
		while ((s = r.readLine()) != null) {
			tor.add(Double.parseDouble(s));
		}
		return tor;
	}
	
	public static List<Integer> readIntegerList(String file) throws IOException {
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		List<Integer> tor = new ArrayList<Integer>();
		while ((s = r.readLine()) != null) {
			tor.add(Integer.parseInt(s));
		}
		return tor;
	}
	
	public static List<Pair<String,String>> chrompairs(String[] args, Genome g) throws IOException {
		String pairfile = Args.parseString(args, "pairfile", "");
		String excludefile = Args.parseString(args, "excludefile", "");
		List<String> chromlist = g.getChromList();
		List<Pair<String,String>> chrompairs = new ArrayList<Pair<String,String>>();
		if (pairfile.equals("")) {
			Set<Pair<String,String>> exset = new HashSet<Pair<String,String>>();
			if (!excludefile.equals("")) {
				BufferedReader r = new BufferedReader(new FileReader(excludefile));
				String s;
				String[] split;
				while ((s = r.readLine()) != null) {
					split = s.split("\t");
					exset.add(new Pair<String,String>(split[0],split[1]));
				}
			}
			if (Args.parseFlags(args).contains("same")) {
				for (int i=0; i<chromlist.size(); i++) {
					String chrom1 = chromlist.get(i);
					Pair<String,String> tmp = new Pair<String,String>(chrom1,chrom1);
					if (!exset.contains(tmp)) {
						chrompairs.add(tmp);
					}
				}
			} else {
				for (int i=0; i<chromlist.size(); i++) {
					String chrom1 = chromlist.get(i);
					for (int j=i; j<chromlist.size(); j++) {
						String chrom2 = chromlist.get(j);
						Pair<String,String> tmp = new Pair<String,String>(chrom1,chrom2);
						if (!exset.contains(tmp)) {
							chrompairs.add(tmp);
						}
					}
				}
			}
		} else {
			BufferedReader r = new BufferedReader(new FileReader(pairfile));
			String s;
			String[] split;
			while ((s = r.readLine()) != null) {
				split = s.split("\t");
				chrompairs.add(new Pair<String,String>(split[0],split[1]));
			}
		}
		return chrompairs;
	}
	
	public static Genome parseGenome(String[] args) throws NotFoundException {
		Genome g;
		String genomefile = Args.parseString(args, "genomefile", "");
		System.err.println("genome: "+genomefile);
		if (!genomefile.equals("")) {
			String[] split = genomefile.split("/");
			if (split[split.length-1].endsWith("_1.txt")) {
				System.err.println("using dbids");
				g = new Genome(genomefile,new File(genomefile),false);
			} else {
				g = new Genome(genomefile,new File(genomefile),true);
			}
		} else {
			g = Args.parseGenome(args).cdr();
		}
		return g;
	}
	
	public static boolean usingDBIDS(String[] args) {
		String genomefile = Args.parseString(args, "genomefile", "");
		if (!genomefile.equals("")) {
			String[] split = genomefile.split("/");
			if (split[split.length-1].endsWith("_1.txt")) {
				return true;
			}
		}
		return false;
	}
	
	public static Region findEnclosingRegion(Point p, SortedSet<Region> set) {
		Region r = p.expand(0);
		SortedSet<Region> subset = set.headSet(r);
		if (!subset.isEmpty()) {
			Region last = subset.last();
			if (last.overlaps(r)) {
				return last;
			}
		}
		subset = set.tailSet(r);
		if (!subset.isEmpty()) {
			Region first = subset.first();
			if (first.overlaps(r)) {
				return first;
			}
		}
		return null;
	}
	
	public static long genomeArea(Genome g) {
		List<String> chromlist = g.getChromList();
		long all = 0;
		for (int i=0; i<chromlist.size(); i++) {
			int isize = g.getChromLength(chromlist.get(i));
			for (int j=i; j<chromlist.size(); j++) {
				int jsize = g.getChromLength(chromlist.get(j));
				all += isize*jsize;
			}
		}
		return all;
	}
	
}
