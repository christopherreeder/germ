package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import org.broad.igv.bbfile.BBFileReader;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.ChromosomeGenerator;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class InteractionSpecificity {

	private static DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		Genome g = SproutUtils.parseGenome(args);
		int binsize = Args.parseInteger(args, "binsize", 0);
		String margfile = Args.parseString(args, "margfile", "");
		String pointfile = Args.parseString(args, "pointfile", "");
		String tssfile = Args.parseString(args, "tssfile", "");
		String outfile = Args.parseString(args, "outfile", "");
		PrintStream out = new PrintStream(outfile);
		String jointsumfile = Args.parseString(args, "jointsumfile", "");
		PrintStream jointsumout = new PrintStream(jointsumfile);
		double interthresh = Args.parseDouble(args, "interthresh", 0);
		double backval = Args.parseDouble(args, "backval", 0);
		int numsamps = Args.parseInteger(args, "numsamps", 0);
		int runnum = Args.parseInteger(args, "runnum", 0);
		int perrun = Args.parseInteger(args, "perrun", 0);
		boolean allchrom = Args.parseFlags(args).contains("allchrom");

		int tstart = (runnum-1)*perrun;
		int tend = runnum*perrun;

		BBFileReader bbfr = new BBFileReader(margfile);
		ApproxInterDeconvolve aid = new ApproxInterDeconvolve(args);
		TreeMap<Point,String> tssmap = new TreeMap<Point,String>();
		TreeMap<Point,String> points = new TreeMap<Point,String>();

		BufferedReader r;
		String s;
		String[] split;

		int tindex = 0;
		r = new BufferedReader(new FileReader(tssfile));
		while ((s = r.readLine()) != null) {
			if (tindex>=tstart && tindex<tend) {
				split = s.split("\t");
				Point tmp = Point.fromString(g, split[0]);
				if (tmp!=null) {
					tssmap.put(tmp, s);
				} else {
					System.err.println(s);
				}
			}
			tindex++;
		}

		if (!pointfile.equals("")) {
			r = new BufferedReader(new FileReader(pointfile));
			while ((s = r.readLine()) != null) {
				split = s.split("\t");
				Point tmp = Point.fromString(g, split[0]);
				if (tmp!=null) {
					points.put(tmp, s);
				} else {
					System.err.println(s);
				}
			}
		}
		int index = 0;
		int inc = (tssmap.size() / 1000) + 1;
		for (Point tss : tssmap.keySet()) {
			Region tssreg = new Region(g, tss.getChrom(), tss.getLocation(), tss.getLocation()+binsize);
			double jointsum = 0;
			try {
				Iterator<Region> chromiter;
				if (allchrom) {
					ChromosomeGenerator chromgen = new ChromosomeGenerator();
					chromiter = chromgen.execute(g);
				} else {
					List<Region> chromreg = new ArrayList<Region>();
					chromreg.add(SproutUtils.chromRegion(tssreg.getChrom(), g));
					chromiter = chromreg.iterator();
				}
				//Region chromreg = SproutUtils.chromRegion(tssreg.getChrom(), g);
				while (chromiter.hasNext()) {
					Region chromreg = chromiter.next();
					double[][] joint = aid.deconvolve(tssreg,chromreg,null);
					int[][] counts = aid.deconvolveCount(tssreg, chromreg, null); //NEED THIS? YES
					int localsamps = 0; //aid.deconvolveNumSamps(tssreg, chromreg, null); //NEED THIS? NO
					//double backprob = 1d / ((double)joint[0].length);
					int start = 0;
					boolean inreg = false;
					int maxcount = 0;
					double density = 0;
					if (joint!=null) {
						int backsize = 0; //NEED THIS?
						for (int i=0; i<joint[0].length; i++) {
							jointsum += joint[0][i];
							/*
						if (joint[0][i]>interthresh) {
							backsize++;
						}
							 */
						}
						//double backprob = 1d / ((double)backsize);
						for (int i=0; i<joint[0].length; i++) {
							if (!inreg && joint[0][i]>interthresh) {
								inreg = true;
								maxcount = counts[0][i];
								start = i*binsize;
								density = joint[0][i];
							} else if (inreg && joint[0][i]<interthresh) {
								inreg = false;
								int end = i*binsize;
								Point startp = new Point(g, tss.getChrom(), start);
								Point endp = new Point(g, tss.getChrom(), end);
								SortedMap<Point,String> subpoint = points.subMap(startp, endp);
								Region tmpreg = new Region(g, chromreg.getChrom(), start, end);
								double tmpback = 0; //((double)((end-start)/binsize))*backprob;
								if (tmpback==1) tmpback = 0.99d;
								//density /= jointsum;
								//double pval = ExtractDomains.binomialDifference(density, tmpback, localsamps);
								out.println(tssmap.get(tss)+"\t"+tmpreg+"\t"+(end-start)+"\t"+subpoint.size()+"\t"+maxcount+"\t"+density+"\t"+jointsum+"\t"+tmpback+"\t"+localsamps);
							} else if (inreg) {
								if (counts[0][i]>maxcount) {
									maxcount = counts[0][i];
								}
								density += joint[0][i];
							}
						}
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
				System.err.println("Bad Chromosome: "+tssreg);
			}
			if (index % inc == 0) {
				System.err.println(index+" interactions considered "+dfm.format(new Date()));
			}
			index++;
			jointsumout.println(tss+"\t"+jointsum);
		}
		out.flush();
		out.close();
		jointsumout.flush();
		jointsumout.close();
	}

}
