package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;

public class DeterminePseudo {

	private static DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		Genome g = SproutUtils.parseGenome(args);
		int numsamps = Args.parseInteger(args, "numsamps", -1);
		double pcutoff = Args.parseDouble(args, "pcutoff", 0.05);
		double fdrcutoff = Args.parseDouble(args, "fdrcutoff", 0.1);
		int binsize = Args.parseInteger(args, "binsize", 0);
		String margfile = Args.parseString(args, "margfile", "");
		String infile = Args.parseString(args, "infile", "");
		String jointsumfile = Args.parseString(args, "jointsumfile", "");
		String auxfile = Args.parseString(args, "auxfile", "");
		double genomebinsize = Args.parseDouble(args, "genomebinsize", 0);
		double minmult = Args.parseDouble(args, "minmult", 0);
		double maxmult = Args.parseDouble(args, "maxmult", 0);
		double multtries = Args.parseDouble(args, "multtries", 0);
		double startinc = Args.parseDouble(args, "startinc", 0);
		double mininc = Args.parseDouble(args, "mininc", 0);
		boolean inhasp = Args.parseFlags(args).contains("inhasp");
		String outfile = Args.parseString(args, "outfile", "");
		PrintStream out = new PrintStream(outfile);
		String specoutfile = Args.parseString(args, "specoutfile", "");
		PrintStream specout = new PrintStream(specoutfile);

		System.err.println(margfile);

		Map<Point,Double> jointsummap = new HashMap<Point,Double>();
		Map<String,String> idmap = new HashMap<String,String>();
		List<String> namelist = new ArrayList<String>();
		List<Double> totaldensitylist = new ArrayList<Double>();
		List<Double> marglist = new ArrayList<Double>();
		BBFileReader bbfr = new BBFileReader(margfile);
		BigWigIterator bwiter;
		double[] marg = new double[0];
		String currentchrom = "";
		BufferedReader r;
		String s;
		String[] split;
		Set<String> seen = new HashSet<String>();
		List<Result> list = new ArrayList<Result>();
		List<String> speclist = new ArrayList<String>();
		int allcount = 0;
		double pseudo = 0;//9.5e-21;
		double norm = 3.03e-11;
		double genomepseudo = genomebinsize*pseudo;
		if (!jointsumfile.equals("")) {
			r = new BufferedReader(new FileReader(jointsumfile));
			while ((s = r.readLine()) != null) {
				split = s.split("\t");
				Point p = Point.fromString(g, split[0]);
				Double d = Double.valueOf(split[1]);
				jointsummap.put(p, d);
			}
		}
		if (auxfile.equals("")) {
			r = new BufferedReader(new FileReader(infile));
			while ((s = r.readLine()) != null) {
				split = s.split("\t");
				String tmps = "";
				int tmpto = inhasp ? split.length-1 : split.length;
				for (int i=0; i<tmpto; i++) {
					tmps += split[i]+"\t";
				}
				speclist.add(tmps);
				Region reg = Region.fromString(g, split[3]);
				Point p = Point.fromString(g, split[0]);
				if (!p.getChrom().equals(currentchrom)) {
					currentchrom = p.getChrom();
					Region chromreg = SproutUtils.chromRegion(currentchrom, g);
					bwiter = bbfr.getBigWigIterator("chr"+chromreg.getChrom(), chromreg.getStart(), "chr"+chromreg.getChrom(), chromreg.getEnd(), false);
					marg = new double[chromreg.getWidth()/binsize];
					while (bwiter.hasNext()) {
						WigItem wi = bwiter.next();
						int loc = wi.getStartBase();
						marg[loc/binsize] = wi.getWigValue();
					}
				}
				Result tmp = new Result();
				if (p.getChrom().equals(reg.getChrom())) {
					tmp.distance = reg.distance(p);
				} else {
					tmp.distance = -1; 
				}
				tmp.point = p;
				double regsum = 0;
				if (reg.getChrom().equals(p.getChrom())) {
					int tmpbinstart = reg.getStart()/binsize;
					int tmpbinend = reg.getEnd()/binsize + 1;
					if (tmpbinend>=marg.length) {
						System.err.println("reg too long: "+s);
					}
					for (int mloc=tmpbinstart; mloc<Math.min(tmpbinend,marg.length); mloc++) {
						regsum += marg[mloc];
					}
				} else {
					bwiter = bbfr.getBigWigIterator("chr"+reg.getChrom(), reg.getStart(), "chr"+reg.getChrom(), reg.getEnd(), false);
					while (bwiter.hasNext()) {
						WigItem wi = bwiter.next();
						regsum += wi.getWigValue();
					}
				}
				tmp.back = regsum;
				tmp.marg = marg[p.getLocation()/binsize];
				tmp.count = Integer.parseInt(split[6]);
				tmp.totalcount = Integer.parseInt(split[10]);
				tmp.density = Double.parseDouble(split[7]);
				if (!jointsumfile.equals("")) {
					tmp.totaldensity = jointsummap.get(p);
				} else {
					tmp.totaldensity = Double.parseDouble(split[8]);
				}
				tmp.binwidth = Double.parseDouble(split[4]) / binsize;
				tmp.norm = tmp.totaldensity / tmp.marg;
				list.add(tmp);
				if (!seen.contains(split[0])) {
					seen.add(split[0]);
					allcount += tmp.totalcount;
					namelist.add(split[0]);
					marglist.add(tmp.marg);
					totaldensitylist.add(tmp.totaldensity);
				}
			}
		} else {
			r = new BufferedReader(new FileReader(infile));
			BufferedReader ra = new BufferedReader(new FileReader(auxfile));
			String sa = "";
			String[] splita = new String[0];
			while ((s = r.readLine()) != null) {
				sa = ra.readLine();
				split = s.split("\t");
				splita = sa.split("\t");
				if (!(idmap.containsKey(split[0]) && !idmap.get(split[0]).equals(split[1]))) {
					if (!idmap.containsKey(split[0])) {
						idmap.put(split[0], split[1]);
					}
					String tmps = "";
					for (int i=0; i<split.length-1; i++) {
						tmps += split[i]+"\t";
					}
					speclist.add(tmps);
					Region reg = Region.fromString(g, split[3]);
					Point p = Point.fromString(g, split[0]);
					Result tmp = new Result();
					tmp.distance = reg.distance(p);
					tmp.point = p;
					int tmpbinstart = reg.getStart()/binsize;
					int tmpbinend = reg.getEnd()/binsize + 1;
					tmp.back = Double.parseDouble(splita[2]);
					tmp.marg = Double.parseDouble(splita[3]);
					tmp.count = Integer.parseInt(split[6]);
					tmp.totalcount = Integer.parseInt(split[10]);
					tmp.density = Double.parseDouble(split[7]);
					tmp.totaldensity = Double.parseDouble(split[8]);
					tmp.binwidth = Double.parseDouble(split[4]) / binsize;
					tmp.norm = tmp.totaldensity / tmp.marg;
					list.add(tmp);
					if (!seen.contains(split[0])) {
						seen.add(split[0]);
						allcount += tmp.totalcount;
						namelist.add(split[0]);
						marglist.add(tmp.marg);
						totaldensitylist.add(tmp.totaldensity);
					}
				}
			}
		}
		System.err.println("files read "+dfm.format(new Date()));

		int maxpos = -1;
		double maxdens = 0;
		List<Double> d = new ArrayList<Double>();
		for (int i=0; i<marglist.size(); i++) {
			d.add(0d);
			if (totaldensitylist.get(i)>maxdens) {
				maxdens = totaldensitylist.get(i);
				maxpos = i;
			}
		}
		double origmaxdens = maxdens;
		double closeness = Double.MAX_VALUE;
		double closemult = -1;
		if (numsamps!=-1) allcount = numsamps;

		double inc = (maxmult-minmult)/multtries;
		if (inc==0) inc = 1;
		if (startinc!=0) inc = startinc;
		for (double mult=minmult; mult<=maxmult; mult += inc) {
			d = new ArrayList<Double>();
			for (int i=0; i<marglist.size(); i++) {
				d.add(0d);
			}
			maxdens = mult*origmaxdens;
			d.set(maxpos, maxdens);

			double oldabsdiff = Double.MAX_VALUE;
			double absdiff = Double.MAX_VALUE/2d;
			while (oldabsdiff>absdiff) {
				double ysum = 0;
				double xdsum = 0;
				for (int i=0; i<marglist.size(); i++) {
					ysum += marglist.get(i);
					xdsum += totaldensitylist.get(i)+d.get(i);
				}
				for (int i=0; i<marglist.size(); i++) {
					xdsum -= totaldensitylist.get(i) + d.get(i);
					double yminusi = ysum - marglist.get(i);
					double tmp = (marglist.get(i)*xdsum/yminusi)-totaldensitylist.get(i);
					if (i==maxpos) {
						d.set(i, Math.max(maxdens, tmp));
					} else {
						d.set(i, Math.max(0, tmp));
					}
					xdsum += totaldensitylist.get(i) + d.get(i);
				}
				oldabsdiff = absdiff;
				absdiff = 0;
				for (int i=0; i<marglist.size(); i++) {
					double xd = (totaldensitylist.get(i) + d.get(i)) / xdsum;
					if (xd<0 || xd>1) System.err.println("improper x prob: "+xd);
					double yhat = marglist.get(i) / ysum;
					if (yhat<0 || yhat>1) System.err.println("improper y prob"+yhat);
					absdiff += Math.abs(xd-yhat);
				}
				//System.err.println(absdiff);
			}


			Map<String,Double> xdmap = new HashMap<String,Double>();
			for (int i=0; i<marglist.size(); i++) {
				xdmap.put(namelist.get(i), totaldensitylist.get(i)+d.get(i));
			}

			int totalsig = 0;
			int onesig = 0;
			for (int i=0; i<list.size(); i++) {
				Result tmp = list.get(i);
				//double signalprob = (tmp.density + tmp.binwidth*pseudo) / (tmp.totaldensity + genomepseudo);
				//double signalprob = tmp.density / norm;
				//double signalprob = (Math.pow(10,0.35*Math.log10(tmp.density)-1.5341)) / tmp.marg;
				double signalnorm = xdmap.get(tmp.point.toString());
				double signalprob = tmp.density / signalnorm;
				//double backprob = tmp.binwidth / genomebinsize;
				double backprob = tmp.back;
				try {
					double pval = ExtractDomains.fastBinomialDifference(signalprob, backprob, allcount);
					//System.err.println(signalprob+"\t"+backprob+"\t"+tmp.totalcount+"\t"+pval);
					//out.println(tmp.count+"\t"+signalprob+"\t"+backprob+"\t"+tmp.marg+"\t"+tmp.density+"\t"+signalnorm+"\t"+pval);
					if (pval<=pcutoff) {
						totalsig++;
						//System.err.println(tmp.count+"\t"+signalprob+"\t"+backprob+"\t"+tmp.marg+"\t"+pval);
						if (tmp.count==1) {

							onesig++;
						}
					}
				} catch (Exception e) {
					//out.println(tmp.count+"\t"+signalprob+"\t"+backprob+"\t"+tmp.marg+"\t"+tmp.density+"\t"+signalnorm+"\t"+(-1));
					//System.err.println(tmp.count+"\t"+signalprob+"\t"+backprob+"\t"+tmp.marg+"\t"+tmp.density+"\t"+tmp.totaldensity+"\t"+tmp.point);
				}
			}
			double fdr = ((double)onesig)/((double)totalsig);
			System.err.println(onesig+"\t"+totalsig+"\t"+fdr+"\t"+mult+"\t"+dfm.format(new Date()));
			double tmp = Math.abs(fdr-fdrcutoff);
			if (tmp<closeness) {
				closeness = tmp;
				closemult = mult;
			}
			if (startinc==0 && fdr<fdrcutoff) break;
			if (startinc!=0 && fdr<fdrcutoff) {
				mult -= inc;
				if (inc>mininc) {
					inc /= 10;
				} else {
					break;
				}
			}
		}

		System.err.println("closemult: "+closemult);

		d = new ArrayList<Double>();
		for (int i=0; i<marglist.size(); i++) {
			d.add(0d);
		}
		maxdens = closemult*origmaxdens;
		d.set(maxpos, maxdens);
		double oldabsdiff = Double.MAX_VALUE;
		double absdiff = Double.MAX_VALUE/2d;
		while (oldabsdiff>absdiff) {
			double ysum = 0;
			double xdsum = 0;
			for (int i=0; i<marglist.size(); i++) {
				ysum += marglist.get(i);
				xdsum += totaldensitylist.get(i)+d.get(i);
			}
			for (int i=0; i<marglist.size(); i++) {
				xdsum -= totaldensitylist.get(i) + d.get(i);
				double yminusi = ysum - marglist.get(i);
				double tmp = (marglist.get(i)*xdsum/yminusi)-totaldensitylist.get(i);
				if (i==maxpos) {
					d.set(i, Math.max(maxdens, tmp));
				} else {
					d.set(i, Math.max(0, tmp));
				}
				xdsum += totaldensitylist.get(i) + d.get(i);
			}
			oldabsdiff = absdiff;
			absdiff = 0;
			for (int i=0; i<marglist.size(); i++) {
				double xd = (totaldensitylist.get(i) + d.get(i)) / xdsum;
				if (xd<0 || xd>1) System.err.println("improper x prob: "+xd);
				double yhat = marglist.get(i) / ysum;
				if (yhat<0 || yhat>1) System.err.println("improper y prob"+yhat);
				absdiff += Math.abs(xd-yhat);
			}
			//System.err.println(absdiff);
		}
		Map<String,Double> xdmap = new HashMap<String,Double>();
		for (int i=0; i<marglist.size(); i++) {
			xdmap.put(namelist.get(i), totaldensitylist.get(i)+d.get(i));
		}

		int totalsig = 0;
		int onesig = 0;
		for (int i=0; i<list.size(); i++) {
			Result tmp = list.get(i);
			//double signalprob = (tmp.density + tmp.binwidth*pseudo) / (tmp.totaldensity + genomepseudo);
			//double signalprob = tmp.density / norm;
			//double signalprob = (Math.pow(10,0.35*Math.log10(tmp.density)-1.5341)) / tmp.marg;
			double signalnorm = xdmap.get(tmp.point.toString());
			double signalprob = tmp.density / signalnorm;
			//double backprob = tmp.binwidth / genomebinsize;
			double backprob = tmp.back;
			try {
				double pval = ExtractDomains.fastBinomialDifference(signalprob, backprob, allcount);
				//System.err.println(signalprob+"\t"+backprob+"\t"+tmp.totalcount+"\t"+pval);
				out.println(tmp.count+"\t"+signalprob+"\t"+backprob+"\t"+tmp.marg+"\t"+tmp.density+"\t"+signalnorm+"\t"+tmp.distance+"\t"+pval);
				specout.println(speclist.get(i)+pval);
				if (pval<=pcutoff) {
					totalsig++;
					//System.err.println(tmp.count+"\t"+signalprob+"\t"+backprob+"\t"+tmp.marg+"\t"+pval);
					if (tmp.count==1) {

						onesig++;
					}
				}
			} catch (Exception e) {
				out.println(tmp.count+"\t"+signalprob+"\t"+backprob+"\t"+tmp.marg+"\t"+tmp.density+"\t"+signalnorm+"\t"+tmp.distance+"\t"+(-1));
				specout.println(speclist.get(i)+(-1));
				//System.err.println(tmp.count+"\t"+signalprob+"\t"+backprob+"\t"+tmp.marg+"\t"+tmp.density+"\t"+tmp.totaldensity+"\t"+tmp.point);
			}
		}

		/*
		for (int i=0; i<marglist.size(); i++) {
			out.println(namelist.get(i)+"\t"+marglist.get(i)+"\t"+totaldensitylist.get(i)+"\t"+d.get(i));
		}
		 */

		out.flush();
		out.close();
		specout.flush();
		specout.close();
	}

	private static class Result {

		Point point;
		int count, totalcount, distance;
		double density, totaldensity, binwidth, marg, back, norm;

	}

}
