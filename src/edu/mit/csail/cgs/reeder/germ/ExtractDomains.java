package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.List;
import java.util.Random;

import cern.jet.random.Binomial;
import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class ExtractDomains {

	private static DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		/*
		Genome g = SproutUtils.parseGenome(args);
		String resultfile = Args.parseString(args, "resultfile", "");
		String outfile = Args.parseString(args, "outfile", "");
		PrintStream out = new PrintStream(outfile);
		List<Result> results = fromFile(resultfile,g);
		addQVals(results);
		for (int i=0; i<results.size(); i++) {
			out.println(results.get(i).toBroadPeakString());
		}
		out.flush();
		out.close();
		 */

		Genome g = SproutUtils.parseGenome(args);
		String margfile = Args.parseString(args, "margfile", "");
		double thresh = Args.parseDouble(args, "thresh", 0);
		int binsize = Args.parseInteger(args, "binsize", 0);
		int numsamps = Args.parseInteger(args, "numsamps", 0);
		double backval = Args.parseDouble(args, "backval", 0);
		String outfile = Args.parseString(args, "outfile", "");
		int numdomains = Args.parseInteger(args, "numdomains", 0);
		boolean slow = Args.parseFlags(args).contains("slow");
		boolean noqvals = Args.parseFlags(args).contains("noqvals");
		PrintStream out = new PrintStream(outfile);

		List<Region> regions = new ArrayList<Region>();
		List<Double> marg = new ArrayList<Double>();
		List<Result> results = new ArrayList<Result>();
		Random rand = new Random();

		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(margfile));
		s = r.readLine();
		split = s.split("\t");
		if (split.length==1) {
			split = s.split(" ");
		}
		split = split[1].split("=");
		String chrom = split[1].substring(3);
		String nextchrom = chrom;
		while ((s = r.readLine()) != null) {
			if (!s.startsWith("fixedStep")) {
				double tmp = Double.parseDouble(s);
				marg.add(tmp);
			} else {
				split = s.split("\t");
				if (split.length==1) {
					split = s.split(" ");
				}
				split = split[1].split("=");
				nextchrom = split[1].substring(3);
				
				regions = ThresholdFinder.extractDomains(marg, thresh, binsize, g, chrom);
				int inc = regions.size() / 100 + 1;
				for (int i=0; i<regions.size(); i++) {
					Region tmp = regions.get(i);
					int binstart = tmp.getStart()/binsize;
					int binend = tmp.getEnd()/binsize;
					double total = 0;
					double maxval = 0;
					int maxpos = -1;
					for (int j=binstart; j<=binend; j++) {
						double tmpval = marg.get(j);
						if (tmpval>maxval) {
							maxval = tmpval;
							maxpos = j*binsize + binsize/2;
						}
						total += tmpval;
					}
					double tmpback = (double)(binend-binstart+1)*backval;
					double pval;
					if (slow) {
						pval = binomialDifference(total, tmpback, numsamps);
					} else {
						pval = fastBinomialDifference(total, tmpback, numsamps);
					}
					//double pval = rand.nextDouble();
					Result result = new Result();
					result.region = tmp;
					result.pval = pval;
					result.tpval = -Math.log10(pval);
					result.name = "peak"+results.size();
					result.maxval = maxval;
					result.maxpos = maxpos;
					result.avgval = total / ((double)(binend-binstart+1));
					result.totalval = total;
					result.strand = '.';
					results.add(result);
					//out.println("chr"+result.region.getChrom()+"\t"+result.region.getStart()+"\t"+result.region.getEnd()+"\tpeak"+i+"\t"+result.maxval+"\t.\t"+result.avgval+"\t"+(-Math.log10(result.pval))+"\t"+(-1));
					if (i % inc == 0) {
						System.err.println("chr"+chrom+" "+i+" domains "+dfm.format(new Date()));
					}
				}
				chrom = nextchrom;
				marg = new ArrayList<Double>();
				if (numdomains>0 && results.size()>numdomains) {
					break;
				}
			}
		}
		
		
		if (!noqvals) {
			addQVals(results);
		} else {
			Collections.sort(results,new PvalComparator());
		}
		for (int i=0; i<results.size(); i++) {
			out.println(results.get(i).toBroadPeakString());
		}
		out.flush();
		out.close();

	}
	
	/*
	 * p1=signal prob
	 * p2=background prob
	 */
	public static double binomialDifference(double p1, double p2, int numsamps) throws Exception {
		try {
			double tor = 0;
			DRand rand = new DRand(new Date());
			Binomial bin1 = new Binomial(numsamps,p1,rand);
			Binomial bin2 = new Binomial(numsamps,p2,rand);
			for (int i=0; i<=numsamps; i++) {
				if (i==0) {
					tor += bin1.pdf(i);
				} else {
					tor += bin1.pdf(i)*(1-bin2.cdf(i-1));
				}
			}
			return tor;
		} catch (Exception e) {
			System.err.println(p1+"\t"+p2+"\t"+numsamps);
			throw e;
		}
	}
	
	/*
	 * p1=signal prob
	 * p2=background prob
	 */
	public static double fastBinomialDifference(double p1, double p2, int numsamps) {
		double tor = 0;
		double mean1 = ((double)numsamps)*p1;
		double mean2 = ((double)numsamps)*p2;
		double sd1 = Math.sqrt(mean1*(1d-p1));
		double sd2 = Math.sqrt(mean2*(1d-p2));
		int start1 = Math.max((int)(mean1-4d*sd1),0);
		int end1 = Math.min((int)(mean1+4d*sd1), numsamps);
		int start2 = Math.max((int)(mean2-4d*sd2), 0);
		int end2 = Math.min((int)(mean2+4d*sd2), numsamps);
		DRand rand = new DRand(new Date());
		Binomial bin1 = new Binomial(numsamps,p1,rand);
		Binomial bin2 = new Binomial(numsamps,p2,rand);
		for (int i=Math.min(start1,start2); i<=Math.max(end1,end2); i++) {
			if (i==0) {
				tor += bin1.pdf(i);
			} else {
				tor += bin1.pdf(i)*(1-bin2.cdf(i-1));
			}
		}
		return tor;
	}

	public static List<Result> fromFile(String resultfile, Genome g) throws IOException {
		List<Result> results = new ArrayList<Result>();
		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(resultfile));
		r.readLine();
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			String chrom = split[0].substring(3);
			int start = Integer.parseInt(split[1]);
			int end = Integer.parseInt(split[2]);
			Region reg = new Region(g, chrom, start, end);
			String name = split[3];
			double maxval = Double.parseDouble(split[4]);
			char strand = split[5].charAt(0);
			double avgval = Double.parseDouble(split[6]);
			double pval = Double.parseDouble(split[7]);
			double qval = Double.parseDouble(split[8]);
			results.add(new Result(reg, pval, qval, maxval, avgval, name, strand));
		}
		return results;
	}

	public static void addQVals(List<Result> results) {
		Collections.sort(results,new PvalComparator());
		double m = results.size();
		int j = 0;
		double[] pi = new double[39];
		int piindex = 0;
		for (double lam = 0.01; lam<=0.4; lam += 0.01) {
			while (j<results.size() && results.get(j).pval<lam) j++;
			double numer = results.size()-j;
			double denom = m*(1d-lam);
			pi[piindex++] = numer/denom;
		}
		double mindiff = Double.MAX_VALUE;
		int mindiffind = -1;
		for (int i=0; i<pi.length-1; i++) {
			double diff = pi[i]-pi[i+1];
			if (diff<mindiff) {
				mindiff = diff;
				mindiffind = i;
			}
		}
		System.err.println("pi: "+pi[mindiffind]);
		double pihat = pi[mindiffind];
		Result tmp = results.get(results.size()-1);
		tmp.setQVal(pihat*tmp.pval);
		for (int i=results.size()-2; i>=0; i--) {
			tmp = results.get(i);
			Result nexttmp = results.get(i+1);
			tmp.setQVal(Math.min(pihat*m*tmp.pval/((double)i),nexttmp.qval));
		}
	}

	public static class PvalComparator implements Comparator<Result> {

		public int compare(Result arg0, Result arg1) {
			return Double.compare(arg0.pval, arg1.pval);
		}

	}

	public static class Result {

		Region region;
		double pval;
		double tpval;
		double qval;
		double tqval = -1;
		double maxval;
		int maxpos;
		double avgval;
		double totalval;
		String name;
		char strand;

		public Result() {

		}

		public Result(Region region, double tpval, double qval, double maxval, double avgval, String name, char strand) {
			this.region = region;
			this.tpval = tpval;
			this.pval = Math.pow(10d, -tpval);
			this.qval = qval;
			this.maxval = maxval;
			this.avgval = avgval;
			this.name = name;
			this.strand = strand;
		}
		
		public Result(Region region, double maxval, double avgval, String name, char strand) {
			this.region = region;
			this.maxval = maxval;
			this.avgval = avgval;
			this.name = name;
			this.strand = strand;
		}

		public void setQVal(double qval) {
			this.qval = qval;
			this.tqval = -Math.log10(qval);
		}
		
		public void setTQVal(double tqval) {
			this.qval = Math.pow(10d, -tqval);
			this.tqval = tqval;
		}
		
		public void setPVal(double pval) {
			this.pval = pval;
			this.tpval = -Math.log10(pval);
		}
		
		public void setTPVal(double tpval) {
			this.pval = Math.pow(10d, -tpval);
			this.tpval = tpval;
		}

		public String toBroadPeakString() {
			return "chr"+region.getChrom()+"\t"+region.getStart()+"\t"+region.getEnd()+"\t"+name+"\t"+maxval+"\t"+strand+"\t"+avgval+"\t"+tpval+"\t"+tqval;
		}
		
		public static Result fromString(Genome g, String s) {
			String[] split = s.split("\t");
			String chrom = split[0].substring(3);
			int start = Integer.parseInt(split[1]);
			int end = Integer.parseInt(split[2]);
			String name = split[3];
			double maxval = Double.parseDouble(split[4]);
			char strand = split[5].charAt(0);
			double avgval = Double.parseDouble(split[6]);
			double tpval = Double.parseDouble(split[7]);
			double tqval = Double.parseDouble(split[8]);
			Result tmp = new Result(new Region(g, chrom, start, end), maxval, avgval, name, strand);
			tmp.setTPVal(tpval);
			tmp.setTQVal(tqval);
			return tmp;
		}
	}

}
