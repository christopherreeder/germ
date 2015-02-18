package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
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

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class ComputeEstimate {

	private static DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	private SproutStorage storage;
	private List<Double> selfliglik;
	private int binsize, h, readdistbinsize;
	private Genome g;
	private List<Double> readDist;

	public ComputeEstimate(Genome g, int h, int binsize, int storagesize, String readfile, String selfliglikfile) throws IOException {
		this.g = g;
		this.h = h;
		this.binsize = binsize;
		storage = SproutStorage.fromFile(g, readfile, storagesize);
		selfliglik = new ArrayList<Double>();
		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(selfliglikfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			selfliglik.add(Double.parseDouble(split[0]));
		}
	}

	public void readDistFromFile(String file, int readdistbinsize) throws IOException {
		this.readdistbinsize = readdistbinsize;
		readDist = new ArrayList<Double>();
		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(file));
		while ((s = r.readLine()) != null) {
			readDist.add(Double.parseDouble(s));
		}
		System.err.println("readDist size: "+readDist.size());
	}

	/**
	 * @param args
	 * @throws NotFoundException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NotFoundException, IOException {
		Genome g = SproutUtils.parseGenome(args);
		int upstream = Args.parseInteger(args, "upstream", 0);
		int downstream = Args.parseInteger(args, "downstream", 0);
		int h = Args.parseInteger(args, "h", 0);
		int binsize = Args.parseInteger(args, "binsize", 0);
		String readfile = Args.parseString(args, "readfile", "");
		String regionfile = Args.parseString(args, "regionfile", "");
		String selfliglikfile = Args.parseString(args, "selfliglikfile", "");
		int size = Args.parseInteger(args, "size", 0);
		String outfile = Args.parseString(args, "outfile", "");
		PrintStream out = new PrintStream(outfile);

		SproutStorage storage = SproutStorage.fromFile(g, readfile, size);
		List<StrandedRegion> regionlist = new ArrayList<StrandedRegion>();
		List<Double> selfliglik = new ArrayList<Double>();

		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(regionfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			StrandedRegion tmp = StrandedRegion.fromString(g, split[0]);
			if (tmp==null) {
				System.err.println(s);
			} else {
				regionlist.add(tmp.expand(upstream,downstream));
			}
		}

		r = new BufferedReader(new FileReader(selfliglikfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			selfliglik.add(Double.parseDouble(split[0]));
		}

		int maxlength = 0;
		int totallength = 0;
		for (StrandedRegion reg : regionlist) {
			int length = reg.getWidth();
			totallength += length;
			if (length>maxlength) {
				maxlength = length;
			}
		}

		double[][] tor = new double[totallength/binsize+1][maxlength/binsize+1];
		int numvbins = maxlength / binsize + 1;

		for (StrandedRegion reg : regionlist) {
			int numbins = reg.getWidth()/binsize + 1;
			StrandedRegion bigreg = reg.expand(3*h, 3*h);
			Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg);
			List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
			for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
				Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
				if (PairedReadDistribution.isSelfLigation(pairn) && bigreg.contains(pairn.cdr())) {
					localreads.add(pairn);
				}
			}

			int inc = Math.max(numbins/100,1);
			for (int i=0; i<numbins; i++) {
				StrandedPoint pi = new StrandedPoint(g,reg.getChrom(),reg.getStart()+i*binsize+binsize/2,'-');
				for (int j=i; j<numvbins; j++) {
					StrandedPoint pj = new StrandedPoint(g,reg.getChrom(),reg.getStart()+j*binsize+binsize/2,'+');
					Pair<StrandedPoint,StrandedPoint> pairi = new Pair<StrandedPoint,StrandedPoint>(pi,pj);
					for (int k=0; k<localreads.size(); k++) {
						Pair<StrandedPoint,StrandedPoint> pairj = localreads.get(k);
						tor[i][j] += selfliglik.get(pairj.car().distance(pairj.cdr()))*CrossValidation.normal(((double)pairi.car().distance(pairj.car()))/((double)h),((double)pairi.cdr().distance(pairj.cdr()))/((double)h),1d)/((double)h);
					}
				}
				if (i % (inc) == 0) {
					System.err.println("processed "+i+" bins for jack "+dfm.format(new Date()));
				}
			}
		}

		for (int i=0; i<tor.length; i++) {
			for (int j=0; j<tor[i].length; j++) {
				out.print(tor[i][j]+"\t");
			}
			out.println();
		}
		out.flush();
		out.close();
	}

	public double[][] computeEstimate(StrandedRegion reg) {
		int length = reg.getWidth();
		double[][] tor = new double[length/binsize+1][length/binsize+1];
		int numvbins = length / binsize + 1;

		int numbins = reg.getWidth()/binsize + 1;
		StrandedRegion bigreg = reg.expand(3*h, 3*h);
		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg);
		List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (PairedReadDistribution.isSelfLigation(pairn) && bigreg.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}

		int inc = Math.max(numbins/100,1);
		for (int i=0; i<numbins; i++) {
			StrandedPoint pi = new StrandedPoint(g,reg.getChrom(),reg.getStart()+i*binsize+binsize/2,'-');
			for (int j=i; j<numvbins; j++) {
				StrandedPoint pj = new StrandedPoint(g,reg.getChrom(),reg.getStart()+j*binsize+binsize/2,'+');
				Pair<StrandedPoint,StrandedPoint> pairi = new Pair<StrandedPoint,StrandedPoint>(pi,pj);
				for (int k=0; k<localreads.size(); k++) {
					Pair<StrandedPoint,StrandedPoint> pairj = localreads.get(k);
					tor[i][j] += selfliglik.get(pairj.car().distance(pairj.cdr()))*CrossValidation.normal(((double)pairi.car().distance(pairj.car()))/((double)h),((double)pairi.cdr().distance(pairj.cdr()))/((double)h),1d)/((double)h);
				}
			}
			if (i % (inc) == 0) {
				System.err.println("processed "+i+" bins for jack "+dfm.format(new Date()));
			}
		}

		return tor;
	}

	public double[][] computeEstimate(StrandedRegion reg, int height) {
		int length = reg.getWidth();
		double[][] tor = new double[length/binsize+1][height/binsize+1];
		int numvbins = tor[0].length;

		int numbins = tor.length;
		StrandedRegion bigreg = reg.expand(3*h, 3*h);
		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg);
		List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (PairedReadDistribution.isSelfLigation(pairn) && bigreg.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}

		int inc = Math.max(numbins/100,1);
		for (int i=0; i<numbins; i++) {
			//StrandedPoint pi = new StrandedPoint(g,reg.getChrom(),reg.getStart()+i*binsize+binsize/2,'-');
			for (int j=0; j<numvbins; j++) {
				//StrandedPoint pj = new StrandedPoint(g,reg.getChrom(),pi.getLocation()+j*binsize+binsize/2,'+');
				int midpoint = reg.getStart()+i*binsize+binsize/2;
				int width = j*binsize+binsize/2;
				StrandedPoint pi = new StrandedPoint(g,reg.getChrom(),midpoint-width/2,'+');
				StrandedPoint pj = new StrandedPoint(g,reg.getChrom(),midpoint+width/2,'+');
				Pair<StrandedPoint,StrandedPoint> pairi = new Pair<StrandedPoint,StrandedPoint>(pi,pj);
				for (int k=0; k<localreads.size(); k++) {
					Pair<StrandedPoint,StrandedPoint> pairj = localreads.get(k);
					tor[i][j] += selfliglik.get(pairj.car().distance(pairj.cdr()))*CrossValidation.normal(((double)pairi.car().distance(pairj.car()))/((double)h),((double)pairi.cdr().distance(pairj.cdr()))/((double)h),1d)/((double)h);
				}
			}
			if (i % (inc) == 0) {
				System.err.println("processed "+i+" bins for jack "+dfm.format(new Date()));
			}
		}

		return tor;
	}

	public double[] computeBandEstimate(StrandedRegion reg, int band) {
		int length = reg.getWidth();
		double[] tor = new double[length/binsize+1];

		int numbins = tor.length;
		StrandedRegion bigreg = reg.expand(3*h, 3*h);
		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg);
		List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (PairedReadDistribution.isSelfLigation(pairn) && bigreg.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}

		int inc = Math.max(numbins/100,1);
		for (int i=0; i<numbins; i++) {
			int midpoint = reg.getStart()+i*binsize+binsize/2;
			int width = band*binsize+binsize/2;
			StrandedPoint pi = new StrandedPoint(g,reg.getChrom(),midpoint-width/2,'+');
			StrandedPoint pj = new StrandedPoint(g,reg.getChrom(),midpoint+width/2,'+');
			Pair<StrandedPoint,StrandedPoint> pairi = new Pair<StrandedPoint,StrandedPoint>(pi,pj);
			for (int k=0; k<localreads.size(); k++) {
				Pair<StrandedPoint,StrandedPoint> pairj = localreads.get(k);
				tor[i] += selfliglik.get(pairj.car().distance(pairj.cdr()))*CrossValidation.normal(((double)pairi.car().distance(pairj.car()))/((double)h),((double)pairi.cdr().distance(pairj.cdr()))/((double)h),1d)/((double)h);
			}
			if (i % (inc) == 0) {
				System.err.println("processed "+i+" bins for jack "+dfm.format(new Date()));
			}
		}

		return tor;
	}

	public double[] computeSortedBandEstimate(StrandedRegion reg, int band) {
		int length = reg.getWidth();
		double[] tor = new double[length/binsize+1];

		int numbins = tor.length;
		StrandedRegion bigreg = reg.expand(3*h, 3*h);
		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg);
		List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (PairedReadDistribution.isSelfLigation(pairn) && bigreg.contains(pairn.cdr()) && width(pairn)<=band*binsize+2000) {
				localreads.add(pairn);
			}
		}
		Collections.sort(localreads,new MidpointComp());
		int start = 0;

		int inc = Math.max(numbins/100,1);
		if (localreads.size()>0) {
			for (int i=0; i<numbins; i++) {
				int midpoint = reg.getStart()+i*binsize+binsize/2;
				int width = band*binsize+binsize/2;
				StrandedPoint pi = new StrandedPoint(g,reg.getChrom(),midpoint-width/2,'+');
				StrandedPoint pj = new StrandedPoint(g,reg.getChrom(),midpoint+width/2,'+');
				Pair<StrandedPoint,StrandedPoint> pairi = new Pair<StrandedPoint,StrandedPoint>(pi,pj);
				while (midpoint-midpoint(localreads.get(start)) > 2000 && start<localreads.size()-1) start++;
				for (int k=start; k<localreads.size() && midpoint(localreads.get(k))-midpoint<=2000; k++) {
					Pair<StrandedPoint,StrandedPoint> pairj = localreads.get(k);
					tor[i] += selfliglik.get(pairj.car().distance(pairj.cdr()))*CrossValidation.normal(((double)pairi.car().distance(pairj.car()))/((double)h),((double)pairi.cdr().distance(pairj.cdr()))/((double)h),1d)/((double)h);
				}
				if (i % (inc) == 0) {
					//System.err.println("processed "+i+" bins for jack "+dfm.format(new Date()));
				}
			}
		}

		return tor;
	}

	/*assumes left reads will fall in reg */
	public int countPairsBetween(StrandedRegion reg, StrandedRegion query) {
		StrandedRegion bigreg = reg.expand(readDist.size()/2, readDist.size()/2);
		StrandedRegion bigquery = query.expand(readDist.size()/2, readDist.size()/2);
		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg);
		List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (bigquery.contains(pairn.cdr())) {
				localreads.add(pairn.swap());
			} else if (bigquery.contains(pairn.car()) && bigreg.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}
		return localreads.size();
	}

	/* Assumes left reads will fall in reg */
	public double[] computeSortedInterEstimate(StrandedRegion reg, StrandedRegion query) {
		int length = reg.getWidth();
		double[] tor = new double[length/binsize+1];
		int torbinstart = reg.getStart() / binsize;
		int torbinend = reg.getEnd() / binsize;
		int querybinstart = query.getStart() / binsize;
		int querybinend = query.getEnd() / binsize;

		int numbins = tor.length;
		StrandedRegion bigreg = reg.expand(readDist.size()/2, readDist.size()/2);
		StrandedRegion bigquery = query.expand(readDist.size()/2, readDist.size()/2);
		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg);
		List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (bigquery.contains(pairn.cdr())) {
				localreads.add(pairn.swap());
			} else if (bigquery.contains(pairn.car()) && bigreg.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}
		Collections.sort(localreads,new RightComp());
		int start = 0;

		int index = 0;
		//System.err.println(localreads.size()+" local reads from "+query+" in "+reg);
		int inc = Math.max(localreads.size()/100,1);
		for (Pair<StrandedPoint,StrandedPoint> p : localreads) {
			int rbinpos = p.cdr().getLocation() / binsize;
			int lbinpos = p.car().getLocation() / binsize;
			int rstart = Math.max(rbinpos - readDist.size()/2, torbinstart);
			int rend = Math.min(rbinpos + readDist.size()/2, torbinend+1);
			int lstart = Math.max(lbinpos - readDist.size()/2, querybinstart);
			int lend = Math.min(lbinpos + readDist.size()/2, querybinend+1);
			//System.err.println(rbinpos+" "+lbinpos+" "+rstart+" "+rend+" "+lstart+" "+lend);
			for (int i=rstart; i<rend; i++) {
				double rprob = 0d;
				if (p.cdr().getStrand()=='+') {
					rprob = readDist.get(i-rbinpos+readDist.size()/2);
				} else {
					rprob = readDist.get(readDist.size()-(i-rbinpos+readDist.size()/2+1));
				}
				for (int j=lstart; j<lend; j++) {
					double lprob = 0d;
					if (p.car().getStrand()=='+') {
						lprob = readDist.get(j-lbinpos+readDist.size()/2);
					} else {
						lprob = readDist.get(readDist.size()-(j-lbinpos+readDist.size()/2+1));
					}
					tor[i-torbinstart] += rprob*lprob;
					//System.err.println(i+" "+rprob+" "+lprob);
				}
			}
			index++;
			if (index % (inc) == 0) {
				//System.err.println("processed "+index+" reads for inter estimate "+dfm.format(new Date()));
			}
		}
		return tor;
	}

	public int width(Pair<StrandedPoint,StrandedPoint> pair) {
		return pair.cdr().getLocation() - pair.car().getLocation();
	}

	public int midpoint(Pair<StrandedPoint,StrandedPoint> pair) {
		return (pair.car().getLocation()+pair.cdr().getLocation()) / 2;
	}

	public double[][] computeEstimate(Region reg1, Region reg2, char leftstrand, char rightstrand) {
		int length1 = reg1.getWidth();
		int length2 = reg2.getWidth();
		double[][] tor = new double[length1/binsize+1][length2/binsize+1];
		int numvbins = length2 / binsize + 1;

		int numbins = length1/binsize + 1;
		//StrandedRegion bigreg = reg.expand(3*h, 3*h);
		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(reg1);
		List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (pairn.car().getStrand()==leftstrand && pairn.cdr().getStrand()==rightstrand && reg2.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}

		int inc = Math.max(numbins/100,1);
		for (int i=0; i<numbins; i++) {
			StrandedPoint pi = new StrandedPoint(g,reg1.getChrom(),reg1.getStart()+i*binsize+binsize/2,leftstrand);
			for (int j=0; j<numvbins; j++) {
				StrandedPoint pj = new StrandedPoint(g,reg2.getChrom(),reg2.getStart()+j*binsize+binsize/2,rightstrand);
				Pair<StrandedPoint,StrandedPoint> pairi = new Pair<StrandedPoint,StrandedPoint>(pi,pj);
				for (int k=0; k<localreads.size(); k++) {
					Pair<StrandedPoint,StrandedPoint> pairj = localreads.get(k);
					double weight = 1;
					if (leftstrand=='-' && rightstrand=='+' && pairj.car().getChrom().equals(pairj.cdr().getChrom()) &&
							pairj.car().distance(pairj.cdr())<selfliglik.size()) {
						weight = 1-selfliglik.get(pairj.car().distance(pairj.cdr()));
					}
					tor[i][j] += weight*CrossValidation.normal(((double)pairi.car().distance(pairj.car()))/((double)h),((double)pairi.cdr().distance(pairj.cdr()))/((double)h),1d)/((double)h);
				}
			}
			if (i % (inc) == 0) {
				System.err.println("processed "+i+" bins for jack "+dfm.format(new Date()));
			}
		}

		return tor;
	}

	public double[][] computeInterEstimate(Region reg1, Region reg2) {
		int length1 = reg1.getWidth();
		int length2 = reg2.getWidth();
		double[][] tor = new double[length1/binsize+1][length2/binsize+1];
		int numvbins = length2 / binsize + 1;

		int numbins = length1/binsize + 1;
		//StrandedRegion bigreg = reg.expand(3*h, 3*h);
		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(reg1);
		List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (reg2.contains(pairn.cdr())) {
				localreads.add(pairn);
				//				System.err.println(pairn);
			}
		}

		int inc = Math.max(numbins/100,1);
		for (int i=0; i<numbins; i++) {
			StrandedPoint pi = new StrandedPoint(g,reg1.getChrom(),reg1.getStart()+i*binsize+binsize/2,'+');
			for (int j=0; j<numvbins; j++) {
				StrandedPoint pj = new StrandedPoint(g,reg2.getChrom(),reg2.getStart()+j*binsize+binsize/2,'+');
				Pair<StrandedPoint,StrandedPoint> pairi = new Pair<StrandedPoint,StrandedPoint>(pi,pj);
				for (int k=0; k<localreads.size(); k++) {
					Pair<StrandedPoint,StrandedPoint> pairj = localreads.get(k);
					double weight = 1;
					//					if (leftstrand=='-' && rightstrand=='+' && pairj.car().getChrom().equals(pairj.cdr().getChrom()) &&
					//							pairj.car().distance(pairj.cdr())<selfliglik.size()) {
					//						weight = 1-selfliglik.get(pairj.car().distance(pairj.cdr()));
					//					}
					double prob = 0;
					int xdist = (pi.getLocation()-pairj.car().getLocation())/readdistbinsize+readDist.size()/2;
					int ydist = (pj.getLocation()-pairj.cdr().getLocation())/readdistbinsize+readDist.size()/2;
					if (pairj.car().getStrand()=='-') {
						if (xdist>readDist.size()/2) {
							xdist = readDist.size()/2 - (xdist-readDist.size()/2);
						} else {
							xdist = readDist.size()/2 + (readDist.size()/2-xdist);
						}
					}
					if (pairj.cdr().getStrand()=='-') {
						if (ydist>readDist.size()/2) {
							ydist = readDist.size()/2 - (ydist-readDist.size()/2);
						} else {
							ydist = readDist.size()/2 + (readDist.size()/2-ydist);
						}
					}
					if (xdist>=0 && xdist<readDist.size() && ydist>=0 && ydist<readDist.size()) {
						prob = readDist.get(xdist)*readDist.get(ydist);
					}
					if (i==450 && j==5 && prob>0) {
						System.err.println(pi+"\t"+pj+"\t"+pairi+"\t"+pairj+"\t"+xdist+"\t"+ydist);
					}
					tor[i][j] += weight*prob;
				}
			}
			if (i % (inc) == 0) {
				System.err.println("processed "+i+" bins for jack "+dfm.format(new Date()));
			}
		}

		return tor;
	}

	public double[][] smooth(double[][] in) {
		int kernelrad = 3*h/binsize;
		double[][] kernel = new double[2*kernelrad][2*kernelrad];
		StrandedPoint pi = new StrandedPoint(g,"1",kernelrad*binsize,'-');
		StrandedPoint pj = new StrandedPoint(g,"1",kernelrad*binsize,'+');
		Pair<StrandedPoint,StrandedPoint> pairi = new Pair<StrandedPoint,StrandedPoint>(pi,pj);
		for (int k=0; k<kernel.length; k++) {
			StrandedPoint pk = new StrandedPoint(g, "1", 0+k*binsize+binsize/2,'-');
			for (int l=0; l<kernel[k].length; l++) {
				StrandedPoint pl = new StrandedPoint(g, "1", 0+l*binsize+binsize/2,'+');
				Pair<StrandedPoint,StrandedPoint> pairj = new Pair<StrandedPoint,StrandedPoint>(pk,pl);
				kernel[k][l] += CrossValidation.normal(((double)pairi.car().distance(pairj.car()))/((double)h),((double)pairi.cdr().distance(pairj.cdr()))/((double)h),1d)/((double)h);
			}
		}
		double kernelsum = 0;
		for (int i=0; i<kernel.length; i++) {
			for (int k=0; k<kernel[i].length; k++) {
				kernelsum += kernel[i][k];
			}
		}
		for (int i=0; i<kernel.length; i++) {
			for (int j=0; j<kernel[i].length; j++) {
				kernel[i][j] /= kernelsum;
			}
		}
		double[][] tor = new double[in.length][in[0].length];
		int numvbins = tor[0].length;

		int numbins = tor.length;
		int inc = Math.max(numbins/100,1);
		for (int i=0; i<numbins; i++) {
			pi = new StrandedPoint(g,"1",0+i*binsize+binsize/2,'-');
			for (int j=0; j<numvbins; j++) {
				pj = new StrandedPoint(g,"1",0+j*binsize+binsize/2,'+');
				pairi = new Pair<StrandedPoint,StrandedPoint>(pi,pj);
				for (int k=i-kernelrad+1; k<i+kernelrad-1; k++) {
					StrandedPoint pk = new StrandedPoint(g, "1", 0+k*binsize+binsize/2,'-');
					for (int l=j-kernelrad+1; l<i+kernelrad-1; l++) {
						StrandedPoint pl = new StrandedPoint(g, "1", 0+l*binsize+binsize/2,'+');
						Pair<StrandedPoint,StrandedPoint> pairj = new Pair<StrandedPoint,StrandedPoint>(pk,pl);
						if (k>=0 && k<in.length && l>=0 && l<in[k].length && k-i+kernelrad>=0 && k-i+kernelrad<kernel.length && l-j+kernelrad>=0 && l-j+kernelrad<kernel.length) {
							tor[i][j] += in[k][l]*kernel[k-i+kernelrad][l-j+kernelrad];
						}
					}
				}
			}
			if (i % (inc) == 0) {
				System.err.println("processed "+i+" bins for jack "+dfm.format(new Date()));
			}
		}

		return tor;
	}

	public int computeCount(Region reg1, Region reg2, char leftstrand, char rightstrand) {
		Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(reg1);
		List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
		for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
			Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
			if (pairn.car().getStrand()==leftstrand && pairn.cdr().getStrand()==rightstrand && reg2.contains(pairn.cdr())) {
				localreads.add(pairn);
			}
		}
		return localreads.size();
	}

	public class MidpointComp implements Comparator<Pair<StrandedPoint,StrandedPoint>> {

		public int compare(Pair<StrandedPoint, StrandedPoint> arg0,
				Pair<StrandedPoint, StrandedPoint> arg1) {
			Integer midpoint0 = (arg0.car().getLocation() + arg0.cdr().getLocation()) / 2;
			Integer midpoint1 = (arg1.car().getLocation() + arg1.cdr().getLocation()) / 2;
			return midpoint0.compareTo(midpoint1);
		}

	}

	public class RightComp implements Comparator<Pair<StrandedPoint,StrandedPoint>> {

		public int compare(Pair<StrandedPoint, StrandedPoint> arg0,
				Pair<StrandedPoint, StrandedPoint> arg1) {
			// TODO Auto-generated method stub
			return arg0.cdr().compareTo(arg1.cdr());
		}

	}

}
