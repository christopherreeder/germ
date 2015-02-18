package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class CrossValidation {
	
	private static DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

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
		double roottwoh = Math.sqrt(2)*((double)h);
		double hsqrd = Math.pow(h, 2);
		double twohsqrd = 2*hsqrd;
		String readfile = Args.parseString(args, "readfile", "");
		String regionfile = Args.parseString(args, "regionfile", "");
		String selfliglikfile = Args.parseString(args, "selfliglikfile", "");
		int size = Args.parseInteger(args, "size", 0);
		Random rand = new Random();
		
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
				regionlist.add(tmp);
			}
		}
		
		r = new BufferedReader(new FileReader(selfliglikfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			selfliglik.add(Double.parseDouble(split[0]));
		}
		
		double fhatsqr = 0;
		double wsum = 0;
		double N = 0;
		int count = 0;
		int inc = regionlist.size() / 100;
		for (StrandedRegion reg : regionlist) {
			StrandedRegion bigreg = reg.expand(upstream+3*h, downstream+3*h);
			Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg);
			List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
			for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
				Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
				if (PairedReadDistribution.isSelfLigation(pairn) && bigreg.contains(pairn.cdr())) {
					localreads.add(pairn);
				}
			}
			for (int i=0; i<localreads.size(); i++) {
				Pair<StrandedPoint,StrandedPoint> pairi = localreads.get(i);
				double wi = selfliglik.get(pairi.car().distance(pairi.cdr()));
				wsum += wi;
				for (int j=0; j<localreads.size(); j++) {
					Pair<StrandedPoint,StrandedPoint> pairj = localreads.get(j);
					fhatsqr += wi*selfliglik.get(pairj.car().distance(pairj.cdr()))*normal(((double)pairi.car().distance(pairj.car()))/roottwoh,((double)pairi.cdr().distance(pairj.cdr()))/roottwoh,1d)/roottwoh;
				}
			}
			count++;
			if (count % (inc) == 0) {
				System.err.println("processed "+count+" regions for fhatsqr "+dfm.format(new Date()));
			}
		}
		
		count = 0;
		double jack = 0;
		for (StrandedRegion reg : regionlist) {
			StrandedRegion bigreg = reg.expand(upstream+3*h, downstream+3*h);
			Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(bigreg);
			List<Pair<StrandedPoint,StrandedPoint>> localreads = new ArrayList<Pair<StrandedPoint,StrandedPoint>>();
			for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
				Pair<StrandedPoint,StrandedPoint> pairn = storage.getLeftPair(n);
				if (PairedReadDistribution.isSelfLigation(pairn) && bigreg.contains(pairn.cdr())) {
					localreads.add(pairn);
				}
			}
			N += localreads.size();
			
			for (int i=0; i<localreads.size(); i++) {
				//compute the pointwise estimate of the integral of the product of f and \hat{f}
				double fhat = 0;
				Pair<StrandedPoint,StrandedPoint> pairi = localreads.get(i);
				for (int j=0; j<localreads.size(); j++) {
					if (i!=j) {
						Pair<StrandedPoint,StrandedPoint> pairj = localreads.get(j);
						fhat += selfliglik.get(pairj.car().distance(pairj.cdr()))*normal(((double)pairi.car().distance(pairj.car()))/((double)h),((double)pairi.cdr().distance(pairj.cdr()))/((double)h),1d)/((double)h);
					}
				}
				jack += fhat / (wsum - selfliglik.get(pairi.car().distance(pairi.cdr())));
			}
			count++;
			if (count % (inc) == 0) {
				System.err.println("processed "+count+" regions for jack "+dfm.format(new Date()));
			}
		}
		
		double ise = fhatsqr - 2d * jack / N;
		System.out.println(h+"\t"+ise);
	}
	
	public static double distance(Pair<StrandedPoint,StrandedPoint> pair1, Pair<StrandedPoint,StrandedPoint> pair2) {
		return Math.sqrt(Math.pow(pair1.car().distance(pair2.car()), 2) + Math.pow(pair1.cdr().distance(pair2.cdr()), 2));
	}
	
	public static double normal(double left, double right, double sigmasqrd) {
		return Math.exp(-1*(Math.pow(left,2)/sigmasqrd + Math.pow(right, 2)/sigmasqrd)/2d) / (2d * Math.PI * sigmasqrd);
	}

}
