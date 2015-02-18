package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class ThresholdFinder {

	private static DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	/**
	 * @param args
	 * @throws IOException 
	 * @throws NotFoundException 
	 */
	public static void main(String[] args) throws IOException, NotFoundException {
		Genome g = SproutUtils.parseGenome(args);
		String margfile = Args.parseString(args, "margfile", "");
		int expmin = Args.parseInteger(args, "expmin", 0);
		int expmax = Args.parseInteger(args, "expmax", 0);
		double expstep = Args.parseDouble(args, "expstep", 1);
		int binsize = Args.parseInteger(args, "binsize", 0);
		boolean printall = Args.parseFlags(args).contains("printall");
		String outfile = Args.parseString(args, "outfile", "");
		PrintStream out = new PrintStream(outfile);
		
		List<List<Region>> regionlists = new ArrayList<List<Region>>();
		List<Region> regions = new ArrayList<Region>();
		List<Double> marg = new ArrayList<Double>();

		double margmin = Double.MAX_VALUE;
		double margmax = 0;
		int index = 0;
		boolean inreg = false;
		int start = -1;
		int end = -1;
		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(margfile));
		s = r.readLine();
		System.err.println(s);
		split = s.split("\t");
		if (split.length==1) {
			split = s.split(" ");
		}
		split = split[1].split("=");
		String chrom = split[1].substring(3);
		int maxsize = Integer.MAX_VALUE/10;
		System.err.println("maxsize: "+maxsize);
		while (((s = r.readLine()) != null) && (marg.size()<maxsize)) {
			if (!s.startsWith("fixedStep")) {
				double tmp = Double.parseDouble(s);
				marg.add(tmp);
				if (tmp>margmax) margmax = tmp;
				if (tmp<margmin && tmp>0) margmin = tmp;
			} else {
				System.err.println(s);
			}
		}
		expmin = (int)Math.floor(Math.log10(margmin));
		expmax = (int)Math.ceil(Math.log10(margmax));
		System.err.println("margmin: "+margmin);
		System.err.println("margmax: "+margmax);
		System.err.println("expmin: "+expmin);
		System.err.println("expmax: "+expmax);
		expmin = expmax-10;
		System.err.println("new expmin: "+expmin);
		int numiters = (expmax-expmin+1)*20;
		for (double exp=expmin; exp<expmax; exp++) {
			int inc = marg.size()/100;
			for (double i=0; i<1; i += expstep) {
				double thresh = Math.pow(10d, exp+i);
				System.err.println("thresh: "+thresh);
				regions = extractDomains(marg, thresh, binsize, g, chrom);
				int mean = meanRegionSize(regions);
				int var = varianceRegionSize(regions, mean);
				if (printall) {
					regionlists.add(regions);
				} else {
					out.println(thresh+"\t"+regions.size()+"\t"+mean+"\t"+var);
				}
			}
		}
		if (printall) {
			int maxnum = 0;
			for (int i=0; i<regionlists.size(); i++) {
				if (regionlists.get(i).size()>maxnum) maxnum = regionlists.get(i).size();
			}
			for (int i=0; i<regionlists.size(); i++) {
				List<Region> tmp = regionlists.get(i);
				for (int j=0; j<maxnum; j++) {
					if (j<tmp.size()) {
						out.print(tmp.get(j).getWidth()+"\t");
					} else {
						out.print("0\t");
					}
				}
				out.println();
			}
		}
		out.flush();
		out.close();
	}

	public static int meanRegionSize(List<Region> list) {
		if (list.size()==0) {
			return 0;
		} else {
			int sum = 0;
			for (int i=0; i<list.size(); i++) {
				sum += list.get(i).getWidth();
			}
			return sum / list.size();
		}
	}

	public static int varianceRegionSize(List<Region> list, int mean) {
		if (list.size()==0) {
			return 0;
		} else {
			int sum = 0;
			for (int i=0; i<list.size(); i++) {
				sum += Math.pow(list.get(i).getWidth()-mean, 2);
			}
			return sum / list.size();
		}
	}

	public static List<Region> extractDomains(List<Double> marg, double thresh, int binsize, Genome g, String chrom) {
		List<Region> regions = new ArrayList<Region>();
		int index = 0;
		boolean inreg = false;
		int end = -1;
		int start = -1;
		for (Double tmpd : marg) {
			double tmp = tmpd;
			if (tmp<thresh) {
				if (inreg) {
					end = (index-1)*binsize + binsize/2;
					regions.add(new Region(g, chrom, start, end));
					inreg = false;
				}
			} else if (!inreg) {
				start = index*binsize + binsize/2;
				inreg = true;
			}
			index++;
		}
		return regions;
	}

}
