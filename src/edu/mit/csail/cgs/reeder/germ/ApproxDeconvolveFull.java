package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class ApproxDeconvolveFull {

	/**
	 * @param args
	 * @throws NotFoundException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NotFoundException, IOException {
		Genome g = SproutUtils.parseGenome(args);
		int h = Args.parseInteger(args, "h", 0);
		int binsize = Args.parseInteger(args, "binsize", 0);
		String readfile = Args.parseString(args, "readfile", "");
		String selfliglikfile = Args.parseString(args, "selfliglikfile", "");
		String chromfile = Args.parseString(args, "chromfile", "");
		int storagesize = Args.parseInteger(args, "storagesize", 0);
		int band = Args.parseInteger(args, "band", 0);
		String foutfile = Args.parseString(args, "foutfile", "");
		PrintStream fout = new PrintStream(foutfile);
		
		Map<String,double[]> chromtoq = new HashMap<String,double[]>();
		List<Region> regions = new ArrayList<Region>();
		ComputeEstimate ce;
		ce = new ComputeEstimate(g, h, binsize, storagesize, readfile, selfliglikfile);
		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(chromfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			regions.add(SproutUtils.chromRegion(split[0],g));
		}
		
		for (Region tmp : regions) {
			StrandedRegion reg = new StrandedRegion(tmp,'+');
			double[] f = ce.computeSortedBandEstimate(reg, band);
			chromtoq.put(tmp.getChrom(), f);
		}
		
		long length = 0;
		double sum = 0;
		for (String chrom : chromtoq.keySet()) {
			double[] q = chromtoq.get(chrom);
			for (int i=0; i<q.length; i++) {
				sum += q[i];
			}
			length += q.length;
		}
		System.err.println("total length: "+length);
		System.err.println("sum: "+sum);
		
		System.err.println("writing");
		for (String chrom : chromtoq.keySet()) {
			fout.println("fixedStep"+"\t"+"chrom=chr"+chrom+"\t"+"start="+1+"\tstep="+binsize);
			double[] tmplist = chromtoq.get(chrom);
			for (int i=0; i<tmplist.length; i++) {
				fout.println((tmplist[i]/sum));
			}
		}
		fout.flush();
		fout.close();
	}

}
