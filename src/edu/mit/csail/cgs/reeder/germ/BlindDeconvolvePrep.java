package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class BlindDeconvolvePrep {

	private static DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		Genome g = SproutUtils.parseGenome(args);
		int upstream = Args.parseInteger(args, "upstream", 0);
		int downstream = Args.parseInteger(args, "downstream", 0);
		int targetsize = Args.parseInteger(args, "targetsize", -1);
		int h = Args.parseInteger(args, "h", 0);
		int binsize = Args.parseInteger(args, "binsize", 0);
		String readfile = Args.parseString(args, "readfile", "");
		String region = Args.parseString(args, "region", "");
		String regionfile = Args.parseString(args, "regionfile", "");
		String selfliglikfile = Args.parseString(args, "selfliglikfile", "");
		int storagesize = Args.parseInteger(args, "storagesize", 0);
		String foutbase = Args.parseString(args, "foutbase", "");
		String coutbase = Args.parseString(args, "coutbase", "");
		String clenoutbase = Args.parseString(args, "clenoutbase", "");
		boolean execute = Args.parseFlags(args).contains("execute");
		String ginfile = Args.parseString(args, "ginfile", "");
		int foffset = 27;
		int fsmooth = 100;
		
		double[][] gmat = new double[0][0];
		if (execute) {
			int index = 0;
			boolean ginit = false;
			BufferedReader r = new BufferedReader(new FileReader(ginfile));
			String s;
			String[] split;
			while ((s = r.readLine()) != null) {
				split = s.split("\t");
				if (!ginit) {
					gmat = new double[split.length][split.length];
					ginit = true;
				}
				for (int i=0; i<split.length; i++) {
					double tmp = Double.parseDouble(split[i]);
					gmat[index][i] = tmp;
				}
				index++;
			}
		}

		ComputeEstimate ce = new ComputeEstimate(g, h, binsize, storagesize, readfile, selfliglikfile);
		List<StrandedRegion> regions = new ArrayList<StrandedRegion>();
		BufferedReader r;
		String s;
		String[] split;
		if (region.equals("")) {
			r = new BufferedReader(new FileReader(regionfile));
			while ((s = r.readLine()) != null) {
				split = s.split("\t");
				Region tmp = Region.fromString(g, split[0]);
				if (null==tmp) {
					System.err.println("no region "+s);
				} else {
					StrandedRegion tmps = new StrandedRegion(tmp,'+');
					if (targetsize==-1) {
						regions.add(tmps.expand(upstream, downstream));
					} else {
						if (tmps.getWidth()<targetsize) {
							int tmpex = (targetsize-tmps.getWidth()) / 2;
							regions.add(tmps.expand(tmpex,tmpex));
						} else {
							int tmpex = (tmps.getWidth()-targetsize);
							regions.add(new StrandedRegion(g,tmps.getChrom(),tmps.getStart(),tmps.getStart()+tmpex,'+'));
						}
					}
				}
			}
		} else {
			Region tmp = Region.fromString(g, region);
			if (null==tmp) {
				System.err.println("no region "+region);
			} else {
				StrandedRegion tmps = new StrandedRegion(tmp,'+');
				if (targetsize==-1) {
					regions.add(tmps.expand(upstream, downstream));
				} else {
					if (tmps.getWidth()<targetsize) {
						int tmpex = (targetsize-tmps.getWidth()) / 2;
						regions.add(tmps.expand(tmpex,tmpex));
					} else {
						int tmpex = (tmps.getWidth()-targetsize);
						regions.add(new StrandedRegion(g,tmps.getChrom(),tmps.getStart(),tmps.getStart()+tmpex,'+'));
					}
				}
			}
		}

		for (int i=0; i<regions.size(); i++) {
			double[][] c = ce.computeEstimate(regions.get(i));
			double csum = 0;
			for (int x=0; x<c.length; x++) {
				for (int y=0; y<c[x].length; y++) {
					csum += c[x][y];
				}
			}
			System.err.println("csum: "+csum);
			for (int x=0; x<c.length; x++) {
				for (int y=0; y<c[x].length; y++) {
					c[x][y] /= csum;
				}
			}
			double[] f = new double[c.length];
			for (int x=foffset; x<f.length-foffset; x++) {
				f[x] = c[x-foffset][x];
			}
			for (int x=0; x<fsmooth; x++) {
				f[x] = f[fsmooth-1]*(x+1)/fsmooth;
			}
			for (int x=f.length-fsmooth; x<f.length; x++) {
				f[x] = f[f.length-fsmooth]*(f.length-x)/fsmooth;
			}
			double fsum = 0;
			for (int x=0; x<f.length; x++) {
				fsum += f[x];
			}
			for (int x=0; x<f.length; x++) {
				f[x] /= fsum;
			}
			if (execute) {
				double[][][] bigc = new double[1][1][1];
				bigc[0] = c;
				BlindDeconvolve bd = new BlindDeconvolve(foffset, true, 1, 20, f, bigc, gmat);
				PrintStream fout = new PrintStream(foutbase+"out"+i+".txt");
				bd.execute(fout);
				fout.flush();
				fout.close();
			} else {
				PrintStream cout = new PrintStream(coutbase+i+".txt");
				PrintStream fout = new PrintStream(foutbase+i+".txt");
				PrintStream clenout = new PrintStream(clenoutbase+i+".txt");
				for (int x=0; x<c.length; x++) {
					for (int y=0; y<c[x].length; y++) {
						cout.print(c[x][y]+"\t");
					}
					cout.println();
				}
				for (int x=0; x<f.length; x++) {
					fout.println(f[x]);
				}
				clenout.println(c.length);
				cout.flush();
				cout.close();
				fout.flush();
				fout.close();
				clenout.flush();
				clenout.close();
			}
			System.err.println("processed region "+i+" "+dfm.format(new Date()));
		}
	}

}
