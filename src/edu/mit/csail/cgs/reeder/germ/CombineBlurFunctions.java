package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import edu.mit.csail.cgs.tools.utils.Args;

public class CombineBlurFunctions {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		int minfile = Args.parseInteger(args, "minfile", 0);
		int maxfile = Args.parseInteger(args, "maxfile", 0);
		int trim = Args.parseInteger(args, "trim", 0);
		String gbase = Args.parseString(args, "gbase", "");
		String outfile = Args.parseString(args, "outfile", "");
		PrintStream out = new PrintStream(outfile);

		double[][] masterg = new double[0][0];
		BufferedReader r;
		String s;
		String[] split;
		int index = 0;
		boolean ginit = false;
		r = new BufferedReader(new FileReader(gbase+minfile+".txt"));
		while ((s = r.readLine()) != null) {
			if (index<trim) {
				index++;
				continue;
			}
			split = s.split("\t");
			if (!ginit) {
				masterg = new double[split.length-2*trim][split.length-2*trim];
				ginit = true;
			}
			if (index>=split.length-trim) {
				break;
			}
			for (int i=trim; i<split.length-trim; i++) {
				masterg[index-trim][i-trim] = Double.parseDouble(split[i]);
			}
			index++;
		}
		double max = 0;
		for (int i=0; i<masterg.length; i++) {
			for (int j=0; j<masterg[i].length; j++) {
				if (masterg[i][j]>max) {
					max = masterg[i][j];
				}
			}
		}
		for (int i=0; i<masterg.length; i++) {
			for (int j=0; j<masterg[i].length; j++) {
				masterg[i][j] /= max;
			}
		}
		for (int file=minfile+1; file<maxfile; file++) {
			double[][] tmpg = new double[masterg.length][masterg.length];
			index = 0;
			boolean exists = false;
			boolean good = true;
			r = new BufferedReader(new FileReader(gbase+file+".txt"));
			while ((s = r.readLine()) != null) {
				exists = true;
				if (index<trim) {
					index++;
					continue;
				}
				split = s.split("\t");
				if (index>=split.length-trim) {
					break;
				}
				for (int i=trim; i<split.length-trim; i++) {
					double tmpd = Double.parseDouble(split[i]);
					if (!Double.isNaN(tmpd)) {
						tmpg[index-trim][i-trim] = tmpd;
					} else {
						good = false;
					}
				}
				index++;
			}
			if (exists && !good) {
				System.err.println("contains NaN: "+file);
			}
			if (exists && good) {
				max = 0;
				for (int i=0; i<tmpg.length; i++) {
					for (int j=0; j<tmpg[i].length; j++) {
						if (tmpg[i][j]>max) {
							max = tmpg[i][j];
						}
					}
				}
				if (max==0) {
					System.err.println("max is 0: "+file);
				} else {
					for (int i=0; i<tmpg.length; i++) {
						for (int j=0; j<tmpg[i].length; j++) {
							tmpg[i][j] /= max;
						}
					}
					for (int i=0; i<tmpg.length; i++) {
						for (int j=0; j<tmpg[i].length; j++) {
							masterg[i][j] += tmpg[i][j];
						}
					}
				}
			}
		}
		double sum = 0;
		for (int i=0; i<masterg.length; i++) {
			for (int j=0; j<masterg[i].length; j++) {
				sum += masterg[i][j];
			}
		}
		for (int i=0; i<masterg.length; i++) {
			for (int j=0; j<masterg[i].length; j++) {
				masterg[i][j] /= sum;
			}
		}
		for (int i=0; i<masterg.length; i++) {
			for (int j=0; j<masterg[i].length; j++) {
				out.print(masterg[i][j]+"\t");
			}
			out.println();
		}
		out.flush();
		out.close();
	}

}
