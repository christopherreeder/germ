package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import edu.mit.csail.cgs.tools.utils.Args;

public class SubmitTSSSpec {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String submitfile = Args.parseString(args, "submitfile", "");
		String genomefile = Args.parseString(args, "genomefile", "");
		int binsize = Args.parseInteger(args, "binsize", 10);
		String margfile = Args.parseString(args, "margfile", "");
		String outbase = Args.parseString(args, "outbase", "");
		int numruns = Args.parseInteger(args, "numruns", 0);
		String pointfile = Args.parseString(args, "peakfile", "");
		double interthresh = Args.parseDouble(args, "interthresh", 0);
		int interthreshlow = Args.parseInteger(args, "interthreshlow", 0);
		int interthreshhigh = Args.parseInteger(args, "interthreshhigh", 0);
		double backval = Args.parseDouble(args, "backval", 0);
		int numsamps = Args.parseInteger(args, "numsamps", 0);
		String readdistfile = Args.parseString(args, "readdistfile", "");
		String readfile = Args.parseString(args, "readfile", "");
		int storagesize = Args.parseInteger(args, "storagesize", 0);
		int perrun = Args.parseInteger(args, "perrun", 0);
		String tssfile = Args.parseString(args, "tssfile", "");
		boolean allchrom = Args.parseFlags(args).contains("allchrom");
		boolean jointsum = Args.parseFlags(args).contains("jointsum");
		String errorfile = Args.parseString(args, "errorfile", "");

		String wd = Args.parseString(args, "wd", "/");
		PrintStream out = new PrintStream(submitfile);
		out.println("#!/bin/bash");
		out.println();

		List<Integer> numlist = new ArrayList<Integer>();
		if (!errorfile.equals("")) {
			BufferedReader r;
			String s;
			String[] split;
			r = new BufferedReader(new FileReader(errorfile));
			while ((s = r.readLine()) != null) {
				split = s.split(" ");
				split = split[1].split("_");
				numlist.add(Integer.parseInt(split[2].substring(0, split[2].length()-5)));
			}
		} else {
			for (int i=1; i<=numruns; i++) {
				numlist.add(i);
			}
		}

		for (int ithresh=interthreshlow; ithresh<=interthreshhigh; ithresh++) {
			interthresh = Math.pow(10d, ((double)ithresh));
			System.err.println(ithresh);
			for (int index=0; index<numlist.size(); index++) {
				int i = numlist.get(index);
				String outfile = outbase+"_"+(Math.abs(ithresh))+"_"+i+".txt";
				String jointsumfile = "";
				if (jointsum) {
					jointsumfile = outbase+"_"+(Math.abs(ithresh))+"_"+i+"jointsum.txt";
				}


				List<String> command = new ArrayList<String>();
				command.add("/usr/bin/qsub");
				command.add("-r");
				command.add("y");
				command.add("-l");
				command.add("mem_free=10G");
				command.add("-wd");
				command.add(wd);
				command.add("-v");
				command.add("'genomefile="+genomefile+",binsize="+binsize+",margfile="+margfile+",pointfile="+pointfile+
						",interthresh="+interthresh+",backval="+backval+",numsamps="+numsamps+",readdistfile="+readdistfile+
						",readfile="+readfile+",storagesize="+storagesize+",runnum="+i+",perrun="+perrun+",tssfile="+tssfile+
						",jointsumfile="+jointsumfile+",outfile="+outfile+"'");
				command.add("-o");
				command.add(outbase+"_"+(Math.abs(ithresh))+"stdout"+i+".txt");
				command.add("-e");
				command.add(outbase+"_"+(Math.abs(ithresh))+"stderr"+i+".txt");
				if (allchrom) {
					command.add("../runtssspecallchrom.sh");
				} else {
					command.add("../runtssspec.sh");
				}

				for (int k=0; k<command.size(); k++) {
					out.print(command.get(k)+" ");
				}
				out.println();
				out.println("sleep 1");
			}
		}
		out.flush();
		out.close();
	}

}
