package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class SubmitBlindDeconvolve {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws NotFoundException 
	 */
	public static void main(String[] args) throws IOException, NotFoundException {
		String submitfile = Args.parseString(args, "submitfile", "");
		int minnum = Args.parseInteger(args, "minnum", 0);
		int maxnum = Args.parseInteger(args, "maxnum", 0);
		int numouterinters = Args.parseInteger(args, "numouterinters", 0);
		int numinnerinters = Args.parseInteger(args, "numinnerinters", 0);
		String postname = Args.parseString(args, "postname", "");
		String finbase = Args.parseString(args, "finbase", "");
		String cinbase = Args.parseString(args, "cinbase", "");
		String ginbase = Args.parseString(args, "ginbase", "");
		String clenbase = Args.parseString(args, "clenbase", "");
		int numc = Args.parseInteger(args, "numc", 0);
		String goutbase = Args.parseString(args, "goutbase", "");
		String foutbase = Args.parseString(args, "foutbase", "");
		String coutbase = Args.parseString(args, "coutbase", "");
		String outbase = Args.parseString(args, "outbase", "");
		String wd = Args.parseString(args, "wd", "/");
		PrintStream out = new PrintStream(submitfile);
		out.println("#!/bin/bash");
		out.println();

		for (int i=minnum; i<maxnum; i++) {
			String finfile = finbase+i+postname;
			String cinfile = cinbase+i+postname;
			String ginfile = ginbase+i+postname;
			String clenfile = clenbase+i+postname;
			String goutfile = goutbase+i+".txt";
			String foutfile = foutbase+i+".txt";
			String coutfile = coutbase+i+".txt";
			
			List<String> command = new ArrayList<String>();
			command.add("/usr/bin/qsub");
			command.add("-l");
			command.add("mem_free=2G");
			command.add("-wd");
			command.add("/cluster/ccr/deconvolve"+wd);
			command.add("-q");
			command.add("batch");
			command.add("-v");
			command.add("'numouterinters="+numouterinters+",numinnerinters="+numinnerinters+",finfile="+finfile+",cinfile="+
							cinfile+",ginfile="+"../emp2kg12k_norm.txt"+",clenfile="+clenfile+",numc="+numc+
							",goutfile="+goutfile+",foutfile="+foutfile+",coutfile="+coutfile+"'");
			command.add("-o");
			command.add(outbase+"stdout"+i+".txt");
			command.add("-e");
			command.add(outbase+"stderr"+i+".txt");
			command.add("../runblinddeconvolve.sh");

			for (int k=0; k<command.size(); k++) {
				out.print(command.get(k)+" ");
			}
			out.println();
		}
		out.flush();
		out.close();
	}

}
