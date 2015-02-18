package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import edu.mit.csail.cgs.tools.utils.Args;

public class SubmitPrep {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String submitfile = Args.parseString(args, "submitfile", "");
		String genomefile = Args.parseString(args, "genomefile", "");
		int targetsize = Args.parseInteger(args, "targetsize", 0);
		int h = Args.parseInteger(args, "h", 0);
		int binsize = Args.parseInteger(args, "binsize", 0);
		String readfile = Args.parseString(args, "readfile", "");
		String regionfile = Args.parseString(args, "regionfile", "");
		String selfliglikfile = Args.parseString(args, "selfliglikfile", "");
		int storagesize = Args.parseInteger(args, "storagesize", 0);
		
		String clenoutbasebase = Args.parseString(args, "clenoutbasebase", "");
		String foutbasebase = Args.parseString(args, "foutbasebase", "");
		String coutbasebase = Args.parseString(args, "coutbasebase", "");
		String outbase = Args.parseString(args, "outbase", "");
		String wd = Args.parseString(args, "wd", "/");
		PrintStream out = new PrintStream(submitfile);
		out.println("#!/bin/bash");
		out.println();
		
		List<String> regionlist = new ArrayList<String>();
		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(regionfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			regionlist.add(split[0]);
		}

		for (int i=0; i<regionlist.size(); i++) {
			String clenoutbase = clenoutbasebase+i+"_";
			String foutbase = foutbasebase+i+"_";
			String coutbase = coutbasebase+i+"_";
			String region = regionlist.get(i);
			
			List<String> command = new ArrayList<String>();
			command.add("/usr/bin/qsub");
			command.add("-l");
			command.add("mem_free=4G");
			command.add("-wd");
			command.add(wd);
			command.add("-q");
			command.add("batch");
			command.add("-v");
			command.add("'genomefile="+genomefile+",targetsize="+targetsize+",h="+h+",binsize="+
							binsize+",readfile="+readfile+",region="+region+",selfliglikfile="+selfliglikfile+
							",storagesize="+storagesize+",foutbase="+foutbase+",clenoutbase="+clenoutbase+",coutbase="+coutbase+"'");
			command.add("-o");
			command.add(outbase+"stdout"+i+".txt");
			command.add("-e");
			command.add(outbase+"stderr"+i+".txt");
			command.add("../runprep.sh");

			for (int k=0; k<command.size(); k++) {
				out.print(command.get(k)+" ");
			}
			out.println();
		}
		out.flush();
		out.close();
	}

}
