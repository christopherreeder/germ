package edu.mit.csail.cgs.reeder.germ;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

import edu.mit.csail.cgs.tools.utils.Args;

public class SubmitPseudo {

	/**
	 * @param args
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException {
		String base = Args.parseString(args, "base", "");
		String margfile = Args.parseString(args, "margfile", "");
		String genomefile = Args.parseString(args, "genomefile", "");
		int numsamps = Args.parseInteger(args, "numsamps", 0);
		double minmult = Args.parseDouble(args, "minmult", 7);
		double maxmult = Args.parseDouble(args, "maxmult", 100);
		double startinc = Args.parseDouble(args, "startinc", 1);
		double mininc = Args.parseDouble(args, "mininc", 0.001);
		String wd = Args.parseString(args, "wd", "");
		int cutofflow = Args.parseInteger(args, "cutofflow", 0);
		int cutoffhigh = Args.parseInteger(args, "cutoffhigh", 0);
		String outfile = Args.parseString(args, "outfile", "");
		PrintStream out = new PrintStream(outfile);
		
		for (int cutoff=cutofflow; cutoff<=cutoffhigh; cutoff++) {
			int abscutoff = Math.abs(cutoff);
			List<String> command = new ArrayList<String>();
			command.add("echo \"java -Xmx10g -cp /cluster/ccr/germ.jar edu.mit.csail.cgs.reeder.germ.DeterminePseudo --genomefile");
			command.add("\\\""+genomefile+"\\\"");
			command.add("--binsize 10");
			command.add("--infile");
			command.add("\\\""+base+"_"+abscutoff+"_all.txt\\\"");
			command.add("--auxfile \\\"\\\"");
			command.add("--margfile");
			command.add("\\\""+margfile+"\\\"");
			command.add("--outfile");
			command.add("\\\""+base+"_"+abscutoff+"_all_pvalsaux.txt\\\"");
			command.add("--specoutfile");
			command.add("\\\""+base+"_"+abscutoff+"_all_pvals.txt\\\"");
			command.add("--minmult");
			command.add(Double.toString(minmult));
			command.add("--maxmult");
			command.add(Double.toString(maxmult));
			command.add("--multtries 1");
			command.add("--startinc");
			command.add(Double.toString(startinc));
			command.add("--mininc");
			command.add(Double.toString(mininc));
			command.add("--pcutoff 0.05 --fdrcutoff 0.1");
			command.add("--numsamps");
			command.add(numsamps+"\" | /usr/bin/qsub -l mem_free=10G");
			command.add("-wd");
			command.add(wd);
			command.add("-o");
			command.add(base+"_"+abscutoff+"out.txt");
			command.add("-e");
			command.add(base+"_"+abscutoff+"err.txt");
			for (int i=0; i<command.size(); i++) {
				out.print(command.get(i)+" ");
			}
			out.println();
		}
		out.flush();
		out.close();
	}

}
