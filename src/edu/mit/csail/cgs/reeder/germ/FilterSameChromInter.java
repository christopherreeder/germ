package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class FilterSameChromInter {

	/**
	 * @param args
	 * @throws NotFoundException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NotFoundException, IOException {
		Genome g = SproutUtils.parseGenome(args);
		boolean self = Args.parseFlags(args).contains("self");
		String infile = Args.parseString(args, "infile", "");
		String outfile = Args.parseString(args, "outfile", "");
		int distcutoff = Args.parseInteger(args, "distcutoff", 0);
		int selfdistcutoff = Args.parseInteger(args, "selfdistcutoff", 10000);
		PrintStream out = new PrintStream(outfile);
		
		BufferedReader r;
		String s;
		String[] split;
		r = new BufferedReader(new FileReader(infile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			StrandedPoint tmp1 = StrandedPoint.fromString(g, split[0]);
			StrandedPoint tmp2 = StrandedPoint.fromString(g, split[1]);
			if (tmp1==null || tmp2==null) {
				System.err.println(s);
			} else {
				boolean samechrom = tmp1.getChrom().equals(tmp2.getChrom());
				int distance = samechrom ? tmp1.distance(tmp2) : Integer.MAX_VALUE;
				if (self && tmp1.getChrom().equals(tmp2.getChrom()) && (tmp1.getStrand()=='-' && tmp2.getStrand()=='+') && distance>distcutoff) {
					out.println(s);
				} else if (!self && tmp1.getChrom().equals(tmp2.getChrom()) && (distance>selfdistcutoff || !(tmp1.getStrand()=='-' && tmp2.getStrand()=='+')) && distance>distcutoff) {
					out.println(s);
				}
			}
		}
		out.flush();
		out.close();
	}

}
