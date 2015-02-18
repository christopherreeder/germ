package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import edu.mit.csail.cgs.tools.utils.Args;

public class BlindDeconvolve {

	private static DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	private int numouterinters, numinnerinters, foffset;
	double[] f;
	double[] oldf;
	double[][][] c;
	double[][][] tmpc;
	double[][] g;
	double[][] oldg;
	boolean nonblind;

	public BlindDeconvolve(int foffset, boolean nonblind, int numouterinters, int numinnerinters, double[] f, double[][][] c, double[][] g) {
		this.f = f;
		this.oldf = new double[f.length];
		for (int i=0; i<f.length; i++) {
			this.oldf[i] = f[i];
		}
		this.c = c;
		this.tmpc = new double[c.length][c[0].length][c[0][0].length];
		for (int i=0; i<c.length; i++) {
			for (int j=0; j<c[i].length; j++) {
				for (int k=0; k<c[i][j].length; k++) {
					this.tmpc[i][j][k] = c[i][j][k];
				}
			}
		}
		this.g = g;
		this.oldg = new double[g.length][g[0].length];
		for (int i=0; i<g.length; i++) {
			for (int j=0; j<g[i].length; j++) {
				this.oldg[i][j] = g[i][j];
			}
		}
		this.numouterinters = numouterinters;
		this.numinnerinters = numinnerinters;
		this.nonblind = nonblind;
		this.foffset = foffset;
	}

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		boolean nonblind = Args.parseFlags(args).contains("nonblind");
		int numouterinters = Args.parseInteger(args, "numouterinters", 0);
		int numinnerinters = Args.parseInteger(args, "numinnerinters", 0);
		String finfile = Args.parseString(args, "finfile", "");
		String cinfile = Args.parseString(args, "cinfile", "");
		String ginfile = Args.parseString(args, "ginfile", "");
		String clenfile = Args.parseString(args, "clenfile", "");
		int numc = Args.parseInteger(args, "numc", 0);
		String goutfile = Args.parseString(args, "goutfile", "");
		PrintStream gout = new PrintStream(goutfile);
		String foutfile = Args.parseString(args, "foutfile", "");
		PrintStream fout = new PrintStream(foutfile);
		String coutfile = Args.parseString(args, "coutfile", "");
		PrintStream cout = new PrintStream(coutfile);
		int foffset = Args.parseInteger(args, "foffset", 0);
		//String tmpoutfile = Args.parseString(args, "tmpoutfile", "");
		//PrintStream tmpout = new PrintStream(tmpoutfile);

		double lastmaxtmp = Double.MAX_VALUE;
		double lastmintmp = Double.MAX_VALUE;
		double[] f = new double[0];
		double[] oldf = new double[0];
		double[][][] c = new double[numc][0][0];
		double[][][] tmpc = new double[numc][0][0];
		double[][] g = new double[0][0];
		double[][] oldg = new double[0][0];

		BufferedReader r;
		String s;
		String[] split;
		int index = 0;
		r = new BufferedReader(new FileReader(clenfile));
		while ((s = r.readLine()) != null) {
			int tmp = Integer.parseInt(s);
			c[index] = new double[tmp][0];
			tmpc[index] = new double[tmp][0];
			f = new double[tmp];
			oldf = new double[tmp];
			index++;
		}

		index = 0;
		r = new BufferedReader(new FileReader(finfile));
		while ((s = r.readLine()) != null) {
			double tmp = Double.parseDouble(s);
			f[index] = tmp;
			oldf[index] = tmp;
			index++;
		}

		int geneindex = 0;
		index = 0;
		r = new BufferedReader(new FileReader(cinfile));
		while ((s = r.readLine()) != null) {
			if (s.equals("")) {
				geneindex++;
				index = 0;
			} else {
				split = s.split("\t");
				c[geneindex][index] = new double[split.length];
				tmpc[geneindex][index] = new double[split.length];
				for (int i=0; i<split.length; i++) {
					Double tmp = Double.parseDouble(split[i]);
					c[geneindex][index][i] = tmp;
					tmpc[geneindex][index][i] = tmp;
				}
				index++;
			}
		}

		index = 0;
		boolean ginit = false;
		r = new BufferedReader(new FileReader(ginfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			if (!ginit) {
				g = new double[split.length][split.length];
				oldg = new double[split.length][split.length];
				ginit = true;
			}
			for (int i=0; i<split.length; i++) {
				double tmp = Double.parseDouble(split[i]);
				g[index][i] = tmp;
				oldg[index][i] = tmp;
			}
			index++;
		}
		/*
		for (int x = 0; x<f.length; x++) {
			fout.print(f[x]+"\t");
		}
		fout.println();

		for (int x=0; x<g.length; x++) {
			for (int y=0; y<g[x].length; y++) {
				gout.print(g[x][y]+"\t");
			}
			gout.println();
		}

		tmpc[0] = new double[tmpc[0].length][tmpc[0][0].length];
		for (int x=0; x<tmpc[0].length; x++) {
			for (int y=0; y<tmpc[0][x].length; y++) {
				int m = Math.max(0,x-(g.length/2)+1);
				int my = Math.max(0,y-(g.length/2)+1);
				m = Math.max(m,my);
				int maxm = Math.min(f.length-1,x+(g.length/2)+1);
				int maxmy = Math.min(f.length-1,y+(g.length/2)+1);
				maxm = Math.min(maxm, maxmy);
				for (;m<maxm; m++) {
					tmpc[0][x][y] += g[x-m+(g.length/2)+1][y-m+(g.length/2)+1]*f[m];
				}
			}
		}

		for (int x=0; x<tmpc[0].length; x++) {
			for (int y=0; y<tmpc[0][x].length; y++) {
				cout.print(tmpc[0][x][y]+"\t");
			}
			cout.println();
		}
		 */
		int totaliters = numouterinters*numinnerinters*2;

		//fix the boundaries of f

		/*		for (int i=0; i<numinnerinters; i++) {
			for (int gene=0; gene<c.length; gene++) {
				double[][] tmpgene = c[gene];
				double[][] tmp1 = new double[tmpgene.length][tmpgene[0].length];
				//for (int x=2*g.length; x<tmpgene.length-(2*g.length); x++) {
				for (int x=0; x<tmpgene.length; x++) {
					double[] tmpgenex = tmpgene[x];
					//for (int y=2*g.length; y<tmpgenex.length-(2*g.length); y++) {
					for (int y=0; y<tmpgenex.length; y++) {	

						int m = Math.max(0,-(g.length/2)+x);
						m = Math.max(m, -(g[0].length/2)+y);
						int maxm = Math.min(oldf.length, (g.length/2)+x+1);
						maxm = Math.min(maxm, (g[0].length)/2+y+1);

						int m = 0;
						int maxm = f.length;
						for (; m<maxm; m++) {
							if (x-m>=-oldg.length/2 && y-m>=-oldg.length/2 && x-m<g.length/2 && y-m<g.length/2) {
								tmp1[x][y] += g[x-m+(g.length/2)][y-m+(g.length/2)]*oldf[m];
							}
						}
//						if (tmpgene[x][y]==0 && tmp1[x][y]==0) {
						if (tmp1[x][y]==0) {
							tmp1[x][y] = 1;
						} else {

							if (tmp1[x][y]==0) {
								throw new Exception(tmpgene[x][y]+" "+x+" "+y);
							}

							tmp1[x][y] = tmpgene[x][y] / tmp1[x][y];
						}
					}
				}
				tmpc[gene] = tmp1;
			}
			int x = 0;
			int maxx = f.length;
			for (; x<maxx; x++) {
				double tmp = 0;
				boolean nan = false;
				//int count = 0;
				for (int gene=0; gene<c.length; gene++) {
					double[][] tmp1 = tmpc[gene];

					int m = -(Math.min(tmp1.length-x-1, g.length/2-1));
					int my = -(Math.min(tmp1[0].length-x-1, g.length/2-1));
					int miny = my;
					//m = Math.max(m, my);
					int maxm = Math.min(x, g.length/2+1);
					int maxmy = Math.min(x, g.length/2+1);
					//maxm = Math.min(maxm,maxmy);
					//if (x==0) System.err.println(m+" "+maxm+" "+my+" "+maxmy);


					int m = Math.max(-f.length+2*g.length+x+1, -g.length/2 +1 +1);
					int my = Math.max(-f.length + 2*g.length + x + 1, -g.length/2 + 1 + 1);
					int miny = my;
					int maxm = Math.min(-2*g.length + x + 1,g.length/2 + 1 + 1);
					int maxmy = Math.min(-2*g.length + x + 1, g.length / 2 + 1 + 1);

					int m = Math.max(-f.length+x+1, -g.length/2 +1 +1);
					int my = Math.max(-f.length + x + 1, -g.length/2 + 1 + 1);
					int miny = my;
					int maxm = Math.min(x + 1,g.length/2 + 1 + 1);
					int maxmy = Math.min(x + 1, g.length / 2 + 1 + 1);
					for (; m<maxm; m++) {
						for (my=miny;my<maxmy; my++) {
							//if (x==0 && g[-m+(g.length/2)+1][-my+(g.length/2)+1]>0) System.err.println(m);
							//if (x==0) System.err.println(m+" "+tmp1[x-m][x-my]+" "+g[-m+(g.length/2)+1][-my+(g.length/2)+1]);
							//if (x-m>=2*g.length && x-my>=2*g.length && -m>=-(g.length/2)-1 && -my>=-(g.length/2)-1 && -m<(g.length/2)-1 && -my<(g.length/2)-1 && x-m<f.length-(2*g.length) && x-my<f.length-(2*g.length)) {
							if (x-m>=0 && x-my>=0 && -m>=-(g.length/2)-1 && -my>=-(g.length/2)-1 && -m<(g.length/2)-1 && -my<(g.length/2)-1 && x-m<f.length && x-my<f.length) {
								tmp += tmp1[x-m][x-my]*g[-m+(g.length/2)+1][-my+(g.length/2)+1];
								//count++;
							}
							if (Double.isNaN(tmp) && !nan) {
								System.err.println(tmp1[x-m][x-my]);
								nan = true;
							}
						}
					}
				}
				//if (tmp==0) System.err.println(x);
				f[x] = tmp*oldf[x];
				//System.err.println("count: "+count);
			}
			for (x = 0; x<f.length; x++) {
				oldf[x] = f[x];
			}
			for (x = 0; x<f.length; x++) {
				fout.print(f[x]+"\t");
			}
			fout.println();
			System.err.println(i+"pre iterations "+dfm.format(new Date()));
		}*/

		//Now start the actual AM procedure
		for (int k=0; k<numouterinters; k++) {
			if (!nonblind) {
				for (int i=0; i<numinnerinters; i++) {
					double maxtmp = 0;
					for (int gene=0; gene<c.length; gene++) {
						double[][] tmpgene = c[gene];
						double[][] tmp1 = new double[tmpgene.length][tmpgene[0].length];
						//for (int x=2*g.length; x<tmpgene.length-(2*g.length); x++) {
						for (int x=0; x<tmpgene.length; x++) {
							double[] tmpgenex = tmpgene[x];
							//for (int y=2*g.length; y<tmpgenex.length-(2*g.length); y++) {
							for (int y=0; y<tmpgenex.length; y++) {
								/*
							int m = Math.max(0,-(oldg.length/2)+x);
							m = Math.max(m, -(oldg[0].length/2)+y);
							int maxm = Math.min(f.length, (oldg.length/2)+x+1);
							maxm = Math.min(maxm, (oldg[0].length)/2+y+1);
								 */
								int m = 0;
								int maxm = f.length;
								for (; m<maxm; m++) {
									if (x-m>=-oldg.length/2 && y-m>=-oldg.length/2 && x-m<g.length/2 && y-m<g.length/2) {
										tmp1[x][y] += oldg[x-m+(oldg.length/2)][y-m+(oldg.length/2)]*f[m];
									}
								}
								//							if (tmpgene[x][y]==0 && tmp1[x][y]==0) {
								if (tmp1[x][y]==0) {
									tmp1[x][y] = 1;
								} else {
									tmp1[x][y] = tmpgene[x][y] / tmp1[x][y];
								}
							}
						}
						tmpc[gene] = tmp1;
					}
					int x = -(g.length/2)-1;
					int maxx = (g.length/2);
					for (;x<maxx; x++) {
						int y = -(g.length/2)-1;
						int ymax = (g.length/2);
						for (; y<ymax; y++) {
							double tmp = 0;
							for (int gene=0; gene<c.length; gene++) {
								double[][] tmp1 = tmpc[gene];
								/*
							int m = -(Math.min(tmp1.length-1, tmp1.length-x-1));
							int my = -(Math.min(tmp1[0].length-1, tmp1[0].length-y-1));
							m = Math.max(m, my);
							if (x<0) {
								m = Math.max(m, -x);
							}
							if (y<0) {
								m = Math.max(m,-y);
							}
							for (;m<=0; m++) {
								try {
									tmp += tmp1[x-m][y-m]*f[-m];
								} catch (Exception e) {
									System.err.println(x+"\t"+y+"\t"+m);
									throw e;
								}
							}
								 */
								for (int m=-f.length; m<f.length; m++) {
									//if ((x-m>=2*g.length) && (y-m>=2*g.length) && (-m>=0) && (x-m<f.length-(2*g.length)) && (y-m<f.length-(2*g.length)) && (-m<f.length)) {
									if ((x-m>=0) && (y-m>=0) && (-m>=0) && (x-m<f.length) && (y-m<f.length) && (-m<f.length)) {
										tmp += tmp1[x-m][y-m]*f[-m];
									}
								}
							}
							g[x+(g.length/2)+1][y+(g.length/2)+1] = tmp*oldg[x+(g.length/2)+1][y+(g.length)/2+1];
						}
					}
					for (x=0; x<g.length; x++) {
						for (int y=0; y<g[x].length; y++) {
							double tmp = g[x][y] / oldg[x][y];
							if (tmp>maxtmp) maxtmp = tmp;
							//tmpout.print(tmp+"\t");
							oldg[x][y] = g[x][y];
						}
						//tmpout.println();
					}
					//tmpout.flush();
					/*
				for (x=0; x<g.length; x++) {
					for (int y=0; y<g[x].length; y++) {
						gout.print(g[x][y]+"\t");
					}
					gout.println();
				}
				gout.flush();
					 */
					System.err.println("g maxtmp: "+maxtmp);
					System.err.println((k*numinnerinters*2)+i+" iterations "+dfm.format(new Date()));
					if (maxtmp<=1.1) break;
				}
			}
			/*
			for (int x=0; x<g.length; x++) {
				for (int y=0; y<g[x].length; y++) {
					gout.print(g[x][y]+"\t");
				}
				gout.println();
			}
			tmpc[0] = new double[tmpc[0].length][tmpc[0][0].length];
			for (int x=0; x<tmpc[0].length; x++) {
				for (int y=0; y<tmpc[0][x].length; y++) {
					int m = Math.max(0,x-(g.length/2)+1);
					int my = Math.max(0,y-(g.length/2)+1);
					m = Math.max(m,my);
					int maxm = Math.min(f.length-1,x+(g.length/2)+1);
					int maxmy = Math.min(f.length-1,y+(g.length/2)+1);
					maxm = Math.min(maxm, maxmy);
					for (;m<maxm; m++) {
						tmpc[0][x][y] += g[x-m+(g.length/2)+1][y-m+(g.length/2)+1]*f[m];
					}
				}
			}
			for (int x=0; x<tmpc[0].length; x++) {
				for (int y=0; y<tmpc[0][x].length; y++) {
					cout.print(tmpc[0][x][y]+"\t");
				}
				cout.println();
			}
			 */
			for (int i=0; i<numinnerinters; i++) {
				double maxtmp = 0;
				double mintmp = Double.MAX_VALUE;
				double tmpsum = 0;
				for (int gene=0; gene<c.length; gene++) {
					double[][] tmpgene = c[gene];
					double[][] tmp1 = new double[tmpgene.length][tmpgene[0].length];
					//for (int x=2*g.length; x<tmpgene.length-(2*g.length); x++) {
					for (int x=0; x<tmpgene.length; x++) {
						double[] tmpgenex = tmpgene[x];
						//for (int y=2*g.length; y<tmpgenex.length-(2*g.length); y++) {
						for (int y=0; y<tmpgenex.length; y++) {	
							/*
							int m = Math.max(0,-(g.length/2)+x);
							m = Math.max(m, -(g[0].length/2)+y);
							int maxm = Math.min(oldf.length, (g.length/2)+x+1);
							maxm = Math.min(maxm, (g[0].length)/2+y+1);
							 */
							int m = 0;
							int maxm = f.length;
							for (; m<maxm; m++) {
								if (x-m>=-oldg.length/2 && y-m>=-oldg.length/2 && x-m<g.length/2 && y-m<g.length/2) {
									tmp1[x][y] += g[x-m+(g.length/2)][y-m+(g.length/2)]*oldf[m];
								}
							}
							//							if (tmpgene[x][y]==0 && tmp1[x][y]==0) {
							if (tmp1[x][y]==0) {
								tmp1[x][y] = 1;
							} else {

								if (tmp1[x][y]==0) {
									throw new Exception(tmpgene[x][y]+" "+x+" "+y);
								}

								tmp1[x][y] = tmpgene[x][y] / tmp1[x][y];
							}
						}
					}
					tmpc[gene] = tmp1;
				}
				int x = foffset;
				int maxx = f.length-foffset;
				for (; x<maxx; x++) {
					double tmp = 0;
					boolean nan = false;
					//int count = 0;
					for (int gene=0; gene<c.length; gene++) {
						double[][] tmp1 = tmpc[gene];
						/*
						int m = -(Math.min(tmp1.length-x-1, g.length/2-1));
						int my = -(Math.min(tmp1[0].length-x-1, g.length/2-1));
						int miny = my;
						//m = Math.max(m, my);
						int maxm = Math.min(x, g.length/2+1);
						int maxmy = Math.min(x, g.length/2+1);
						//maxm = Math.min(maxm,maxmy);
						//if (x==0) System.err.println(m+" "+maxm+" "+my+" "+maxmy);
						 */
						/*
						int m = Math.max(-f.length+2*g.length+x+1, -g.length/2 +1 +1);
						int my = Math.max(-f.length + 2*g.length + x + 1, -g.length/2 + 1 + 1);
						int miny = my;
						int maxm = Math.min(-2*g.length + x + 1,g.length/2 + 1 + 1);
						int maxmy = Math.min(-2*g.length + x + 1, g.length / 2 + 1 + 1);
						 */
						int m = Math.max(-f.length+x+1, -g.length/2 +1 +1);
						int my = Math.max(-f.length + x + 1, -g.length/2 + 1 + 1);
						int miny = my;
						int maxm = Math.min(x + 1,g.length/2 + 1 + 1);
						int maxmy = Math.min(x + 1, g.length / 2 + 1 + 1);
						for (; m<maxm; m++) {
							for (my=miny;my<maxmy; my++) {
								//if (x==0 && g[-m+(g.length/2)+1][-my+(g.length/2)+1]>0) System.err.println(m);
								//if (x==0) System.err.println(m+" "+tmp1[x-m][x-my]+" "+g[-m+(g.length/2)+1][-my+(g.length/2)+1]);
								//if (x-m>=2*g.length && x-my>=2*g.length && -m>=-(g.length/2)-1 && -my>=-(g.length/2)-1 && -m<(g.length/2)-1 && -my<(g.length/2)-1 && x-m<f.length-(2*g.length) && x-my<f.length-(2*g.length)) {
								if (x-m>=0 && x-my>=0 && -m>=-(g.length/2)-1 && -my>=-(g.length/2)-1 && -m<(g.length/2)-1 && -my<(g.length/2)-1 && x-m<f.length && x-my<f.length) {
									tmp += tmp1[x-m][x-my]*g[-m+(g.length/2)+1][-my+(g.length/2)+1];
									//count++;
								}
								if (Double.isNaN(tmp) && !nan) {
									System.err.println(tmp1[x-m][x-my]);
									nan = true;
								}
							}
						}
					}
					//if (tmp==0) System.err.println(x);
					f[x] = tmp*oldf[x];
					//System.err.println("count: "+count);
				}
				for (x = 0; x<f.length; x++) {
					double tmp = f[x]/oldf[x];
					if (tmp>maxtmp) maxtmp = tmp;
					if (tmp<mintmp) mintmp = tmp;
					tmpsum += Math.abs(Math.log(tmp));
					oldf[x] = f[x];
				}

				for (x = 0; x<f.length; x++) {
					fout.print(f[x]+"\t");
				}
				fout.println();

				System.err.println("f maxtmp: "+maxtmp+" f mintmp: "+mintmp+" f tmpsum: "+tmpsum);
				System.err.println((k*numinnerinters*2+numinnerinters)+i+" iterations "+dfm.format(new Date()));
				if (Math.abs(lastmaxtmp-maxtmp)<=.001 && Math.abs(lastmintmp-mintmp)<=.001) {
					break;
				} else {
					lastmaxtmp = maxtmp;
					lastmintmp = mintmp;
				}
			}
			/*
			tmpc[0] = new double[tmpc[0].length][tmpc[0][0].length];
			for (int x=0; x<tmpc[0].length; x++) {
				for (int y=0; y<tmpc[0][x].length; y++) {
					int m = Math.max(0,x-(g.length/2)+1);
					int my = Math.max(0,y-(g.length/2)+1);
					m = Math.max(m,my);
					int maxm = Math.min(f.length-1,x+(g.length/2)+1);
					int maxmy = Math.min(f.length-1,y+(g.length/2)+1);
					maxm = Math.min(maxm, maxmy);
					for (;m<maxm; m++) {
						tmpc[0][x][y] += g[x-m+(g.length/2)+1][y-m+(g.length/2)+1]*f[m];
					}
				}
			}
			for (int x=0; x<tmpc[0].length; x++) {
				for (int y=0; y<tmpc[0][x].length; y++) {
					cout.print(tmpc[0][x][y]+"\t");
				}
				cout.println();
			}
			 */
		}
		if (!nonblind) {
			for (int x=0; x<g.length; x++) {
				for (int y=0; y<g[x].length; y++) {
					gout.print(g[x][y]+"\t");
				}
				gout.println();
			}
		}
		for (int x = 0; x<f.length; x++) {
			fout.print(f[x]+"\t");
		}
		fout.println();
		gout.flush();
		gout.close();
		fout.flush();
		fout.close();
		cout.flush();
		cout.close();
		//tmpout.close();
	}

	public void execute(PrintStream fout) throws Exception {
		double lastmaxtmp = Double.MAX_VALUE;
		double lastmintmp = Double.MAX_VALUE;
		for (int k=0; k<numouterinters; k++) {
			if (!nonblind) {
				for (int i=0; i<numinnerinters; i++) {
					double maxtmp = 0;
					double mintmp = Double.MAX_VALUE;
					double tmpsum = 0;
					for (int gene=0; gene<c.length; gene++) {
						double[][] tmpgene = c[gene];
						double[][] tmp1 = new double[tmpgene.length][tmpgene[0].length];
						//for (int x=2*g.length; x<tmpgene.length-(2*g.length); x++) {
						for (int x=0; x<tmpgene.length; x++) {
							double[] tmpgenex = tmpgene[x];
							//for (int y=2*g.length; y<tmpgenex.length-(2*g.length); y++) {
							for (int y=0; y<tmpgenex.length; y++) {
								/*
							int m = Math.max(0,-(oldg.length/2)+x);
							m = Math.max(m, -(oldg[0].length/2)+y);
							int maxm = Math.min(f.length, (oldg.length/2)+x+1);
							maxm = Math.min(maxm, (oldg[0].length)/2+y+1);
								 */
								int m = 0;
								int maxm = f.length;
								for (; m<maxm; m++) {
									if (x-m>=-oldg.length/2 && y-m>=-oldg.length/2 && x-m<g.length/2 && y-m<g.length/2) {
										tmp1[x][y] += oldg[x-m+(oldg.length/2)][y-m+(oldg.length/2)]*f[m];
									}
								}
								//							if (tmpgene[x][y]==0 && tmp1[x][y]==0) {
								if (tmp1[x][y]==0) {
									tmp1[x][y] = 1;
								} else {
									tmp1[x][y] = tmpgene[x][y] / tmp1[x][y];
								}
							}
						}
						tmpc[gene] = tmp1;
					}
					int x = -(g.length/2)-1;
					int maxx = (g.length/2);
					for (;x<maxx; x++) {
						int y = -(g.length/2)-1;
						int ymax = (g.length/2);
						for (; y<ymax; y++) {
							double tmp = 0;
							for (int gene=0; gene<c.length; gene++) {
								double[][] tmp1 = tmpc[gene];
								/*
							int m = -(Math.min(tmp1.length-1, tmp1.length-x-1));
							int my = -(Math.min(tmp1[0].length-1, tmp1[0].length-y-1));
							m = Math.max(m, my);
							if (x<0) {
								m = Math.max(m, -x);
							}
							if (y<0) {
								m = Math.max(m,-y);
							}
							for (;m<=0; m++) {
								try {
									tmp += tmp1[x-m][y-m]*f[-m];
								} catch (Exception e) {
									System.err.println(x+"\t"+y+"\t"+m);
									throw e;
								}
							}
								 */
								for (int m=-f.length; m<f.length; m++) {
									//if ((x-m>=2*g.length) && (y-m>=2*g.length) && (-m>=0) && (x-m<f.length-(2*g.length)) && (y-m<f.length-(2*g.length)) && (-m<f.length)) {
									if ((x-m>=0) && (y-m>=0) && (-m>=0) && (x-m<f.length) && (y-m<f.length) && (-m<f.length)) {
										tmp += tmp1[x-m][y-m]*f[-m];
									}
								}
							}
							g[x+(g.length/2)+1][y+(g.length/2)+1] = tmp*oldg[x+(g.length/2)+1][y+(g.length)/2+1];
						}
					}
					for (x=0; x<g.length; x++) {
						for (int y=0; y<g[x].length; y++) {
							double tmp = g[x][y] / oldg[x][y];
							if (tmp>maxtmp) maxtmp = tmp;
							//tmpout.print(tmp+"\t");
							oldg[x][y] = g[x][y];
						}
						//tmpout.println();
					}
					//tmpout.flush();
					/*
				for (x=0; x<g.length; x++) {
					for (int y=0; y<g[x].length; y++) {
						gout.print(g[x][y]+"\t");
					}
					gout.println();
				}
				gout.flush();
					 */
					System.err.println("g maxtmp: "+maxtmp);
					System.err.println((k*numinnerinters*2)+i+" iterations "+dfm.format(new Date()));
					if (maxtmp<=1.1) break;
				}
			}
			/*
			for (int x=0; x<g.length; x++) {
				for (int y=0; y<g[x].length; y++) {
					gout.print(g[x][y]+"\t");
				}
				gout.println();
			}
			tmpc[0] = new double[tmpc[0].length][tmpc[0][0].length];
			for (int x=0; x<tmpc[0].length; x++) {
				for (int y=0; y<tmpc[0][x].length; y++) {
					int m = Math.max(0,x-(g.length/2)+1);
					int my = Math.max(0,y-(g.length/2)+1);
					m = Math.max(m,my);
					int maxm = Math.min(f.length-1,x+(g.length/2)+1);
					int maxmy = Math.min(f.length-1,y+(g.length/2)+1);
					maxm = Math.min(maxm, maxmy);
					for (;m<maxm; m++) {
						tmpc[0][x][y] += g[x-m+(g.length/2)+1][y-m+(g.length/2)+1]*f[m];
					}
				}
			}
			for (int x=0; x<tmpc[0].length; x++) {
				for (int y=0; y<tmpc[0][x].length; y++) {
					cout.print(tmpc[0][x][y]+"\t");
				}
				cout.println();
			}
			 */
			for (int i=0; i<numinnerinters; i++) {
				double maxtmp = 0;
				double mintmp = Double.MAX_VALUE;
				double tmpsum = 0;
				for (int gene=0; gene<c.length; gene++) {
					double[][] tmpgene = c[gene];
					double[][] tmp1 = new double[tmpgene.length][tmpgene[0].length];
					//for (int x=2*g.length; x<tmpgene.length-(2*g.length); x++) {
					for (int x=0; x<tmpgene.length; x++) {
						double[] tmpgenex = tmpgene[x];
						//for (int y=2*g.length; y<tmpgenex.length-(2*g.length); y++) {
						for (int y=0; y<tmpgenex.length; y++) {	
							/*
							int m = Math.max(0,-(g.length/2)+x);
							m = Math.max(m, -(g[0].length/2)+y);
							int maxm = Math.min(oldf.length, (g.length/2)+x+1);
							maxm = Math.min(maxm, (g[0].length)/2+y+1);
							 */
							int m = 0;
							int maxm = f.length;
							for (; m<maxm; m++) {
								if (x-m>=-oldg.length/2 && y-m>=-oldg.length/2 && x-m<g.length/2 && y-m<g.length/2) {
									tmp1[x][y] += g[x-m+(g.length/2)][y-m+(g.length/2)]*oldf[m];
								}
							}
							//							if (tmpgene[x][y]==0 && tmp1[x][y]==0) {
							if (tmp1[x][y]==0) {
								tmp1[x][y] = 1;
							} else {

								if (tmp1[x][y]==0) {
									throw new Exception(tmpgene[x][y]+" "+x+" "+y);
								}

								tmp1[x][y] = tmpgene[x][y] / tmp1[x][y];
							}
						}
					}
					tmpc[gene] = tmp1;
				}
				int x = foffset;
				int maxx = f.length-foffset;
				for (; x<maxx; x++) {
					double tmp = 0;
					boolean nan = false;
					//int count = 0;
					for (int gene=0; gene<c.length; gene++) {
						double[][] tmp1 = tmpc[gene];
						/*
						int m = -(Math.min(tmp1.length-x-1, g.length/2-1));
						int my = -(Math.min(tmp1[0].length-x-1, g.length/2-1));
						int miny = my;
						//m = Math.max(m, my);
						int maxm = Math.min(x, g.length/2+1);
						int maxmy = Math.min(x, g.length/2+1);
						//maxm = Math.min(maxm,maxmy);
						//if (x==0) System.err.println(m+" "+maxm+" "+my+" "+maxmy);
						 */
						/*
						int m = Math.max(-f.length+2*g.length+x+1, -g.length/2 +1 +1);
						int my = Math.max(-f.length + 2*g.length + x + 1, -g.length/2 + 1 + 1);
						int miny = my;
						int maxm = Math.min(-2*g.length + x + 1,g.length/2 + 1 + 1);
						int maxmy = Math.min(-2*g.length + x + 1, g.length / 2 + 1 + 1);
						 */
						int m = Math.max(-f.length+x+1, -g.length/2 +1 +1);
						int my = Math.max(-f.length + x + 1, -g.length/2 + 1 + 1);
						int miny = my;
						int maxm = Math.min(x + 1,g.length/2 + 1 + 1);
						int maxmy = Math.min(x + 1, g.length / 2 + 1 + 1);
						for (; m<maxm; m++) {
							for (my=miny;my<maxmy; my++) {
								//if (x==0 && g[-m+(g.length/2)+1][-my+(g.length/2)+1]>0) System.err.println(m);
								//if (x==0) System.err.println(m+" "+tmp1[x-m][x-my]+" "+g[-m+(g.length/2)+1][-my+(g.length/2)+1]);
								//if (x-m>=2*g.length && x-my>=2*g.length && -m>=-(g.length/2)-1 && -my>=-(g.length/2)-1 && -m<(g.length/2)-1 && -my<(g.length/2)-1 && x-m<f.length-(2*g.length) && x-my<f.length-(2*g.length)) {
								if (x-m>=0 && x-my>=0 && -m>=-(g.length/2)-1 && -my>=-(g.length/2)-1 && -m<(g.length/2)-1 && -my<(g.length/2)-1 && x-m<f.length && x-my<f.length) {
									tmp += tmp1[x-m][x-my]*g[-m+(g.length/2)+1][-my+(g.length/2)+1];
									//count++;
								}
								if (Double.isNaN(tmp) && !nan) {
									System.err.println(tmp1[x-m][x-my]);
									nan = true;
								}
							}
						}
					}
					//if (tmp==0) System.err.println(x);
					f[x] = tmp*oldf[x];
					//System.err.println("count: "+count);
				}
				for (x = 0; x<f.length; x++) {
					double tmp = f[x]/oldf[x];
					if (tmp>maxtmp) maxtmp = tmp;
					if (tmp<mintmp) mintmp = tmp;
					tmpsum += Math.abs(Math.log(tmp));
					oldf[x] = f[x];
				}

				for (x = 0; x<f.length; x++) {
					fout.print(f[x]+"\t");
				}
				fout.println();

				System.err.println("f maxtmp: "+maxtmp+" f mintmp: "+mintmp+" f tmpsum: "+tmpsum);
				System.err.println((k*numinnerinters*2+numinnerinters)+i+" iterations "+dfm.format(new Date()));
				if (Math.abs(lastmaxtmp-maxtmp)<=.001 && Math.abs(lastmintmp-mintmp)<=.001) {
					break;
				} else {
					lastmaxtmp = maxtmp;
					lastmintmp = mintmp;
				}
			}
			/*
			tmpc[0] = new double[tmpc[0].length][tmpc[0][0].length];
			for (int x=0; x<tmpc[0].length; x++) {
				for (int y=0; y<tmpc[0][x].length; y++) {
					int m = Math.max(0,x-(g.length/2)+1);
					int my = Math.max(0,y-(g.length/2)+1);
					m = Math.max(m,my);
					int maxm = Math.min(f.length-1,x+(g.length/2)+1);
					int maxmy = Math.min(f.length-1,y+(g.length/2)+1);
					maxm = Math.min(maxm, maxmy);
					for (;m<maxm; m++) {
						tmpc[0][x][y] += g[x-m+(g.length/2)+1][y-m+(g.length/2)+1]*f[m];
					}
				}
			}
			for (int x=0; x<tmpc[0].length; x++) {
				for (int y=0; y<tmpc[0][x].length; y++) {
					cout.print(tmpc[0][x][y]+"\t");
				}
				cout.println();
			}
			 */
		}
	}

}


