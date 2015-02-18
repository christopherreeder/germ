package edu.mit.csail.cgs.reeder.germ;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class PairedReadDistribution {

	private DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	public static final int SELF_LIGATION_DISTANCE = 10000;
	public static final int MINUS = 0;
	public static final int PLUS = 1;
	public static final int MINUSMINUS = 0;
	public static final int MINUSPLUS = 1;
	public static final int PLUSMINUS = 2;
	public static final int PLUSPLUS = 3;
	public static final int[] ORIENTATIONS = {MINUSMINUS, MINUSPLUS, PLUSMINUS, PLUSPLUS};
	public static final int NUM_EVENT_TYPES = 3;
	public static final int NUM_ORIENTATIONS = 4;
	public static final int SELF_INDICATOR = 0;
	public static final int PROX_INDICATOR = 1;
	public static final int DISTAL_INDICATOR = 2;
	public static final int POSTSS = 0;
	public static final int NEGTSS = 1;
	public static final int NONTSS = 2;
	public static final int[] TYPES = {POSTSS, NEGTSS, NONTSS};
	public int radius;
	private float[][][] selfdist; // event type x position x position
	private float[][][][] proxinterdist; // event type x strand orientation x position x position
	private float[][][] distalinterdist; // event type x strand orientation x position

	public static void main(String[] args) throws NotFoundException, IOException {
		/*
		Genome g = Args.parseGenome(args).cdr();
		int window = Args.parseInteger(args, "window", 0);
		String annofile = Args.parseString(args, "annofile", "");
		String peakfile = Args.parseString(args, "peakfile", "");
		String directory = Args.parseString(args, "directory", "");
		String readfile = Args.parseString(args, "readfile", "");
		int size = Args.parseInteger(args, "size", 0);
		SproutStorage storage = SproutStorage.fromFile(g, readfile, size);
		PairedReadDistribution dist = new PairedReadDistribution(Args.parseInteger(args, "radius", 0));
		dist.initializeFromAnnotation(annofile, peakfile, storage, g, window);
		dist.writeToDirectory(directory);
		*/
		/*
		Genome g = Args.parseGenome(args).cdr();
		int window = Args.parseInteger(args, "window", 0);
		int numTSSs = Args.parseInteger(args, "numTSSs", 0);
		String annofile = Args.parseString(args, "annofile", "");
		String peakfile = Args.parseString(args, "peakfile", "");
		String directory = Args.parseString(args, "directory", "");
		String readfile = Args.parseString(args, "readfile", "");
		int size = Args.parseInteger(args, "size", 0);
		SproutStorage storage = SproutStorage.fromFile(g, readfile, size);
		PairedReadDistribution dist = new PairedReadDistribution(Args.parseInteger(args, "radius", 0));
		dist.initializeFromAnnotationStrong(annofile, peakfile, storage, g, window, numTSSs);
		dist.writeToDirectory(directory);
		*/
		/*
		PairedReadDistribution dist = new PairedReadDistribution(Args.parseInteger(args, "radius", 0));
		dist.initializeFromDirectory(Args.parseString(args, "directory", ""));
		dist.smooth(Args.parseInteger(args, "window", 0));
		dist.writeToDirectory(Args.parseString(args, "newdirectory", ""));
		*/
		
		/*
		PairedReadDistribution dist = new PairedReadDistribution(Args.parseInteger(args, "radius", 0));
		dist.initializeFromDirectory(Args.parseString(args, "directory", ""));
		dist.normalize();
		dist.writeToDirectory(Args.parseString(args, "newdirectory", ""));
		*/
		
		/*
		Genome g = Args.parseGenome(args).cdr();
		String annofile = Args.parseString(args, "annofile", "");
		String peakfile = Args.parseString(args, "peakfile", "");
		String directory = Args.parseString(args, "directory", "");
		String readfile = Args.parseString(args, "readfile", "");
		int window = Args.parseInteger(args, "window", 0);
		int size = Args.parseInteger(args, "size", 0);
		SproutStorage storage = SproutStorage.fromFile(g, readfile, size);
		PairedReadDistribution dist = new PairedReadDistribution(Args.parseInteger(args, "radius", 0));
		dist.initializeFromAnnotationStage1(annofile, peakfile, storage, g, window);
		dist.writeToDirectoryStage1(directory);
		*/
		
		Genome g = Args.parseGenome(args).cdr();
		String peakfile = Args.parseString(args, "peakfile", "");
		String directory = Args.parseString(args, "directory", "");
		String readfile = Args.parseString(args, "readfile", "");
		int window = Args.parseInteger(args, "window", 0);
		int size = Args.parseInteger(args, "size", 0);
		SproutStorage storage = SproutStorage.fromFileEven(g, readfile);
		PairedReadDistribution dist = new PairedReadDistribution(Args.parseInteger(args, "radius", 0));
		dist.initializeFromPeaks(peakfile, storage, g, window);
		dist.writeToDirectoryStage1(directory);
		
	}
	
	public PairedReadDistribution(int radius) {
		this.radius = radius;
	}
	
	public void initializeFromPeaks(String peakfile, SproutStorage storage, Genome g, int window) throws IOException {
		Set<StrandedPoint> peaks = new HashSet<StrandedPoint>();
		int diameter = 2*radius;
		selfdist = new float[1][diameter][diameter];
		proxinterdist = new float[1][4][diameter][diameter];
		distalinterdist = new float[1][4][diameter];
		BufferedReader r = new BufferedReader(new FileReader(peakfile));
		String s;
		String[] split;
		r.readLine();
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			Point p = Point.fromString(g, split[0]);
			StrandedPoint tmp = new StrandedPoint(p,'+');
			peaks.add(tmp);
		}
		r.close();
		computeSelfDist(peaks, selfdist[0], storage);
		selfdist[0] = smooth(selfdist[0],window);
		normalize(selfdist[0]);
		for (int i=0; i<ORIENTATIONS.length; i++) {
			computeProxDist(peaks, proxinterdist[0][i], storage, i);
			proxinterdist[0][i] = smooth(proxinterdist[0][i],window);
			normalize(proxinterdist[0][i]);
		}
		computeDistalDist(peaks, distalinterdist[0][MINUS], storage, false);
		computeDistalDist(peaks, distalinterdist[0][PLUS], storage, true);
		distalinterdist[0][MINUS] = smooth(distalinterdist[0][MINUS],window);
		normalize(distalinterdist[0][MINUS]);
		distalinterdist[0][PLUS] = smooth(distalinterdist[0][PLUS],window);
		normalize(distalinterdist[0][PLUS]);
	}
	
	public void initializeFromAnnotationStage1(String annofile, String peakfile, SproutStorage storage, Genome g, int window) throws IOException {

		Set<String> uniqueTSSs = new HashSet<String>();
		Set<String> nonUniqueTSSs = new HashSet<String>();
		Set<StrandedPoint> posTSSs = new HashSet<StrandedPoint>();
		Set<StrandedPoint> negTSSs = new HashSet<StrandedPoint>();
		Set<StrandedPoint> nonTSSs = new HashSet<StrandedPoint>();
		SortedSet<StrandedPoint> sorted = new TreeSet<StrandedPoint>();
		BufferedReader r = new BufferedReader(new FileReader(annofile));
		int diameter = 2*radius;
		selfdist = new float[1][diameter][diameter];
		proxinterdist = new float[1][4][diameter][diameter];
		distalinterdist = new float[1][4][diameter];
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			if (nonUniqueTSSs.contains(split[1])) {
				continue;
			} else if (uniqueTSSs.contains(split[1])) {
				uniqueTSSs.remove(split[1]);
				nonUniqueTSSs.add(split[1]);
			} else {
				uniqueTSSs.add(split[1]);
			}
			StrandedPoint tmp = StrandedPoint.fromString(g, split[0]);
			if (tmp==null) {
				System.err.println(s);
			} else {
				sorted.add(tmp);
			}
		}
		r.close();
		System.err.println("unique TSSs found "+dfm.format(new Date()));
		r = new BufferedReader(new FileReader(annofile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			if (uniqueTSSs.contains(split[1])) {
				StrandedPoint tss = StrandedPoint.fromString(g, split[0]);
				if (tss==null) {
					System.err.println(s);
				} else {
					if (tss.getStrand()=='+') {
						posTSSs.add(tss);
					} else {
						negTSSs.add(tss);
					}
				}
			}
		}
		r.close();
		System.err.println("unique TSSs parsed "+dfm.format(new Date()));
		r = new BufferedReader(new FileReader(peakfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			Point p = Point.fromString(g, split[2]);
			if (p==null) {
				System.err.println(s);
			} else {
				StrandedPoint tmp = new StrandedPoint(p,'+');
				StrandedPoint min = new StrandedPoint(g, tmp.getChrom(), tmp.getLocation()-diameter, '+');
				StrandedPoint max = new StrandedPoint(g, tmp.getChrom(), tmp.getLocation()+diameter, '+');
				SortedSet<StrandedPoint> subset = sorted.subSet(min, max);
				if (subset.isEmpty()) {
					nonTSSs.add(tmp);
				}
			}
		}
		r.close();
		System.err.println("nonTSSs found "+dfm.format(new Date()));

		computeSelfDist(posTSSs, selfdist[0], storage);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeSelfDist(negTSSs, selfdist[0], storage);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeSelfDist(nonTSSs, selfdist[0], storage);
		System.err.println("dist computed "+dfm.format(new Date()));
		selfdist[0] = smooth(selfdist[0],window);
		normalize(selfdist[0]);
		for (int i=0; i<ORIENTATIONS.length; i++) {
			computeProxDist(posTSSs, proxinterdist[0][i], storage, i);
			System.err.println("dist computed "+dfm.format(new Date()));
			computeProxDist(negTSSs, proxinterdist[0][i], storage, i);
			System.err.println("dist computed "+dfm.format(new Date()));
			computeProxDist(nonTSSs, proxinterdist[0][i], storage, i);
			System.err.println("dist computed "+dfm.format(new Date()));
			proxinterdist[0][i] = smooth(proxinterdist[0][i],window);
			normalize(proxinterdist[0][i]);
		}
		computeDistalDist(posTSSs, distalinterdist[0][MINUS], storage, false);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(posTSSs, distalinterdist[0][PLUS], storage, true);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(negTSSs, distalinterdist[0][MINUS], storage, false);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(negTSSs, distalinterdist[0][PLUS], storage, true);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(nonTSSs, distalinterdist[0][MINUS], storage, false);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(nonTSSs, distalinterdist[0][PLUS], storage, true);
		System.err.println("dist computed "+dfm.format(new Date()));
		distalinterdist[0][MINUS] = smooth(distalinterdist[0][MINUS],window);
		normalize(distalinterdist[0][MINUS]);
		distalinterdist[0][PLUS] = smooth(distalinterdist[0][PLUS],window);
		normalize(distalinterdist[0][PLUS]);
	}

	public void initializeFromAnnotation(String annofile, String peakfile, SproutStorage storage, Genome g, int window) throws IOException {

		Set<String> uniqueTSSs = new HashSet<String>();
		Set<String> nonUniqueTSSs = new HashSet<String>();
		Set<StrandedPoint> posTSSs = new HashSet<StrandedPoint>();
		Set<StrandedPoint> negTSSs = new HashSet<StrandedPoint>();
		Set<StrandedPoint> nonTSSs = new HashSet<StrandedPoint>();
		SortedSet<StrandedPoint> sorted = new TreeSet<StrandedPoint>();
		BufferedReader r = new BufferedReader(new FileReader(annofile));
		int diameter = 2*radius;
		selfdist = new float[3][diameter][diameter];
		proxinterdist = new float[3][4][diameter][diameter];
		distalinterdist = new float[3][4][diameter];
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			if (nonUniqueTSSs.contains(split[1])) {
				continue;
			} else if (uniqueTSSs.contains(split[1])) {
				uniqueTSSs.remove(split[1]);
				nonUniqueTSSs.add(split[1]);
			} else {
				uniqueTSSs.add(split[1]);
			}
			StrandedPoint tmp = StrandedPoint.fromString(g, split[0]);
			if (tmp==null) {
				System.err.println(s);
			} else {
				sorted.add(tmp);
			}
		}
		r.close();
		System.err.println(uniqueTSSs.size()+" unique TSSs found "+dfm.format(new Date()));
		r = new BufferedReader(new FileReader(annofile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			if (uniqueTSSs.contains(split[1])) {
				StrandedPoint tss = StrandedPoint.fromString(g, split[0]);
				if (tss==null) {
					System.err.println(s);
				} else {
					if (tss.getStrand()=='+') {
						posTSSs.add(tss);
					} else {
						negTSSs.add(tss);
					}
				}
			}
		}
		r.close();
		System.err.println("unique TSSs parsed "+dfm.format(new Date()));
		r = new BufferedReader(new FileReader(peakfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			Point p = Point.fromString(g, split[2]);
			if (p==null) {
				System.err.println(s);
			} else {
				StrandedPoint tmp = new StrandedPoint(p,'+');
				StrandedPoint min = new StrandedPoint(g, tmp.getChrom(), tmp.getLocation()-diameter, '+');
				StrandedPoint max = new StrandedPoint(g, tmp.getChrom(), tmp.getLocation()+diameter, '+');
				SortedSet<StrandedPoint> subset = sorted.subSet(min, max);
				if (subset.isEmpty()) {
					nonTSSs.add(tmp);
				}
			}
		}
		r.close();
		System.err.println("nonTSSs found "+dfm.format(new Date()));

		computeSelfDist(posTSSs, selfdist[POSTSS], storage);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeSelfDist(negTSSs, selfdist[NEGTSS], storage);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeSelfDist(nonTSSs, selfdist[NONTSS], storage);
		System.err.println("dist computed "+dfm.format(new Date()));
		selfdist[POSTSS] = smooth(selfdist[POSTSS],window);
		normalize(selfdist[POSTSS]);
		selfdist[POSTSS] = smooth(selfdist[NEGTSS],window);
		normalize(selfdist[NEGTSS]);
		selfdist[POSTSS] = smooth(selfdist[NONTSS],window);
		normalize(selfdist[NONTSS]);
		for (int i=0; i<ORIENTATIONS.length; i++) {
			computeProxDist(posTSSs, proxinterdist[POSTSS][i], storage, i);
			System.err.println("dist computed "+dfm.format(new Date()));
			computeProxDist(negTSSs, proxinterdist[NEGTSS][i], storage, i);
			System.err.println("dist computed "+dfm.format(new Date()));
			computeProxDist(nonTSSs, proxinterdist[NONTSS][i], storage, i);
			System.err.println("dist computed "+dfm.format(new Date()));
			proxinterdist[POSTSS][i] = smooth(proxinterdist[POSTSS][i],window);
			normalize(proxinterdist[POSTSS][i]);
			proxinterdist[NEGTSS][i] = smooth(proxinterdist[NEGTSS][i],window);
			normalize(proxinterdist[NEGTSS][i]);
			proxinterdist[NONTSS][i] = smooth(proxinterdist[NONTSS][i],window);
			normalize(proxinterdist[NONTSS][i]);
		}
		computeDistalDist(posTSSs, distalinterdist[POSTSS][MINUS], storage, false);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(posTSSs, distalinterdist[POSTSS][PLUS], storage, true);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(negTSSs, distalinterdist[NEGTSS][MINUS], storage, false);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(negTSSs, distalinterdist[NEGTSS][PLUS], storage, true);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(nonTSSs, distalinterdist[NONTSS][MINUS], storage, false);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(nonTSSs, distalinterdist[NONTSS][PLUS], storage, true);
		System.err.println("dist computed "+dfm.format(new Date()));
		distalinterdist[POSTSS][MINUS] = smooth(distalinterdist[POSTSS][MINUS],window);
		normalize(distalinterdist[POSTSS][MINUS]);
		distalinterdist[POSTSS][PLUS] = smooth(distalinterdist[POSTSS][PLUS],window);
		normalize(distalinterdist[POSTSS][PLUS]);
		distalinterdist[NEGTSS][MINUS] = smooth(distalinterdist[NEGTSS][MINUS],window);
		normalize(distalinterdist[NEGTSS][MINUS]);
		distalinterdist[NEGTSS][PLUS] = smooth(distalinterdist[NEGTSS][PLUS],window);
		normalize(distalinterdist[NEGTSS][PLUS]);
		distalinterdist[NONTSS][MINUS] = smooth(distalinterdist[NONTSS][MINUS],window);
		normalize(distalinterdist[NONTSS][MINUS]);
		distalinterdist[NONTSS][PLUS] = smooth(distalinterdist[NONTSS][PLUS],window);
		normalize(distalinterdist[NONTSS][PLUS]);
	}
	
	public void initializeFromAnnotationStrong(String annofile, String peakfile, SproutStorage storage, Genome g, int window, int numTSSs) throws IOException {

		Set<String> uniqueTSSs = new HashSet<String>();
		Set<String> nonUniqueTSSs = new HashSet<String>();
		SortedSet<Pair<StrandedPoint,Integer>> posTSSs = new TreeSet<Pair<StrandedPoint,Integer>>(new PeakComparator());
		SortedSet<Pair<StrandedPoint,Integer>> negTSSs = new TreeSet<Pair<StrandedPoint,Integer>>(new PeakComparator());
		Set<StrandedPoint> nonTSSs = new HashSet<StrandedPoint>();
		SortedSet<StrandedPoint> sorted = new TreeSet<StrandedPoint>();
		BufferedReader r = new BufferedReader(new FileReader(annofile));
		int diameter = 2*radius;
		selfdist = new float[3][diameter][diameter];
		proxinterdist = new float[3][4][diameter][diameter];
		distalinterdist = new float[3][4][diameter];
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			if (nonUniqueTSSs.contains(split[1])) {
				continue;
			} else if (uniqueTSSs.contains(split[1])) {
				uniqueTSSs.remove(split[1]);
				nonUniqueTSSs.add(split[1]);
			} else {
				uniqueTSSs.add(split[1]);
			}
			StrandedPoint tmp = StrandedPoint.fromString(g, split[0]);
			if (tmp==null) {
				System.err.println(s);
			} else {
				sorted.add(tmp);
			}
		}
		r.close();
		System.err.println(uniqueTSSs.size()+" unique TSSs found "+dfm.format(new Date()));
		r = new BufferedReader(new FileReader(annofile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			if (uniqueTSSs.contains(split[1])) {
				StrandedPoint tss = StrandedPoint.fromString(g, split[0]);
				if (tss==null) {
					System.err.println(s);
				} else {
					Pair<Integer,Integer> minmax = storage.getLeftMinMaxn(tss.expand(window));
					if (tss.getStrand()=='+') {
						posTSSs.add(new Pair<StrandedPoint,Integer>(tss,minmax.cdr()-minmax.car()));
						if (posTSSs.size()>numTSSs) {
							posTSSs.remove(posTSSs.first());
						}
					} else {
						negTSSs.add(new Pair<StrandedPoint,Integer>(tss,minmax.cdr()-minmax.car()));
						if (negTSSs.size()>numTSSs) {
							negTSSs.remove(negTSSs.first());
						}
					}
				}
			}
		}
		r.close();
		System.err.println("unique TSSs parsed "+dfm.format(new Date()));
		r = new BufferedReader(new FileReader(peakfile));
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			Point p = Point.fromString(g, split[2]);
			if (p==null) {
				System.err.println(s);
			} else {
				StrandedPoint tmp = new StrandedPoint(p,'+');
				StrandedPoint min = new StrandedPoint(g, tmp.getChrom(), tmp.getLocation()-diameter, '+');
				StrandedPoint max = new StrandedPoint(g, tmp.getChrom(), tmp.getLocation()+diameter, '+');
				SortedSet<StrandedPoint> subset = sorted.subSet(min, max);
				if (subset.isEmpty()) {
					nonTSSs.add(tmp);
				}
			}
		}
		r.close();
		System.err.println("nonTSSs found "+dfm.format(new Date()));
		
		Set<StrandedPoint> posTSSSet = new HashSet<StrandedPoint>();
		for (Pair<StrandedPoint,Integer> pair : posTSSs) {
			posTSSSet.add(pair.car());
		}
		Set<StrandedPoint> negTSSSet = new HashSet<StrandedPoint>();
		for (Pair<StrandedPoint,Integer> pair : negTSSs) {
			negTSSSet.add(pair.car());
		}

		computeSelfDist(posTSSSet, selfdist[POSTSS], storage);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeSelfDist(negTSSSet, selfdist[NEGTSS], storage);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeSelfDist(nonTSSs, selfdist[NONTSS], storage);
		System.err.println("dist computed "+dfm.format(new Date()));
		selfdist[POSTSS] = smooth(selfdist[POSTSS],window);
		normalize(selfdist[POSTSS]);
		selfdist[NEGTSS] = smooth(selfdist[NEGTSS],window);
		normalize(selfdist[NEGTSS]);
		selfdist[NONTSS] = smooth(selfdist[NONTSS],window);
		normalize(selfdist[NONTSS]);
		for (int i=0; i<ORIENTATIONS.length; i++) {
			computeProxDist(posTSSSet, proxinterdist[POSTSS][i], storage, i);
			System.err.println("dist computed "+dfm.format(new Date()));
			computeProxDist(negTSSSet, proxinterdist[NEGTSS][i], storage, i);
			System.err.println("dist computed "+dfm.format(new Date()));
			computeProxDist(nonTSSs, proxinterdist[NONTSS][i], storage, i);
			System.err.println("dist computed "+dfm.format(new Date()));
			proxinterdist[POSTSS][i] = smooth(proxinterdist[POSTSS][i],window);
			normalize(proxinterdist[POSTSS][i]);
			proxinterdist[NEGTSS][i] = smooth(proxinterdist[NEGTSS][i],window);
			normalize(proxinterdist[NEGTSS][i]);
			proxinterdist[NONTSS][i] = smooth(proxinterdist[NONTSS][i],window);
			normalize(proxinterdist[NONTSS][i]);
		}
		computeDistalDist(posTSSSet, distalinterdist[POSTSS][MINUS], storage, false);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(posTSSSet, distalinterdist[POSTSS][PLUS], storage, true);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(negTSSSet, distalinterdist[NEGTSS][MINUS], storage, false);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(negTSSSet, distalinterdist[NEGTSS][PLUS], storage, true);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(nonTSSs, distalinterdist[NONTSS][MINUS], storage, false);
		System.err.println("dist computed "+dfm.format(new Date()));
		computeDistalDist(nonTSSs, distalinterdist[NONTSS][PLUS], storage, true);
		System.err.println("dist computed "+dfm.format(new Date()));
		distalinterdist[POSTSS][MINUS] = smooth(distalinterdist[POSTSS][MINUS],window);
		normalize(distalinterdist[POSTSS][MINUS]);
		distalinterdist[POSTSS][PLUS] = smooth(distalinterdist[POSTSS][PLUS],window);
		normalize(distalinterdist[POSTSS][PLUS]);
		distalinterdist[NEGTSS][MINUS] = smooth(distalinterdist[NEGTSS][MINUS],window);
		normalize(distalinterdist[NEGTSS][MINUS]);
		distalinterdist[NEGTSS][PLUS] = smooth(distalinterdist[NEGTSS][PLUS],window);
		normalize(distalinterdist[NEGTSS][PLUS]);
		distalinterdist[NONTSS][MINUS] = smooth(distalinterdist[NONTSS][MINUS],window);
		normalize(distalinterdist[NONTSS][MINUS]);
		distalinterdist[NONTSS][PLUS] = smooth(distalinterdist[NONTSS][PLUS],window);
		normalize(distalinterdist[NONTSS][PLUS]);
	}
	
	public void writeToDirectoryStage1(String directory) throws FileNotFoundException {
		write2DDist(directory+"selfpos.txt", selfdist[0]);
		for (int i=0; i<ORIENTATIONS.length; i++) {
			write2DDist(directory+"proxpos"+i+".txt", proxinterdist[0][i]);
		}
		write1DDist(directory+"distalposminus.txt", distalinterdist[0][MINUS]);
		write1DDist(directory+"distalposplus.txt", distalinterdist[0][PLUS]);
	}

	public void writeToDirectory(String directory) throws FileNotFoundException {
		write2DDist(directory+"selfpos.txt", selfdist[POSTSS]);
		write2DDist(directory+"selfneg.txt", selfdist[NEGTSS]);
		write2DDist(directory+"selfnon.txt", selfdist[NONTSS]);
		for (int i=0; i<ORIENTATIONS.length; i++) {
			write2DDist(directory+"proxpos"+i+".txt", proxinterdist[POSTSS][i]);
			write2DDist(directory+"proxneg"+i+".txt", proxinterdist[NEGTSS][i]);
			write2DDist(directory+"proxnon"+i+".txt", proxinterdist[NONTSS][i]);
		}
		write1DDist(directory+"distalposminus.txt", distalinterdist[POSTSS][MINUS]);
		write1DDist(directory+"distalposplus.txt", distalinterdist[POSTSS][PLUS]);
		write1DDist(directory+"distalnegminus.txt", distalinterdist[NEGTSS][MINUS]);
		write1DDist(directory+"distalnegplus.txt", distalinterdist[NEGTSS][PLUS]);
		write1DDist(directory+"distalnonminus.txt", distalinterdist[NONTSS][MINUS]);
		write1DDist(directory+"distalnonplus.txt", distalinterdist[NONTSS][PLUS]);
	}

	public void write2DDist(String file, float[][] dist) throws FileNotFoundException {
		PrintStream out = new PrintStream(file);
		for (int i=0; i<dist.length; i++) {
			for (int j=0; j<dist.length; j++) {
				out.print(dist[i][j]+"\t");
			}
			out.println();
		}
		out.flush();
		out.close();
	}
	
	public float[][] read2DDist(String file) throws IOException {
		float[][] tor;
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		s = r.readLine();
		split = s.split("\t");
		tor = new float[split.length][split.length];
		for (int i=0; i<split.length; i++) {
			tor[0][i] = Float.valueOf(split[i]);
		}
		int index = 0;
		while ((s = r.readLine()) != null) {
			index++;
			split = s.split("\t");
			for (int i=0; i<split.length; i++) {
				tor[index][i] = Float.valueOf(split[i]);
			}
		}
		r.close();
		return tor;
	}

	public void write1DDist(String file, float[] dist) throws FileNotFoundException {
		PrintStream out = new PrintStream(file);
		for (int i=0; i<dist.length; i++) {
			out.println(dist[i]);
		}
		out.flush();
		out.close();
	}
	
	public float[] read1DDist(String file) throws IOException {
		float[] tor;
		ArrayList<Float> torlist = new ArrayList<Float>();
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			torlist.add(Float.valueOf(s));
		}
		tor = new float[torlist.size()];
		for (int i=0; i<torlist.size(); i++) {
			tor[i] = torlist.get(i);
		}
		r.close();
		return tor;
	}

	private void computeSelfDist(Set<StrandedPoint> points, float[][] dist, SproutStorage storage) {
		for (StrandedPoint p : points) {
			Region pregion = p.expand(radius);
			Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(pregion);
			for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
				Pair<StrandedPoint,StrandedPoint> pair = storage.getLeftPair(n);
				if (isSelfLigation(pair)) {
					int x = pair.car().getLocation()-p.getLocation();
					int y = pair.cdr().getLocation()-p.getLocation();
					if (Math.abs(x)<radius && Math.abs(y)<radius) {
						dist[x+radius][y+radius]++;
					}
				}
			}
		}
		/*
		float sum = 0;
		for (int i=0; i<dist.length; i++) {
			for (int j=0; j<dist[i].length; j++) {
				sum += dist[i][j];
			}
		}
		for (int i=0; i<dist.length; i++) {
			for (int j=0; j<dist[i].length; j++) {
				dist[i][j] /= sum;
			}
		}
		*/
	}

	private void computeProxDist(Set<StrandedPoint> points, float[][] dist, SproutStorage storage, int orientation) {
		for (StrandedPoint p : points) {
			Region pregion = p.expand(radius);
			Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(pregion);
			for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
				Pair<StrandedPoint,StrandedPoint> pair = storage.getLeftPair(n);
				if (pair.car().getChrom().equals(pair.cdr().getChrom()) && hasOrientation(pair, orientation)) {
					int x = pair.car().getLocation()-p.getLocation();
					int y = pair.cdr().getLocation()-p.getLocation();
					if (Math.abs(x)<radius && Math.abs(y)<radius) {
						dist[x+radius][y+radius]++;
					}
				}
			}
		}
		/*
		float sum = 0;
		for (int i=0; i<dist.length; i++) {
			for (int j=0; j<dist[i].length; j++) {
				sum += dist[i][j];
			}
		}
		for (int i=0; i<dist.length; i++) {
			for (int j=0; j<dist[i].length; j++) {
				dist[i][j] /= sum;
			}
		}
		*/
	}

	private void computeDistalDist(Set<StrandedPoint> points, float[] dist, SproutStorage storage, boolean orientation) {
		for (StrandedPoint p : points) {
			Set<Pair<StrandedPoint,StrandedPoint>> seen = new HashSet<Pair<StrandedPoint,StrandedPoint>>();
			Region pregion = p.expand(radius);
			Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(pregion);
			for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
				Pair<StrandedPoint,StrandedPoint> pair = storage.getLeftPair(n);
				seen.add(pair);
				if (orientation ? pair.car().getStrand()=='+' : pair.car().getStrand()=='-') {
					int x = pair.car().getLocation()-p.getLocation();
					if (Math.abs(x)<radius) {
						dist[x+radius]++;
					}
				}
			}
			minmaxn = storage.getRightMinMaxn(pregion);
			for (int n=minmaxn.car(); n<minmaxn.cdr(); n++) {
				Pair<StrandedPoint,StrandedPoint> pair = storage.getRightPair(n);
				if (!seen.contains(pair) && (orientation ? pair.cdr().getStrand()=='+' : pair.cdr().getStrand()=='-')) {
					int y = pair.cdr().getLocation()-p.getLocation();
					if (Math.abs(y)<radius) {
						dist[y+radius]++;
					}
				}
			}
		}
		/*
		float sum = 0;
		for (int i=0; i<dist.length; i++) {
			sum += dist[i];
		}
		for (int i=0; i<dist.length; i++) {
			dist[i] /= sum;
		}
		*/
	}
	
	public void smooth(int window) {
		float[][][] newselfdist = new float[NUM_EVENT_TYPES][selfdist[0].length-window][selfdist[0][0].length-window];
		newselfdist[0] = smooth(selfdist[0],window);
		newselfdist[1] = smooth(selfdist[1],window);
		newselfdist[2] = smooth(selfdist[2],window);
		selfdist = newselfdist;
		float[][][][] newproxinterdist = new float[proxinterdist.length][proxinterdist[0].length][proxinterdist[0][0].length][proxinterdist[0][0][0].length];
		for (int i=0; i<NUM_EVENT_TYPES; i++) {
			for (int j=0; j<NUM_ORIENTATIONS; j++) {
				newproxinterdist[i][j] = smooth(proxinterdist[i][j],window);
			}
		}
		proxinterdist = newproxinterdist;
		float[][][] newdistalinterdist = new float[distalinterdist.length][distalinterdist[0].length][distalinterdist[0][0].length];
		for (int i=0; i<NUM_EVENT_TYPES; i++) {
			for (int j=0; j<2; j++) {
				newdistalinterdist[i][j] = smooth(distalinterdist[i][j],window);
			}
		}
		distalinterdist = newdistalinterdist;
		radius = selfdist[0].length/2;
	}
	
	public void normalize() {
		normalize(selfdist[0]);
		normalize(selfdist[1]);
		normalize(selfdist[2]);
		for (int i=0; i<NUM_EVENT_TYPES; i++) {
			for (int j=0; j<NUM_ORIENTATIONS; j++) {
				normalize(proxinterdist[i][j]);
			}
		}
		for (int i=0; i<NUM_EVENT_TYPES; i++) {
			for (int j=0; j<2; j++) {
				normalize(distalinterdist[i][j]);
			}
		}
	}
	
	public static void normalize(float[] arr) {
		float sum = 0;
		for (int i=0; i<arr.length; i++) {
			sum += arr[i];
		}
		for (int i=0; i<arr.length; i++) {
			arr[i] /= sum;
		}
	}
	
	public static void normalize(float[][] arr) {
		float sum = 0;
		for (int i=0; i<arr.length; i++) {
			for (int j=0; j<arr[i].length; j++) {
				sum += arr[i][j];
			}
		}
		for (int i=0; i<arr.length; i++) {
			for (int j=0; j<arr[i].length; j++) {
				arr[i][j] /= sum;
			}
		}
	}
	
	public static float[] smooth(float[] arr, int window) {
		float[] tor = new float[arr.length-window];
		for (int i=0; i<tor.length; i++) {
			float sum = 0;
			for (int j=i; j<i+window; j++) {
				sum += arr[j];
			}
			tor[i] = sum / ((float)window);
		}
		return tor;
	}
	
	public static float[][] smooth(float[][] arr, int window) {
		float[][] tor = new float[arr.length-window][arr[0].length-window];
		for (int i=0; i<tor.length; i++) {
			for (int j=0; j<tor[i].length; j++) {
				float sum = 0;
				for (int k=i; k<i+window; k++) {
					for (int l=j; l<j+window; l++) {
						sum += arr[k][l];
					}
				}
				tor[i][j] = sum/((float)(window*window));
			}
		}
		return tor;
	}

	public float probability(Pair<StrandedPoint,StrandedPoint> pet, Point mu1, Point mu2, int tau1, int tau2) {
		if (mu2==null) {
			//self-ligation
			if (!isSelfLigation(pet)) {
				return 0;
			} else {
				return prob(pet.car().getLocation()-mu1.getLocation(),pet.cdr().getLocation()-mu1.getLocation(),selfdist[tau1]);
			}
		} else if (mu1.equals(mu2)) {
			//proximal inter-ligation
			return prob(pet.car().getLocation()-mu1.getLocation(),pet.cdr().getLocation()-mu1.getLocation(),proxinterdist[tau1][strandOrientation(pet)]);
		} else {
			//distal inter-ligation
			return prob(pet.car().getLocation()-mu1.getLocation(),distalinterdist[tau1][strandValue(pet.car().getStrand())])*
					prob(pet.cdr().getLocation()-mu2.getLocation(),distalinterdist[tau2][strandValue(pet.cdr().getStrand())]);
		}
	}
	
	public float probability(StrandedPoint read, Point mu, int tau) {
		return prob(read.getLocation()-mu.getLocation(),distalinterdist[tau][strandValue(read.getStrand())]);
	}
	
	public void initializeFromDirectoryStage1(String directory) throws IOException {
		float[][] selfpos = read2DDist(directory+"selfpos.txt");
		radius = selfpos.length/2;
		selfdist = new float[1][selfpos.length][selfpos.length];
		proxinterdist = new float[1][NUM_ORIENTATIONS][selfpos.length][selfpos.length];
		distalinterdist = new float[1][2][selfpos.length];
		selfdist[0] = selfpos;
		for (int i=0; i<ORIENTATIONS.length; i++) {
			proxinterdist[0][i] = read2DDist(directory+"proxpos"+i+".txt");
		}
		distalinterdist[0][MINUS] = read1DDist(directory+"distalposminus.txt");
		distalinterdist[0][PLUS] = read1DDist(directory+"distalposplus.txt");
	}

	public void initializeFromDirectory(String directory) throws IOException {
		float[][] selfpos = read2DDist(directory+"selfpos.txt");
		radius = selfpos.length/2;
		selfdist = new float[NUM_EVENT_TYPES][selfpos.length][selfpos.length];
		proxinterdist = new float[NUM_EVENT_TYPES][NUM_ORIENTATIONS][selfpos.length][selfpos.length];
		distalinterdist = new float[NUM_EVENT_TYPES][2][selfpos.length];
		selfdist[POSTSS] = selfpos;
		selfdist[NEGTSS] = read2DDist(directory+"selfneg.txt");
		selfdist[NONTSS] = read2DDist(directory+"selfnon.txt");
		for (int i=0; i<ORIENTATIONS.length; i++) {
			proxinterdist[POSTSS][i] = read2DDist(directory+"proxpos"+i+".txt");
			proxinterdist[NEGTSS][i] = read2DDist(directory+"proxneg"+i+".txt");
			proxinterdist[NONTSS][i] = read2DDist(directory+"proxnon"+i+".txt");
		}
		distalinterdist[POSTSS][MINUS] = read1DDist(directory+"distalposminus.txt");
		distalinterdist[POSTSS][PLUS] = read1DDist(directory+"distalposplus.txt");
		distalinterdist[NEGTSS][MINUS] = read1DDist(directory+"distalnegminus.txt");
		distalinterdist[NEGTSS][PLUS] = read1DDist(directory+"distalnegplus.txt");
		distalinterdist[NONTSS][MINUS] = read1DDist(directory+"distalnonminus.txt");
		distalinterdist[NONTSS][PLUS] = read1DDist(directory+"distalnonplus.txt");
	}
	
	public void initializeFromFile(String file) throws IOException {
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		s = r.readLine();
		split = s.split("\t");
		int diameter = split.length;
		radius = diameter/2;
		selfdist = new float[NUM_EVENT_TYPES][diameter][diameter];
		proxinterdist = new float[NUM_EVENT_TYPES][NUM_ORIENTATIONS][diameter][diameter];
		distalinterdist = new float[NUM_EVENT_TYPES][NUM_ORIENTATIONS][diameter];
		int indicator = SELF_INDICATOR;
		int currentOrientation = 0;
		int currentType = 0;
		int currentRow = 0;
		float[][] currentDist = selfdist[currentType];
		for (int i=0; i<diameter; i++) {
			currentDist[currentRow][i] = Float.valueOf(split[i]);
		}
		while ((s = r.readLine()) != null) {
			currentRow++;
			if (currentRow==diameter) {
				if (!s.equals("")) {
					System.err.println("not at end of matrix");
				}
				currentRow = 0;
				if (indicator==PROX_INDICATOR) {
					currentOrientation++;
					if (currentOrientation==NUM_ORIENTATIONS) {
						currentOrientation = 0;
						currentType++;
						if (currentType==NUM_EVENT_TYPES) {
							currentType = 0;
							indicator = DISTAL_INDICATOR;
						}
					}
				} else if (indicator==DISTAL_INDICATOR) {
					currentOrientation++;
					if (currentOrientation==NUM_ORIENTATIONS) {
						currentOrientation = 0;
						currentType++;
						if (currentType==NUM_EVENT_TYPES) {
							break;
						}
					}
				} else {
					currentType++;
					if (currentType==NUM_EVENT_TYPES) {
						currentType = 0;
						indicator = PROX_INDICATOR;
					}
				}
			}
			split = s.split("\t");
			if (indicator==SELF_INDICATOR) {
				for (int i=0; i<diameter; i++) {
					selfdist[currentType][currentRow][i] = Float.valueOf(split[i]);
				}
			} else if (indicator==PROX_INDICATOR) {
				for (int i=0; i<diameter; i++) {
					proxinterdist[currentType][currentOrientation][currentRow][i] = Float.valueOf(split[i]);
				}
			} else {
				distalinterdist[currentType][currentOrientation][currentRow] = Float.valueOf(split[0]);
			}
		}
		r.close();
	}

	private int strandOrientation(Pair<StrandedPoint,StrandedPoint> pet) {
		return strandValue(pet.car().getStrand())*1 + strandValue(pet.cdr().getStrand())*2;
	}

	private int strandValue(char strand) {
		return strand=='-' ? 0 : 1;
	}

	public static boolean isSelfLigation(Pair<StrandedPoint,StrandedPoint> pet) {
		return pet.car().getChrom().equals(pet.cdr().getChrom()) && pet.cdr().getLocation()-pet.car().getLocation() < SELF_LIGATION_DISTANCE 
				&& pet.car().getStrand()=='-' && pet.cdr().getStrand()=='+';
	}

	private boolean hasOrientation(Pair<StrandedPoint,StrandedPoint> pet, int orientation) {
		if (orientation==MINUSMINUS) {
			return pet.car().getStrand()=='-' && pet.cdr().getStrand()=='-';
		}
		if (orientation==MINUSPLUS) {
			return pet.car().getStrand()=='-' && pet.cdr().getStrand()=='+';
		}
		if (orientation==PLUSMINUS) {
			return pet.car().getStrand()=='+' && pet.cdr().getStrand()=='-';
		}
		if (orientation==PLUSPLUS) {
			return pet.car().getStrand()=='+' && pet.cdr().getStrand()=='+';
		}
		return false;
	}

	private float prob(int relativedist1, int relativedist2, float[][] dist) {
		if (Math.abs(relativedist1)>radius-1 || Math.abs(relativedist2)>radius-1) {
			return 0;
		} else {
			return dist[relativedist1+radius][relativedist2+radius];
		}
	}

	private float prob(int relativedist1, float[] dist) {
		if (Math.abs(relativedist1)>radius-1) {
			return 0;
		} else {
			return dist[relativedist1+radius];
		}
	}
	
	private static class PeakComparator implements Comparator<Pair<StrandedPoint,Integer>> {

		public int compare(Pair<StrandedPoint, Integer> arg0,
				Pair<StrandedPoint, Integer> arg1) {
			return arg0.cdr().compareTo(arg1.cdr());
		}
		
	}

}
