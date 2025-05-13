package ru.ac.phyche.gcms.svekla.javafxgui;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import ru.ac.phyche.gcms.svekla.Chemoinformatics;

public class SurgeIsomers {

	public static String concatNewLine(String[] s) {
		StringBuilder result = new StringBuilder();
		for (int i = 0; i < s.length; i++) {
			result = result.append(s[i] + "\n");
		}
		return result.toString();
	}

	private static String[] surgeParams(String formula, boolean noTriple, boolean noSmallRings, boolean noXXX,
			boolean custom, String customLine, boolean countOnly) {
		String execName = null;
		String os = System.getProperty("os.name").toLowerCase();
		if (os.contains("nix") || os.contains("nux")) {
			execName = "./surge-linux-v1.0";
		}
		if (os.contains("win")) {
			execName = "./surge-windows-v1.0.exe";
		}
		if (custom) {
			String x = execName;
			if (countOnly) {
				x = x + " -u";
			} else {
				x = x + " -S";
			}
			x = x + " " + customLine.trim() + " " + formula.trim();
			return x.split("\\s+");
		}
		String result = execName;
		if (countOnly) {
			result = result + " -u";
		} else {
			result = result + " -S";
		}
		if (noTriple) {
			result = result + " -T";
		}
		if (noSmallRings) {
			result = result + " -t0 -f0";
		}

		result = result + " -P -B1,2,3,4,";
		if (noXXX) {
			result = result + "5,";
		}
		result = result + "6,7,8,9";
		result = result + " " + formula.trim();
		return result.split("\\s+");
	}

	public static int surgeCount(String formula, boolean noTriple, boolean noSmallRings, boolean noXXX, boolean custom,
			String customLine) throws IOException {
		String[] surgeparams = surgeParams(formula, noTriple, noSmallRings, noXXX, custom, customLine, true);
		ProcessBuilder b = new ProcessBuilder(surgeparams);
		b.directory(new File("./"));
		Process p = b.start();
		Scanner sc = new Scanner(p.getErrorStream());
		int result = -1;
		while (sc.hasNextLine()) {
			String s = sc.nextLine();
			System.out.println(s);
			if (s.contains(">Z generated")) {
				String[] x = s.split("\\s+");
				if (x[5].equals("->")) {
					if (x[7].equals("in")) {
						result = Integer.parseInt(x[6]);
					}
				}
			}
		}
		sc.close();
		return result;
	}

	public static String[] surgeGenerate(String formula, boolean noTriple, boolean noSmallRings, boolean noXXX,
			boolean custom, String customLine) throws IOException {
		String[] surgeparams = surgeParams(formula, noTriple, noSmallRings, noXXX, custom, customLine, false);
		ProcessBuilder b = new ProcessBuilder(surgeparams);
		b.directory(new File("./"));
		Process p = b.start();
		Scanner sc = new Scanner(p.getInputStream());
		ArrayList<String> result = new ArrayList<String>();
		int i=0;
		while (sc.hasNextLine()) {
			String s = sc.nextLine();
			//System.out.println(i+" "+s);
			if (s.trim().split("\\s+").length == 1) {
				result.add(s);
			}
			i=i+1;
		}
		sc.close();
		return result.toArray(new String[result.size()]);
	}

}
