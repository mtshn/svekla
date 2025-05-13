package ru.ac.phyche.gcms.ei2fp_java;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import org.openscience.cdk.exception.CDKException;

import ru.ac.phyche.gcms.svekla.Chemoinformatics;

public class CreateTrainSet {

	public static float[] fp(String smiles) throws CDKException {
		String s = Chemoinformatics.canonical(smiles.trim(), false);
		return Chemoinformatics.fingerprints(s, Chemoinformatics.FingerprintsType.CIRCULAR_6_1024);
	}

	public static String fpToSparseString(float[] fp) throws CDKException {
		String result = "";
		for (int i = 0; i < fp.length; i++) {
			if (fp[i] != 0) {
				result = result + i + "  " + fp[i] + " ";
			}
		}
		if(result.trim().equals("")) {
			throw new CDKException("Zero fingerprint");
		}
		return result.trim();
	}

	public static String[] convertString(String s) throws CDKException {
		String[] splt = s.trim().split("\\|");
		String fpS = fpToSparseString(fp(splt[1].trim()));
		return new String[] { fpS, splt[2].trim() };
	}

	private static int[] intsrnd(int n) {
		int[] a = new int[n];
		for (int i = 0; i < n; i++) {
			a[i] = i;
		}
		Random rnd = new Random();
		for (int i = 0; i < n; i++) {
			int x = a[i];
			int j = rnd.nextInt(n);
			int y = a[j];
			a[i] = y;
			a[j] = x;
		}
		return a;
	}

	public static String[][] convertManyStrings(String[] s) {
		String[][] result = new String[s.length][];
		Arrays.stream(intsrnd(s.length)).parallel().forEach(i -> {
			try {
				result[i] = convertString(s[i]);
			} catch (Exception e) {
				e.printStackTrace();
				result[i] = new String[] { s[i] + " " + e.getMessage(), null };
			}
		});
		return result;
	}

	public static boolean convertManyStrings(int n, FileWriter fp, FileWriter spectrum, FileWriter bad,
			BufferedReader input) throws IOException {
		ArrayList<String> al = new ArrayList<String>();
		String s = input.readLine();
		int i = 0;
		while ((s != null) && (i < n - 1)) {
			al.add(s);
			s = input.readLine();
			i++;
		}
		if (s != null) {
			al.add(s);
		}
		System.out.println(al.size());
		String[][] result = convertManyStrings(al.toArray(new String[al.size()]));
		for (int j = 0; j < result.length; j++) {
			if (result[j][1] == null) {
				bad.write(result[j][0] + "\n");
			} else {
				fp.write(result[j][0] + "\n");
				spectrum.write(result[j][1] + "\n");
			}
		}
		return s != null;
	}

	public static void main(String args[]) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader("./train_ei2fp/mainlib.ms"));
		FileWriter fw1 = new FileWriter("./train_ei2fp/fp.txt");
		FileWriter fw2 = new FileWriter("./train_ei2fp/spectra.txt");
		FileWriter fw3 = new FileWriter("./train_ei2fp/bad.txt");
		int n = 10000;
		boolean t = convertManyStrings(n, fw1, fw2, fw3, br);
		while (t) {
			t = convertManyStrings(n, fw1, fw2, fw3, br);
		}
		fw1.close();
		fw2.close();
		fw3.close();
		br.close();
	}
}
