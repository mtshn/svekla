package ru.ac.phyche.gcms.svekla.javafxgui;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Scanner;

import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;

public class MSpectra {
	private static final float THRESHOLD = 2.0f;

	public static class SpectrumMS {
		public float[] intensities = new float[1000];
		public int minMZ = 0;
		public int maxMZ = 999;
	}

	private static int mzToPixels(int width, int minMZ, int maxMZ, int mz) {
		float p = (1.0f * mz - 1.0f * minMZ) / (1.0f * maxMZ - 1.0f * minMZ);
		int pixel = (int) (16 + p * (width - 34));
		return pixel;
	}

	public static void displaySpectrumOnCanvas(Canvas canvas, SpectrumMS spectrum, int minMZ, int maxMZ,
			SpectrumMS spectrumForComparison, boolean squareRoot, boolean markAllPeaks, String caption) {
		GraphicsContext g = canvas.getGraphicsContext2D();
		int width = (int) canvas.getWidth();
		int height = (int) canvas.getHeight();
		g.setFill(Color.valueOf("#FFFFFF"));
		g.fillRect(1, 1, width - 1, height - 1);
		g.setLineWidth(0.85);
		g.setStroke(Color.valueOf("#000000"));
		g.setFont(Font.font(14));
		g.strokeText(caption, 15, 15);
		g.setLineWidth(2);
		g.setStroke(Color.valueOf("#0044ee"));
		g.strokeRect(2, 2, width - 4, height - 4);
		g.setStroke(Color.valueOf("#006622"));
		g.strokeLine(2, height - 30, width - 4, height - 30);
		int newMinMZ = (int) (10 * Math.floor(minMZ / 10.0));
		int newMaxMZ = (int) (10 * Math.floor(maxMZ / 10.0)) + 10;
		int sign = -100;
		for (int mz = newMinMZ; mz <= newMaxMZ; mz = mz + 10) {
			g.setLineWidth(2);
			g.setStroke(Color.valueOf("#110011"));
			int pix = mzToPixels(width, newMinMZ, newMaxMZ, mz);
			if (pix >= sign + 52) {
				g.strokeLine(pix, height - 30, pix, height - 20);
				g.setLineWidth(0.5);
				g.setStroke(Color.valueOf("#000000"));
				g.setFont(Font.font(12));
				g.strokeText("" + mz, pix - 10, height - 8);
				sign = pix;
			} else {
				g.strokeLine(pix, height - 30, pix, height - 26);
			}
		}
		for (int mz = newMinMZ; mz <= newMaxMZ; mz = mz + 1) {
			int pix = mzToPixels(width, newMinMZ, newMaxMZ, mz);
			float intens = spectrum.intensities[mz] / 999.0f;
			if (squareRoot) {
				intens = (float) Math.sqrt(intens);
			}
			int maxpeakLength = height - 60;
			int peakLength = (int) (maxpeakLength * intens);
			g.setLineWidth(2);
			g.setStroke(Color.valueOf("#BB0022"));
			if (spectrumForComparison.intensities[mz] > 0) {
				g.setStroke(Color.valueOf("#337700"));
				g.setLineWidth(1.5);
			}
			g.strokeLine(pix, height - 30, pix, height - 30 - peakLength);
			if (markAllPeaks && (intens > 0)) {
				g.setLineWidth(2.5);
				g.setStroke(Color.valueOf("#00EE11"));
				g.strokeLine(pix, height - 30, pix, height - 34);
			}
		}

	}

	public static SpectrumMS load(String filename) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String s = br.readLine();
			String linesL = "";
			while (s != null) {
				if (!s.trim().equals("")) {
					linesL += s.trim() + "\n";
				}
				s = br.readLine();
			}
			br.close();
			String[] lines = linesL.split("\\r?\\n|\\)\\s+\\(|\\)\\(");
			return load(lines);
		} catch (Exception e) {
			return null;
		}
	}

	public static SpectrumMS load(String[] lines) {
		try {
			SpectrumMS result = new SpectrumMS();

			for (int l = 0; l < lines.length; l++) {
				String[] sp = lines[l].trim().split("\\s+");
				float mz = -100;
				float intens = -100;
				try {
					mz = Float.parseFloat(sp[0].replace(',', '.'));
					intens = Float.parseFloat(sp[1].replace(',', '.'));
				} catch (Throwable e) {
				}
				if ((mz > 0) && (intens > 0) && (mz < result.intensities.length)) {
					result.intensities[Math.round(mz)] += intens;
				}
			}
			float maxint = -1000;
			int minmz = 0;

			for (int i = 0; i < result.intensities.length; i++) {
				if (result.intensities[i] > maxint) {
					maxint = result.intensities[i];
				}
				if (result.intensities[i] != 0) {
					if (minmz == 0) {
						minmz = i;
					}
				}
			}
			for (int i = 0; i < result.intensities.length; i++) {
				result.intensities[i] = result.intensities[i] * 999.0f / maxint;
			}
			for (int i = 0; i < result.intensities.length; i++) {
				result.intensities[i] = result.intensities[i] >= THRESHOLD ? result.intensities[i] : 0.0f;
			}
			int maxmz = 0;

			for (int i = 0; i < result.intensities.length; i++) {
				if (result.intensities[i] != 0) {
					maxmz = i;
				}
			}

			result.minMZ = minmz;
			result.maxMZ = maxmz;
			return result;
		} catch (Exception e) {
			return null;
		}
	}

	public static String spectralSimilarity(SpectrumMS s1, SpectrumMS s2) {
		int min = Math.max(s2.minMZ, s1.minMZ);
		int max = Math.max(s2.maxMZ, s1.maxMZ);
		float a1 = 0;
		float a2 = 0;
		float a12 = 0;
		int n1 = 0;
		int n2 = 0;
		int n12 = 0;
		float i1 = 0;
		float i2 = 0;
		float i12_1 = 0;
		float i12_2 = 0;
		for (int i = min; i < max + 1; i++) {
			float x1 = s1.intensities[i];
			float x2 = s2.intensities[i];
			a1 += Math.sqrt(x1 * x1) * i * i;
			a2 += Math.sqrt(x2 * x2) * i * i;
			a12 += Math.sqrt(x1 * x2) * i * i;
			if (x1 > 0) {
				n1++;
				if (x2 > 0) {
					n12++;
					i12_1 += x1;
					i12_2 += x2;
				}
			}
			if (x2 > 0) {
				n2++;
			}
			i1 += x1;
			i2 += x2;
		}
		float dp = 0;
		if (a12 != 0) {
			dp = a12 * a12 / (a1 * a2);
		}
		float jaccard = 0;
		if (n1 + n2 - n12 != 0) {
			jaccard = n12 * 1.0f / (n1 * 1.0f + n2 - n12);
		}
		float wp = i12_1 / i1;
		float wr = i12_2 / i2;
		dp = Math.round(dp * 1000) / 1000f;
		jaccard = Math.round(jaccard * 1000) / 1000f;
		wr = Math.round(wr * 1000) / 1000f;
		wp = Math.round(wp * 1000) / 1000f;
		String result = "Dot product " + dp + " Jaccard: " + jaccard + " Weighted recall/precision: " + wr + " / " + wp;
		return result;
	}

	public static SpectrumMS sp(String smiles) {
		try {
			String os = System.getProperty("os.name").toLowerCase();
			Scanner sc = null;
			SpectrumMS result = new SpectrumMS();

			if (os.contains("nix") || os.contains("nux")) {
				FileWriter sh = new FileWriter("cfm.sh");
				sh.write("export LD_LIBRARY_PATH=./rdkit;./cfm-predict \"");
				sh.write(smiles + "\"");
				sh.close();
				ProcessBuilder b = new ProcessBuilder("sh", "cfm.sh");
				b.directory(new File("./"));
				Process p = b.start();
				sc = new Scanner(p.getInputStream());
			}
			if (os.contains("win")) {
				ProcessBuilder b = new ProcessBuilder("cfm-predict.exe", smiles);
				b.directory(new File("./"));
				Process p = b.start();
				sc = new Scanner(p.getInputStream());
			}

			boolean t = false;
			while (sc.hasNextLine()) {
				String s = sc.nextLine();
				System.out.println(s);
				String[] sp = s.split("\\s+");
				if (sp.length == 1) {
					if (sp[0].trim().equals("energy0")) {
						t = true;
					}
				}
				if (t) {
					if (sp.length == 2) {
						int mz = Math.round(Float.parseFloat(sp[0]));
						float intens = Float.parseFloat(sp[1]);
						result.intensities[mz] += intens;
					}
				}
			}
			sc.close();
			if (!t) {
				throw (new Exception("CFM ERROR"));
			}
			float maxint = -1000;
			int minmz = 0;

			for (int i = 0; i < result.intensities.length; i++) {
				if (result.intensities[i] > maxint) {
					maxint = result.intensities[i];
				}
				if (result.intensities[i] != 0) {
					if (minmz == 0) {
						minmz = i;
					}
				}
			}
			for (int i = 0; i < result.intensities.length; i++) {
				result.intensities[i] = result.intensities[i] * 999.0f / maxint;
			}
			for (int i = 0; i < result.intensities.length; i++) {
				result.intensities[i] = result.intensities[i] >= THRESHOLD ? result.intensities[i] : 0.0f;
			}
			int maxmz = 0;

			for (int i = 0; i < result.intensities.length; i++) {
				if (result.intensities[i] != 0) {
					maxmz = i;
				}
			}

			result.minMZ = minmz;
			result.maxMZ = maxmz;
			 Files.deleteIfExists(Paths.get("./cfm.sh"));
			return result;
		} catch (Exception e) {
			e.printStackTrace();
			SpectrumMS result = new SpectrumMS();
			result.minMZ = 40;
			result.maxMZ = 400;
			return result;
		}
	}
}
