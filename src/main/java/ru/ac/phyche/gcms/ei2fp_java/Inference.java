package ru.ac.phyche.gcms.ei2fp_java;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;

import ru.ac.phyche.gcms.svekla.Chemoinformatics;

public class Inference {

	public static class ModelFP {

		public static final int SPECTRUM_LENGTH = 800;
		public static final int FP_LENGTH = 1024;

		private float[][] dense1W;
		private float[][] dense2W;
		private float[][] dense3W;
		private float[][] dense4W;
		private float[][] dense5W;
		private float[] dense1Bias;
		private float[] dense2Bias;
		private float[] dense3Bias;
		private float[] dense4Bias;
		private float[] dense5Bias;
		private float[] batchnormMu;
		private float[] batchnormSigma2;
		private float[] batchnormBeta;
		private float[] batchnormGamma;

		public boolean loaded() {
			if ((dense1W != null) && (batchnormGamma != null)) {
				return true;
			}
			return false;
		}

		public float[] predict(float[] spectrum) {
			float[] x = dense(spectrum, dense1W, dense1Bias);
			x = silu(x);
			x = silu(dense(x, dense2W, dense2Bias));
			x = silu(dense(x, dense3W, dense3Bias));
			x = silu(dense(x, dense4W, dense4Bias));
			x = batchnorm(x, batchnormGamma, batchnormBeta, batchnormSigma2, batchnormMu);
			x = sigmoid(dense(x, dense5W, dense5Bias));
			return x;
		}

		public void loadFromFolder(String folderName) throws IOException {
			dense1W = loadWeight(folderName + "/dense1_weight.txt");
			dense2W = loadWeight(folderName + "/dense2_weight.txt");
			dense3W = loadWeight(folderName + "/dense3_weight.txt");
			dense4W = loadWeight(folderName + "/dense4_weight.txt");
			dense5W = loadWeight(folderName + "/dense5_weight.txt");
			dense1Bias = loadBias(folderName + "/dense1_bias.txt");
			dense2Bias = loadBias(folderName + "/dense2_bias.txt");
			dense3Bias = loadBias(folderName + "/dense3_bias.txt");
			dense4Bias = loadBias(folderName + "/dense4_bias.txt");
			dense5Bias = loadBias(folderName + "/dense5_bias.txt");
			batchnormMu = loadBias(folderName + "/batchnorm_mu.txt");
			batchnormSigma2 = loadBias(folderName + "/batchnorm_sigma2.txt");
			batchnormBeta = loadBias(folderName + "/batchnorm_beta.txt");
			batchnormGamma = loadBias(folderName + "/batchnorm_gamma.txt");
		}
	}

	public static float[] relu(float[] x) {
		float[] r = new float[x.length];
		for (int i = 0; i < x.length; i++) {
			if (x[i] > 0) {
				r[i] = x[i];
			} else {
				r[i] = 0f;
			}
		}
		return r;
	}

	public static float[] dense(float[] input, float[][] weight, float[] bias) {
		float[] r = new float[bias.length];
		if (r.length != weight.length) {
			throw new RuntimeException("Dimensions mismatch");
		}
		if (input.length != weight[0].length) {
			throw new RuntimeException("Dimensions mismatch");
		}
		for (int i = 0; i < r.length; i++) {
			float x = bias[i];
			for (int j = 0; j < input.length; j++) {
				x = x + input[j] * weight[i][j];
			}
			r[i] = x;
		}
		return r;
	}

	public static float[] loadBias(String filename) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(filename));
		ArrayList<Float> f = new ArrayList<Float>();
		String s = br.readLine();
		while (s != null) {
			if (!s.trim().equals("")) {
				f.add(Float.parseFloat(s.trim()));
			}
			s = br.readLine();
		}
		br.close();
		float[] r = new float[f.size()];
		for (int i = 0; i < r.length; i++) {
			r[i] = f.get(i);
		}
		return r;
	}

	public static float[][] loadWeight(String filename) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(filename));
		ArrayList<float[]> f = new ArrayList<float[]>();
		String s = br.readLine();
		while (s != null) {
			if (!s.trim().equals("")) {
				String[] sp = s.trim().split("\\s+");
				float[] f1 = new float[sp.length];
				for (int i = 0; i < f1.length; i++) {
					f1[i] = Float.parseFloat(sp[i]);
				}
				f.add(f1);
			}
			s = br.readLine();
		}
		br.close();
		float[][] r = new float[f.size()][];
		for (int i = 0; i < r.length; i++) {
			r[i] = f.get(i);
		}
		return r;
	}

	public static float sigmoid(float x) {
		double y = Math.exp(0 - x);
		y = 1 / (1 + y);
		return (float) y;
	}

	public static float[] sigmoid(float[] x) {
		float[] result = new float[x.length];
		for (int i = 0; i < x.length; i++) {
			result[i] = sigmoid(x[i]);
		}
		return result;
	}

	public static float silu(float x) {
		return x * sigmoid(x);
	}

	public static float[] silu(float[] x) {
		float[] result = new float[x.length];
		for (int i = 0; i < x.length; i++) {
			result[i] = silu(x[i]);
		}
		return result;
	}

	public static float[] batchnorm(float[] input, float[] gamma, float[] beta, float[] sigma2, float[] mu) {
		if ((input.length != gamma.length) || (input.length != beta.length) || (input.length != sigma2.length)
				|| (input.length != mu.length)) {
			throw new RuntimeException("Dimensions mismatch");
		}
		float[] result = new float[input.length];
		for (int i = 0; i < input.length; i++) {
			float x = input[i];
			float a = (float) (Math.sqrt(sigma2[i] + (1E-5)));
			x = (x - mu[i]) / a;
			x = x * gamma[i] + beta[i];
			result[i] = x;
		}
		return result;
	}

	public static float crossentropy(float[] pred, float ref[]) {
		if (pred.length != ref.length) {
			throw new RuntimeException("Dimensions mismatch");
		}
		double x = 0;
		for (int i = 0; i < ref.length; i++) {
			x = x + ref[i] * Math.log(pred[i] + 1E-7);
			x = x + (1 - ref[i]) * Math.log(1 - pred[i] + 1E-7);
		}
		x = 0 - x / pred.length;
		return (float) x;
	}

	public static float[] crossentropy(float[][] pred, float ref[][]) {
		if (pred.length != ref.length) {
			throw new RuntimeException("Dimensions mismatch");
		}
		float[] result = new float[pred.length];
		for (int i = 0; i < pred.length; i++) {
			result[i] = crossentropy(pred[i], ref[i]);
		}
		return result;
	}

	public static float[] crossentropy(float[] predictiedFingerprints, String[] smiles) {
		int n = smiles.length;
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

		float[] fp0 = predictiedFingerprints;
		float[] result = new float[smiles.length];
		AtomicInteger ai = new AtomicInteger();
		Arrays.stream(a).parallel().forEach(i -> {
			try {
				String s = Chemoinformatics.canonical(smiles[i].trim(), false);
				float[] fp1 = Chemoinformatics.fingerprints(s, Chemoinformatics.FingerprintsType.CIRCULAR_6_1024);
				result[i] = crossentropy(fp0, fp1);
				if(ai.incrementAndGet()%50000==0) {
					System.out.println(ai.get()+" of "+result.length);
				}
			} catch (Exception e) {
				e.printStackTrace();
				result[i] = Float.POSITIVE_INFINITY;
			}

		});
		return result;
	}

	public static void main(String[] args) throws IOException {
		ModelFP m = new ModelFP();
		m.loadFromFolder("./fingerprints_model");
		BufferedReader br = new BufferedReader(new FileReader("./train_ei2fp/spectra.txt"));
		String s = null;
		for (int i = 0; i < 5; i++) {
			s = br.readLine();
		}
		float[] sp = new float[800];
		String[] splt = s.split("\\s+");
		for (int i = 0; i < splt.length; i = i + 2) {
			sp[Integer.parseInt(splt[i])] = Float.parseFloat(splt[i + 1]) / 1000;
		}
		float[] z = m.predict(sp);
		System.out.println(z.length);
		for (int i = 0; i < z.length; i++) {
			if (z[i] > 0.5) {
				System.out.print(i + " ");
			}
		}
		System.out.println();
		for (int i = 0; i < z.length; i++) {
			System.out.print(z[i] + " ");
		}
		br.close();
	}

}
