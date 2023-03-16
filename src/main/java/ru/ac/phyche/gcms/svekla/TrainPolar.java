package ru.ac.phyche.gcms.svekla;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;

import org.apache.commons.lang3.tuple.Pair;
import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.graph.MergeVertex;
import org.deeplearning4j.nn.conf.layers.Convolution1DLayer;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.GlobalPoolingLayer;
import org.deeplearning4j.nn.conf.layers.LossLayer;
import org.deeplearning4j.nn.conf.layers.PoolingType;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.activations.Activation;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.learning.config.Adam;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.openscience.cdk.exception.CDKException;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;

/**
 * 
 * Utils for training neural networks (MLP and CNN) for Standard polar columns
 * (e.g. DB-WAX). Trainsfer learning is used: at first stage models are trained
 * for non-polar columns and than trained for polar columns.
 */
public class TrainPolar {
	public static final int MIN_POLAR_COLUMN_NUMBER = 15; // This number should be equal Columns.nonPolarMain.length+1
	public static final int MAX_COLUMN_NUMBER = 35;
	private static final String[] polarMainColums = new String[] { "DB-Wax", "Carbowax_20M", "Supelcowax-10", "OV-351",
			"HP-Innowax", "PEG-20M", "BP-20", "FFAP", "CP-Wax_52CB", "HP-Innowax_FSC", "Innowax", "Innowax_FSC",
			"RTX-Wax", "PEG_4000", "DB-FFAP", "Carbowax", "ZB-Wax", "HP-Wax", "AT-Wax", "Stabilwax" };

	private static int polarColumnNum(String column) {
		for (int i = 0; i < polarMainColums.length; i++) {
			if (polarMainColums[i].equals(column)) {
				return (i + MIN_POLAR_COLUMN_NUMBER);
			}
		}
		return polarMainColums.length + MIN_POLAR_COLUMN_NUMBER;
	}

	@SuppressWarnings("unused")
	private static void convertNISTTmpFileToDataSetPolarColumns(String inp, String outp) throws IOException {
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outp))));
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inp))));
		String s = reader.readLine();

		while (s != null) {
			ArrayList<String> block = new ArrayList<String>();
			while (!s.equals("***END***")) {
				block.add(s);
				s = reader.readLine();
			}
			String[] blockArr = new String[block.size()];
			for (int k = 0; k < block.size(); k++) {
				blockArr[k] = block.get(k);
			}
			for (int j = 4; j < blockArr.length; j++) {
				String[] split = blockArr[j].split("\\s+");
				String columnName = split[6].trim();
				float rI = Float.parseFloat(split[3]);

				if (blockArr[j].contains("Standard_polar")) {
					if (rI != 0) {
						writer.write(blockArr[3].trim() + " " + rI + " " + polarColumnNum(columnName.trim()) + "\n");
					}
				}
			}
			s = reader.readLine();
		}

		reader.close();
		writer.close();
	}

	// An alias to MLPFromDescriptorsAndFingerprints
	/**
	 * Almost the same as MLPFromDescriptorsAndFingerprints Parallel computation of
	 * fingerprints is employed while training
	 */
	public static class MLP extends MLPFromDescriptorsAndFingerprints {
		@Override
		public INDArray[] nextBatchInput(RetentionsDataset dataSet, int[] indices) throws CDKException {
			float[][] featuresDescriptors = new float[indices.length][];
			final float[][] fingerprints = new float[indices.length][];
			int[] ints = new int[indices.length];
			for (int i = 0; i < indices.length; i++) {
				ints[i] = i;
			}
			Arrays.stream(ints).parallel().forEach(i -> {
				try {
					fingerprints[i] = dataSet.fingerprints(
							Chemoinformatics.FingerprintsType.ADDITIVE_CIRCULAR_4_1024_NO_SCALE, indices[i]);
					featuresDescriptors[i] = descriptorsFeatures(dataSet, indices[i]);
				} catch (CDKException e) {
					throw (new RuntimeException(e.getMessage()));
				}
			});
			;
			INDArray descriptorsINDArray = Nd4j.create(featuresDescriptors);
			INDArray fingerprintsINDArray = Nd4j.create(fingerprints);
			return new INDArray[] { descriptorsINDArray, fingerprintsINDArray };
		}
	}

	// An alias to CNN1DFromSMILESModel
	/**
	 * Almost the same as CNN1DFromSMILESModel
	 * cudnnAlgoMode(Convolution1DLayer.AlgoMode.NO_WORKSPACE) is used to avoid some
	 * internal bugs in Deepleaning4j
	 */
	public static class CNN extends CNN1DFromSMILESModel {
		@Override
		public void initNN() {
			ComputationGraphConfiguration conf = new NeuralNetConfiguration.Builder().weightInit(WeightInit.RELU)
					.updater(new Adam(0.0003)).cudnnAlgoMode(Convolution1DLayer.AlgoMode.NO_WORKSPACE).graphBuilder()
					.addLayer("CNN0",
							new Convolution1DLayer.Builder(6, 1).nIn(36).nOut(300).activation(Activation.RELU).build(),
							"INPUT")
					.addLayer("CNN1",
							new Convolution1DLayer.Builder(6, 1).nIn(300).nOut(300).activation(Activation.RELU).build(),
							"CNN0")
					.addLayer("POOLING", new GlobalPoolingLayer.Builder(PoolingType.AVG).build(), "CNN1")
					.addVertex("MERGE", new MergeVertex(), "POOLING", "COLUMN_INFO")
					.addLayer("DENSE0", new DenseLayer.Builder().nIn(338).nOut(600).activation(Activation.RELU).build(),
							"MERGE")
					.addInputs("INPUT", "COLUMN_INFO")
					.addLayer("DENSE1",
							new DenseLayer.Builder().nIn(600).nOut(1).activation(Activation.IDENTITY).build(), "DENSE0")
					.addLayer("OUT", new LossLayer.Builder()
							.lossFunction(LossFunctions.LossFunction.MEAN_ABSOLUTE_ERROR).build(), "DENSE1")
					.setOutputs("OUT").build();
			ComputationGraph nn = new ComputationGraph(conf);
			nn.init();
			this.setNn(nn);
		}
	}

	/**
	 * Split for cross-validation
	 * 
	 * @param n    n-fold
	 * @param data dataset (will not modified)
	 * @return n data sets
	 * @throws CDKException cdk
	 */
	public static RetentionsDataset[] cvSplit(int n, RetentionsDataset data) throws CDKException {
		RetentionsDataset cv0 = data.copy();
		RetentionsDataset testSets[] = new RetentionsDataset[n];
		int compoundsInSubset = (Math.round(((float) cv0.compounds().size()) / n));
		for (int i = 0; i < n - 1; i++) {
			testSets[i] = cv0.compoundsBasedSplitAndShuffle(compoundsInSubset).copy();
		}
		testSets[n - 1] = cv0.copy();
		return testSets;
	}

	/**
	 * Simple model that averages reuslts of other models
	 *
	 */
	public static class SimpleAverageModel extends Model {

		private ArrayList<Model> models = new ArrayList<Model>();

		/**
		 * Add model to average
		 * 
		 * @param m
		 */
		public void addModel(Model m) {
			models.add(m);
		}

		@Override
		public float predictRI(RetentionsDataset dataSet, int entry) throws CDKException {
			float sum = 0;
			for (Model m : models) {
				sum += m.predictRI(dataSet, entry);
			}
			return sum / models.size();
		}

		@Override
		public void save(String filename) throws IOException {
			throw new UnsupportedOperationException("This model can not be saved");
		}

		@Override
		public void load(String filename) throws IOException {
			throw new UnsupportedOperationException("This model can not be loaded");
		}

	}

	private static void trainMLPandCNNForPolarAndNonPolar(int batchSize, int maxIterationsNonPolar,
			int maxIterationsPolar, String outMLPFileNonpolar, String outCNNFileNonpolar, String outMLPFilePolar,
			String outCNNFilePolar, String descriptorsInfoFile, RetentionsDataset nonpolar, RetentionsDataset polar,
			RetentionsDataset[] excludeFromTrain, String outputFolder, FileWriter fw) throws CDKException, IOException {
		Pair<float[], float[]> minmaxArray = Descriptors.readFromFile(descriptorsInfoFile).getMinMaxArray();
		Descriptors d = Descriptors.instance(Descriptors.descriptors2DBut_nAtomLAC_And_MolIP, minmaxArray.getLeft(),
				minmaxArray.getRight(), true);
		RetentionsDataset trainNonpolar = nonpolar.copy();
		RetentionsDataset trainPolar = polar.copy();
		trainNonpolar.makeCanoncalAll(false);
		trainPolar.makeCanoncalAll(false);

		fw.write((new Date()).toString() + "\n");
		fw.write("Polar data set (before exclusions): " + trainPolar.size() + "  Individual compounds: "
				+ trainPolar.compounds().size() + "\n");
		fw.write("Non-polar data set (before exclusions): " + trainNonpolar.size() + "  Individual compounds: "
				+ trainNonpolar.compounds().size() + "\n");
		fw.flush();
		for (int i = 0; i < excludeFromTrain.length; i++) {
			excludeFromTrain[i].makeCanoncalAll(false);
			trainNonpolar.filterIdentical(excludeFromTrain[i]);
			trainNonpolar.filterIdenticalByInchi(excludeFromTrain[i]);
			trainPolar.filterIdentical(excludeFromTrain[i]);
			trainPolar.filterIdenticalByInchi(excludeFromTrain[i]);
		}
		fw.write("Polar data set (after exclusions): " + trainPolar.size() + "  Individual compounds: "
				+ trainPolar.compounds().size() + "\n");
		fw.write("Non-polar data set (after exclusions): " + trainNonpolar.size() + "  Individual compounds: "
				+ trainNonpolar.compounds().size() + "\n");
		fw.write("\n\n");
		fw.flush();
		fw.write((new Date()).toString() + "\n");
		d.precompute(trainNonpolar.compounds(), false);
		d.precompute(trainPolar.compounds(), false);
		RetentionsDataset validationNonpolar = trainNonpolar.compoundsBasedSplitAndShuffle(0.1F);
		RetentionsDataset validationPolar = trainPolar.compoundsBasedSplitAndShuffle(0.1F);
		validationPolar.saveToFile(outputFolder + "/validationPolar.ri");
		validationNonpolar.saveToFile(outputFolder + "/validationNonpolar.ri");

		MLP mlp = new MLP();
		mlp.init(trainNonpolar, validationNonpolar, d);
		fw.write((new Date()).toString() + "\n");
		fw.write("Training of MLP for non-polar columns...\n");
		int bestiter = mlp.trainMultipleIterations(maxIterationsNonPolar, 500, 2500, batchSize, outMLPFileNonpolar,
				false, outputFolder + "/mlpnonpolar.log");
		fw.write("Best iteration at " + bestiter + "\n");
		fw.write("Accuracy for validation set: "
				+ mlp.validate(validationNonpolar, outputFolder + "/mlpnonpolar_val.txt") + "\n\n");
		mlp.setTrainSet(trainPolar);
		mlp.setValidationSet(validationPolar);
		fw.flush();
		fw.write((new Date()).toString() + "\n");
		fw.write("Training of MLP for polar columns...\n");
		bestiter = mlp.trainMultipleIterations(maxIterationsPolar, 100, 500, batchSize, outMLPFilePolar, false,
				outputFolder + "/mlppolar.log");
		fw.write("Best iteration at " + bestiter + "\n");
		fw.write("Accuracy for validation set: " + mlp.validate(validationPolar, outputFolder + "/mlppolar_val.txt")
				+ "\n\n");
		mlp = null;

		CNN cnn = new CNN();
		cnn.init(trainNonpolar, validationNonpolar, d);
		fw.flush();
		fw.write((new Date()).toString() + "\n");
		fw.write("Training of CNN for non-polar columns...\n");
		bestiter = cnn.trainMultipleIterations(maxIterationsNonPolar, 500, 2500, batchSize, outCNNFileNonpolar, false,
				outputFolder + "/cnnnonpolar.log");
		fw.write("Best iteration at " + bestiter + "\n");
		fw.write("Accuracy for validation set: "
				+ cnn.validate(validationNonpolar, outputFolder + "/cnnnonpolar_val.txt") + "\n\n");
		cnn.setTrainSet(trainPolar);
		cnn.setValidationSet(validationPolar);
		fw.flush();
		fw.write((new Date()).toString() + "\n");
		fw.write("Training of CNN for polar columns...\n");
		bestiter = cnn.trainMultipleIterations(maxIterationsPolar, 100, 500, batchSize, outCNNFilePolar, false,
				outputFolder + "/cnnpolar.log");
		fw.write("Best iteration at " + bestiter + "\n");
		fw.write("Accuracy for validation set: " + cnn.validate(validationPolar, outputFolder + "/cnnpolar_val.txt")
				+ "\n\n");
		cnn = null;
		fw.write((new Date()).toString() + "\n");
		fw.flush();
	}

	private static String validateMLPandCNN(String mlpFile, String cnnFile, String descriptorsFile,
			RetentionsDataset testSet, String mlpPredictions, String cnnPredictions, String averagePredictions)
			throws IOException, CDKException {
		Pair<float[], float[]> minmaxArray = Descriptors.readFromFile(descriptorsFile).getMinMaxArray();
		Descriptors d = Descriptors.instance(Descriptors.descriptors2DBut_nAtomLAC_And_MolIP, minmaxArray.getLeft(),
				minmaxArray.getRight(), false);

		MLP mlp = new MLP();
		mlp.load(mlpFile);
		mlp.setDescriptorsGenerator(d);
		String result = "MLP: " + mlp.validate(testSet, mlpPredictions);
		CNN cnn = new CNN();
		cnn.load(cnnFile);
		result += "\nCNN: " + cnn.validate(testSet, cnnPredictions);
		SimpleAverageModel a = new SimpleAverageModel();
		a.models.add(cnn);
		a.models.add(mlp);
		a.setDescriptorsGenerator(mlp.getDescriptorsGenerator());
		result += "\nAverage: " + a.validate(testSet, averagePredictions);
		return result;
	}

	private static void trainMLPandCNNForPolarAndNonPolar(int batchSize, int maxIterationsNonPolar,
			int maxIterationsPolar, String outMLPFileNonpolar, String outCNNFileNonpolar, String outMLPFilePolar,
			String outCNNFilePolar, String descriptorsInfoFile, RetentionsDataset nonpolar, RetentionsDataset polar,
			String excludeFromTrainFiles, String outputFolder, FileWriter fw) throws CDKException, IOException {
		BufferedReader reader = new BufferedReader(
				new InputStreamReader(new FileInputStream(new File(excludeFromTrainFiles))));
		String s = reader.readLine();
		ArrayList<String> list = new ArrayList<String>();
		while (s != null) {
			if (!s.trim().equals("")) {
				list.add(s);
			}
			s = reader.readLine();
		}
		reader.close();
		RetentionsDataset[] excludeFromTrain = new RetentionsDataset[list.size()];
		RetentionsDataset nonpolarCopy = nonpolar.copy();
		RetentionsDataset polarCopy = polar.copy();
		nonpolarCopy.makeCanoncalAll(false);
		polarCopy.makeCanoncalAll(false);
		fw.write((new Date()).toString() + "\n");
		for (int i = 0; i < list.size(); i++) {
			System.out.println(i);
			excludeFromTrain[i] = RetentionsDataset.loadFromFile(list.get(i));
			excludeFromTrain[i].makeCanoncalAll(false);
			fw.write(list.get(i) + " " + i + " size (records, compounds):" + excludeFromTrain[i].size() + " "
					+ excludeFromTrain[i].compounds().size() + "\n");
			fw.write("overlap " + excludeFromTrain[i].countIdenticalByInchi(nonpolarCopy) + " "
					+ excludeFromTrain[i].countIdenticalByInchi(polarCopy) + "\n");
			fw.flush();
		}
		fw.write("Overlaps:\n");
		for (int i = 0; i < list.size(); i++) {
			for (int j = 0; j < list.size(); j++) {
				fw.write(excludeFromTrain[i].countIdenticalByCanonicalSmiles(excludeFromTrain[j]) + " ");
			}
			fw.write("\n");
			fw.flush();
		}
		fw.flush();
		fw.write("\n\n");
		trainMLPandCNNForPolarAndNonPolar(batchSize, maxIterationsNonPolar, maxIterationsPolar, outMLPFileNonpolar,
				outCNNFileNonpolar, outMLPFilePolar, outCNNFilePolar, descriptorsInfoFile, nonpolarCopy, polarCopy,
				excludeFromTrain, outputFolder, fw);
	}

	private static void printUsageInfo() {
		System.out.println("Usage:");
		System.out.println("For training of MLP and CNN neural networks for both non-polar and polar columns:");
		System.out.println(" .... Train <batch size> <max training iterations nonpolar> <max training iterations polar>"
				+ " <folder, where trained neural networks should be saved> <descriptors info file>"
				+ " <non polar dataset file> <polar data set file>"
				+ " <file that contains list of data set files (one file name with path per line), that should be excluded from train set>"
				+ " <output folder>" + "<log file name>");
		System.out.println();
		System.out.println("For validation of MLP and CNN neural networks for a data set:");
		System.out.println(".... Validate <MLP neural network> <CNN neural network> <descriptors info file>"
				+ "	<data set> <MLP predictions file (output)>"
				+ " <CNN predictions file (output)> <average predictions file (output)> ");
		System.out.println();
		System.out.println("Convert neural networks from txt to binary file format:");
		System.out.println(".... txt2nn <MLP neural network text file> <CNN neural network text file>"
				+ " <MLP neural network binary file> <CNN neural network binary file>");
		System.out.println("Convert neural networks from txt to binary file format:");
		System.out.println(".... nn2txt <MLP neural network binary file> <CNN neural network binary file>"
				+ " <MLP neural network text file> <CNN neural network text file>");
		System.out.println();
	}

	public static void main(String[] args) {
		try {
			Nd4j.getMemoryManager().togglePeriodicGc(true);
			Nd4j.getMemoryManager().setAutoGcWindow(12000);
			if (args.length == 0) {
				printUsageInfo();
				System.exit(1);
			}

			if (args[0].equals("Train")) {
				if (args.length != 11) {
					printUsageInfo();
					System.exit(1);
				}
				FileWriter fw = new FileWriter(args[10]);
				RetentionsDataset nonpolar = RetentionsDataset.loadFromFile(args[6]);
				RetentionsDataset polar = RetentionsDataset.loadFromFile(args[7]);

				trainMLPandCNNForPolarAndNonPolar(Integer.parseInt(args[1]), Integer.parseInt(args[2]),
						Integer.parseInt(args[3]), args[4] + "/mlp.nn", args[4] + "/cnn.nn", args[4] + "/mlpPolar.nn",
						args[4] + "/cnnPolar.nn", args[5], nonpolar, polar, args[8], args[9], fw);
				System.exit(0);
			}

			if (args[0].equals("Validate")) {
				if (args.length != 8) {
					printUsageInfo();
					System.exit(1);
				}
				RetentionsDataset data = RetentionsDataset.loadFromFile(args[4]);
				System.out.println(validateMLPandCNN(args[1], args[2], args[3], data, args[5], args[6], args[7]));
				System.exit(0);
			}

			if (args[0].equals("nn2txt")) {
				if (args.length != 5) {
					printUsageInfo();
					System.exit(1);
				}
				ComputationGraph mlp = ComputationGraph.load(new File(args[1]), false);
				ComputationGraph cnn = ComputationGraph.load(new File(args[2]), false);
				Nd4j.writeTxt(mlp.params(), args[3]);
				Nd4j.writeTxt(cnn.params(), args[4]);
				System.exit(0);
			}

			if (args[0].equals("txt2nn")) {
				if (args.length != 5) {
					printUsageInfo();
					System.exit(1);
				}
				MLP mlpModel = new MLP();
				mlpModel.initNN();
				INDArray params = Nd4j.readTxt(args[1]);
				mlpModel.getNn().setParams(params);
				mlpModel.getNn().save(new File(args[3]), false);

				CNN cnnModel = new CNN();
				cnnModel.initNN();
				params = Nd4j.readTxt(args[2]);
				cnnModel.getNn().setParams(params);
				cnnModel.getNn().save(new File(args[4]), false);
				System.exit(0);
			}

			printUsageInfo();
			System.exit(1);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
