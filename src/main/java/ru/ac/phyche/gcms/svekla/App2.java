package ru.ac.phyche.gcms.svekla;

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

import javax.imageio.ImageIO;

import org.apache.commons.lang3.tuple.Pair;
import org.nd4j.linalg.factory.Nd4j;
import org.openscience.cdk.exception.CDKException;

public class App2 {

	private static int DB_WAX_NUM = TrainPolar.MIN_POLAR_COLUMN_NUMBER;

	private static void predictForSqualaneDB1DB5DB624WAX(String smiles_, Model polar, Model nonpolar, Model db624)
			throws CDKException {
		String smiles = Chemoinformatics.canonical(smiles_, false);
		float riSqualane = nonpolar.predictRI(smiles, Columns.columnNum("Squalane"));
		float riDB1 = nonpolar.predictRI(smiles, Columns.columnNum("DB-1"));
		float riDB5 = nonpolar.predictRI(smiles, Columns.columnNum("DB-5"));
		float riWAX = polar.predictRI(smiles, DB_WAX_NUM);
		float riDB624 = db624.predictRI(smiles, 0);
		System.out.println(smiles + " " + riSqualane + " " + riDB1 + " " + riDB5 + " " + riDB624 + " " + riWAX);
	}

	private static void printUsageInfo() {
		System.out.println("Usage:");
		System.out.println("1) Base prediction for 5 common stationary phases:");
		System.out.println(".... Predict " + " <SMILES string or file with SMILES strings, one per line>");
		System.out.println("When this option is used models and descriptors info"
				+ " should be located in the ./models_polar folder");
		System.out.println();
		System.out.println("2) Base prediction for 5 common stationary phases:");
		System.out.println(".... Predict5 <folder with neural networks: for polar and nonpolar phases, CNN and MLP>"
				+ " <descriptors info file>" + " <second level model file (SVR for DB624)>"
				+ " <a SMILES string or a file with SMILES strings, one per line>");
		System.out.println();
		System.out.println("3) Predict retention index using CNN and MLP:");
		System.out.println(".... PredictNN <MLP neural network text file> <CNN neural network text file>"
				+ " <descriptors info file> <column number (Integer)>"
				+ " <a SMILES string or a file with SMILES strings, one per line>");
		System.out.println();
		System.out.println("4) Predict retention index using CNN, MLP, second level model:");
		System.out
				.println(".... PredictSL <folder with four neural networks: for polar and nonpolar phases, CNN and MLP>"
						+ " <descriptors info file>" + " <second level model file>" + " <second level model type>"
						+ " <a SMILES string or a file with SMILES strings, one per line>");
		System.out.println();
		System.out.println("5) Make cross-validation for the second-level model:");
		System.out.println(".... CV <folder with four neural networks: for polar and nonpolar phases, CNN and MLP>"
				+ " <descriptors info file>" + " <second level model file>" + " <second level model type>"
				+ " <data set file> <log file (cross-validation)>");
		System.out.println("6) Train second-level model for a data set:");
		System.out.println(".... Train <folder with four neural networks: for polar and nonpolar phases, CNN and MLP>"
				+ " <descriptors info file>" + " <second level model file>" + " <second level model type>"
				+ " <data set file>");
		System.out.println("7) Validate second-level model for a data set:");
		System.out
				.println(".... Validate <folder with four neural networks: for polar and nonpolar phases, CNN and MLP>"
						+ " <descriptors info file>" + " <second level model file>" + " <second level model type>"
						+ " <data set file> <predictions file>");
		System.out.println("8) Remove from the data set compounds that are contained in other data set:");
		System.out.println(
				".... RemoveOverlap <data set 1 file> <data set 2 file>" + " <out file (data set 1 after removing)>");
		System.out.println("The folder with neural networks should contain following files:");
		System.out.println("mlp.nn, cnn.nn, mlpPolar.nn, cnnPolar.nn");
		System.out.println("Second level models: 0 - SVR, 1 - SVR linear, 2 - Random forest, 3 - Gradient boosting");
	}

	private static void predict5(String nnFolder, String descriptorsFile, String secondLevelModelFile,
			String smilesOrSmilesFile) throws IOException {
		TrainPolar.MLP mlp = new TrainPolar.MLP();
		TrainPolar.CNN cnn = new TrainPolar.CNN();
		TrainPolar.MLP mlpPolar = new TrainPolar.MLP();
		TrainPolar.CNN cnnPolar = new TrainPolar.CNN();
		mlp.load(nnFolder + "/mlp.nn");
		cnn.load(nnFolder + "/cnn.nn");
		mlpPolar.load(nnFolder + "/mlpPolar.nn");
		cnnPolar.load(nnFolder + "/cnnPolar.nn");
		TrainPolar.SimpleAverageModel averageNonpolar = new TrainPolar.SimpleAverageModel();
		TrainPolar.SimpleAverageModel averagePolar = new TrainPolar.SimpleAverageModel();
		averageNonpolar.addModel(mlp);
		averageNonpolar.addModel(cnn);
		averagePolar.addModel(mlpPolar);
		averagePolar.addModel(cnnPolar);

		SecondLevelModelForPolar model = new SecondLevelModelForPolar.SecondLevelSVRModel();
		Pair<float[], float[]> minmaxArray = Descriptors.readFromFile(descriptorsFile).getMinMaxArray();
		Descriptors d = Descriptors.instance(Descriptors.descriptors2DBut_nAtomLAC_And_MolIP, minmaxArray.getLeft(),
				minmaxArray.getRight(), false);
		mlp.setDescriptorsGenerator(d);
		mlpPolar.setDescriptorsGenerator(d);
		model.setDescriptorsGenerator(d);
		model.setNeuralNetworks(cnn.getNn(), mlp.getNn(), cnnPolar.getNn(), mlpPolar.getNn());
		model.load(secondLevelModelFile);

		try {
			String smiles = Chemoinformatics.canonical(smilesOrSmilesFile, false);
			System.out.println("***************************");
			System.out.println("SMILES Squalane DB-1 DB-5 DB-624 DB-WAX");
			predictForSqualaneDB1DB5DB624WAX(smiles, averagePolar, averageNonpolar, model);
			System.out.println("***************************");
		} catch (CDKException e) {
			BufferedReader reader = new BufferedReader(
					new InputStreamReader(new FileInputStream(new File(smilesOrSmilesFile))));
			System.out.println("***************************");
			System.out.println("Predictions:");
			System.out.println("SMILES Squalane DB-1 DB-5 DB-624 DB-WAX");
			String s = reader.readLine();
			while (s != null) {
				if (!s.trim().equals("")) {
					String smiles = s.trim().split("\\s+")[0];
					try {
						smiles = Chemoinformatics.canonical(smiles, false);
						predictForSqualaneDB1DB5DB624WAX(smiles, averagePolar, averageNonpolar, model);
					} catch (CDKException e1) {
						System.out.println(smiles + " ERROR");
					}
				}
				s = reader.readLine();
			}
			reader.close();
		}
	}

	private static void predictNN(String mlpFile, String cnnFile, String descriptorsFile, String columnType,
			String smilesOrSmilesFile) throws IOException {
		TrainPolar.MLP mlp = new TrainPolar.MLP();
		TrainPolar.CNN cnn = new TrainPolar.CNN();
		TrainPolar.SimpleAverageModel a = new TrainPolar.SimpleAverageModel();
		mlp.load(mlpFile);
		cnn.load(cnnFile);
		Pair<float[], float[]> minmaxArray = Descriptors.readFromFile(descriptorsFile).getMinMaxArray();
		Descriptors d = Descriptors.instance(Descriptors.descriptors2DBut_nAtomLAC_And_MolIP, minmaxArray.getLeft(),
				minmaxArray.getRight(), false);
		mlp.setDescriptorsGenerator(d);
		a.setDescriptorsGenerator(d);
		a.addModel(cnn);
		a.addModel(mlp);
		try {
			String smiles = Chemoinformatics.canonical(smilesOrSmilesFile, false);
			System.out.println("***************************");
			System.out.println("Prediction:");
			System.out.println("Prediction CNN: " + cnn.predictRI(smiles, Integer.parseInt(columnType)));
			System.out.println("Prediction MLP: " + mlp.predictRI(smiles, Integer.parseInt(columnType)));
			System.out.println("Prediction Average: " + a.predictRI(smiles, Integer.parseInt(columnType)));
			System.out.println("***************************");
		} catch (CDKException e) {
			BufferedReader reader = new BufferedReader(
					new InputStreamReader(new FileInputStream(new File(smilesOrSmilesFile))));
			System.out.println("***************************");
			System.out.println("Predictions:");
			System.out.println("SMILES CNN MLP Average");
			String s = reader.readLine();
			while (s != null) {
				if (!s.trim().equals("")) {
					String smiles = s.trim().split("\\s+")[0];
					try {
						smiles = Chemoinformatics.canonical(smiles, false);
						float cnnRI = cnn.predictRI(smiles, Integer.parseInt(columnType));
						float mlpRI = mlp.predictRI(smiles, Integer.parseInt(columnType));
						float averageRI = a.predictRI(smiles, Integer.parseInt(columnType));
						System.out.println(smiles + " " + cnnRI + " " + mlpRI + " " + averageRI);
					} catch (CDKException e1) {
						System.out.println(smiles + " ERROR");
					}
				}
				s = reader.readLine();
			}
			reader.close();
		}
	}

	private static SecondLevelModelForPolar createModelAndLoadNNsAndDescriptors(String nnFolder, String descriptorsFile,
			int secondLevelModelType, boolean precomputeDescriptors) throws IOException {
		TrainPolar.MLP mlp = new TrainPolar.MLP();
		TrainPolar.CNN cnn = new TrainPolar.CNN();
		TrainPolar.MLP mlpPolar = new TrainPolar.MLP();
		TrainPolar.CNN cnnPolar = new TrainPolar.CNN();

		mlp.load(nnFolder + "/mlp.nn");
		cnn.load(nnFolder + "/cnn.nn");
		mlpPolar.load(nnFolder + "/mlpPolar.nn");
		cnnPolar.load(nnFolder + "/cnnPolar.nn");

		SecondLevelModelForPolar model = null;
		if (secondLevelModelType == 0) {
			model = new SecondLevelModelForPolar.SecondLevelSVRModel();
		}
		if (secondLevelModelType == 1) {
			model = new SecondLevelModelForPolar.SecondLevelLinearSVRModel();
		}
		if (secondLevelModelType == 2) {
			model = new SecondLevelModelForPolar.SecondLevelRandomForestModel();
		}
		if (secondLevelModelType == 3) {
			model = new SecondLevelModelForPolar.SecondLevelGBMModel();
		}

		Pair<float[], float[]> minmaxArray = Descriptors.readFromFile(descriptorsFile).getMinMaxArray();
		Descriptors d = Descriptors.instance(Descriptors.descriptors2DBut_nAtomLAC_And_MolIP, minmaxArray.getLeft(),
				minmaxArray.getRight(), precomputeDescriptors);
		mlp.setDescriptorsGenerator(d);
		mlpPolar.setDescriptorsGenerator(d);
		model.setDescriptorsGenerator(d);
		model.setNeuralNetworks(cnn.getNn(), mlp.getNn(), cnnPolar.getNn(), mlpPolar.getNn());
		return model;
	}

	private static void predictSecondLevel(String nnFolder, String descriptorsFile, String secondLevelModelFile,
			int secondLevelModelType, String smilesOrSmilesFile) throws IOException {

		SecondLevelModelForPolar model = createModelAndLoadNNsAndDescriptors(nnFolder, descriptorsFile,
				secondLevelModelType, false);
		model.load(secondLevelModelFile);

		try {
			String smiles = Chemoinformatics.canonical(smilesOrSmilesFile, false);
			System.out.println("***************************");
			System.out.println("Prediction: " + model.predictRI(smiles, 0));
			System.out.println("***************************");
		} catch (CDKException e) {
			BufferedReader reader = new BufferedReader(
					new InputStreamReader(new FileInputStream(new File(smilesOrSmilesFile))));
			System.out.println("***************************");
			System.out.println("Predictions:");
			String s = reader.readLine();
			while (s != null) {
				if (!s.trim().equals("")) {
					String smiles = s.trim().split("\\s+")[0];
					try {
						smiles = Chemoinformatics.canonical(smiles, false);
						System.out.println(smiles + " " + model.predictRI(smiles, 0));
					} catch (CDKException e1) {
						System.out.println(smiles + " ERROR");
					}
				}
				s = reader.readLine();
			}
			reader.close();
		}
	}

	private static void crossValidationSecondLevel(String nnFolder, String descriptorsFile, String secondLevelModelFile,
			int secondLevelModelType, RetentionsDataset dataset, String crossValidationFile)
			throws IOException, CDKException {
		SecondLevelModelForPolar model = createModelAndLoadNNsAndDescriptors(nnFolder, descriptorsFile,
				secondLevelModelType, true);
		RetentionsDataset train = dataset.copy();
		train.makeCanoncalAll(false);
		model.getDescriptorsGenerator().precompute(train.compounds(), false);
		// validation set - one compound. zero is not allowed while really we DO NOT
		// NEED a validation set
		RetentionsDataset val = train.compoundsBasedSplitAndShuffle(1);
		model.init(train, val, model.getDescriptorsGenerator());
		model.crossValidation(model.defaultParameters(), 10, train, crossValidationFile);
		model.save(secondLevelModelFile);
	}

	private static void trainSecondLevel(String nnFolder, String descriptorsFile, String secondLevelModelFile,
			int secondLevelModelType, RetentionsDataset dataset) throws IOException, CDKException {
		SecondLevelModelForPolar model = createModelAndLoadNNsAndDescriptors(nnFolder, descriptorsFile,
				secondLevelModelType, true);
		RetentionsDataset train = dataset.copy();
		train.makeCanoncalAll(false);
		model.getDescriptorsGenerator().precompute(train.compounds(), false);
		// validation set - one compound. zero is not allowed while really we DO NOT
		// NEED a validation set
		RetentionsDataset val = train.compoundsBasedSplitAndShuffle(1);
		model.init(train, val, model.getDescriptorsGenerator());
		model.train(model.defaultParameters());
		model.save(secondLevelModelFile);
	}

	private static void validationSecondLevel(String nnFolder, String descriptorsFile, String secondLevelModelFile,
			int secondLevelModelType, RetentionsDataset dataset, String predictionsFile)
			throws IOException, CDKException {
		SecondLevelModelForPolar model = createModelAndLoadNNsAndDescriptors(nnFolder, descriptorsFile,
				secondLevelModelType, false);
		model.load(secondLevelModelFile);
		RetentionsDataset test = dataset.copy();
		test.makeCanoncalAll(false);
		String result = model.validate(test, predictionsFile);
		System.out.println(result);
	}

	public static void main(String[] args) {
		try {
			Nd4j.getMemoryManager().togglePeriodicGc(true);
			Nd4j.getMemoryManager().setAutoGcWindow(12000);

			if (args.length == 0) {
				printUsageInfo();
				System.exit(1);
			}

			if (args[0].equals("Predict")) {
				if (args.length != 2) {
					printUsageInfo();
					System.exit(1);
				}
				predict5("./models_polar", "./models_polar/descriptors_info.txt", "./models_polar/db624.svr", args[1]);
				System.exit(0);
			}
			if (args[0].equals("Predict5")) {
				if (args.length != 5) {
					printUsageInfo();
					System.exit(1);
				}
				predict5(args[1], args[2], args[3], args[4]);
				System.exit(0);
			}
			if (args[0].equals("PredictNN")) {
				if (args.length != 6) {
					printUsageInfo();
					System.exit(1);
				}
				predictNN(args[1], args[2], args[3], args[4], args[5]);
				System.exit(0);
			}
			if (args[0].equals("PredictSL")) {
				if (args.length != 6) {
					printUsageInfo();
					System.exit(1);
				}
				predictSecondLevel(args[1], args[2], args[3], Integer.parseInt(args[4]), args[5]);
				System.exit(0);
			}
			if (args[0].equals("CV")) {
				if (args.length != 7) {
					printUsageInfo();
					System.exit(1);
				}
				RetentionsDataset data = RetentionsDataset.loadFromFile(args[5]);
				crossValidationSecondLevel(args[1], args[2], args[3], Integer.parseInt(args[4]), data, args[6]);
				System.exit(0);
			}
			if (args[0].equals("Train")) {
				if (args.length != 6) {
					printUsageInfo();
					System.exit(1);
				}
				RetentionsDataset data = RetentionsDataset.loadFromFile(args[5]);
				trainSecondLevel(args[1], args[2], args[3], Integer.parseInt(args[4]), data);
				System.exit(0);
			}
			if (args[0].equals("Validate")) {
				if (args.length != 7) {
					printUsageInfo();
					System.exit(1);
				}
				RetentionsDataset data = RetentionsDataset.loadFromFile(args[5]);
				validationSecondLevel(args[1], args[2], args[3], Integer.parseInt(args[4]), data, args[6]);
				System.exit(0);
			}
			if (args[0].equals("RemoveOverlap")) {
				if (args.length != 4) {
					printUsageInfo();
					System.exit(1);
				}
				RetentionsDataset data1 = RetentionsDataset.loadFromFile(args[1]);
				RetentionsDataset data2 = RetentionsDataset.loadFromFile(args[2]);
				data1.makeCanoncalAll(false);
				data2.makeCanoncalAll(false);
				System.out.println("Overlap " + data1.countIdenticalByCanonicalSmiles(data2));
				data1.filterIdentical(data2);
				data1.filterIdenticalByInchi(data2);
				data1.saveToFile(args[3]);
				System.exit(0);
			}
			if (args[0].equals("SMILESToDepiction")) {
				if (args.length != 3) {
					printUsageInfo();
					System.exit(1);
				}
				String smiles = Chemoinformatics.canonical(args[1], false);
				BufferedImage img = Chemoinformatics.smilesToImageAsBufferedImage(smiles);
				ImageIO.write(img, "png", new File(args[2]));
				System.exit(0);
			}
			printUsageInfo();
			System.exit(1);
		} catch (Exception e) {
			System.out.println("Exception " + e.getMessage());
			e.printStackTrace();
		}
	}
}
