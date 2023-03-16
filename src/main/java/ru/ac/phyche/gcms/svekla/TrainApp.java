package ru.ac.phyche.gcms.svekla;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import org.nd4j.linalg.factory.Nd4j;
import org.openscience.cdk.exception.CDKException;

/**
 * 
 * Command-line application for training of the model.
 *
 */
public class TrainApp {

	private static final float validationFraction = 0.05F;
	private static final float forMetaModelFraction = 0.1F;
	private static final int iterations1DCNN = 200001;
	private static final int iterations2DCNN = 100001;
	private static final int iterationsMLP = 120001;
	private static final int iterationsLinear = 30001;
	private static final int boostingTuneHyperparametersSteps = 50;
	private static final int boostingNEstimators = 800;
	private static final int batchSize = 16;

	private static HashMap<String, Object> xgboostParameters() {
		return XGBoostModel.xgboostDefaultParameters();
	}

	private static void createFolder(String folder) throws IOException {
		File file = new File(folder);
		if (file.exists()) {
			if (!file.isDirectory()) {
				throw new IOException("Can not create folder!");
			}
		} else {
			boolean t0 = file.mkdir();
			if (!t0) {
				throw new IOException("Can not create folder!");
			}
		}
	}

	private static void createFolders(String folder) throws IOException {
		createFolder(folder);
		createFolder(folder + "/logs");
		createFolder(folder + "/models");
		createFolder(folder + "/predictions");
	}

	private static void calculateOverlapping(FileWriter fw, RetentionsDataset first, RetentionsDataset second)
			throws IOException, CDKException {
		fw.write(first.countIdenticalByCanonicalSmiles(second) + " " + first.countIdenticalByInchi(second) + " "
				+ first.countIdenticalByInchikeys(second) + " " + second.countIdenticalByCanonicalSmiles(first) + " "
				+ second.countIdenticalByInchi(first) + " " + second.countIdenticalByInchikeys(first) + "\n");
	}

	private static void validationAndTrainingSetSized(FileWriter fw, Model model) throws IOException {
		fw.write("Training and validation set sizes (data set records): " + model.getTrainSet().size() + " "
				+ model.getValidationSet().size() + "\n");
		fw.write("Training and validation set sizes (compounds): " + model.getTrainSet().compounds().size() + " "
				+ model.getValidationSet().compounds().size() + "\n");
	}

	private static Model trainValidateNeuralNetwork(FileWriter fw, RetentionsDataset trainSet,
			RetentionsDataset validationSet, Descriptors descriptorGenerator, String folder, NeuralNetModel model,
			int iterations, String modelName) throws IOException, CDKException {
		model.init(trainSet, validationSet, descriptorGenerator);
		fw.write("\nTraining " + modelName
				+ ". Checking if no overlappng between training set and validation set. Here should be only zeros... ");
		calculateOverlapping(fw, model.getTrainSet(), model.getValidationSet());
		validationAndTrainingSetSized(fw, model);
		fw.flush();
		int bestIteration = model.trainMultipleIterations(iterations, 500, 2500, batchSize,
				folder + "/models/" + modelName + ".nn", true, folder + "/logs/" + modelName + ".log");
		fw.write(modelName + "model succesfully trained. The best performance on iteration " + bestIteration + "\n");
		fw.write("Accuracy for validation set: "
				+ model.validate(model.getValidationSet(), folder + "/predictions/" + modelName + "validationSet.txt")
				+ "\n");
		return model;
	}

	private static Model trainValidateBoosting(FileWriter fw, RetentionsDataset trainSet,
			RetentionsDataset validationSet, Descriptors descriptorGenerator, String folder, XGBoostModel model,
			int estimators, String modelName, HashMap<String, Object> params) throws IOException, CDKException {
		model.init(trainSet, validationSet, descriptorGenerator);
		fw.write("\nTraining " + modelName
				+ ". Checking if no overlappng between training set and validation set. Here should be only zeros... ");
		calculateOverlapping(fw, model.getTrainSet(), model.getValidationSet());
		validationAndTrainingSetSized(fw, model);
		fw.flush();
		if (params == null) {
			HashMap<String, Object> bestParams = model.hyperParametersTuning(boostingTuneHyperparametersSteps,
					estimators, folder + "/models/" + modelName + ".xgboost", folder + "/logs/" + modelName + ".log");
			fw.write(modelName + "model succesfully trained. The best performance was achieved for parameters "
					+ model.paramsToString(bestParams) + "\n");
		} else {
			model.train(xgboostParameters(), boostingNEstimators);
			model.save(folder + "/models/" + modelName + ".xgboost");
			fw.write(modelName + "model succesfully trained. \n");
		}
		fw.write("Accuracy for validation set: "
				+ model.validate(model.getValidationSet(), folder + "/predictions/" + modelName + "validationSet.txt")
				+ "\n");
		fw.flush();
		return model;
	}

	private static Model[] trainNeuralNetworkAndBoostingModels(FileWriter fw, RetentionsDataset trainSet,
			RetentionsDataset validationSet, Descriptors descriptorGenerator, String folder)
			throws IOException, CDKException {
		Model[] models = new Model[4];
		CNN1DFromSMILESModel modelCNN1D = new CNN1DFromSMILESModel();
		models[0] = trainValidateNeuralNetwork(fw, trainSet, validationSet, null, folder, modelCNN1D, iterations1DCNN,
				"CNN1D");
		CNN2DFromDepictionModel modelCNN2D = new CNN2DFromDepictionModel();
		models[1] = trainValidateNeuralNetwork(fw, trainSet, validationSet, null, folder, modelCNN2D, iterations2DCNN,
				"CNN2D");
		MLPFromDescriptorsAndFingerprints modelMLP = new MLPFromDescriptorsAndFingerprints();
		models[2] = trainValidateNeuralNetwork(fw, trainSet, validationSet, descriptorGenerator, folder, modelMLP,
				iterationsMLP, "MLP2");
		XGBoostModel modelXGBoost = new XGBoostModel();
		models[3] = trainValidateBoosting(fw, trainSet, validationSet, descriptorGenerator, folder, modelXGBoost,
				boostingNEstimators, "XGBoost", xgboostParameters());
		return models;
	}

	private static StackingLinearMetaLearnerModel trainStackingMetamodel(FileWriter fw,
			RetentionsDataset trainSetForMetamodel, RetentionsDataset validationSet, String folder, Model[] models,
			Descriptors descriptorGenerator) throws CDKException, IOException {
		StackingLinearMetaLearnerModel model = new StackingLinearMetaLearnerModel();
		model.init(trainSetForMetamodel, validationSet, descriptorGenerator, models);
		fw.write("\nTraining linear meta-model...\n");
		validationAndTrainingSetSized(fw, model);
		for (int i = 0; i < models.length; i++) {
			fw.write("Model " + i + " . Accuracy for validation set "
					+ models[i].validate(model.getValidationSet(), null) + "\n");
		}
		fw.flush();
		model = model.trainMultipleIterations(iterationsLinear, 500, 30002);
		model.save(folder + "/models/linearMetaModel.nn");
		fw.write("Linear meta model succesfully trained.\n");
		fw.write("Accuracy for validation set "
				+ model.validate(model.getValidationSet(), folder + "/predictions/linearMetaModelValidationSet.txt")
				+ "\n");
		fw.flush();
		return model;
	}

	private static void trainAllModels(FileWriter fw, RetentionsDataset trainSet, String folder)
			throws CDKException, IOException {
		RetentionsDataset trainSetCopy = trainSet.copy();
		RetentionsDataset forMetaModel = trainSetCopy.compoundsBasedSplitAndShuffle(forMetaModelFraction);
		RetentionsDataset validationSet = trainSetCopy.compoundsBasedSplitAndShuffle(validationFraction);
		calculateOverlapping(fw, trainSetCopy, forMetaModel);
		calculateOverlapping(fw, trainSetCopy, validationSet);
		calculateOverlapping(fw, forMetaModel, validationSet);

		fw.write("Dataset sizes for training and validation and for meta model (data set records): "
				+ trainSetCopy.size() + " " + validationSet.size() + " " + forMetaModel.size() + "\n");
		fw.write("Dataset sizes for training and validation and for meta model (compounds): "
				+ trainSetCopy.compounds().size() + " " + validationSet.compounds().size() + " "
				+ forMetaModel.compounds().size() + "\n");
		forMetaModel.saveToFile(folder + "/forMetaModel.ri");
		validationSet.saveToFile(folder + "/validationSet.ri");

		Descriptors d = Descriptors.instance(Descriptors.descriptors2DBut_nAtomLAC_And_MolIP);
		d.precompute(trainSet.compounds(), true);
		Model[] models = trainNeuralNetworkAndBoostingModels(fw, trainSetCopy, validationSet, d, folder);
		trainStackingMetamodel(fw, forMetaModel, validationSet, folder, models, d);
		Descriptors noPrecomputed = Descriptors.instance(Descriptors.descriptors2DBut_nAtomLAC_And_MolIP,
				d.getMinMaxArray().getLeft(), d.getMinMaxArray().getRight(), true);
		noPrecomputed.saveToFile(folder + "/models/descriptors_info.txt");
		d.saveToFile(folder + "/logs/descriptorsTrain.txt");
	}

	private static void testPreTrainedModels(FileWriter fw, RetentionsDataset holdOutTestSet, String folder)
			throws IOException, CDKException {
		String fileCNN1D = folder + "/models/CNN1D.nn";
		String fileCNN2D = folder + "/models/CNN2D.nn";
		String fileMLP = folder + "/models/MLP2.nn";
		String fileXGBoost = folder + "/models/XGBoost.xgboost";
		String fileLinearMeta = folder + "/models/linearMetaModel.nn";
		String fileDescriptors = folder + "/models/descriptors_info.txt";
		StackingLinearMetaLearnerModel twoLevelModelLinear = new StackingLinearMetaLearnerModel();
		twoLevelModelLinear.loadFiveModelsFromFiles(fileCNN1D, fileCNN2D, fileMLP, fileXGBoost, fileLinearMeta,
				fileDescriptors);
		for (int i = 0; i < twoLevelModelLinear.getModels().length; i++) {
			fw.write("Model " + i + " . Accuracy for test set: " + twoLevelModelLinear.getModels()[i]
					.validate(holdOutTestSet, folder + "/predictions/model" + i + "testSet.txt") + "\n");
		}
		fw.write("Linear stacked model. Accuracy for test set: "
				+ twoLevelModelLinear.validate(holdOutTestSet,
						folder + "/predictions/predictionForTestSetLinearStackingModel.txt")
				+ "\n\nLinear model parameters:\n");
		float linearParams[] = twoLevelModelLinear.getParams();
		for (int i = 0; i < linearParams.length; i++) {
			fw.write(linearParams[i] + " ");
		}
	}

	private static void trainAndTest(FileWriter fw, RetentionsDataset trainAndValidationSet,
			RetentionsDataset holdOutTestSet, String folder) throws CDKException, IOException {
		createFolders(folder);
		RetentionsDataset trainAndValidationSetCopy = trainAndValidationSet.copy();
		fw.write("Dataset sizes before filtering (data set records): " + trainAndValidationSetCopy.size() + " "
				+ holdOutTestSet.size() + "\n");
		fw.write("Dataset sizes before filtering (compounds): " + trainAndValidationSetCopy.compounds().size() + " "
				+ holdOutTestSet.compounds().size() + "\n");
		fw.write(
				"Overlapping between test set and training+validation set BEFORE filtering (overlapping compounds will be removed from train set): ");
		calculateOverlapping(fw, trainAndValidationSetCopy, holdOutTestSet);
		trainAndValidationSetCopy.filterIdentical(holdOutTestSet);
		trainAndValidationSetCopy.filterIdenticalByInchi(holdOutTestSet);
		fw.write("Overlapping between test set and training+validation set after filtering (should be zeros here): ");
		calculateOverlapping(fw, trainAndValidationSetCopy, holdOutTestSet);
		fw.write("Dataset sizes after filtering (data set records). Full train-validation set "
				+ trainAndValidationSetCopy.size() + " Hold out test set: " + holdOutTestSet.size() + "\n");
		fw.write("Dataset sizes after filtering (compounds). Full train-validation set "
				+ trainAndValidationSetCopy.compounds().size() + " Hold out test set: "
				+ holdOutTestSet.compounds().size() + "\n");
		fw.flush();
		fw.write("\n\n\nTraining started\n\n");
		trainAllModels(fw, trainAndValidationSetCopy, folder);
		fw.write("\n\n\nTraining finished\n\nTesing\n\n");
		testPreTrainedModels(fw, holdOutTestSet, folder);
	}

	public static void main(String[] args) {
		if (args.length != 4) {
			System.out.println("This script needs 4 space separated command-line parameters:");
			System.out.println("1) Training data set file");
			System.out.println("2) Test data set file");
			System.out.println("3) Output file (results of validation will be logged to this file).");
			System.out.println("4) Output folder.  Full logs and trained models will be saved here.");
			System.out.println("Dataset file format: multiple lines, one retention index record per line.");
			System.out
					.println("<SMILES molecule representation> <Kovac retention index> <Column type (integer value)>");
			System.out.println("Example (butane, DB-5 column):");
			System.out.println("CCCC 400 16");
			System.exit(1);
		}
		try {
			String trainSetFile = args[0];
			String testSetFile = args[1];
			String outFile = args[2];
			String outFolder = args[3];

			Nd4j.getMemoryManager().togglePeriodicGc(true);
			Nd4j.getMemoryManager().setAutoGcWindow(12000);
			RetentionsDataset nist = RetentionsDataset.loadFromFile(trainSetFile);
			RetentionsDataset test = RetentionsDataset.loadFromFile(testSetFile);
			nist.makeCanoncalAll(false);
			test.makeCanoncalAll(false);
			FileWriter fw = new FileWriter(outFile);
			trainAndTest(fw, nist, test, outFolder);
			fw.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(e.getMessage());
			System.out.println("Training falied!!!");
		}

	}

}
