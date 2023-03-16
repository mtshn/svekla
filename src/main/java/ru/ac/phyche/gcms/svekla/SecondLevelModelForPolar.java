package ru.ac.phyche.gcms.svekla;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.openscience.cdk.exception.CDKException;

import com.thoughtworks.xstream.XStream;

import smile.base.cart.Loss;
import smile.data.DataFrame;
import smile.data.formula.Formula;
import smile.math.kernel.GaussianKernel;
import smile.math.kernel.LinearKernel;
import smile.regression.GradientTreeBoost;
import smile.regression.KernelMachine;
import smile.regression.RandomForest;
import smile.regression.SVR;

/**
 * It is machine learning model that predicts gas chromatographic retention for
 * smaller data sets. It uses retention indices for common standard polar and
 * standard non-polar stationary phases and molecular descriptors as input
 * features.
 * 
 *
 */
public abstract class SecondLevelModelForPolar extends Model {
	private static final int MIN_POLAR_COLUMN_NUMBER = TrainPolar.MIN_POLAR_COLUMN_NUMBER; // This number should be
																							// equal
																							// Columns.nonPolarMain.length+1
	private static final int MAX_COLUMN_NUMBER = TrainPolar.MAX_COLUMN_NUMBER;
	private ComputationGraph polarcnn = null;
	private ComputationGraph polarmlp = null;
	private ComputationGraph nonpolarcnn = null;
	private ComputationGraph nonpolarmlp = null;
	private Pair<double[][], double[]> trainSetDouble;
	private boolean useOnlyDescriptors = false;
	private boolean useOnlyRetentionIndices = false;

	/**
	 * Implementation of SecondLevelModelForPolar using random forest as
	 * second-level model. Hyperparameters are not optimal, maybe...
	 * 
	 *
	 */
	public static class SecondLevelRandomForestModel extends SecondLevelModelForPolar {
		private RandomForest forest;

		private static final Random rnd = new Random();

		@Override
		public HashMap<String, Object> randomParameters() {
			HashMap<String, Object> params = new HashMap<String, Object>();
			params.put("mtry", 10 + rnd.nextInt(300));
			params.put("ntrees", 100 + rnd.nextInt(10000));
			params.put("maxDepth", 2 + rnd.nextInt(30));
			params.put("maxNodes", 2 + rnd.nextInt(50));
			params.put("nodeSize", 2 + rnd.nextInt(30));
			params.put("subsample", 0.7 + 0.3 * Math.random());
			return params;
		}

		@Override
		public HashMap<String, Object> defaultParameters() {
			HashMap<String, Object> params = new HashMap<String, Object>();
			params.put("mtry", 116);
			params.put("ntrees", 2000);
			params.put("maxDepth", 25);
			params.put("maxNodes", 50);
			params.put("nodeSize", 4);
			params.put("subsample", 1.0);
			return params;
		}

		@Override
		public float predictRI(RetentionsDataset dataSet, int entry) throws CDKException {
			double[] features = featuresDouble(dataSet, entry);
			return (float) forest.predict(doubleFeaturesToDataFrame(features))[0];
		}

		@Override
		public void save(String filename) throws IOException {
			XStream xstream = new XStream();
			String xml = xstream.toXML(forest);
			FileWriter fw = new FileWriter(filename);
			fw.write(xml);
			fw.close();
		}

		@Override
		public void load(String filename) throws IOException {
			XStream xstream = new XStream();
			forest = (RandomForest) xstream.fromXML(new File(filename));
		}

		@Override
		public String train(HashMap<String, Object> params) throws CDKException, IOException {
			System.out.println("Random Forest training started...");
			forest = RandomForest.fit(Formula.lhs("label"), trainSetAsDataFrame(), (int) params.get("ntrees"),
					(int) params.get("mtry"), (int) params.get("maxDepth"), (int) params.get("maxNodes"),
					(int) params.get("nodeSize"), (double) params.get("subsample"));
			String val = this.validate(getValidationSet(), null);
			System.out.println(val);
			return val;
		}

	}

	/**
	 * Implementation of SecondLevelModelForPolar using gradient boosting as
	 * second-level model. Hyperparameters are not optimal, maybe...
	 * 
	 *
	 */
	public static class SecondLevelGBMModel extends SecondLevelModelForPolar {
		private GradientTreeBoost gbm;

		private static final Random rnd = new Random();

		@Override
		public HashMap<String, Object> randomParameters() {
			HashMap<String, Object> params = new HashMap<String, Object>();
			params.put("p", 0.001 + 0.3 * Math.random());
			params.put("ntrees", 100 + rnd.nextInt(10000));
			params.put("maxDepth", 10 + rnd.nextInt(10));
			params.put("maxNodes", 2 + rnd.nextInt(3));
			params.put("nodeSize", 2 + rnd.nextInt(3));
			params.put("shrinkage", 0.001 + 0.1 * Math.random());
			params.put("subsample", 0.1 + 0.8 * Math.random());
			return params;
		}

		@Override
		public HashMap<String, Object> defaultParameters() {
			HashMap<String, Object> params = new HashMap<String, Object>();
			params.put("p", 0.07);
			params.put("ntrees", 8000);
			params.put("maxDepth", 17);
			params.put("maxNodes", 3);
			params.put("nodeSize", 3);
			params.put("shrinkage", 0.027);
			params.put("subsample", 0.53);
			return params;
		}

		@Override
		public float predictRI(RetentionsDataset dataSet, int entry) throws CDKException {
			double[] features = featuresDouble(dataSet, entry);
			return (float) gbm.predict(doubleFeaturesToDataFrame(features))[0];
		}

		@Override
		public void save(String filename) throws IOException {
			XStream xstream = new XStream();
			String xml = xstream.toXML(gbm);
			FileWriter fw = new FileWriter(filename);
			fw.write(xml);
			fw.close();
		}

		@Override
		public void load(String filename) throws IOException {
			XStream xstream = new XStream();
			gbm = (GradientTreeBoost) xstream.fromXML(new File(filename));
		}

		@Override
		public String train(HashMap<String, Object> params) throws CDKException, IOException {
			System.out.println("GBM training started...");
			gbm = GradientTreeBoost.fit(Formula.lhs("label"), trainSetAsDataFrame(),
					Loss.huber((double) params.get("p")), (int) params.get("ntrees"), (int) params.get("maxDepth"),
					(int) params.get("maxNodes"), (int) params.get("nodeSize"), (double) params.get("shrinkage"),
					(double) params.get("subsample"));
			String val = this.validate(getValidationSet(), null);
			System.out.println(val);
			return val;
		}

	}

	/**
	 * Implementation of SecondLevelModelForPolar using support vector regression
	 * with gaussian kernel as second-level model.
	 * 
	 *
	 */
	public static class SecondLevelSVRModel extends SecondLevelModelForPolar {
		private KernelMachine<double[]> svr;

		@Override
		public HashMap<String, Object> randomParameters() {
			HashMap<String, Object> params = new HashMap<String, Object>();
			params.put("sigma", Math.random() * 50);
			params.put("eps", Math.random() * 5);
			params.put("C", Math.random() * 30000);
			params.put("tol", Math.random() * 3);
			return params;
		}

		@Override
		public HashMap<String, Object> defaultParameters() {
			HashMap<String, Object> params = new HashMap<String, Object>();
			params.put("sigma", 31.0);
			params.put("eps", 0.82);
			params.put("C", 31000.0);
			params.put("tol", 2.7);

			return params;
		}

		@Override
		public float predictRI(RetentionsDataset dataSet, int entry) throws CDKException {
			double[] features = featuresDouble(dataSet, entry);
			return (float) svr.predict(features);
		}

		@Override
		public void save(String filename) throws IOException {
			XStream xstream = new XStream();
			String xml = xstream.toXML(svr);
			FileWriter fw = new FileWriter(filename);
			fw.write(xml);
			fw.close();
		}

		@SuppressWarnings("unchecked")
		@Override
		public void load(String filename) throws IOException {
			XStream xstream = new XStream();
			svr = (KernelMachine<double[]>) xstream.fromXML(new File(filename));
		}

		@Override
		public String train(HashMap<String, Object> params) throws CDKException, IOException {
			System.out.println("SVR training started...");
			svr = SVR.fit(getTrainSetDouble().getLeft(), getTrainSetDouble().getRight(),
					new GaussianKernel((double) params.get("sigma")), (double) params.get("eps"),
					(double) params.get("C"), (double) params.get("tol"));
			String val = this.validate(getValidationSet(), null);
			System.out.println(val);
			return val;
		}
	}

	/**
	 * Implementation of SecondLevelModelForPolar using support vector regression
	 * with linear kernel as second-level model.
	 * 
	 *
	 */
	public static class SecondLevelLinearSVRModel extends SecondLevelModelForPolar {
		private KernelMachine<double[]> svr;

		public HashMap<String, Object> randomParameters() {
			HashMap<String, Object> params = new HashMap<String, Object>();
			params.put("tol", Math.random() * 10);
			params.put("C", Math.random() * 1000);
			params.put("eps", Math.random() * 10);
			return params;
		}

		@Override
		public HashMap<String, Object> defaultParameters() {
			HashMap<String, Object> params = new HashMap<String, Object>();
			params.put("eps", 1.0);
			params.put("C", 7.5);
			params.put("tol", 1.0);

			return params;
		}

		@Override
		public float predictRI(RetentionsDataset dataSet, int entry) throws CDKException {
			double[] features = featuresDouble(dataSet, entry);
			return (float) svr.predict(features);
		}

		@Override
		public void save(String filename) throws IOException {
			XStream xstream = new XStream();
			String xml = xstream.toXML(svr);
			FileWriter fw = new FileWriter(filename);
			fw.write(xml);
			fw.close();
		}

		@SuppressWarnings("unchecked")
		@Override
		public void load(String filename) throws IOException {
			XStream xstream = new XStream();
			svr = (KernelMachine<double[]>) xstream.fromXML(new File(filename));
		}

		@Override
		public String train(HashMap<String, Object> params) throws CDKException, IOException {
			System.out.println("SVR training started...");
			svr = SVR.fit(getTrainSetDouble().getLeft(), getTrainSetDouble().getRight(), new LinearKernel(),
					(double) params.get("eps"), (double) params.get("C"), (double) params.get("tol"));
			String val = this.validate(getValidationSet(), null);
			System.out.println(val);
			return val;
		}
	}

	/**
	 *
	 * @param useOnlyRetentionIndices_ if true only retention indices will be used
	 *                                 as input features for second-level model (no
	 *                                 descriptors)
	 */
	public void setUseOnlyRetentionIndices(boolean useOnlyRetentionIndices_) {
		useOnlyRetentionIndices = useOnlyRetentionIndices_;
	}

	/**
	 *
	 * @param useOnlyDescriptors_ if true only molecular descriptors will be used as
	 *                            input features for second-level model (no
	 *                            retention indices)
	 */
	public void setUseOnlyDescriptors(boolean useOnlyDescriptors_) {
		useOnlyDescriptors = useOnlyDescriptors_;
	}

	@Override
	public void init(RetentionsDataset trainSet_, float validationFraction) throws CDKException {
		throw new UnsupportedOperationException(
				"Descriptors generator is required for this model. Initialization without descriptors generator is not supported");
	}

	@Override
	public void init(RetentionsDataset trainSet_, RetentionsDataset validationSet_, Descriptors descriptorsGenerator_)
			throws CDKException {
		super.init(trainSet_, validationSet_, descriptorsGenerator_);
		trainSetDouble = dataSetToDouble(this.getTrainSet());
	}

	@Override
	public void init(RetentionsDataset trainSet_, float validationFraction_, Descriptors descriptorsGenerator_)
			throws CDKException {
		super.init(trainSet_, validationFraction_, descriptorsGenerator_);
		trainSetDouble = dataSetToDouble(this.getTrainSet());
	}

	/**
	 * Set neural networks that will be used for creation of input features
	 * (retention indices for polar and non-polar columns)
	 * 
	 * @param cnn1 CNN trained to predict retention indices for non-polar columns
	 * @param mlp1 MLP trained to predict retention indices for non-polar columns
	 * @param cnn2 CNN trained to predict retention indices for polar columns
	 * @param mlp2 MLP trained to predict retention indices for polar columns
	 */
	public void setNeuralNetworks(ComputationGraph cnn1, ComputationGraph mlp1, ComputationGraph cnn2,
			ComputationGraph mlp2) {
		polarcnn = cnn1;
		polarmlp = mlp1;
		nonpolarcnn = cnn2;
		nonpolarmlp = mlp2;
	}

	/**
	 * Set neural networks that will be used for creation of input features
	 * (retention indices for polar and non-polar columns) from the other second
	 * level model
	 * 
	 * @param a other second level model
	 */
	public void setPretrainedModelsFrom(SecondLevelModelForPolar a) {
		this.polarcnn = a.polarcnn;
		this.polarmlp = a.polarmlp;
		this.nonpolarcnn = a.nonpolarcnn;
		this.nonpolarmlp = a.nonpolarmlp;
	}

	/**
	 * 
	 * @return random hyperparameters (for hyperparameters tuning)
	 */
	public abstract HashMap<String, Object> randomParameters();

	/**
	 * 
	 * @return default (optimized) hyperparameters. Recommended to use with train()
	 *         method
	 */
	public abstract HashMap<String, Object> defaultParameters();

	/**
	 * 
	 * @param params hyperparameters (use defaultParameters() method)
	 * @return accuracy measures for validation set
	 * @throws CDKException cdk
	 * @throws IOException  io
	 */
	public abstract String train(HashMap<String, Object> params) throws CDKException, IOException;

	private float[] retentionForAllColumns(float[] fingerprintsArray, float[] descriptorsArray, float[][] smilesArray,
			int startColumnNum, ComputationGraph mlp, ComputationGraph cnn) {
		float[][] fingerprints = new float[MAX_COLUMN_NUMBER + 1 - startColumnNum][];
		float[][] descriptorsColumn = new float[MAX_COLUMN_NUMBER + 1 - startColumnNum][];
		float[][][] smiles = new float[MAX_COLUMN_NUMBER + 1 - startColumnNum][][];
		float[][] columns = new float[MAX_COLUMN_NUMBER + 1 - startColumnNum][];
		for (int i = startColumnNum; i < MAX_COLUMN_NUMBER + 1; i++) {
			fingerprints[i - startColumnNum] = fingerprintsArray;
			smiles[i - startColumnNum] = smilesArray;
			columns[i - startColumnNum] = Columns.columnAndColumnTypeOneHot(i);
			descriptorsColumn[i - startColumnNum] = RetentionsDataset.mergeArrays(columns[i - startColumnNum],
					descriptorsArray);
		}
		float[] retentionsMLP = mlp
				.output(new INDArray[] { Nd4j.create(descriptorsColumn), Nd4j.create(fingerprints) })[0]
						.toFloatVector();
		float[] retentionsCNN = cnn.output(new INDArray[] { Nd4j.create(smiles), Nd4j.create(columns) })[0]
				.toFloatVector();
		return RetentionsDataset.mergeArrays(retentionsMLP, retentionsCNN);
	}

	private float[] features(RetentionsDataset dataSet, int entry) throws CDKException {
		float[] descriptors = RetentionsDataset.mergeArrays(dataSet.descriptorsNoNaNs(entry, getDescriptorsGenerator()),
				dataSet.funcGroups(entry));

		if (useOnlyDescriptors) {
			return descriptors;
		}
		float[] fingerprints = dataSet.fingerprints(Chemoinformatics.FingerprintsType.ADDITIVE_CIRCULAR_4_1024_NO_SCALE,
				entry);
		float[][] smiles = new float[Chemoinformatics.SMILES_TOKENS][Chemoinformatics.SMILES_LEN];
		int[] tokens = Chemoinformatics.tokenize(dataSet.getSmiles(entry));
		for (int j = 0; j < tokens.length; j++) {
			smiles[tokens[j]][j] = 1.0F;
		}

		float[] polarOut = retentionForAllColumns(fingerprints, descriptors, smiles, MIN_POLAR_COLUMN_NUMBER, polarmlp,
				polarcnn);
		float[] nonpolarOut = retentionForAllColumns(fingerprints, descriptors, smiles, 0, nonpolarmlp, nonpolarcnn);
		if (!useOnlyRetentionIndices) {
			return RetentionsDataset.mergeArrays(descriptors, RetentionsDataset.mergeArrays(polarOut, nonpolarOut));
		}
		return RetentionsDataset.mergeArrays(polarOut, nonpolarOut);
	}

	/**
	 * 
	 * @param dataSet data set
	 * @param entry   number of entry in the data set
	 * @return features (double array) for this entry
	 * @throws CDKException
	 */
	public double[] featuresDouble(RetentionsDataset dataSet, int entry) throws CDKException {
		float[] f = features(dataSet, entry);
		double[] result = new double[f.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = f[i];
		}
		return result;
	}

	private Pair<double[][], double[]> dataSetToDouble(RetentionsDataset dataset) throws CDKException {
		double[][] resultFeatures = new double[dataset.size()][];
		double[] resultLabels = new double[dataset.size()];
		for (int i = 0; i < resultFeatures.length; i++) {
			resultFeatures[i] = featuresDouble(dataset, i);
			resultLabels[i] = dataset.getRetention(i);
		}
		return Pair.of(resultFeatures, resultLabels);
	}

	/**
	 * 
	 * @return data sets as double arrays. First entry - features, second entry -
	 *         labels
	 */
	public Pair<double[][], double[]> getTrainSetDouble() {
		return trainSetDouble;
	}

	private DataFrame doubleDatasetToDataFrame(Pair<double[][], double[]> data) {
		double[][] features = data.getLeft();
		double[] retentionsFlat = data.getRight();
		double[][] retentions = new double[features.length][1];

		String[] names = new String[features[0].length];
		for (int i = 0; i < names.length; i++) {
			names[i] = "" + i;
		}
		for (int i = 0; i < features.length; i++) {
			retentions[i] = new double[] { retentionsFlat[i] };
		}
		DataFrame dataFrame = DataFrame.of(retentions, "label").merge(DataFrame.of(features, names));
		return dataFrame;
	}

	private static RetentionsDataset[] cvSplit(int n, RetentionsDataset data) throws CDKException {
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
	 * 
	 * @return train set as DataFrame for Smile framework
	 */
	public DataFrame trainSetAsDataFrame() {
		return doubleDatasetToDataFrame(trainSetDouble);
	}

	/**
	 * 
	 * @param features double array
	 * @return features double array converted to DataFrame object for Smile
	 *         framework
	 */
	public static DataFrame doubleFeaturesToDataFrame(double[] features) {
		String[] names = new String[features.length];
		for (int i = 0; i < names.length; i++) {
			names[i] = "" + i;
		}
		DataFrame dataFrame = DataFrame.of(new double[][] { features }, names);
		return dataFrame;
	}

	/**
	 * Random hyperparameters tuning
	 * @param nAttempts number of tries
	 * @param bestModelFileName save best model to this file
	 * @param paramTuningLogFile log file name
	 * @return best hyperparameters
	 * @throws IOException io
	 * @throws CDKException cdk
	 */
	public HashMap<String, Object> hyperParametersTuning(int nAttempts, String bestModelFileName,
			String paramTuningLogFile) throws IOException, CDKException {
		float bestResult = Float.POSITIVE_INFINITY;
		HashMap<String, Object> bestParameters = null;
		FileWriter log = null;
		if (paramTuningLogFile != null) {
			log = new FileWriter(paramTuningLogFile);
		}
		for (int i = 0; i < nAttempts; i++) {
			HashMap<String, Object> parameters = randomParameters();
			String validationString = this.train(parameters);
			float mae = this.mae(validationString);
			if (mae < bestResult) {
				bestResult = mae;
				bestParameters = parameters;
				if (bestModelFileName != null) {
					this.save(bestModelFileName);
				}
			}
			if (log != null) {
				log.write(paramsToString(parameters) + " " + validationString + "\n");
				log.flush();
			}
		}
		if (bestModelFileName != null) {
			this.load(bestModelFileName);
		}
		if (log != null) {
			log.close();
		}
		return bestParameters;
	}

	/**
	 * Cross-validation
	 * @param params hyperparameters
	 * @param n number of folds
	 * @param data data set for cross-validation
	 * @param predictionsFile output file (for creation of plot)
	 * @return overall accuracy
	 * @throws CDKException cdk
	 * @throws IOException io
	 */
	public String crossValidation(HashMap<String, Object> params, int n, RetentionsDataset data, String predictionsFile)
			throws CDKException, IOException {
		boolean makePred = false;
		FileWriter pred = null;
		if (predictionsFile != null) {
			if (!predictionsFile.trim().equals("")) {
				makePred = true;
				pred = new FileWriter(predictionsFile);
			}
		}

		ArrayList<Float> deviations = new ArrayList<Float>();
		ArrayList<Float> percentageErrors = new ArrayList<Float>();
		RetentionsDataset testSets[] = cvSplit(n, data);
		String[] metrics = new String[n];
		RetentionsDataset oldTrain = this.getTrainSet();
		RetentionsDataset oldVal = this.getValidationSet();

		for (int i = 0; i < n; i++) {
			RetentionsDataset train = data.copy();
			RetentionsDataset test = testSets[i];
			train.filterIdentical(test);
			train.filterIdenticalByInchi(test);
			// validation set - one compound. zero is not allowed while really we DO NOT
			// NEED a validation set
			this.init(train, train.compoundsBasedSplitAndShuffle(1), this.getDescriptorsGenerator());
			this.train(params);
			int j = 0;
			float[] deviationsSubset = new float[test.size()];
			float[] percentageErrorsSubset = new float[test.size()];
			while (j < test.size()) {
				float testretention = test.getRetention(j);
				float retentionOut = predictRI(test, j);
				deviationsSubset[j] = Math.abs(testretention - retentionOut);
				percentageErrorsSubset[j] = Math.abs(100 * (testretention - retentionOut) / testretention);
				if (makePred) {
					pred.write(test.getSmiles(j) + " " + testretention + " " + retentionOut + "\n");
				}
				j++;
			}
			float rmse = 0.0F;
			for (int qj = 0; qj < test.size(); qj++) {
				rmse += deviationsSubset[qj] * deviationsSubset[qj];
				deviations.add(deviationsSubset[qj]);
				percentageErrors.add(percentageErrorsSubset[qj]);
			}
			rmse = (float) Math.sqrt(rmse / test.size());

			metrics[i] = "RMSE: " + rmse + " MAE: " + RetentionsDataset.mean(deviationsSubset) + " MPE: "
					+ RetentionsDataset.mean(percentageErrorsSubset) + " MdAE: "
					+ RetentionsDataset.median(deviationsSubset) + " MdPE: "
					+ RetentionsDataset.median(percentageErrorsSubset);
		}

		float[] deviationsFloat = new float[deviations.size()];
		float[] percentageErrorsFloat = new float[percentageErrors.size()];
		float rmse = 0.0F;
		for (int i = 0; i < deviations.size(); i++) {
			deviationsFloat[i] = deviations.get(i);
			percentageErrorsFloat[i] = percentageErrors.get(i);
			rmse += deviationsFloat[i] * deviationsFloat[i];
		}
		rmse = (float) Math.sqrt(rmse / deviations.size());
		String result = "RMSE: " + rmse + " MAE: " + RetentionsDataset.mean(deviationsFloat) + " MPE: "
				+ RetentionsDataset.mean(percentageErrorsFloat) + " MdAE: " + RetentionsDataset.median(deviationsFloat)
				+ " MdPE: " + RetentionsDataset.median(percentageErrorsFloat);
		if (makePred) {
			for (int i = 0; i < n; i++) {
				pred.write("Subset" + i + " " + metrics[i] + "\n");
			}
			pred.write("\nCV " + result);
			pred.close();
		}
		for (int i = 0; i < n; i++) {
			System.out.print("Subset" + i + " " + metrics[i] + "\n");
		}
		System.out.println("\nCV " + result);
		this.setTrainSet(oldTrain);
		this.setValidationSet(oldVal);
		return (result);
	}

	private String paramsToString(HashMap<String, Object> params) {
		String result = "";
		for (Map.Entry<String, Object> e : params.entrySet()) {
			result += e.getKey() + " " + e.getValue() + " ";
		}
		return result;
	}

}
