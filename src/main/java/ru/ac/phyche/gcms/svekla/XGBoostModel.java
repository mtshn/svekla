package ru.ac.phyche.gcms.svekla;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.openscience.cdk.exception.CDKException;

import ml.dmlc.xgboost4j.java.Booster;
import ml.dmlc.xgboost4j.java.DMatrix;
import ml.dmlc.xgboost4j.java.XGBoost;
import ml.dmlc.xgboost4j.java.XGBoostError;

/**
 * XGBoost for retention index prediction.
 *
 */
public class XGBoostModel extends Model {

	private Booster x;
	private DMatrix validationSetDMatrix = null;
	private DMatrix trainSetDMatrix = null;
	private Pair<String, String> temporaryFileName = Pair.of("./XGBoost.train.ixt", "./XGBoost.val.ixt");
	private static final Random rnd = new Random();

	/**
	 * Random XGBoost parameters for random hyperparameters tuning! See source code!
	 * 
	 * @return random XGBoost hyperparameters
	 */
	public HashMap<String, Object> randomParameters() {
		HashMap<String, Object> params = new HashMap<String, Object>();
		params.put("eta", randomFloat(new float[] { 0.01F, 0.05F, 0.1F, 0.2F, 0.3F }));
		params.put("gamma", randomFloat(new float[] { 0.0F, 0.05F, 0.1F, 0.5F, 1F, 2F }));
		params.put("lambda", randomFloat(new float[] { 0.0F, 0.01F, 0.05F, 0.1F, 0.5F, 1F }));
		params.put("max_depth", randomInt(new int[] { 3, 5, 8, 12, 15, 18, 21 }));
		params.put("min_child_weight", randomInt(new int[] { 9, 12, 15, 18, 21, 24 }));
		params.put("subsample", randomFloat(new float[] { 0.3F, 0.5F, 1F }));
		params.put("colsample_bytree", randomFloat(new float[] { 0.4F, 0.5F, 0.6F, 0.7F }));
		params.put("objective", "reg:squarederror");
		return params;
	}

	/**
	 * Default XGBoost parameters! See XGBoost4j docs!!! ("eta", 0.05F), ("gamma",
	 * 0.05F), ("lambda", 0.05F), ("max_depth", 21), ("min_child_weight", 21),
	 * ("subsample", 0.5F), ("colsample_bytree", 0.5F), ("objective",
	 * "reg:squarederror");
	 * 
	 * @return default XGBoost hyperparameters
	 */
	public static HashMap<String, Object> xgboostDefaultParameters() {
		HashMap<String, Object> params = new HashMap<String, Object>();
		params.put("eta", 0.05F);
		params.put("gamma", 0.05F);
		params.put("lambda", 0.05F);
		params.put("max_depth", 21);
		params.put("min_child_weight", 21);
		params.put("subsample", 0.5F);
		params.put("colsample_bytree", 0.5F);
		params.put("objective", "reg:squarederror");
		return params;
	}

	/**
	 * File names for temporary files
	 * 
	 * @param train      training set (XGBoost file format)
	 * @param validation validation set (XGBoost file format)
	 */
	public void setTemporaryFileNames(String train, String validation) {
		this.temporaryFileName = Pair.of(train, validation);
	}

	/**
	 * Creates input features for entry-th entry of dataSet data set. Input
	 * features: molecular descriptors, functional groups, one-hot encoded polarity
	 * of the column and one-hot encoded column name. See
	 * Chemoinformatics.funcGroups method, Descriptors and Columns classes.
	 * 
	 * @param dataSet a data set
	 * @param entry   number of entry
	 * @return input features for XGBoost
	 * @throws CDKException CDK
	 */
	public float[] features(RetentionsDataset dataSet, int entry) throws CDKException {
		float[] d = dataSet.descriptorsNoNaNs(entry, getDescriptorsGenerator());
		float[] columns = Columns.columnAndColumnTypeOneHot(dataSet.getColumn(entry));
		float[] g = dataSet.funcGroups(entry);
		return RetentionsDataset.mergeArrays(RetentionsDataset.mergeArrays(columns, d), g);
	}

	@Override
	public void init(RetentionsDataset trainSet_, float validationFraction) throws CDKException {
		throw new UnsupportedOperationException(
				"Descriptors generator is required for this model. Initialization without descriptors generator is not supported");
	}

	@Override
	public void init(RetentionsDataset trainSet_, float validationFraction, Descriptors descriptorsGenerator_)
			throws CDKException {
		super.init(trainSet_, validationFraction, descriptorsGenerator_);
		trainSetDMatrix = dataSetToDMatrix(this.getTrainSet(), temporaryFileName.getLeft());
		validationSetDMatrix = dataSetToDMatrix(this.getValidationSet(), temporaryFileName.getRight());
	}

	@Override
	public void init(RetentionsDataset trainSet_, RetentionsDataset validationSet_, Descriptors descriptorsGenerator_)
			throws CDKException {
		super.init(trainSet_, validationSet_, descriptorsGenerator_);
		trainSetDMatrix = dataSetToDMatrix(this.getTrainSet(), temporaryFileName.getLeft());
		validationSetDMatrix = dataSetToDMatrix(this.getValidationSet(), temporaryFileName.getRight());
	}

	/**
	 * Train XGBoost model.
	 * 
	 * @param params XGBoost parameters. See XGBoost4j docs!
	 * @param n      number of estimators. See XGBoost4j docs!
	 * @return accuracy string (see Model.validate(...) method) for trained model
	 *         and validation set.
	 * @throws CDKException wrapped XGBoost exception
	 * @throws IOException  IO (temporary files are used)
	 */
	public String train(HashMap<String, Object> params, int n) throws CDKException, IOException {
		HashMap<String, DMatrix> watches = new HashMap<String, DMatrix>();
		watches.put("val", validationSetDMatrix);
		watches.put("train", trainSetDMatrix);
		try {
			System.out.println("XGBoost training started...");
			x = XGBoost.train(trainSetDMatrix, params, n, watches, null, null);
		} catch (XGBoostError e) {
			throw new CDKException(e.getMessage());
		}
		String val = this.validate(getValidationSet(), null);
		System.out.println(val);
		return val;
	}

	/**
	 * 
	 * @param params XGBoost parameters. See XGBoost4j docs!
	 * @return parameters as String
	 */
	public String paramsToString(HashMap<String, Object> params) {
		String result = "";
		for (Map.Entry<String, Object> e : params.entrySet()) {
			result += e.getKey() + " " + e.getValue() + " ";
		}
		return result;
	}

	/**
	 * Random hyperparameters tuning! Select the best hyperparameters set using
	 * validation set.
	 * 
	 * @param nAttempts          number of random sets of hyperparameters.
	 * @param estimators         number of XGBoost trees
	 * @param bestModelFileName  name of file to which parameters of the best
	 *                           XGBoost model will be saved. Can be null.
	 * @param paramTuningLogFile log file name
	 * @return the best set of hyperparameters
	 * @throws IOException  IO (temporary files are used)
	 * @throws CDKException wrapped XGBoost exception
	 */
	public HashMap<String, Object> hyperParametersTuning(int nAttempts, int estimators, String bestModelFileName,
			String paramTuningLogFile) throws IOException, CDKException {
		float bestResult = Float.POSITIVE_INFINITY;
		HashMap<String, Object> bestParameters = null;
		FileWriter log = null;
		if (paramTuningLogFile != null) {
			log = new FileWriter(paramTuningLogFile);
		}
		for (int i = 0; i < nAttempts; i++) {
			HashMap<String, Object> parameters = randomParameters();
			String validationString = this.train(parameters, estimators);
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

	@Override
	public float predictRI(RetentionsDataset dataSet, int entry) throws CDKException {
		float[] floatF = features(dataSet, entry);
		try {
			DMatrix f = new DMatrix(floatF, 1, floatF.length);
			return x.predict(f)[0][0];
		} catch (XGBoostError e) {
			throw new CDKException(e.getMessage());
		}
	}

	@Override
	public void save(String filename) throws IOException {
		try {
			x.saveModel(filename);
		} catch (XGBoostError e) {
			throw (new IOException(e.getMessage()));
		}
	}

	@Override
	public void load(String filename) throws IOException {
		try {
			x = XGBoost.loadModel(filename);
		} catch (XGBoostError e) {
			throw (new IOException(e.getMessage()));
		}
	}

	private DMatrix dataSetToDMatrix(RetentionsDataset dataset, String tmpFilename) throws CDKException {
		try {
			DMatrix result = dataSetToXGBooostDMatrixFile(dataset, tmpFilename);
			return result;
		} catch (IOException e) {
			throw new CDKException(e.getMessage());
		}
	}

	private DMatrix dataSetToXGBooostDMatrixFile(RetentionsDataset dataset, String filename)
			throws IOException, CDKException {
		FileWriter f = new FileWriter("./XGBoost.tmp");
		for (int i = 0; i < dataset.size(); i++) {
			f.write(dataset.getRetention(i) + " ");
			float[] features = features(dataset, i);
			for (int j = 0; j < features.length; j++) {
				f.write(j + ":" + features[j] + " ");
			}
			f.write("\n");
		}
		f.close();
		try {
			DMatrix mat = new DMatrix("./XGBoost.tmp");
			return mat;
		} catch (XGBoostError e) {
			throw new IOException(e.getMessage());
		}
	}

	private static float randomFloat(float[] a) {
		return a[rnd.nextInt(a.length)];
	}

	private static int randomInt(int[] a) {
		return a[rnd.nextInt(a.length)];
	}
}
