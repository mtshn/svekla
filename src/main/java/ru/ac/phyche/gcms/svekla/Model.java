package ru.ac.phyche.gcms.svekla;

import java.io.FileWriter;
import java.io.IOException;

import org.openscience.cdk.exception.CDKException;

/**
 * 
 * Machine learning model for the retention index prediction. It is an abstract
 * class. StackingLinearMetaLearnerModel, NeuralNetModel, XGBoostModel classes
 * implement this class. Internally it stores data sets for validation and
 * training: RetentionsDataset validationSet and RetentionsDataset trainSet and
 * generator of descriptors. NOTE! MAE (mae) = mean absolute error, RMSE (rmse)
 * = root mean square error, MPE (mpe) = mean percentage error, MdAE (mdae) =
 * median absolute error, MdPE (mdpe) = median percentage error.
 *
 */
public abstract class Model {

	private RetentionsDataset trainSet;
	private RetentionsDataset validationSet;
	private Descriptors descriptorsGenerator = null;

	/**
	 * Predict retention index value for one entry from a data set. This method uses
	 * SMILES string and chromatographic column number. It (of course!) should not
	 * use retention index value that is stored with the entry.
	 * 
	 * @param dataSet data set that contains the entry for which retention index
	 *                should be predicted.
	 * @param entry   number of the entry in the data set
	 * @return predicted retention index value for entry-th entry.
	 * @throws CDKException CDK
	 */
	public abstract float predictRI(RetentionsDataset dataSet, int entry) throws CDKException;

	/**
	 * Save the trained model to file. It saves only parameters of model (in
	 * Deeplearning4j or XGBoost file format). It doesn't save any metadata,
	 * descriptor generator, validation set and so on.
	 * 
	 * @param filename file name
	 * @throws IOException IO
	 */
	public abstract void save(String filename) throws IOException;

	/**
	 * Load the trained model from file. It loads only parameters of model (in
	 * Deeplearning4j or XGBoost file format). It doesn't load any metadata,
	 * descriptor generator, validation set and so on. For neural networks: a state
	 * of the updater isn't saved by save(filename) and can't be loaded with this
	 * method.
	 * 
	 * @param filename file name
	 * @throws IOException IO
	 */
	public abstract void load(String filename) throws IOException;

	/**
	 * 
	 * @return generator (manager) of descriptors, which is used in this instance
	 */
	public Descriptors getDescriptorsGenerator() {
		return descriptorsGenerator;
	}

	/**
	 * 
	 * @param descriptorsGenerator generator (manager) of molecular descriptors
	 */
	public void setDescriptorsGenerator(Descriptors descriptorsGenerator) {
		this.descriptorsGenerator = descriptorsGenerator;
	}

	/**
	 * Simple getter.
	 * 
	 * @return training set
	 */
	public RetentionsDataset getTrainSet() {
		return trainSet;
	}

	/**
	 * Simple setter.
	 * 
	 * @param trainSet training set
	 */
	public void setTrainSet(RetentionsDataset trainSet) {
		this.trainSet = trainSet;
	}

	/**
	 * Simple getter.
	 * 
	 * @return validation set
	 */
	public RetentionsDataset getValidationSet() {
		return validationSet;
	}

	/**
	 * Simple setter
	 * 
	 * @param validationSet validation set
	 */
	public void setValidationSet(RetentionsDataset validationSet) {
		this.validationSet = validationSet;
	}

	/**
	 * This method doesn't alter trainSet_ instance. It creates deep copy of it
	 * internally. It sets up training and validation sets and initializes model.
	 * Note! This version can not be used for models which use descriptors. For them
	 * it will throw UnsupportedOperationException. Split is compounds based i.e. no
	 * identical compounds will be in both data sets. See
	 * RetentionsDataset.compoundsBasedSplitAndShuffle(fraction) method.
	 * 
	 * @param trainSet_          training and validation set (will be randomly
	 *                           split)
	 * @param validationFraction fraction of COMPOUNDS that will be used for
	 *                           validation (all entries for each compound will be
	 *                           used either for validation or training).
	 * @throws CDKException CDK. Also it throws UnsupportedOperationException if it
	 *                      is used for models that use descriptors.
	 */
	public void init(RetentionsDataset trainSet_, float validationFraction) throws CDKException {
		init(trainSet_, validationFraction, null);
	}

	/**
	 * This method doesn't alter trainSet_ instance. It creates deep copy of it
	 * internally. It sets up training and validation sets and initializes model.
	 * Note! This version can not be used for StackingLinearMetaLearnerModel. For it
	 * it will throw UnsupportedOperationException. Split is compounds based i.e. no
	 * identical compounds will be in both data sets. See
	 * RetentionsDataset.compoundsBasedSplitAndShuffle(fraction) method.
	 * 
	 * @param trainSet_             training and validation set (will be randomly
	 *                              split)
	 * @param validationFraction    fraction of COMPOUNDS that will be used for
	 *                              validation (all entries for each compound will
	 *                              be used either for validation or training).
	 * @param descriptorsGenerator_ molecular descriptors manager
	 * @throws CDKException CDK. Also it throws UnsupportedOperationException if it
	 *                      is used for StackingLinearMetaLearnerModel.
	 */
	public void init(RetentionsDataset trainSet_, float validationFraction, Descriptors descriptorsGenerator_)
			throws CDKException {
		this.descriptorsGenerator = descriptorsGenerator_;
		RetentionsDataset trainSetCopy = trainSet_.copy();
		trainSetCopy.shuffle();
		this.setValidationSet(trainSetCopy.compoundsBasedSplitAndShuffle(validationFraction));
		this.setTrainSet(trainSetCopy);
		if ((getTrainSet().countIdenticalByInchi(getValidationSet()) != 0)
				|| (getValidationSet().countIdenticalByInchi(getTrainSet()) != 0)) {
			throw new RuntimeException("Overlap between training set and validation set. Aborted!");
		}
	}

	/**
	 * This method doesn't alter trainSet_ instance. It creates deep copy of both
	 * trainingSet_ and validationSet_ internally. It sets up training and
	 * validation sets and initializes model. Note! This version can not be used for
	 * StackingLinearMetaLearnerModel. For it it will throw
	 * UnsupportedOperationException.
	 * 
	 * @param trainSet_             training set
	 * @param validationSet_        validation set. If there are any compounds which
	 *                              are contained in both training set and
	 *                              validation set the method will throws
	 *                              RuntimeException.
	 * @param descriptorsGenerator_ descriptors manager
	 * @throws CDKException CDK. Also it throws UnsupportedOperationException if it
	 *                      is used for StackingLinearMetaLearnerModel. If there are
	 *                      any compounds which are contained in both training set
	 *                      and validation set the method will throws
	 *                      RuntimeException.
	 * 
	 */
	public void init(RetentionsDataset trainSet_, RetentionsDataset validationSet_, Descriptors descriptorsGenerator_)
			throws CDKException {
		this.descriptorsGenerator = descriptorsGenerator_;
		RetentionsDataset trainSetCopy = trainSet_.copy();
		trainSetCopy.shuffle();
		this.setValidationSet(validationSet_.copy());
		this.setTrainSet(trainSetCopy);
		if ((getTrainSet().countIdenticalByInchi(getValidationSet()) != 0)
				|| (getValidationSet().countIdenticalByInchi(getTrainSet()) != 0)) {
			throw new RuntimeException("Overlap between training set and validation set. Aborted!");
		}
	}

	/**
	 * Predict retention index for specified SMILES string and type of
	 * chromatographic column. Model have to be trained or parameters have to be
	 * loaded.
	 * 
	 * @param smiles     SMILES string
	 * @param columnType type of chromatographic column (see Columns class).
	 * @return Kovac retention index value
	 * @throws CDKException CDK
	 */
	public float predictRI(String smiles, int columnType) throws CDKException {
		RetentionsDataset oneItem = RetentionsDataset
				.create(new RetentionsEntry[] { RetentionsEntry.instance(smiles, 0, columnType) });
		return predictRI(oneItem, 0);
	}

	/**
	 * Make prediction using this model for all entries from test data set.
	 * Calculate model accuracy (MAE, RMSE, MPE, MdAE, MdPE) for this model.
	 * 
	 * @param test            test data set
	 * @param predictionsFile name of file to which predictions will be saved. Can
	 *                        be null!
	 * @return A string with accuracy measures. E.g. "RMSE: 192.06718 MAE: 34.59384
	 *         MPE: 2.159631 MdAE: 17.081543 MdPE: 1.1565123"
	 * @throws CDKException Chemoinformatics exception
	 * @throws IOException  IO
	 */
	public String validate(RetentionsDataset test, String predictionsFile) throws CDKException, IOException {
		boolean makePred = false;
		FileWriter pred = null;
		if (predictionsFile != null) {
			if (!predictionsFile.trim().equals("")) {
				makePred = true;
				pred = new FileWriter(predictionsFile);
			}
		}

		int j = 0;
		float[] deviations = new float[test.size()];
		float[] percentageErrors = new float[test.size()];
		while (j < test.size()) {
			float testretention = test.getRetention(j);
			float retentionOut = predictRI(test, j);
			deviations[j] = Math.abs(testretention - retentionOut);
			percentageErrors[j] = Math.abs(100 * (testretention - retentionOut) / testretention);
			if (makePred) {
				pred.write(test.getSmiles(j) + " " + testretention + " " + retentionOut + "\n");
			}
			j++;
		}
		if (makePred) {
			pred.close();
		}
		float rmse = 0.0F;
		for (int qj = 0; qj < test.size(); qj++) {
			rmse += deviations[qj] * deviations[qj];
		}
		rmse = (float) Math.sqrt(rmse / test.size());

		String metrics = "RMSE: " + rmse + " MAE: " + RetentionsDataset.mean(deviations) + " MPE: "
				+ RetentionsDataset.mean(percentageErrors) + " MdAE: " + RetentionsDataset.median(deviations)
				+ " MdPE: " + RetentionsDataset.median(percentageErrors);
		return (metrics);
	}

	/**
	 * 
	 * @param test test data set
	 * @return MAE value for this data set
	 * @throws CDKException CDK
	 * @throws IOException  IO
	 */
	public float mae(RetentionsDataset test) throws CDKException, IOException {
		return Float.parseFloat(validate(test, null).split("\\s+")[3]);
	}

	/**
	 * 
	 * @param test test data set
	 * @return MAE value for this data set
	 * @throws CDKException CDK
	 * @throws IOException  IO
	 */
	public float mdae(RetentionsDataset test) throws CDKException, IOException {
		return Float.parseFloat(validate(test, null).split("\\s+")[7]);
	}

	/**
	 * 
	 * @param test test data set
	 * @return RMSE value for this data set
	 * @throws CDKException CDK
	 * @throws IOException  IO
	 */
	public float rmse(RetentionsDataset test) throws CDKException, IOException {
		return Float.parseFloat(validate(test, null).split("\\s+")[1]);
	}

	/**
	 * 
	 * @param test test data set
	 * @return MPE value for this data set
	 * @throws CDKException CDK
	 * @throws IOException  IO
	 */
	public float mpe(RetentionsDataset test) throws CDKException, IOException {
		return Float.parseFloat(validate(test, null).split("\\s+")[5]);
	}

	/**
	 * 
	 * @param test test data set
	 * @return MDPE value for this data set
	 * @throws CDKException CDK
	 * @throws IOException  IO
	 */
	public float mdpe(RetentionsDataset test) throws CDKException, IOException {
		return Float.parseFloat(validate(test, null).split("\\s+")[9]);
	}

	/**
	 * 
	 * @param validationString output string of validation(...) method
	 * @return MAE value
	 * @throws CDKException CDK
	 * @throws IOException  IO
	 */
	public float mae(String validationString) throws CDKException, IOException {
		return Float.parseFloat(validationString.split("\\s+")[3]);
	}

	/**
	 * 
	 * @param validationString output string of validation(...) method
	 * @return MdAE value
	 * @throws CDKException CDK
	 * @throws IOException  IO
	 */
	public float mdae(String validationString) throws CDKException, IOException {
		return Float.parseFloat(validationString.split("\\s+")[7]);
	}

	/**
	 * 
	 * @param validationString output string of validation(...) method
	 * @return RMSE value
	 * @throws CDKException CDK
	 * @throws IOException  IO
	 */
	public float rmse(String validationString) throws CDKException, IOException {
		return Float.parseFloat(validationString.split("\\s+")[1]);
	}

	/**
	 * 
	 * @param validationString output string of validation(...) method
	 * @return MPE value
	 * @throws CDKException CDK
	 * @throws IOException  IO
	 */
	public float mpe(String validationString) throws CDKException, IOException {
		return Float.parseFloat(validationString.split("\\s+")[5]);
	}

	/**
	 * 
	 * @param validationString output string of validation(...) method
	 * @return MdPE value
	 * @throws CDKException CDK
	 * @throws IOException  IO
	 */
	public float mdpe(String validationString) throws CDKException, IOException {
		return Float.parseFloat(validationString.split("\\s+")[9]);
	}
}
