package ru.ac.phyche.gcms.svekla;

import java.io.File;
import java.io.IOException;

import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.LossLayer;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.activations.Activation;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.learning.config.Adam;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.openscience.cdk.exception.CDKException;

/**
 * 
 * Linear model which combines predictions from other models and allows to
 * obtain more accurate retention index predictions.
 *
 */
public class StackingLinearMetaLearnerModel extends Model {
	private Model[] models = null;
	private ComputationGraph nn = null;
	private INDArray trainFeaturesINDArray = null;
	private INDArray retentionsINDArray = null;

	/**
	 * The instance stores multiple instances of Model class, that will be used as
	 * input for the linear model.
	 * 
	 * @return array of 1-st level models.
	 */
	public Model[] getModels() {
		return models;
	}

	/**
	 * @param models array of 1-st level models.
	 */
	public void setModels(Model[] models) {
		this.models = models;
	}

	/**
	 * 
	 * @return parameters of linear regression.
	 */
	public float[] getParams() {
		INDArray params = nn.params();
		return params.toFloatVector();
	}

	private ComputationGraphConfiguration createComputationGraphConfiguration(int modelsLength) {
		ComputationGraphConfiguration conf = new NeuralNetConfiguration.Builder().weightInit(WeightInit.XAVIER)
				.updater(new Adam(0.001)).graphBuilder()
				.addLayer("DENSE0",
						new DenseLayer.Builder().nIn(modelsLength).nOut(1).activation(Activation.IDENTITY).build(),
						"INPUT")
				.addLayer("OUT",
						new LossLayer.Builder().lossFunction(LossFunctions.LossFunction.MEAN_ABSOLUTE_ERROR).build(),
						"DENSE0")
				.addInputs("INPUT").setOutputs("OUT").build();
		return conf;
	}

	@Override
	public void init(RetentionsDataset trainSet_, float validationFraction) throws CDKException {
		throw new UnsupportedOperationException(
				"Descriptors generator is required for this model. Initialization without descriptors generator is not supported");
	}

	@Override
	public void init(RetentionsDataset trainSet_, float validationFraction, Descriptors descriptorsGenerator_)
			throws CDKException {
		throw new UnsupportedOperationException(
				"Pre-trained models are required for this model. Initialization without pre-trained models are not supported");
	}

	@Override
	public void init(RetentionsDataset trainSet_, RetentionsDataset validationSet, Descriptors descriptorsGenerator_)
			throws CDKException {
		throw new UnsupportedOperationException(
				"Pre-trained models are required for this model. Initialization without pre-trained models are not supported");
	}

	/**
	 * This method doesn't alter trainSet_ instance. It creates deep copy of both
	 * trainingSet_ and validationSet_ internally. It sets up training and
	 * validation sets and initializes model.
	 * 
	 * @param trainSet_             training set
	 * @param validationSet_        validation set. If there are any compounds which
	 *                              are contained in both training set and
	 *                              validation set the method will throws
	 *                              RuntimeException.
	 * @param models_               array of 1-st level models. The instance stores
	 *                              multiple instances of Model class, that will be
	 *                              used as input for the linear model.
	 * @param descriptorsGenerator_ descriptors manager
	 * @throws CDKException CDK. Also it throws UnsupportedOperationException if it
	 *                      is used for StackingLinearMetaLearnerModel. If there are
	 *                      any compounds which are contained in both training set
	 *                      and validation set the method will throws
	 *                      RuntimeException.
	 * 
	 */
	public void init(RetentionsDataset trainSet_, RetentionsDataset validationSet_, Descriptors descriptorsGenerator_,
			Model[] models_) throws CDKException {
		super.init(trainSet_, validationSet_, descriptorsGenerator_);
		this.models = models_;
		float[][] inputFeatures = new float[getTrainSet().size()][models.length];
		float[][] retentions = new float[getTrainSet().size()][];
		for (int i = 0; i < getTrainSet().size(); i++) {
			retentions[i] = new float[] { getTrainSet().getRetention(i) / 1000F };
			for (int j = 0; j < models.length; j++) {
				inputFeatures[i][j] = models[j].predictRI(getTrainSet(), i) / 1000F;
			}
		}
		this.trainFeaturesINDArray = Nd4j.create(inputFeatures);
		this.retentionsINDArray = Nd4j.create(retentions);
		ComputationGraphConfiguration conf = createComputationGraphConfiguration(models.length);
		this.nn = new ComputationGraph(conf);
		this.nn.init();
	}

	/**
	 * Train linear model.
	 * 
	 * @param nIters                     number of training iteration
	 * @param printScoreEveryKIterations k. Print score every k iterations
	 * @param validateEveryMIterations   m. Validate and print accuracy every m
	 *                                   iterations.
	 * @return trained model
	 * @throws IOException  IO
	 * @throws CDKException CDK (while feature creation)
	 */
	public StackingLinearMetaLearnerModel trainMultipleIterations(int nIters, int printScoreEveryKIterations,
			int validateEveryMIterations) throws IOException, CDKException {
		int n = nIters;
		int k = printScoreEveryKIterations;
		int m = validateEveryMIterations;
		for (int i = 0; i < n; i++) {
			DataSet batch0 = new DataSet(trainFeaturesINDArray, retentionsINDArray);
			nn.fit(batch0);
			if (i % k == 0) {
				System.out.println("Score: " + nn.score(batch0));
			}
			if (i % m == 0) {
				System.gc();
				String validationString = this.validate(getValidationSet(), null);
				System.out.println("Training iteration " + i + " ; Accuracy: " + validationString);
			}
		}
		return this;
	}

	/**
	 * Load this linear combining (stacking) model and all supported 1-st level
	 * models. Order of models: CNN1D, CNN2D, MLP, XGBoost.
	 * 
	 * @param fileCNN1D       name of file with CNN1D parameters
	 * @param fileCNN2D       name of file with CNN2D parameters
	 * @param fileMLP         name of file with MLP parameters
	 * @param fileXGBoost     name of file with XGBoost parameters
	 * @param fileLinearMeta  name of file with linear model parameters
	 *                        (deeplearning4j file format)
	 * @param fileDescriptors descriptor generator
	 * @throws IOException IO
	 */
	public void loadFiveModelsFromFiles(String fileCNN1D, String fileCNN2D, String fileMLP, String fileXGBoost,
			String fileLinearMeta, String fileDescriptors) throws IOException {
		Descriptors descriptorsGenerator = Descriptors.readFromFile(fileDescriptors);
		Model[] models = new Model[4];
		CNN1DFromSMILESModel modelCNN1D = new CNN1DFromSMILESModel();
		CNN2DFromDepictionModel modelCNN2D = new CNN2DFromDepictionModel();
		MLPFromDescriptorsAndFingerprints modelMLP = new MLPFromDescriptorsAndFingerprints();
		XGBoostModel modelXGBoost = new XGBoostModel();
		modelCNN1D.load(fileCNN1D);
		modelCNN2D.load(fileCNN2D);
		modelMLP.load(fileMLP);
		modelMLP.setDescriptorsGenerator(descriptorsGenerator);
		modelXGBoost.load(fileXGBoost);
		modelXGBoost.setDescriptorsGenerator(descriptorsGenerator);
		models[0] = modelCNN1D;
		models[1] = modelCNN2D;
		models[2] = modelMLP;
		models[3] = modelXGBoost;
		this.setModels(models);
		this.load(fileLinearMeta);
	}

	@Override
	public float predictRI(RetentionsDataset dataSet, int entry) throws CDKException {
		float[][] input = new float[1][models.length];
		for (int i = 0; i < models.length; i++) {
			input[0][i] = models[i].predictRI(dataSet, entry) / 1000F;
		}
		return 1000F * nn.output(Nd4j.create(input))[0].toFloatVector()[0];
	}

	@Override
	public void save(String filename) throws IOException {
		nn.save(new File(filename), false);
	}

	@Override
	public void load(String filename) throws IOException {
		nn = ComputationGraph.load(new File(filename), false);
	}

}
