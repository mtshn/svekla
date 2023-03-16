package ru.ac.phyche.gcms.svekla;

import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.graph.ElementWiseVertex;
import org.deeplearning4j.nn.conf.graph.MergeVertex;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.LossLayer;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.activations.Activation;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.learning.config.Adam;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.openscience.cdk.exception.CDKException;

/**
 * 
 * Deep residual multilayer perceptron with two inputs for retention index
 * prediction. First input features set: molecular fingerprints, second input
 * features set: information about column and molecular descriptors. See .pdf
 * files with the article and the supplementary material and source code of this
 * class.
 *
 */
public class MLPFromDescriptorsAndFingerprints extends NeuralNetModel {
	/**
	 * Creates input features (second input) for entry-th entry of dataSet data set.
	 * Input features: molecular descriptors, functional groups, one-hot encoded
	 * polarity of the column and one-hot encoded column name. See
	 * Chemoinformatics.funcGroups method, Descriptors and Columns classes.
	 * 
	 * @param dataSet a data set
	 * @param entry   number of entry
	 * @return input features for XGBoost
	 * @throws CDKException CDK
	 */
	public float[] descriptorsFeatures(RetentionsDataset dataSet, int entry) throws CDKException {
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
	public void initNN() {
		ComputationGraphConfiguration conf = new NeuralNetConfiguration.Builder().weightInit(WeightInit.RELU)
				.updater(new Adam(0.0003)).graphBuilder()
				.addLayer("DESCRIPTORS0",
						new DenseLayer.Builder().nIn(366).nOut(300).activation(Activation.TANH).build(), "INPUT0")
				.addLayer("DESCRIPTORS1",
						new DenseLayer.Builder().nIn(300).nOut(300).activation(Activation.RELU).build(), "DESCRIPTORS0")
				.addLayer("FINGERPRINTS_IN",
						new DenseLayer.Builder().nIn(1024).nOut(1200).activation(Activation.RELU).build(), "INPUT1")
				.addLayer("FINGERPRINTS0",
						new DenseLayer.Builder().nIn(1200).nOut(1200).dropOut(0.95).activation(Activation.RELU).build(),
						"FINGERPRINTS_IN")
				.addLayer("FINGERPRINTS1",
						new DenseLayer.Builder().nIn(1200).nOut(1200).dropOut(0.95).activation(Activation.RELU).build(),
						"FINGERPRINTS0")
				.addVertex("ADD0", new ElementWiseVertex(ElementWiseVertex.Op.Add), "FINGERPRINTS_IN", "FINGERPRINTS1")
				.addLayer("FINGERPRINTS2",
						new DenseLayer.Builder().nIn(1200).nOut(1200).dropOut(0.95).activation(Activation.RELU).build(),
						"ADD0")
				.addLayer("FINGERPRINTS3",
						new DenseLayer.Builder().nIn(1200).nOut(1200).dropOut(0.95).activation(Activation.RELU).build(),
						"FINGERPRINTS2")
				.addVertex("ADD1", new ElementWiseVertex(ElementWiseVertex.Op.Add), "ADD0", "FINGERPRINTS3")

				.addVertex("MERGE", new MergeVertex(), "ADD1", "DESCRIPTORS1")
				.addLayer("DENSE0", new DenseLayer.Builder().nIn(1500).nOut(600).activation(Activation.RELU).build(),
						"MERGE")
				.addInputs("INPUT0", "INPUT1")
				.addLayer("DENSE1", new DenseLayer.Builder().nIn(600).nOut(1).activation(Activation.IDENTITY).build(),
						"DENSE0")
				.addLayer("OUT",
						new LossLayer.Builder().lossFunction(LossFunctions.LossFunction.MEAN_ABSOLUTE_ERROR).build(),
						"DENSE1")
				.setOutputs("OUT").build();
		ComputationGraph nn = new ComputationGraph(conf);
		nn.init();
		this.setNn(nn);
	}

	@Override
	public INDArray[] nextBatchInput(RetentionsDataset dataSet, int[] indices) throws CDKException {
		float[][] featuresDescriptors = new float[indices.length][];
		float[][] fingerprints = new float[indices.length][];
		for (int i = 0; i < indices.length; i++) {
			featuresDescriptors[i] = descriptorsFeatures(dataSet, indices[i]);
			fingerprints[i] = dataSet.fingerprints(Chemoinformatics.FingerprintsType.ADDITIVE_CIRCULAR_4_1024_NO_SCALE,
					indices[i]);
		}
		INDArray descriptorsINDArray = Nd4j.create(featuresDescriptors);
		INDArray fingerprintsINDArray = Nd4j.create(fingerprints);
		return new INDArray[] { descriptorsINDArray, fingerprintsINDArray };
	}
}
