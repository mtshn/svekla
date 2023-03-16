package ru.ac.phyche.gcms.svekla;

import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.graph.MergeVertex;
import org.deeplearning4j.nn.conf.layers.Convolution2D;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.GlobalPoolingLayer;
import org.deeplearning4j.nn.conf.layers.LossLayer;
import org.deeplearning4j.nn.conf.layers.PoolingType;
import org.deeplearning4j.nn.conf.layers.SubsamplingLayer;
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
 * Deep 2D deep convolutional neural network with two inputs for the retention
 * index prediction. First input features set: the one-hot encoded
 * 2D-representation of molecule (structure sketch), second input features set:
 * information about chromatographic column. See .pdf files with the article and
 * the supplementary material and source code of this class. See also
 * Chemoinformatics.representation2d(...) method.
 *
 */
public class CNN2DFromDepictionModel extends NeuralNetModel {

	@Override
	public void initNN() {
		ComputationGraphConfiguration conf = new NeuralNetConfiguration.Builder().weightInit(WeightInit.RELU)
				.updater(new Adam(0.0003)).graphBuilder()
				.addLayer("CNN4",
						new Convolution2D.Builder(new int[] { 4, 4 }, new int[] { 1, 1 }).nIn(29).nOut(50)
								.activation(Activation.RELU).build(),
						"INPUT")
				.addLayer("SUBSAMPLING0",
						new SubsamplingLayer.Builder(SubsamplingLayer.PoolingType.MAX, new int[] { 2, 2 },
								new int[] { 2, 2 }).build(),
						"CNN4")
				.addLayer("CNN5",
						new Convolution2D.Builder(new int[] { 4, 4 }, new int[] { 1, 1 }).nIn(50).nOut(300)
								.activation(Activation.RELU).build(),
						"SUBSAMPLING0")
				.addLayer("SUBSAMPLING1",
						new SubsamplingLayer.Builder(SubsamplingLayer.PoolingType.MAX, new int[] { 2, 2 },
								new int[] { 2, 2 }).build(),
						"CNN5")
				.addLayer("CNN6",
						new Convolution2D.Builder(new int[] { 4, 4 }, new int[] { 1, 1 }).nIn(300).nOut(300)
								.activation(Activation.RELU).build(),
						"SUBSAMPLING1")
				.addLayer("POOLING", new GlobalPoolingLayer.Builder(PoolingType.AVG).build(), "CNN6")
				.addVertex("MERGE", new MergeVertex(), "POOLING", "COLUMN_INFO")
				.addLayer("DENSE0", new DenseLayer.Builder().nIn(338).nOut(600).activation(Activation.RELU).build(),
						"MERGE")
				.addInputs("INPUT", "COLUMN_INFO")
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
		float[][][][] features = new float[indices.length][][][];
		float[][] columnInfo = new float[indices.length][];
		for (int i = 0; i < indices.length; i++) {
			features[i] = dataSet.representation2d(indices[i]);
			columnInfo[i] = Columns.columnAndColumnTypeOneHot(dataSet.getColumn(indices[i]));
		}
		INDArray columnInfoINDArray = Nd4j.create(columnInfo);
		INDArray featuresINDArray = Nd4j.create(features);
		return new INDArray[] { featuresINDArray, columnInfoINDArray };
	}
}
