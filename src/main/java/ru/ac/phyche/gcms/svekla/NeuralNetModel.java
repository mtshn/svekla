package ru.ac.phyche.gcms.svekla;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.MultiDataSet;
import org.nd4j.linalg.factory.Nd4j;
import org.openscience.cdk.exception.CDKException;

/**
 * 
 * Abstract class for neural net-based retention index prediction model.
 *
 */
public abstract class NeuralNetModel extends Model {

	private int exampleNow = 0;
	private ComputationGraph nn = null;

	/**
	 * Initializes neural network (create ComputationGraph)
	 */
	public abstract void initNN();

	/**
	 * Input features for some data set entries grouped to the batch. If indices[]
	 * is {3,4,5} - input batch for 3-th, 4-th and 5-th entry from array will be
	 * created.
	 * 
	 * @param dataSet a data set
	 * @param indices indices of entries in the data set for which the batch should
	 *                be created.
	 * @return batch input for ComputationGraph with multiple inputs
	 * @throws CDKException CDK
	 */
	public abstract INDArray[] nextBatchInput(RetentionsDataset dataSet, int[] indices) throws CDKException;

	/**
	 * getter
	 * 
	 * @return computation graph (neural network)
	 */
	public ComputationGraph getNn() {
		return nn;
	}

	/**
	 * Setter
	 * 
	 * @param nn computation graph (neural network)
	 */
	public void setNn(ComputationGraph nn) {
		this.nn = nn;
	}

	@Override
	public void init(RetentionsDataset trainSet_, float validationFraction) throws CDKException {
		super.init(trainSet_, validationFraction);
		this.initNN();
	}

	@Override
	public void init(RetentionsDataset trainSet_, float validationFraction, Descriptors descriptorsGenerator_)
			throws CDKException {
		super.init(trainSet_, validationFraction, descriptorsGenerator_);
		this.initNN();
	}

	@Override
	public void init(RetentionsDataset trainSet_, RetentionsDataset validationSet_, Descriptors descriptorsGenerator_)
			throws CDKException {
		super.init(trainSet_, validationSet_, descriptorsGenerator_);
		this.initNN();
	}

	@Override
	public float predictRI(RetentionsDataset dataSet, int entry) throws CDKException {
		return 1000F * nn.output(this.nextBatchInput(dataSet, new int[] { entry }))[0].toFloatVector()[0];
	}

	private int nextTrainSample() {
		if (this.exampleNow < this.getTrainSet().size() - 1) {
			this.exampleNow = this.exampleNow + 1;
			return this.exampleNow - 1;
		} else {
			this.exampleNow = 0;
			return (getTrainSet().size() - 1);
		}
	}

	private int[] nextTrainBatch(int batchSize) {
		int[] result = new int[batchSize];
		for (int i = 0; i < batchSize; i++) {
			result[i] = nextTrainSample();
		}
		return result;
	}

	private float[][] retentions(int[] indices) {
		float[][] result = new float[indices.length][];
		for (int i = 0; i < indices.length; i++) {
			result[i] = new float[] { getTrainSet().getRetention(indices[i]) / 1000F };
		}
		return result;
	}

	/**
	 * Train one iteration (one batch)
	 * 
	 * @param printScore if TRUE it will print score for the batch via
	 *                   system.out.println()
	 * @param batchSize  batch size for neural network
	 * @throws CDKException e
	 */
	public void trainOneIteration(boolean printScore, int batchSize) throws CDKException {
		int[] batch = nextTrainBatch(batchSize);
		INDArray[] input = this.nextBatchInput(getTrainSet(), batch);
		INDArray[] retentions = new INDArray[] { Nd4j.create(retentions(batch)) };

		MultiDataSet batch0 = new MultiDataSet(input, retentions);
		nn.fit(batch0);
		if (printScore) {
			System.out.println("Score: " + nn.score(batch0));
		}
	}

	/**
	 * Train neural network.
	 * 
	 * @param nIters                     number of training iteration
	 * @param printScoreEveryKIterations k. Print score every k iterations
	 * @param validateEveryMIterations   m. Validate and print accuracy every m
	 *                                   iterations.
	 * @param batchSize                  batch size
	 * @param bestModelFileName          name of file to which model parameters will
	 *                                   be saved for that iteration at which there
	 *                                   were best score for validation set. Can be
	 *                                   null. Don't use bestModelFileName == null
	 *                                   and loadBest == TRUE.
	 * @param loadBest                   if TRUE parameters for that iteration at
	 *                                   which there were best score for validation
	 *                                   set will be loaded. If TRUE - no further
	 *                                   training will be possible!
	 * @param trainingLogFile            log file name
	 * @return number of that iteration at which there were best score for
	 *         validation set.
	 * @throws IOException  IO
	 * @throws CDKException CDK (while feature creation)
	 */
	public int trainMultipleIterations(int nIters, int printScoreEveryKIterations, int validateEveryMIterations,
			int batchSize, String bestModelFileName, boolean loadBest, String trainingLogFile)
			throws IOException, CDKException {
		int n = nIters;
		int k = printScoreEveryKIterations;
		int m = validateEveryMIterations;
		float bestResult = Float.POSITIVE_INFINITY;
		int bestI = 0;
		FileWriter log = null;
		if (trainingLogFile != null) {
			log = new FileWriter(trainingLogFile);
		}
		for (int i = 0; i < n; i++) {
			trainOneIteration(i % k == 0, batchSize);
			if (i % m == 0) {
				System.gc();
				String validationString = this.validate(getValidationSet(), null);
				float mae = this.mae(validationString);
				if (mae < bestResult) {
					bestResult = mae;
					bestI = i;
					if (bestModelFileName != null) {
						nn.save(new File(bestModelFileName), false);
					}
				}
				System.out.println("Training iteration " + i + " ; Accuracy: " + validationString);
				if (log != null) {
					log.write("Training iteration " + i + " ; Accuracy: " + validationString + "\n");
				}
			}
		}
		if (bestModelFileName != null) {
			if (loadBest) {
				this.load(bestModelFileName);
			}
		}
		if (log != null) {
			log.close();
		}
		return bestI;
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
