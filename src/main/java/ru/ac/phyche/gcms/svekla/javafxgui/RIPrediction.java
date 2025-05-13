package ru.ac.phyche.gcms.svekla.javafxgui;

import org.apache.commons.lang3.tuple.Pair;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import javafx.scene.control.TextField;
import javafx.scene.layout.Background;
import javafx.scene.paint.Color;
import javafx.scene.paint.Paint;
import ru.ac.phyche.gcms.svekla.Chemoinformatics;
import ru.ac.phyche.gcms.svekla.Descriptors;
import ru.ac.phyche.gcms.svekla.Model;
import ru.ac.phyche.gcms.svekla.RetentionsDataset;
import ru.ac.phyche.gcms.svekla.StackingLinearMetaLearnerModel;
import ru.ac.phyche.gcms.svekla.TrainPolar;

public class RIPrediction {

	private static Paint color(float deltaRI) {
		if (deltaRI > 200) {
			return Paint.valueOf("#990000");
		}
		if (deltaRI > 120) {
			return Paint.valueOf("#BB0000");
		}
		if (deltaRI > 100) {
			return Paint.valueOf("#FF0000");
		}
		if (deltaRI > 90) {
			return Paint.valueOf("#FF7700");
		}
		if (deltaRI > 80) {
			return Paint.valueOf("#FFFF00");
		}
		if (deltaRI > 70) {
			return Paint.valueOf("#77FF00");
		}
		if (deltaRI > 50) {
			return Paint.valueOf("#33FF22");
		}
		if (deltaRI > 30) {
			return Paint.valueOf("#11FF11");
		}
		return Paint.valueOf("#11AA11");
	}

	private static StackingLinearMetaLearnerModel loadModel() {
		StackingLinearMetaLearnerModel result = new StackingLinearMetaLearnerModel();
		String fileCNN1D = "./models/CNN1D.nn";
		String fileCNN2D = "./models/CNN2D.nn";
		String fileMLP = "./models/MLP2.nn";
		String fileXGBoost = "./models/XGBoost.xgboost";
		String fileLinearMeta = "./models/linearMetaModel.nn";
		String fileDescriptors = "./models/descriptors_info.txt";
		try {
			result.loadFiveModelsFromFiles(fileCNN1D, fileCNN2D, fileMLP, fileXGBoost, fileLinearMeta, fileDescriptors);
			if (Math.abs(result.predictRI("CCCCCC", 15) - 600) > 25) {
				throw new RuntimeException("Incorrect model");
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(e.getMessage());
			System.out.println("Error while loading and testing model.");
			System.out.println("All model related files must be located in folder \"models\" in current directory.");
			System.out.println("Six files are required.");
			System.out.println("./models/CNN1D.nn");
			System.out.println("./models/CNN2D.nn");
			System.out.println("./models/MLP2.nn");
			System.out.println("./models/XGBoost.xgboost");
			System.out.println("./models/linearMetaModel.nn");
			System.out.println("./models/descriptors_info.txt");
			System.exit(1);
		}
		return result;
	}

	private static Model[] loadModelsSvekla() {

		try {
			Model[] result = new Model[4];
			String nnFolder = "./models_polar";
			TrainPolar.MLP mlp = new TrainPolar.MLP();
			TrainPolar.CNN cnn = new TrainPolar.CNN();
			TrainPolar.MLP mlpPolar = new TrainPolar.MLP();
			TrainPolar.CNN cnnPolar = new TrainPolar.CNN();
			mlp.load(nnFolder + "/mlp.nn");
			cnn.load(nnFolder + "/cnn.nn");
			mlpPolar.load(nnFolder + "/mlpPolar.nn");
			cnnPolar.load(nnFolder + "/cnnPolar.nn");
			Pair<float[], float[]> minmaxArray = Descriptors.readFromFile("./models_polar/descriptors_info.txt")
					.getMinMaxArray();
			Descriptors d = Descriptors.instance(Descriptors.descriptors2DBut_nAtomLAC_And_MolIP, minmaxArray.getLeft(),
					minmaxArray.getRight(), false);
			mlp.setDescriptorsGenerator(d);
			mlpPolar.setDescriptorsGenerator(d);
			result[0] = cnn;
			result[1] = cnnPolar;
			result[2] = mlp;
			result[3] = mlpPolar;
			for (int i = 0; i < 4; i++) {
				if (Math.abs(result[0].predictRI("CCCCCC", 15) - 600) > 25) {
					throw new RuntimeException("Incorrect model");
				}
			}
			return result;
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(e.getMessage());
			System.out.println("Error while loading and testing model.");
			System.exit(1);
			return null;
		}
	}

	public StackingLinearMetaLearnerModel model = loadModel();
	public Model[] modelsSvekla = loadModelsSvekla();

	private static boolean checkIfAlkane(String smiles) {
		String q = smiles.trim().replace('C', ' ').trim();
		if (q.trim().length() == 0) {
			return true;
		}
		return false;
	}

	private static String smilesToMolecularFormula(String s) throws CDKException {
		SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IAtomContainer mol = parser.parseSmiles(s.trim());
		Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		arom.apply(mol);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
		adder.addImplicitHydrogens(mol);
		String molForm = MolecularFormulaManipulator.getString(MolecularFormulaManipulator.getMolecularFormula(mol));
		molForm = molForm + " (" + MolecularFormulaManipulator.getMass(
				MolecularFormulaManipulator.getMolecularFormula(mol), AtomContainerManipulator.MostAbundant) + ")";
		return molForm;
	}

	public float[] predictNonpolarPolar(String smiles) {
		try {
			String canonicalSmiles = Chemoinformatics.canonical(smiles, false);
			int column = 15;
			float stackingIndex = model.predictRI(canonicalSmiles, column);
			float sveklaCNNPolar = modelsSvekla[1].predictRI(canonicalSmiles, column);
			float sveklaMLPPolar = modelsSvekla[3].predictRI(canonicalSmiles, column);
			float polar = (float) (0.5 * sveklaMLPPolar + 0.5 * sveklaCNNPolar);
			return new float[] { stackingIndex, polar };
		} catch (CDKException e) {
			return new float[] { 0, 0 };
		}
	}

	public void predictAndWrite(String smiles, TextField accessResults, TextField nonpolarResult,
			TextField sveklaResults, TextField polarResult, TextField nonpolarDev, TextField polarDev,
			TextField nonpolarTarget, TextField polarTarget, TextField molecularFormula) {
		String canonicalSmiles;
		try {
			canonicalSmiles = Chemoinformatics.canonical(smiles, false);
			if (Chemoinformatics.weight(canonicalSmiles) < 60) {
				throw new CDKException("");
			}
			if (canonicalSmiles.length() < 5) {
				throw new CDKException("");
			}

			int column = 15;
			float cnn1dIndex = (float) 100.0 * canonicalSmiles.trim().length();
			float cnn2dIndex = (float) 100.0 * canonicalSmiles.trim().length();
			float mlpIndex = (float) 100.0 * canonicalSmiles.trim().length();
			float xgboostIndex = (float) 100.0 * canonicalSmiles.trim().length();
			float stackingIndex = (float) 100.0 * canonicalSmiles.trim().length();
			float sveklaCNN = (float) 100.0 * canonicalSmiles.trim().length();
			float sveklaMLP = (float) 100.0 * canonicalSmiles.trim().length();
			float sveklaCNNPolar = (float) 100.0 * canonicalSmiles.trim().length();
			float sveklaMLPPolar = (float) 100.0 * canonicalSmiles.trim().length();

			if (!checkIfAlkane(smiles)) {
				cnn1dIndex = model.getModels()[0].predictRI(canonicalSmiles, column);
				cnn2dIndex = model.getModels()[1].predictRI(canonicalSmiles, column);
				mlpIndex = model.getModels()[2].predictRI(canonicalSmiles, column);
				xgboostIndex = model.getModels()[3].predictRI(canonicalSmiles, column);
				stackingIndex = model.predictRI(canonicalSmiles, column);
				sveklaCNN = modelsSvekla[0].predictRI(canonicalSmiles, column);
				sveklaMLP = modelsSvekla[2].predictRI(canonicalSmiles, column);
				sveklaCNNPolar = modelsSvekla[1].predictRI(canonicalSmiles, column);
				sveklaMLPPolar = modelsSvekla[3].predictRI(canonicalSmiles, column);
			}
			accessResults.setText("CNN1D: " + (int) cnn1dIndex + " CNN2D: " + (int) cnn2dIndex + " MLP: "
					+ (int) mlpIndex + " XGBoost: " + (int) xgboostIndex + " Weighted average: " + (int) stackingIndex
					+ " Median 4 models: "
					+ (int) RetentionsDataset.median(new float[] { cnn1dIndex, cnn2dIndex, mlpIndex, xgboostIndex }));
			nonpolarResult.setText((int) stackingIndex + "");
			sveklaResults.setText("CNN non-polar: " + (int) sveklaCNN + " MLP non-polar: " + (int) sveklaMLP
					+ " CNN polar: " + (int) sveklaCNNPolar + " MLP polar: " + (int) sveklaMLPPolar + " Average polar: "
					+ (int) (0.5 * sveklaMLPPolar + 0.5 * sveklaCNNPolar));
			polarResult.setText("" + (int) (0.5 * sveklaMLPPolar + 0.5 * sveklaCNNPolar));
			float targetPolar = 0;
			float targetNonpolar = 0;
			try {
				targetPolar = Float.parseFloat(polarTarget.getText());
			} catch (Throwable e) {
			}
			try {
				targetNonpolar = Float.parseFloat(nonpolarTarget.getText());
			} catch (Throwable e) {
			}
			if (targetPolar > 0) {
				int diff = (int) Math.abs((0.5 * sveklaMLPPolar + 0.5 * sveklaCNNPolar) - targetPolar);
				polarDev.setText("" + diff);
				polarDev.setBackground(Background.fill(color(diff)));
			} else {
				polarDev.setText("");
				polarDev.setBackground(Background.fill(Paint.valueOf("#FFFFFF")));
			}
			if (targetNonpolar > 0) {
				int diff = (int) Math.abs(stackingIndex - targetNonpolar);
				nonpolarDev.setText("" + diff);
				nonpolarDev.setBackground(Background.fill(color(diff)));

			} else {
				nonpolarDev.setText("");
				nonpolarDev.setBackground(Background.fill(Paint.valueOf("#FFFFFF")));
			}
			molecularFormula.setText(smilesToMolecularFormula(canonicalSmiles));
		} catch (CDKException e) {
			e.printStackTrace();
			polarDev.setText("");
			polarDev.setBackground(Background.fill(Paint.valueOf("#FFFFFF")));
			nonpolarDev.setText("");
			nonpolarDev.setBackground(Background.fill(Paint.valueOf("#FFFFFF")));
			accessResults.setText("");
			nonpolarResult.setText("");
			sveklaResults.setText("");
			polarResult.setText("");
		}
	}

}
