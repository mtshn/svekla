package ru.ac.phyche.gcms.svekla;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import org.openscience.cdk.exception.CDKException;

/**
 * 
 * Command-line application for the retention index prediction using the
 * pretrained model.
 *
 */
public class App {

	private static boolean isSmilesValid(String s) {
		try {
			String canonical = Chemoinformatics.canonical(s, false);
			Chemoinformatics.canonical(canonical, false);
			Chemoinformatics.canonical(canonical, true);
			Chemoinformatics.tokenize(canonical);
			Chemoinformatics.smilesToImage(canonical);
			Chemoinformatics.representation2d(canonical);
			Chemoinformatics.fingerprints(canonical,
					Chemoinformatics.FingerprintsType.ADDITIVE_CIRCULAR_4_1024_NO_SCALE);
			return true;
		} catch (CDKException e) {
			e.printStackTrace();
			return false;
		}
	}

	private static void printUsageInfo() {
		System.out.println("Usage:");
		System.out.println("For prediction of retention index for just one compound use SMILES string as"
				+ " commpand line parameter (default 5% Phenylmethyl siloxane column)");
		System.out.println("or SMILES string and column number, for example");
		System.out.println("CCCCCC 0");
		System.out.println("Other options:");
		System.out.println(
				"-RemoveUnsupportedCompounds <input file name> <output file name> <file name for unsupported compounds>");
		System.out.println("Remove lines from file which contain unsupported compounds. Third parameter is optional");
		System.out.println("-ShuffleAndSplit <Input data set file> <Separated part of the data set (file name)> "
				+ "<Remained part of the data set (file name)> <fraction to split (float number)>");
		System.out.println("Shuffle and split valid data set.");
		System.out.println("Dataset file format: multiple lines, one retention index record per line.");
		System.out.println("<SMILES molecule representation> <Kovac retention index> <Column type (integer value)>");
		System.out.println("Example (butane, DB-5 column):");
		System.out.println("CCCC 400 16");
		System.out.println("-ValidateForDataset <Data set file>");
		System.out.println("Print accuracy of the model for the data set.");
		System.out.println("-CountOverlap <Data set 1 file> <Data set 2 file>");
		System.out.println("Count compounds which are contained in both data sets.");
		System.out.println("-file <file name> <outpt file>");
		System.out.println("Predict retention index for all SMILES strings from file");
		System.out.println(
				"File format: one SMILES per line and optional number which denotes column. Examples of the line:");
		System.out.println("CCCCCC");
		System.out.println("or");
		System.out.println("CCCCCC 0");
		System.out.println("-file -allModels  <file name>  <outpt file>");
		System.out.println("Predict retention index for all SMILES strings from file using all models.");
		System.out.println("-columns");
		System.out.println("Show supported columns and respective numbers.");
	}

	private static void removeUnsupportedCompounds(String inputFile, String outputFile, String unsupportedFile) {
		BufferedReader inp;
		try {
			inp = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			FileWriter fw = new FileWriter(outputFile);
			ArrayList<String> unsupported = new ArrayList<String>();
			String s = inp.readLine();
			while (s != null) {
				if (!s.trim().equals("")) {
					String smiles = s.split("\\s+")[0];
					if (isSmilesValid(smiles)) {
						fw.write(s + "\n");
					} else {
						unsupported.add(s);
					}
				}
				s = inp.readLine();
			}
			fw.close();
			inp.close();
			if (unsupportedFile != null) {
				fw = new FileWriter(unsupportedFile);
				for (String str : unsupported) {
					fw.write(str + "\n");
				}
				fw.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.out.println(e.getMessage());
		}
	}

	private static void shuffleAndSplit(String inputFile, String output1File, String output2File, float fraction) {
		try {
			RetentionsDataset data = RetentionsDataset.loadFromFile(inputFile);
			data.makeCanoncalAll(false);
			RetentionsDataset data1 = data.compoundsBasedSplitAndShuffle(fraction);
			data1.saveToFile(output2File);
			data.saveToFile(output1File);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(e.getMessage());
		}
	}

	private static void validateForDataset(String inputFile) {
		StackingLinearMetaLearnerModel model = loadModel();
		try {
			RetentionsDataset data = RetentionsDataset.loadFromFile(inputFile);
			data.makeCanoncalAll(false);
			System.out.println("Linear_stacking_meta-model: " + model.validate(data, null));
			System.out.println("1DCNN: " + model.getModels()[0].validate(data, null));
			System.out.println("2DCNN: " + model.getModels()[1].validate(data, null));
			System.out.println("MLP: " + model.getModels()[2].validate(data, null));
			System.out.println("XGBoost: " + model.getModels()[3].validate(data, null));
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(e.getMessage());
		}
	}

	private static void predictForFile(String inputFile, String outputFile) {
		try {
			BufferedReader inp = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));
			FileWriter fw = new FileWriter(outputFile);
			String s = inp.readLine();
			fw.write("SMILES SMILES_(canonical_form) Retention_index\n");
			System.out.println("SMILES SMILES_(canonical_form) Retention_index\n");
			StackingLinearMetaLearnerModel model = loadModel();
			while (s != null) {
				if (!s.trim().equals("")) {
					String smiles = s.split("\\s+")[0];
					int column = s.split("\\s+").length == 1 ? 15 : Integer.parseInt(s.split("\\s+")[1]);
					String canonical = Chemoinformatics.canonical(smiles, false);
					String out = smiles + " " + canonical + " " + Columns.column(column) + " "
							+ model.predictRI(canonical, column) + " ";
					System.out.println(out);
					fw.write(out + "\n");
				}
				s = inp.readLine();
			}
			fw.close();
			inp.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(e.getMessage());
		}
	}

	private static void predictForFileAllModels(String inputFile, String outputFile) {
		try {
			BufferedReader inp = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));
			FileWriter fw = new FileWriter(outputFile);
			String s = inp.readLine();
			fw.write("SMILES SMILES_(canonical_form) Column Linear_stacking_meta-model 1DCNN 2DCNN MLP XGBoost\n");
			System.out.println(
					"SMILES SMILES_(canonical_form) Column Linear_stacking_meta-model 1DCNN 2DCNN MLP XGBoost\n");
			StackingLinearMetaLearnerModel model = loadModel();
			while (s != null) {
				if (!s.trim().equals("")) {
					String smiles = s.split("\\s+")[0];
					int column = s.split("\\s+").length == 1 ? 15 : Integer.parseInt(s.split("\\s+")[1]);
					String canonical = Chemoinformatics.canonical(smiles, false);
					String out = smiles + " " + canonical + " " + Columns.column(column) + " "
							+ model.predictRI(canonical, column) + " ";
					out = out + model.getModels()[0].predictRI(canonical, column) + " ";
					out = out + model.getModels()[1].predictRI(canonical, column) + " ";
					out = out + model.getModels()[2].predictRI(canonical, column) + " ";
					out = out + model.getModels()[3].predictRI(canonical, column) + " ";
					System.out.println(out);
					fw.write(out + "\n");
				}
				s = inp.readLine();
			}
			fw.close();
			inp.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(e.getMessage());
		}
	}

	private static void predictForOneCompound(String smiles, int column) {
		StackingLinearMetaLearnerModel model = loadModel();
		try {
			String canonicalSmiles = Chemoinformatics.canonical(smiles, false);
			System.out.println("Predicting retention index for compound " + smiles);
			System.out.println("Canonical form of SMILES " + canonicalSmiles);
			System.out.println("Column type number " + column);
			System.out.println("Column type " + Columns.column(column));
			System.out.println(
					Columns.isNonPolar(column) ? "Standard non-polar column" : "Semi-sandard non-polar column");
			System.out.println("============================================");
			System.out.println("**** Predicted retention index " + model.predictRI(canonicalSmiles, column) + "****");
			System.out.println("============================================");
			System.out.println("Predictions with various models:");
			System.out.println("CNN1D: " + model.getModels()[0].predictRI(canonicalSmiles, column));
			System.out.println("CNN2D: " + model.getModels()[1].predictRI(canonicalSmiles, column));
			System.out.println("MLP: " + model.getModels()[2].predictRI(canonicalSmiles, column));
			System.out.println("XGBoost: " + model.getModels()[3].predictRI(canonicalSmiles, column));
			System.out.println("Linear stacking meta-model: " + model.predictRI(canonicalSmiles, column));
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(e.getMessage());
		}

	}

	private static void countOverlap(String dataset1, String dataset2) {
		try {
			RetentionsDataset d1 = RetentionsDataset.loadFromFile(dataset1);
			RetentionsDataset d2 = RetentionsDataset.loadFromFile(dataset2);
			d1.makeCanoncalAll(false);
			d2.makeCanoncalAll(false);
			System.out.println("Note!!! Stereoisomers (cis/trans and optical) are considered as identical componds");
			System.out.println(d1.countIdenticalByCanonicalSmiles(d2) + " " + d2.countIdenticalByCanonicalSmiles(d1)
					+ " " + d1.countIdenticalByInchi(d2) + " " + d2.countIdenticalByInchi(d1));
		} catch (Exception e) {
			e.printStackTrace();
		}

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

	public static void main(String[] args) {
		if (args.length == 0) {
			printUsageInfo();
			System.exit(1);
		}
		boolean completed = false;
		if (args[0].equals("-RemoveUnsupportedCompounds")) {
			String unsupportedFile = null;
			if ((args.length < 3) || (args.length > 4)) {
				printUsageInfo();
				System.exit(1);
			} else {
				if (args.length == 4) {
					unsupportedFile = args[3];
				}
			}
			removeUnsupportedCompounds(args[1], args[2], unsupportedFile);
			completed = true;
		}
		if (args[0].equals("-ShuffleAndSplit")) {
			if (args.length != 5) {
				printUsageInfo();
				System.exit(1);
			} else {
				shuffleAndSplit(args[1], args[2], args[3], Float.parseFloat(args[4]));
			}
			completed = true;
		}
		if (args[0].equals("-ValidateForDataset")) {
			if (args.length != 2) {
				printUsageInfo();
				System.exit(1);
			} else {
				validateForDataset(args[1]);
			}
			completed = true;
		}
		if (args[0].equals("-CountOverlap")) {
			if (args.length != 3) {
				printUsageInfo();
				System.exit(1);
			} else {
				countOverlap(args[1], args[2]);
			}
			completed = true;
		}
		if (args[0].equals("-file")) {
			if (args.length == 2) {
				printUsageInfo();
				System.exit(1);
			}
			if (args.length == 3) {
				predictForFile(args[1], args[2]);
			}
			if (args.length == 4) {
				if (!args[1].equals("-allModels")) {
					printUsageInfo();
					System.exit(1);
				}
				predictForFileAllModels(args[2], args[3]);
			}
			completed = true;
		}
		if (args[0].equals("-columns")) {
			for (int i = 0; i < 36; i++) {
				System.out.println(i+" "+Columns.column(i));
			}
			completed = true;
		}
		if (completed) {
			System.exit(0);
		}
		if (isSmilesValid(args[0])) {
			if (args.length == 1) {
				predictForOneCompound(args[0], 15);
			} else {
				predictForOneCompound(args[0], Integer.parseInt(args[1]));
			}
			System.exit(0);
		} else {
			System.out.println(args[0] + " - is not correct SMILES representation of a molecule");
			printUsageInfo();
		}
	}
}
