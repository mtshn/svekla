package ru.ac.phyche.gcms.svekla;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class NISTRetentionConverter {

	private static final boolean useStereoChemistry = false;

	private static boolean compareInChIkeys(String inchiKey1, String inchiKey2) {
		if (useStereoChemistry) {
			return (inchiKey1.trim().equals(inchiKey2.trim()));
		} else {
			return (inchiKey1.trim().split("\\-")[0].equals(inchiKey2.trim().split("\\-")[0]));
		}
	}

	private static void convertTmpFileToDataSet(String inp, String outp) throws IOException {
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outp))));
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inp))));
		String s = reader.readLine();

		while (s != null) {
			ArrayList<String> block = new ArrayList<String>();
			while (!s.equals("***END***")) {
				block.add(s);
				s = reader.readLine();
			}
			String[] blockArr = new String[block.size()];
			for (int k = 0; k < block.size(); k++) {
				blockArr[k] = block.get(k);
			}
			for (int j = 4; j < blockArr.length; j++) {
				String[] split = blockArr[j].split("\\s+");
				String columnName = split[6].trim();
				float rI = Float.parseFloat(split[3]);
				if (blockArr[j].contains("Lee_RI")) {
					rI = 127.7F + 4.5269F * rI + 2.6193F * 0.001F * rI * rI
							+ rI * rI * rI * 0.01F * 0.01F * 0.001F * 5F;
					// Lee to Kovac
				}
				if ((blockArr[j].contains("Standard_non-polar")) || (blockArr[j].contains("Semi-standard_non-polar"))) {
					if (rI != 0) {
						writer.write(blockArr[3].trim() + " " + rI + " " + Columns.columnNum(columnName.trim()) + "\n");
					}
				}
			}
			s = reader.readLine();
		}

		reader.close();
		writer.close();
	}

	private static HashMap<String, String> mapRIDatabaseIDName(String smilesFile) throws IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(smilesFile))));
		String s = reader.readLine();
		HashMap<String, String> result = new HashMap<String, String>();
		while (s != null) {
			String name = s.trim();
			name = name.toUpperCase().replace(" ", "_");
			while (!s.contains("Library ID =")) {
				s = reader.readLine();
				if (s.contains("M  END")) {
					reader.close();
					throw (new RuntimeException("Incorrect DB file (no id)"));
				}
			}
			String id = s.split("\\=")[s.split("\\=").length - 1].trim();
			result.put(id, name);
			while (!s.contains("$$$$")) {
				s = reader.readLine();
			}
			s = reader.readLine();
		}
		reader.close();
		return result;
	}

	private static void convertMSPFile(String msp, String outp, String smilesFile, String failedLog, String warningFile,
			String structuresFile) throws IOException {
		HashMap<String, String> mapRIDatabaseIDName = mapRIDatabaseIDName(structuresFile);
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outp))));
		BufferedWriter writerWarnings = new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(new File(warningFile))));

		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(msp))));
		BufferedReader readerSmiles = new BufferedReader(
				new InputStreamReader(new FileInputStream(new File(smilesFile))));
		HashMap<String, String> namesAndSmiles = new HashMap<String, String>();
		BufferedWriter writerFailed = new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(new File(failedLog))));
		String s = readerSmiles.readLine();
		while (s != null) {
			String[] split = s.split("\\s+");
			namesAndSmiles.put(split[0], split[1]);
			s = readerSmiles.readLine();
		}
		readerSmiles.close();
		s = reader.readLine();
		while (s != null) {
			if (!s.contains("Name:")) {
				reader.close();
				writer.close();
				writerFailed.close();
				writerWarnings.close();
				throw (new RuntimeException("Incorrect msp file (no name)"));
			}
			String name = s.split("Name\\:")[1].trim();
			name = name.toUpperCase().replace(" ", "_");
			if (name.length() > 80) {
				name = name.substring(0, 80);
			}
			s = reader.readLine();
			if (!s.contains("InChIKey:")) {
				reader.close();
				writer.close();
				writerFailed.close();
				writerWarnings.close();
				throw (new RuntimeException("Incorrect msp file (no inchikey)"));
			}
			String inchiKey = s.split("InChIKey\\:")[1].trim();
			while ((!s.contains("NIST#:")) && (!s.contains("CAS#:"))) {
				if (s.contains("Name:")) {
					reader.close();
					writer.close();
					writerFailed.close();
					writerWarnings.close();
					throw (new RuntimeException("Incorrect msp file (no id)"));
				}
				s = reader.readLine();
			}
			String id = s.split("\\#\\:")[1].trim();
			if (s.contains("CAS")) {
				id = "C" + id.split("\\-")[0] + id.split("\\-")[1] + id.split("\\-")[2];
			} else {
				int idNum = Integer.parseInt(id);
				if (idNum > 2000000) {
					id = "R" + (idNum - 2000000);
				} else {
					id = "U" + idNum;
				}
			}
			s = reader.readLine();
			if (!s.contains("DB#:")) {
				reader.close();
				writer.close();
				writerFailed.close();
				writerWarnings.close();
				throw (new RuntimeException("Incorrect msp file (no id db)"));
			}
			String idDb = s.split("\\#\\:")[1].trim();
			String smiles = namesAndSmiles.get(mapRIDatabaseIDName.get(idDb));
			if (!name.equals(mapRIDatabaseIDName.get(idDb))) {
				writerWarnings.write("different names in MSP and DBU files" + name + " " + mapRIDatabaseIDName.get(idDb)
						+ " " + idDb + "\n");
			}
			if (smiles == null) {
				writerFailed.write("Can_not_find_structure_by_full_name " + id + " " + name + " " + inchiKey + "\n");
			} else {
				if (!smiles.contains("PARSING_FAILED")) {
					try {
						if (!compareInChIkeys(Chemoinformatics.smilesToInchiKey(smiles), inchiKey.trim())) {
							writerFailed.write("NIST_inchikey_does_not_equal_CDK_inchikey " + id + " " + name + " "
									+ inchiKey + "\n");
						} else {
							writer.write(id + " " + smiles.trim() + " " + name + " " + inchiKey + "\n");
						}
					} catch (CDKException e) {
						writerFailed
								.write("Can_not_convert_smiles_to_inchikey " + id + " " + name + " " + inchiKey + "\n");
					}
				}
			}
			while ((s != null) && (!s.contains("Name"))) {
				s = reader.readLine();
			}
		}
		writerFailed.close();
		writer.close();
		reader.close();
		writerWarnings.close();
	}

	private static void convertRIFile(String ri, String smilesFile, String outp, String failedLog, String warningFile)
			throws IOException {
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outp))));
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(smilesFile))));
		BufferedWriter writerWarnings = new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(new File(warningFile))));
		HashMap<String, String> smiles = new HashMap<String, String>();
		HashMap<String, String> inchiKeys = new HashMap<String, String>();
		HashMap<String, String> names = new HashMap<String, String>();

		BufferedWriter writerFailed = new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(new File(failedLog))));
		String s = reader.readLine();
		while (s != null) {
			String[] split = s.split("\\s+");
			smiles.put(split[0], split[1]);
			names.put(split[0], split[2]);
			inchiKeys.put(split[0], split[3]);
			s = reader.readLine();
		}

		reader.close();
		reader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(ri))));
		s = reader.readLine();
		String casNow = "";
		String smiNow = "NONAME";
		ArrayList<String> compoundInfo = new ArrayList<String>();
		compoundInfo.clear();
		int j = 0;
		HashSet<String> synonims = new HashSet<String>();
		while (s != null) {
			String[] split = s.split("\\s+");
			if (!split[0].trim().equals(casNow)) {
				if (!casNow.equals("")) {
					if (!smiNow.equals("NONAME")) {
						writer.write("***BEGIN***\n");
						writer.write(j + "\n");
						j++;
						writer.write(casNow + "\n");
						writer.write(smiNow + "\n");

						for (int i = 0; i < compoundInfo.size(); i++) {
							writer.write(compoundInfo.get(i) + "\n");
						}
						writer.write("***END***\n");
						if (!synonims.contains(names.get(casNow))) {
							writerWarnings.write("Different names in RI and structures files " + casNow + "\n");
						}
					} else {
						writerFailed.write("No SMILES was found " + casNow + "\n");
					}

				}
				casNow = split[0];
				compoundInfo.clear();
				synonims.clear();
				smiNow = "NONAME";
			}
			if (smiNow.equals("NONAME")) {
				String name = removeNonSymbols(split[1].trim().toUpperCase().replace(" ", "_"));
				if (name.length() > 80) {
					name = name.substring(0, 80);
				}
				String smiGet = smiles.get(casNow);
				if (smiGet != null) {
					smiNow = smiGet;
					try {
						if (!compareInChIkeys(Chemoinformatics.smilesToInchiKey(smiNow), inchiKeys.get(casNow))) {
							throw new CDKException("");
						}
					} catch (CDKException e) {
						writerFailed.close();
						writerWarnings.close();
						writer.close();
						reader.close();
						throw new RuntimeException();
					}
					synonims.add(name);
				}
			}
			compoundInfo.add(s);
			s = reader.readLine();
			if (s != null) {
				if (!s.trim().equals("")) {
					try {
						Double.parseDouble(s.split("\\s+")[3].trim());
					} catch (NumberFormatException e) {
						writerFailed.write("skipped RI record: " + s + "\n");
						s = reader.readLine();
					}
				}
			}
		}
		writerFailed.close();
		writerWarnings.close();
		writer.close();
		reader.close();
	}

	private static boolean isSymbol(char c) {
		return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || (c == '*') || (c == '.')
				|| (c == ',') || (c == ';') || (c == '"') || (c == '\'') || (c == '\\') || (c == '/') || (c == ':')
				|| (c == '_') || (c == '#') || (c == '^') || (c == '%') || (c == '&') || (c == '{') || (c == '}')
				|| (c == '[') || (c == ']') || (c == '(') || (c == ')') || (c == '+') || (c == '-') || (c == '|')
				|| (c == '=') || (c == '@') || (c == '#') || (c == '!') || (c == '$');
	}

	private static String removeNonSymbols(String s) {
		StringBuffer s1 = new StringBuffer(s);
		for (int i = 0; i < s1.length(); i++) {
			if (!isSymbol(s1.charAt(i))) {
				s1.deleteCharAt(i);
			}
		}
		return (s1.toString());
	}

	private static String atomContainerlToSmiles(IAtomContainer mol, boolean stereochemistry) throws CDKException {
		Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		AtomContainerManipulator.suppressHydrogens(mol);
		arom.apply(mol);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
		adder.addImplicitHydrogens(mol);
		SmilesGenerator sg = null;
		if (!stereochemistry) {
			sg = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols ^ SmiFlavor.AtomicMass);
		} else {
			sg = new SmilesGenerator(SmiFlavor.Absolute | SmiFlavor.UseAromaticSymbols ^ SmiFlavor.AtomicMass);
		}
		String smiles = sg.create(mol);
		return smiles;
	}

	private static void sdfToSmiles(String inp, String outp, String failedLog) throws IOException {

		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outp))));
		BufferedWriter writerFailed = new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(new File(failedLog))));

		File sdfFile = new File(inp);
		IteratingSDFReader reader = new IteratingSDFReader(new FileInputStream(sdfFile),
				DefaultChemObjectBuilder.getInstance());

		int failed = 0;
		String title = "";
		while (reader.hasNext()) {
			try {
				title = "No title loaded";
				IAtomContainer mol = (IAtomContainer) reader.next();
				title = mol.getProperty(CDKConstants.TITLE);
				title = title.trim().toUpperCase().replace(" ", "_");
				String canon = Chemoinformatics.canonical(atomContainerlToSmiles(mol, useStereoChemistry),
						useStereoChemistry);
				if (!Chemoinformatics.canonical(canon, useStereoChemistry).equals(canon)) {
					throw (new Exception("Canonical smiles incorrect (stereoisomers)"));
				}
				String canonNoStereo = Chemoinformatics.canonical(canon, false);
				if (!Chemoinformatics.canonical(canonNoStereo, false).equals(canonNoStereo)) {
					throw (new Exception("Canonical smiles incorrect (no stereoisomers)"));
				}
				Chemoinformatics.tokenize(canon);
				Chemoinformatics.tokenize(canonNoStereo);
				Chemoinformatics.tokenize(Chemoinformatics.canonical(canonNoStereo, true));
				Chemoinformatics.representation2d(canon);
				Chemoinformatics.representation2d(canonNoStereo);
				Chemoinformatics.representation2d(Chemoinformatics.canonical(canonNoStereo, true));
				writer.write(title + " " + canon + "\n");
			} catch (Exception e) {
				failed++;
				writer.write(title + " ");
				writer.write("PARSING_FAILED" + "\n");
				writerFailed.write(("Failed to process " + failed + " " + title + " Error: " + e.getMessage() + "\n"));

			}
		}
		writer.close();
		reader.close();
		writerFailed.close();
	}

	private static void extractTableNistRI(String inp, String outp, String tmpFile) throws IOException {
		InputStreamReader nistReader = new InputStreamReader(new FileInputStream(new File(inp)));
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(tmpFile))));
		Integer c = nistReader.read();
		while (c != -1) {
			c = nistReader.read();
			if (c == 0x09 || c == 0x00) {
				writer.write("\n");
			} else if (c == (int) ' ') {
				writer.write("_");
			} else if (isSymbol((char) (int) c)) {
				writer.write(c);
			}
		}
		nistReader.close();
		writer.close();
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(tmpFile))));
		writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outp))));
		String s = reader.readLine();
		int tokensCount = 0;
		String[] line = new String[50];
		while (s != null) {
			if (s.trim().equals("Packed") || s.trim().equals("Capillary") || s.trim().equals("Other")) {
				for (int i = 0; i < tokensCount - 4; i++) {
					writer.write(line[i].trim() + " ");
				}
				writer.write("\n");
				line[0] = line[tokensCount - 4];
				line[1] = line[tokensCount - 3];
				line[2] = line[tokensCount - 2];
				line[3] = line[tokensCount - 1];
				line[4] = s.trim();
				tokensCount = 5;
			} else {
				if (tokensCount < 30) {
					tokensCount++;
					line[tokensCount - 1] = s.trim();
				}
			}
			s = "";
			while (s != null && s.trim().equals("")) {
				s = reader.readLine();
			}

		}
		writer.close();
		reader.close();
	}

	private static void columnList(String inp, String outp) throws IOException {
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outp))));
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inp))));
		String s = reader.readLine();
		HashMap<String, Integer> columns = new HashMap<String, Integer>();
		HashMap<String, String> columnsPolarity = new HashMap<String, String>();

		while (s != null) {
			ArrayList<String> block = new ArrayList<String>();
			while (!s.equals("***END***")) {
				block.add(s);
				s = reader.readLine();
			}
			String[] blockArr = new String[block.size()];
			for (int k = 0; k < block.size(); k++) {
				blockArr[k] = block.get(k);
			}
			for (int j = 4; j < blockArr.length; j++) {
				String[] split = blockArr[j].split("\\s+");
				String columnName = split[6].trim();
				if (columns.containsKey(columnName)) {
					Integer o = columns.get(columnName.trim());
					columns.replace(columnName, o + 1);
				} else {
					columns.put(columnName, 1);
					if (blockArr[j].contains("Standard_non-polar")) {
						columnsPolarity.put(columnName, "Nonpolar");
					}
					if (blockArr[j].contains("Semi-standard_non-polar")) {
						columnsPolarity.put(columnName, "Seminonpolar");
					}
					if (blockArr[j].contains("Standard_polar")) {
						columnsPolarity.put(columnName, "Polar");
					}
				}
			}
			s = reader.readLine();
		}

		for (String col : columns.keySet()) {
			writer.write(col + " " + columns.get(col) + " " + columnsPolarity.get(col) + "\n");
		}

		reader.close();
		writer.close();
	}

	public static void main(String[] args) {
		try {
			extractTableNistRI("./NIST_RI/ri.dat", "./dataSets/NIST_preprocessing/tmp/stage1_preprocessing.tmp",
					"./dataSets/NIST_preprocessing/tmp/stage0_preprocessing.tmp");
			sdfToSmiles("./NIST_RI/USRSTRUC.DB", "./dataSets/NIST_preprocessing/tmp/structures.smi",
					"./dataSets/NIST_preprocessing/logs/sdfParsing.log");
			convertMSPFile("./NIST_RI/NISTRI.MSP", "./dataSets/NIST_preprocessing/tmp/stage2_preprocessing.tmp",
					"./dataSets/NIST_preprocessing/tmp/structures.smi",
					"./dataSets/NIST_preprocessing/logs/stage2_preprocessing.log",
					"./dataSets/NIST_preprocessing/logs/stage2_warnings.txt", "./NIST_RI/USRSTRUC.DB");
			convertRIFile("./dataSets/NIST_preprocessing/tmp/stage1_preprocessing.tmp",
					"./dataSets/NIST_preprocessing/tmp/stage2_preprocessing.tmp",
					"./dataSets/NIST_preprocessing/tmp/stage3_preprocessing.tmp",
					"./dataSets/NIST_preprocessing/logs/stage3_preprocessing.log",
					"./dataSets/NIST_preprocessing/logs/stage3_warnings.txt");
			columnList("./dataSets/NIST_preprocessing/tmp/stage3_preprocessing.tmp",
					"./dataSets/NIST_preprocessing/logs/columns.txt");
			convertTmpFileToDataSet("./dataSets/NIST_preprocessing/tmp/stage3_preprocessing.tmp", "./dataSets/nist.ri");
			RetentionsDataset nist = RetentionsDataset.loadFromFile("./dataSets/nist.ri");
			nist.makeCanoncalAll(false);
			RetentionsDataset[] splitCV = new RetentionsDataset[10];
			splitCV[0] = nist.compoundsBasedSplitAndShuffle(1F / 10F);
			splitCV[1] = nist.compoundsBasedSplitAndShuffle(1F / 9F);
			splitCV[2] = nist.compoundsBasedSplitAndShuffle(1F / 8F);
			splitCV[3] = nist.compoundsBasedSplitAndShuffle(1F / 7F);
			splitCV[4] = nist.compoundsBasedSplitAndShuffle(1F / 6F);
			splitCV[5] = nist.compoundsBasedSplitAndShuffle(1F / 5F);
			splitCV[6] = nist.compoundsBasedSplitAndShuffle(1F / 4F);
			splitCV[7] = nist.compoundsBasedSplitAndShuffle(1F / 3F);
			splitCV[8] = nist.compoundsBasedSplitAndShuffle(1F / 2F);
			splitCV[9] = nist;
			for (int i = 0; i < 10; i++) {
				System.out.println("Split " + i + " " + splitCV[i].size() + " " + splitCV[i].compounds().size());
				splitCV[i].saveToFile("./dataSets/nist_CV_split" + i + ".ri");
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

}
