package ru.ac.phyche.gcms.svekla;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

import org.openscience.cdk.exception.CDKException;

public class RetentionsOfMetabolites {
	private static final boolean useStereoChemistry = false;

	private static boolean compareInChIkeys(String inchiKey1, String inchiKey2) {
		if (useStereoChemistry) {
			return (inchiKey1.trim().equals(inchiKey2.trim()));
		} else {
			return (inchiKey1.trim().split("\\-")[0].equals(inchiKey2.trim().split("\\-")[0]));
		}
	}

	private static RetentionsEntry convertFiehnOrGMDEntry(String inchikey, String inchi, String column,
			String derivatizationType, String retentionIndex, FileWriter log) throws IOException {
		if (inchikey == null) {
			log.write("Information about InChI-key not found. Skip;" + inchikey + " " + inchi + "\n");
			return null;
		}
		if (inchi == null) {
			log.write("Information about InChI not found. Skip;" + inchikey + " " + inchi + "\n");
			return null;
		}
		if (column == null) {
			log.write("Information about column not found. Skip;" + inchikey + " " + inchi + "\n");
			return null;
		}
		if (retentionIndex == null) {
			log.write("Information about retention index not found. Skip;" + inchikey + " " + inchi + "\n");
			return null;
		}
		if (derivatizationType == null) {
			derivatizationType = "0 TMS";
		}
		if (derivatizationType.contains("MEOX")) {
			log.write("MEOX derivatization not supported. Skip;" + inchikey + " " + inchi + "\n");
			return null;
		}
		if ((derivatizationType.split("\\s+").length != 2) || (!derivatizationType.split("\\s+")[1].equals("TMS"))) {
			log.write("Unsupported derivatization. Skip;" + inchikey + " " + inchi + "\n");
			return null;
		}
		int nTMS = Integer.parseInt(derivatizationType.split("\\s+")[0]);
		String smiles = null;
		try {
			smiles = Chemoinformatics.canonical(Chemoinformatics.inchiToSmiles(inchi, useStereoChemistry),
					useStereoChemistry);
			if (!compareInChIkeys(Chemoinformatics.smilesToInchiKey(smiles), inchikey)) {
				throw new CDKException("InChI key computed from InChI doesn't match InChIkey given in file");
			}
			if (nTMS != 0) {
				if (Chemoinformatics.countOH(smiles) != nTMS) {
					throw new CDKException("Ambiguous derivatization");
				} else {
					smiles = Chemoinformatics.replaceOHWithOSiCH33OH(smiles, useStereoChemistry);
				}
				Chemoinformatics.tokenize(smiles);
				Chemoinformatics.representation2d(smiles);
			}
		} catch (CDKException e) {
			log.write(e.getMessage() + inchikey + " " + inchi + "\n");
			return null;
		}
		System.out.println(retentionIndex + " " + column + " " + derivatizationType);
		if (column.equals(
				"Restek Corporation Rtx-5Sil MS (30 m length x 0.25 mm internal diameter with 0.25 µm film made of 95% dimethyl/5%diphenylpolysiloxane)")) {
			// This column is used in FiehnLib library
			float fiehnRI = Float.parseFloat(retentionIndex);
			float x = fiehnRI;
			float kovacRI = (float) (2.333910E-21 * x * x * x * x - 6.623524E-15 * x * x * x + 7.889699E-09 * x * x
					- 1.878315E-03 * x + 1.181852E+3);

			return RetentionsEntry.instance(smiles, kovacRI, 35);
		}
		if (column.equals("5_%_Phenyl_methyl_siloxane")) {
			// GMD library parser store ""5_%_Phenyl_methyl_siloxane" as column name if all
			// are right
			float kovacRI = Float.parseFloat(retentionIndex);
			return RetentionsEntry.instance(smiles, kovacRI, 15);
		}
		log.write("Unsupported column" + inchikey + " " + inchi + "\n");
		return null;
	}

	private static RetentionsDataset convertFiehnLib(String fileName, String logFilename) throws IOException {
		BufferedReader r = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileName))));
		FileWriter log = new FileWriter(logFilename);
		String s = r.readLine();
		ArrayList<RetentionsEntry> result = new ArrayList<RetentionsEntry>();
		while (s != null) {
			String inchikey = null;
			String retentionIndex = null;
			while ((s != null) && (!s.contains("Comments:"))) {
				if (s.contains("InChIKey:")) {
					inchikey = s.split("\\s+")[1];
				}
				if (s.contains("Retention_index:")) {
					retentionIndex = s.split("\\s+")[1];
				}
				s = r.readLine();
			}
			String comments[] = s.split("\\\"\\s\\\"");
			String inchi = null;
			String derivatizationType = null;
			String column = null;

			for (int i = 0; i < comments.length; i++) {
				s = comments[i];
				if (s.split("\\=")[0].equals("InChI")) {
					inchi = s.split("\\=")[1] + '=' + s.split("\\=")[2];
				}
				if (s.split("\\=")[0].equals("derivatization type")) {
					derivatizationType = s.split("\\=")[1];
				}
				if (s.split("\\=")[0].equals("column")) {
					column = s.split("\\=")[1];
				}
			}
			RetentionsEntry entry = convertFiehnOrGMDEntry(inchikey, inchi, column, derivatizationType, retentionIndex,
					log);
			if (entry != null) {
				result.add(entry);
			}
			while ((s != null) && (!s.contains("Name:"))) {
				s = r.readLine();
			}
		}
		r.close();
		log.close();
		return RetentionsDataset.create(result);
	}

	private static RetentionsDataset convertGMDLib(String fileName, String logFilename) throws IOException {
		BufferedReader r = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileName))));
		FileWriter log = new FileWriter(logFilename);
		String s = r.readLine();
		ArrayList<RetentionsEntry> result = new ArrayList<RetentionsEntry>();
		while (s != null) {
			if (!s.contains("Name:")) {
				r.close();
				throw new IOException("No name fonud");
			}
			String name = s.toUpperCase();
			String inchikey = null;
			String retentionIndex = null;
			String inchi = null;
			String derivatizationType = null;
			String column = null;
			while ((s != null) && (!s.contains("Num Peaks:"))) {
				if (s.contains("MET_INCHIKEY:")) {
					inchikey = s.split("\\s+")[1];
				}
				if (s.contains("MET_INCHI:")) {
					inchi = s.split("\\s+")[1];
				}
				if (s.contains("RI:")) {
					retentionIndex = s.split("\\s+")[1];
				}
				if (s.contains("RI VAR5 ALK: TRUE")) {
					column = "5_%_Phenyl_methyl_siloxane";
				}
				s = r.readLine();
			}

			if (name.contains("MEOX")) {
				derivatizationType = "MEOX";
			} else {
				int nTms = 0;
				for (int q = 1; q < 51; q++) {
					if (name.contains(q + "TMS")) {
						nTms = q;
					}
				}
				derivatizationType = nTms + " " + "TMS";
			}

			RetentionsEntry entry = convertFiehnOrGMDEntry(inchikey, inchi, column, derivatizationType, retentionIndex,
					log);
			if (entry != null) {
				result.add(entry);
			}
			while ((s != null) && (!s.contains("Name:"))) {
				s = r.readLine();
			}
		}
		r.close();
		log.close();
		return RetentionsDataset.create(result);
	}

	public static void main(String[] args) {
		try {
			convertFiehnLib("./dataSets/metabolites/MoNA-export-RTX5_Fiehnlib.msp",
					"./dataSets/metabolites/logs/fiehn.log").saveToFile("./dataSets/fiehnlib.ri");
			convertGMDLib("./dataSets/metabolites/GMD_20111121_VAR5_ALK_MSL.txt", "./dataSets/metabolites/logs/gmd.log")
					.saveToFile("./dataSets/gmd.ri");
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
