package ru.ac.phyche.gcms.svekla.javafxgui;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Scanner;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import ru.ac.phyche.gcms.svekla.javafxgui.MSpectra.SpectrumMS;

public class MSMSCFM {
	public static String sp(String smiles, boolean pos, boolean largest) {
		try {
			String os = System.getProperty("os.name").toLowerCase();

			if (os.contains("nix") || os.contains("nux")) {
				FileWriter sh = new FileWriter("cfm.sh");
				sh.write("export LD_LIBRARY_PATH=./rdkit;./cfm-predict \"");
				sh.write(smiles + "\" ");
				sh.write(" 0.001 ./models_msms/" + (pos ? "pos" : "neg") + "_se/param_output0.log ./models_msms/"
						+ (pos ? "pos" : "neg") + "_se/param_config.txt 1 tmp.txt ");
				sh.write((largest ? "1" : "0") + "\n");
				sh.close();
				ProcessBuilder b = new ProcessBuilder("sh", "cfm.sh");
				b.directory(new File("./"));
				Process p = b.start();
				p.waitFor();
			}
			if (os.contains("win")) {
				ProcessBuilder b = new ProcessBuilder("cfm-predict.exe", smiles, "0.001",
						"./models_msms/" + (pos ? "pos" : "neg") + "_se/param_output0.log",
						"./models_msms/" + (pos ? "pos" : "neg") + "_se/param_config.txt", "1", "tmp.txt",
						largest ? "1" : "0");
				b.directory(new File("./"));
				Process p = b.start();
				p.waitFor();
			}

			boolean t = false;
			BufferedReader br = new BufferedReader(new FileReader("tmp.txt"));
			String s = br.readLine();
			String result = pos ? "positive ions\n\n" : "negative ions\n\n";
			while (s != null) {
				System.out.println(s);
				String[] sp = s.trim().split("\\s+");
				if (sp.length == 1) {
					if (sp[0].trim().equals("energy0")) {
						t = true;
					}
				}
				if (t) {
					result = result + s;
					if (sp.length == 3) {
						String molForm = "";
						try {
							SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
							IAtomContainer mol = parser.parseSmiles(sp[2].trim());
							molForm = MolecularFormulaManipulator
									.getString(MolecularFormulaManipulator.getMolecularFormula(mol));
						} catch (Exception e) {

						}
						result = result + " " + molForm;
					}
					result = result + "\n";
				}
				s = br.readLine();
			}
			result = result + "\n";
			br.close();
			if (!t) {
				throw (new Exception("CFM ERROR"));
			}
			 Files.deleteIfExists(Paths.get("./tmp.txt"));
			 Files.deleteIfExists(Paths.get("./cfm.sh"));
			return result;
		} catch (Exception e) {
			e.printStackTrace();
			return "";
		}
	}

}
