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
import javafx.scene.text.Text;
import ru.ac.phyche.gcms.svekla.Chemoinformatics;
import ru.ac.phyche.gcms.svekla.Descriptors;
import ru.ac.phyche.gcms.svekla.Model;
import ru.ac.phyche.gcms.svekla.StackingLinearMetaLearnerModel;
import ru.ac.phyche.gcms.svekla.TrainPolar;
import ru.ac.phyche.gcms.svekla.TrainPolar.CNN;
import ru.ac.phyche.gcms.svekla.TrainPolar.MLP;

public class JSJavaCall {

	public String smiles = "X";
	public String oldSmiles = "X";

	public void jsjavacall(String smiles) {
		System.out.println(smiles);
		String canonicalSmiles;
		try {
			canonicalSmiles = Chemoinformatics.canonical(smiles, false);
			if (Chemoinformatics.weight(canonicalSmiles) < 60) {
				throw new CDKException("");
			}
			if (canonicalSmiles.length() < 5) {
				throw new CDKException("");
			}
			this.smiles = smiles;
		} catch (CDKException e) {
			e.printStackTrace();
		}

	}
}
