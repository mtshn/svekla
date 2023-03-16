package ru.ac.phyche.gcms.svekla;

import java.awt.Font;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import javax.vecmath.Point2d;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.IImplementationSpecification;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.KlekotaRothFingerprinter;
import org.openscience.cdk.fingerprint.LingoFingerprinter;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.geometry.GeometryUtil;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.qsar.DescriptorEngine;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.IDescriptor;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import net.sf.jniinchi.INCHI_RET;

/**
 * 
 * A wrapper class of Chemistry Development Kit (CDK). Basic chemoinformatics
 * routines. CDK 2.3 is supported only; Note, that pretrained models are
 * extremely sensitive to version of CDK and any changes in this class (because
 * "canonical" order of atoms in a SMILES string depends on CDK version, as well
 * as 2D-layout and descriptor values in some cases).
 *
 */
public class Chemoinformatics {
	/**
	 * Timeout for computation of a molecular descriptor. Seconds. 3600 = 1 hour =
	 * effectively NO timeout at all.
	 */
	public static final int DESCRIPTOR_TIMEOUT = 3600; // SECONDS
	/**
	 * Max length of SMILES string. Longer strings are not supported by 1D CNN model
	 */
	public static final int SMILES_LEN = 250;
	/**
	 * A depiction - square image with dimensions DEPICTION_SIZE*DEPICTION_SIZE
	 */
	public static final int DEPICTION_SIZE = 400;
	/**
	 * Possible letters in SMILES string (important for 1D CNN, other symbols are
	 * not supported). Space must be first. Any change here will require
	 * (immediately) retraining of 1D-CNN model.
	 */
	public static final char[] TOKENS = { ' ', 'C', 'c', 'N', 'n', 'H', 'O', 'o', 'F', 'B', 'l', 'r', 'S', 'i', '+',
			'(', ')', '[', ']', '-', '=', '#', '1', '2', '3', '4', '5', '6', '7', '8', '9', 's', 'P', '%', 'I', 's' };
	// The space always comes first! The space has index 0!

	/**
	 * =TOKENS.length;
	 */
	public static final int SMILES_TOKENS = TOKENS.length;

	/**
	 * Supported types of molecular fingerprints; ADDITIVE_CIRCULAR fingerprints can
	 * be not only 0 or 1 but also bigger integers. These fingerprints are similar
	 * to usual circular fingerprints but enumerate occurrences of substructures
	 * rather than denote if a substucture present or doesn't present. See work
	 * 10.1021/acscentsci.9b00085 for explanations on additive fingerprints.
	 */
	public static enum FingerprintsType {
		NONE, MACCS, CIRCULAR_4_1024, CIRCULAR_6_1024, LINGO, PUBCHEM, KLEKOTA, CIRCULAR_6_4096,
		ADDITIVE_CIRCULAR_4_1024_NO_SCALE, ADDITIVE_CIRCULAR_6_1024_NO_SCALE, CIRCULAR_4_512;
	}

	/**
	 * 
	 * @param smiles SMILES string for a molecule
	 * @return molecular weight of the molecule
	 * @throws CDKException CDK exception
	 */
	public static float weight(String smiles) throws CDKException {
		WeightDescriptor molarMass = new WeightDescriptor();
		return ((float) ((DoubleResult) molarMass.calculate(smilesToAtomContainer(smiles)).getValue()).doubleValue());
	}

	/**
	 * Converts SMILES string into integers for 1D CNN
	 * 
	 * @param s SMILES string for a molecule
	 * @return int[SMILES_LEN] - integer array, where each symbol of SMILES string
	 *         is represented as integer according to index of this symbol in TOKENS
	 *         array.
	 * @throws CDKException This method throws an exception if SMILES string of
	 *                      molecule is longer than SMILES_LEN or if the string
	 *                      contains incompatible symbols (i.e. symbols which are
	 *                      not contained in TOKENS array).
	 */
	public static int[] tokenize(String s) throws CDKException {
		int[] result = new int[SMILES_LEN];
		if (s.length() > SMILES_LEN) {
			throw (new CDKException(" Too long SMILES string " + s.length()));
		}
		for (int i = 0; i < s.length(); i++) {
			boolean found = false;
			for (int j = 0; j < TOKENS.length; j++) {
				if (TOKENS[j] == s.charAt(i)) {
					result[i] = j;
					found = true;
				}
			}
			if (!found) {
				throw (new CDKException(" Unknown symbol " + s.charAt(i)));

			}
		}
		return result;
	}

	/**
	 * Backward transform comparing with tokenize().
	 * intsToSmiles(tokenize(s)).trim().equals(s) is true
	 * 
	 * @param q int[SMILES_LEN] array which represents SMILES string. Symbols
	 *          encoded as integers according to TOKENS array.
	 * @return SMILES string (with spaces at end). Use "trim()" always with this
	 *         method.
	 */
	public static String intsToSmiles(int[] q) {
		String s = "";
		for (int i = 0; i < q.length; i++) {
			s = s + TOKENS[q[i]];
		}
		return s;
	}

	/**
	 * Compute InChI-key from SMILES string for a molecule. For some molecules it
	 * often works wrong...
	 * 
	 * @param smiles SMILES string
	 * @return InChI-key
	 * @throws CDKException if something goes wrong...
	 */
	public static String smilesToInchiKey(String smiles) throws CDKException {
		IAtomContainer mol = smilesToAtomContainer(smiles);
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIGenerator gen = factory.getInChIGenerator(mol);
		String inchi = gen.getInchiKey();
		return inchi;
	}

	/**
	 * Compute InChI from SMILES string for a molecule. For some molecules it often
	 * works wrong..
	 * 
	 * @param smiles SMILES string
	 * @return InChI
	 * @throws CDKException CDK exception
	 */
	public static String smilesToInchi(String smiles) throws CDKException {
		IAtomContainer mol = smilesToAtomContainer(smiles);
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIGenerator gen = factory.getInChIGenerator(mol);
		String inchi = gen.getInchi();
		return inchi;
	}

	/**
	 * Create SMILES string from InChI string. It creates SMILES strings with
	 * aromatics symbols, with implicit hydrogens.
	 * 
	 * @param inchi           InChI string
	 * @param stereochemistry if TRUE - cis/trans isomers and optical isomers will
	 *                        be denoted using special symbols. If false - cis/trans
	 *                        and optical isomers will be identical.
	 * @return SMILES string
	 * @throws CDKException CDK internal errors.
	 */
	public static String inchiToSmiles(String inchi, boolean stereochemistry) throws CDKException {
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIToStructure intostruct = factory.getInChIToStructure(inchi, DefaultChemObjectBuilder.getInstance());

		@SuppressWarnings("deprecation")
		INCHI_RET ret = intostruct.getReturnStatus();
		if ((ret != INCHI_RET.OKAY) && (ret != INCHI_RET.WARNING)) {
			throw (new CDKException("INCHI_RET status failed!"));
		}
		IAtomContainer mol = intostruct.getAtomContainer();
		return atomContainerlToSmiles(mol, stereochemistry);
	}

	/**
	 * Convert a SMILES string to a "canonical" SMILES string. For one molecule may
	 * be multiple (almost infinity) SMILES representations. Various atom orders,
	 * and various settings (e.g. use aromatics symbols or not). This methods allow
	 * to ensure that identical molecules will have identical SMILES! Internally it
	 * is inchiToSmiles(smilesToInchi(smiles), stereochemistry); It works because
	 * InChI is unequivocal.
	 * 
	 * 
	 * @param smiles          SMILES string
	 * @param stereochemistry if TRUE - cis/trans isomers and optical isomers will
	 *                        be denoted using special symbols. If false - cis/trans
	 *                        and optical isomers will be identical.
	 * @return canonical SMILES string, output will be identical for identical
	 *         molecules even if inputs are different.
	 * @throws CDKException CDK internal errors.
	 */
	public static String canonical(String smiles, boolean stereochemistry) throws CDKException {
		if ((smiles != null) && (!smiles.trim().equals(""))) {
			return inchiToSmiles(smilesToInchi(smiles), stereochemistry);
		} else {
			return null;
		}
	}

	/**
	 * Molecular fingerprints for a molecule.
	 * 
	 * @param smiles molecule (SMILES string)
	 * @param type   types of fingerprints
	 * @return float array. Typically zeros (0.0F) or ones (1.0F) or values which
	 *         are identical to integers.
	 * @throws CDKException CDK internal errors.
	 */
	public static float[] fingerprints(String smiles, FingerprintsType type) throws CDKException {
		if (type == FingerprintsType.NONE) {
			return new float[] {};
		}
		IAtomContainer mol = smilesToAtomContainer(smiles.trim());
		BitSet fp = null;
		if (type == FingerprintsType.ADDITIVE_CIRCULAR_4_1024_NO_SCALE) {
			return (circularAdditiveFingerPrints(mol, CircularFingerprinter.CLASS_ECFP4, 1024, false));
		}
		if (type == FingerprintsType.ADDITIVE_CIRCULAR_6_1024_NO_SCALE) {
			return (circularAdditiveFingerPrints(mol, CircularFingerprinter.CLASS_ECFP6, 1024, false));
		}
		if (type == FingerprintsType.MACCS) {
			MACCSFingerprinter fpGen = new MACCSFingerprinter();
			fp = fpGen.getBitFingerprint(mol).asBitSet();
		}
		if (type == FingerprintsType.CIRCULAR_4_1024) {
			CircularFingerprinter fpGen = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4, 1024);
			fpGen.calculate(mol);
			fp = fpGen.getBitFingerprint(mol).asBitSet();
		}
		if (type == FingerprintsType.CIRCULAR_4_512) {
			CircularFingerprinter fpGen = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4, 512);
			fpGen.calculate(mol);
			fp = fpGen.getBitFingerprint(mol).asBitSet();
		}
		if (type == FingerprintsType.CIRCULAR_6_1024) {
			CircularFingerprinter fpGen = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP6, 1024);
			fpGen.calculate(mol);
			fp = fpGen.getBitFingerprint(mol).asBitSet();
		}
		if (type == FingerprintsType.CIRCULAR_6_4096) {
			CircularFingerprinter fpGen = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP6, 4096);
			fpGen.calculate(mol);
			fp = fpGen.getBitFingerprint(mol).asBitSet();
		}
		if (type == FingerprintsType.LINGO) {
			LingoFingerprinter fpGen = new LingoFingerprinter(6);
			fp = fpGen.getBitFingerprint(mol).asBitSet();
		}

		if (type == FingerprintsType.PUBCHEM) {
			PubchemFingerprinter fpGen = new PubchemFingerprinter(mol.getBuilder());
			fp = fpGen.getBitFingerprint(mol).asBitSet();
		}

		if (type == FingerprintsType.KLEKOTA) {
			KlekotaRothFingerprinter fpGen = new KlekotaRothFingerprinter();
			fp = fpGen.getBitFingerprint(mol).asBitSet();
		}

		float[] result = new float[fp.size()];
		for (int i = 0; i < fp.size(); i++) {
			if (fp.get(i)) {
				result[i] = 1;
			} else {
				result[i] = 0;
			}
		}
		return result;
	}

	/**
	 * Input: 3 arrays of floats with identical length. Output[i] = (descriptors[j]
	 * - min[j]) / (max[j] - min[j])
	 * 
	 * @param descriptors array of descriptors for a molecule
	 * @param min         minimal values of each descriptor
	 * @param max         max values of each descriptor
	 * @return descriptors scaled to [0,1]
	 */
	public static float[] scaleMinMax(float[] descriptors, float[] min, float[] max) {
		float[] d = new float[descriptors.length];
		for (int j = 0; j < d.length; j++) {
			if ((((max[j]) - (min[j])) == 0) || ((Float) ((max[j]) - (min[j]))).isInfinite()) {
				d[j] = 0;
			} else {
				d[j] = (descriptors[j] - min[j]) / (max[j] - min[j]);
			}
		}
		return d;
	}

	/**
	 * Generates CDK descriptors (2D, 1D only, 3D descriptors are not supported)
	 * scaled to range [0,1].
	 * 
	 * @param smiles          molecule (SMILES string)
	 * @param descriptorNames CDK descriptors names (such as "fragC", "C1SP1",
	 *                        "C2SP1", "C1SP2" etc...).
	 * @param min             min values of each descriptor
	 * @param max             max values of each descriptor
	 * @return descriptors scaled to [0,1]
	 * @throws CDKException CDK internal errors, incorrect SMILES etc.
	 */
	public static float[] descriptors(String smiles, String[] descriptorNames, float[] min, float[] max)
			throws CDKException {
		return scaleMinMax(descriptors(smiles, descriptorNames), min, max);
	}

	/**
	 * Returns all names of descriptors which are supported by CDK (including 3D
	 * descriptors).
	 * 
	 * @return names of descriptors
	 */
	public static String[] allDescriptorNames() {
		List<IDescriptor> descriptors = descriptorList;
		ArrayList<String> result = new ArrayList<String>();
		for (IDescriptor d : descriptors) {
			IMolecularDescriptor md = (IMolecularDescriptor) d;
			result.addAll(Arrays.asList(md.getDescriptorNames()));
		}
		return result.toArray(new String[result.size()]);
	}

	/**
	 * Generates CDK descriptors. 3D descriptors (i.e. descriptors which require 3D
	 * coordinates) are not supported. You can set timeout because some descriptors
	 * for some molecules are calculated too slow.
	 * 
	 * @param smiles          molecule (SMILES string)
	 * @param descriptorNames CDK descriptors names (such as "fragC", "C1SP1",
	 *                        "C2SP1", "C1SP2" etc...).
	 * @return float[descriptorNames.length] array with molecular descriptors
	 * @throws CDKException CDK internal errors, incorrect SMILES etc.
	 */
	public static float[] descriptors(String smiles, String[] descriptorNames) throws CDKException {
		HashSet<String> descriptorNamesSet = new HashSet<String>(Arrays.asList(descriptorNames));
		IAtomContainer mol = smilesToAtomContainer(smiles);

		HashMap<String, Float> resultMap = new HashMap<String, Float>();
		List<IDescriptor> descriptors = descriptorList;

		for (IDescriptor d : descriptors) {
			IMolecularDescriptor md = (IMolecularDescriptor) d;
			HashSet<String> thisDescriptorNamesSet = new HashSet<String>(Arrays.asList(md.getDescriptorNames()));
			thisDescriptorNamesSet.retainAll(descriptorNamesSet);
			if (thisDescriptorNamesSet.size() > 0) {
				boolean failedToCalculate = false;
				DescriptorValue dv = null;
				try {
					try {
						dv = computeDescriptor(mol, md);
					} catch (StackOverflowError | ArrayIndexOutOfBoundsException e) {
						failedToCalculate = true;
					}
				} catch (NullPointerException e) {
					failedToCalculate = true;
					if ((md == null) || (mol == null)) {
						throw e;
					}
				}
				if (dv == null) {
					failedToCalculate = true;
				}
				String[] values = null;
				if (!failedToCalculate) {
					IDescriptorResult dr = dv.getValue();
					values = dr.toString().split(",");
				}
				String[] names = md.getDescriptorNames();

				for (int j = 0; j < names.length; j++) {
					float value = Float.NaN;
					if (!failedToCalculate) {
						value = Float.parseFloat(values[j]);
					}
					resultMap.put(names[j], value);
				}
			}
		}
		float[] result = new float[descriptorNames.length];

		for (int k = 0; k < descriptorNames.length; k++) {
			if ((resultMap.get(descriptorNames[k]) == null)) {
				result[k] = Float.NaN;
			} else {
				float res = resultMap.get(descriptorNames[k]);
				result[k] = Float.isInfinite(res) ? Float.NaN : res;
			}
		}
		return result;
	}

	/**
	 * Create image (depiction of chemical structure) of molecule. Creates 2-color
	 * image with size DEPICTION_SIZE*DEPICTION_SIZE, each pixel is 0.0F and 1.0F.
	 * Blank space (background) 0.0F.
	 * 
	 * @param smiles molecule (SMILES string)
	 * @return float[DEPICTION_SIZE][DEPICTION_SIZE] - depiction
	 * @throws CDKException CDK internal errors
	 */
	public static float[][] smilesToImage(String smiles) throws CDKException {
		BufferedImage img = smilesToImageAsBufferedImage(smiles);
		float a[][] = new float[DEPICTION_SIZE][DEPICTION_SIZE];
		final byte[] pixels = ((DataBufferByte) img.getRaster().getDataBuffer()).getData();

		int x = 0;
		int y = 0;
		for (int i = 1; i + 4 < pixels.length; i += 4) {
			a[x][y] = (pixels[i] == -1) ? 0.0F : 1.0F;
			x++;
			if (x == DEPICTION_SIZE) {
				x = 0;
				y++;
			}
		}
		return a;
	}

	/**
	 * Create image (depiction of chemical structure) of molecule. Creates grayscale
	 * image.
	 * 
	 * @param smiles molecule (SMILES string)
	 * @return image with dimensions DEPICTION_SIZE*DEPICTION_SIZE
	 * @throws CDKException CDK internal errors
	 */
	public static BufferedImage smilesToImageAsBufferedImage(String smiles) throws CDKException {
		IAtomContainer mol = smilesToAtomContainer(smiles);
		DepictionGenerator dg = new DepictionGenerator(new Font(Font.MONOSPACED, Font.PLAIN, 14))
				.withSize(DEPICTION_SIZE, DEPICTION_SIZE).withAromaticDisplay();
		BufferedImage img0 = dg.depict(mol).toImg();
		return img0;
	}

	/**
	 * Count number of functional groups. Groups enumerated according with work
	 * 10.1021/ci600548y, Table 2. Each element of returned array means number of
	 * functional groups according to Table 2 of that work. 0.0F - the group
	 * absents, 1 - presents once, 2 - presents twice etc. NOTE!!! Our understanding
	 * of "groups" doesn't imply that each atom is part of only one group. For
	 * example -C(=O)-NH- fragment gives 1 to group 33, 27, 39. Internally we just
	 * enumerate how many atoms corresponds to each of set of SMARTS queries. List
	 * of SAMRT queries: { "[CX4H3]", "[CX4H2!R]", "[CX4H2R]", "[CX4H1!R]",
	 * "[CX4H1R]", "[CX4H0!R]", "[CX4H0R]", "[CX3H2]", "[CX3H1!R]", "[CX3H1R]",
	 * "[CX3H0!R]", "[CX3H0R]", "*:[cX3H1]:*", "*:[cX3H0](-*):*", "*:[cX3H0](:*):*",
	 * "[CX2H1]", "[CX2H0]", "[OX2H1]", "[CX4H2]-[OX2H1]", "[CX4H1]-[OX2H1]",
	 * "[CX4H0]-[OX2H1]", "a-[OX2H1]", "[OX2H0!R]", "[OX2H0R,oX2H0R]", "O-[OX2H1]",
	 * "[CX3H1!R]=O", "[CX3H0!R]=O", "[CX3H0R]=O", "[O!H1]-[CX3H0!R]=O",
	 * "[O]-[CX3H0R]=O", "[OH1]-[CX3H0!R]=O", "[NH2]-[CX3H0!R]=O",
	 * "[NX3H1]-[CX3H0!R]=O", "[NX3H1]-[CX3H0R]=O", "[NX3H0]-[CX3H0!R]=O",
	 * "[NX3H0]-[CX3H0R]=O", "A-[NH2]", "a-[NH2]", "[NX3H1!R]", "[NX3H1R,nX3H1R]",
	 * "[NX3H0!R]", "[NX3H0R,nX3H0R]", "[NX3H0!R]-[NX2]=[OX1]",
	 * "[NX3H0R,nX3H0R]-[NX2]=[OX1]", "*:[nX2H0]:*", "[NX2H1]", "[NX2H0!R]",
	 * "[NX2H0R]", "[NX2H0R]-[NX3H1R]", "[NX2H0R]=[CX3H0R]-[NX3H0R]",
	 * "[NX2H0R]=[CX3H0R]-[NX3H1R]", "[NX2H0R]=[CX3H1R]-[NX3H0R]",
	 * "[NX2H0R]=[CX3H1R]-[NX3H1R]", "[NX2H0]=[NX2H0]", "[NX2]=C=[SX1]",
	 * "[NX2]=[OX1]", "[OX1][NX3][OX1]", "A-C#[NX1]", "a-C#[NX1]", "A-[F]", "a-[F]",
	 * "[Cl]", "[CX4H2]-[Cl]", "[CX4H1]-[Cl]", "[CX4H0]-[Cl]", "a-[Cl]", "A-[Br]",
	 * "a-[Br]", "A-[I]", "a-[I]", "[SX2H1]", "a-[SX2H1]", "[SX2H0!R]", "[SX2H0R]",
	 * "[SX3][OX1]", "[OX1][SX4][OX1]", "[CX3!R]=[SX1]", "[CX3R]=[SX1]", "[SiX4H1]",
	 * "[SiX4!R]", "[SiX4H0R]", "[PX3]", "[PX4]=O", "[PX4]=S" }
	 * 
	 * @param smiles molecule (SMILES string)
	 * @return an array with numbers of occurrences of each of listed above SMARTS
	 *         queries in the molecule.
	 * @throws CDKException internal CDK
	 */
	public static float[] funcGroups(String smiles) throws CDKException {
		IAtomContainer mol = smilesToAtomContainer(smiles);
		String[] patterns = smartsPatterns();
		float[] result = new float[patterns.length + 1];

		for (int i = 0; i < patterns.length; i++) {
			result[i] = smartsCount(mol, patterns[i]);
		}
		return result;
	}

	/**
	 * 2D representation of molecule for 2D-CNN. At fist stage 2D-coordinates of
	 * atoms are created (the same as for common depiction). Coordinates of centers
	 * of all bonds are calculated too. All this coordinates are rounded to 0.5
	 * units. After it we have a table with size 130*130, where each cell can
	 * contain atom, bond or nothing. Each cell corresponds square 0.5 unit * 0.5
	 * unit on "depiction" plane. Then each cell is one-hot encoded. 26 types of
	 * atoms are considered: "B", "C" (other), "F", "H" (explicit, rare case), "Cl",
	 * "I", "S.planar3", "C.sp2", "C.sp3", "N"(other), "O"(other), "P", "N.nitro",
	 * "S"(other), "X", "S.3", "O.planar3", "N.planar3", "Br", "N.amide", "N.sp2",
	 * "Si", "N.sp3", "C.sp", "O.sp3", "O.sp2". And three types of bonds: single,
	 * double, triple. Finally we have array [26+3][130][130]. 130*130 - spatial
	 * dimensions (65 units*65 units). 29 - number of channels for CNN. Molecule in
	 * CENTERED in the area (center in cell 65:65).
	 * 
	 * @param smiles molecule (SMILES string)
	 * @return [29][130][130] input for 2D-CNN
	 * @throws CDKException CDK internal error, oversized molecules (2D coordinates
	 *                      don't fit square 65*65 units (molecules with chain
	 *                      longer than 50, usually).
	 */
	public static float[][][] representation2d(String smiles) throws CDKException {
		IAtomContainer mol = calculate2dcoordinates(smilesToAtomContainer(smiles));
		float[][][] result = new float[29][130][130];
		for (IAtom atom : mol.atoms()) {
			int[] coordinatesAndType = atomTypeAndCoordinatesToInts(atom);
			try {
				result[coordinatesAndType[0]][coordinatesAndType[1]][coordinatesAndType[2]] += 1;
			} catch (ArrayIndexOutOfBoundsException e) {
				throw new CDKException("Oversized 2D-representation of molecule " + smiles + " " + e.getMessage());
			}
		}
		for (IBond bond : mol.bonds()) {
			int[] coordinatesAndType = bondTypeAndCoordinatesToInts(bond);
			try {
				result[26 + coordinatesAndType[0]][coordinatesAndType[1]][coordinatesAndType[2]] += 1;
			} catch (ArrayIndexOutOfBoundsException e) {
				throw new CDKException("Oversized 2D-representation of molecule " + smiles + " " + e.getMessage());
			}
		}

		return result;
	}

	/**
	 * 
	 * @param smiles molecule (SMILES string)
	 * @return number of -OH groups (including -O-OH, COOH and so on)
	 * @throws CDKException CDK internal error.
	 */
	public static int countOH(String smiles) throws CDKException {
		IAtomContainer mol = smilesToAtomContainer(smiles);
		int count = 0;
		for (IAtom a : mol.atoms()) {
			if (isAtomOH(mol, a)) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Replace all -OH groups (including -O-OH, COOH and so on) with -OSi(CH3)3
	 * groups (i.e. in silico derivatization).
	 * 
	 * @param smiles          molecule (SMILES string)
	 * @param stereochemistry use or not symbols for cis/trans isomers and optical
	 *                        isomers (see canonical(String, boolean) method.
	 * @return derivatized molecule
	 * @throws CDKException CDK internal error.
	 */
	public static String replaceOHWithOSiCH33OH(String smiles, boolean stereochemistry) throws CDKException {
		IAtomContainer mol = smilesToAtomContainer(smiles);
		for (IAtom a : mol.atoms()) {
			if (isAtomOH(mol, a)) {
				replaceOHWithOSiCH33(mol, a);
			}
		}
		return atomContainerlToSmiles(mol, stereochemistry);
	}

	// private methods
	private static final List<IDescriptor> descriptorList = getDescriptorList();

	private static List<IDescriptor> getDescriptorList() {
		List<String> classes = DescriptorEngine
				.getDescriptorClassNameByPackage("org.openscience.cdk.qsar.descriptors.molecular", null);
		DescriptorEngine engine = new DescriptorEngine(classes, null);
		List<IDescriptor> inst = engine.instantiateDescriptors(classes);
		List<IImplementationSpecification> specs = engine.initializeSpecifications(inst);
		engine.setDescriptorInstances(inst);
		engine.setDescriptorSpecifications(specs);
		return engine.getDescriptorInstances();
	}

	private static IAtomContainer smilesToAtomContainer(String s) throws CDKException {
		SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IAtomContainer mol = parser.parseSmiles(s.trim());
		Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		arom.apply(mol);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
		adder.addImplicitHydrogens(mol);
		return mol;
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

	private static float[] circularAdditiveFingerPrints(IAtomContainer mol, int type, int len, boolean scale)
			throws CDKException {
		CircularFingerprinter cf = new CircularFingerprinter(type, len);
		cf.calculate(mol);
		float[] result = new float[len];
		for (int n = 0; n < len; n++) {
			result[n] = 0;
		}
		for (int n = 0; n < cf.getFPCount(); n++) {
			int i = cf.getFP(n).hashCode;
			long b = i >= 0 ? i : ((i & 0x7FFFFFFF) | (1L << 31));
			result[(int) (b % len)] += 1.0;
		}
		if (scale) {
			float max = -1;
			for (int n = 0; n < len; n++) {
				if (result[n] > max) {
					max = result[n];
				}
			}
			if (max != 0) {
				for (int n = 0; n < len; n++) {
					result[n] = result[n] / max;
				}
			} else
				System.out.println("Warning! All fingerprints are zero... " + atomContainerlToSmiles(mol, true));
		}
		return result;
	}

	private static int smartsCount(IAtomContainer mol, String smarts) throws CDKException {
		Pattern querytool = SmartsPattern.create(smarts);
		if (querytool.matches(mol)) {
			int n = querytool.matchAll(mol).uniqueAtoms().countUnique();
			return n;
		} else {
			return 0;
		}

	}

	private static String[] smartsPatterns() {
		return (new String[] { "[CX4H3]", "[CX4H2!R]", "[CX4H2R]", "[CX4H1!R]", "[CX4H1R]", "[CX4H0!R]", "[CX4H0R]",
				"[CX3H2]", "[CX3H1!R]", "[CX3H1R]", "[CX3H0!R]", "[CX3H0R]", "*:[cX3H1]:*", "*:[cX3H0](-*):*",
				"*:[cX3H0](:*):*", "[CX2H1]", "[CX2H0]", "[OX2H1]", "[CX4H2]-[OX2H1]", "[CX4H1]-[OX2H1]",
				"[CX4H0]-[OX2H1]", "a-[OX2H1]", "[OX2H0!R]", "[OX2H0R,oX2H0R]", "O-[OX2H1]", "[CX3H1!R]=O",
				"[CX3H0!R]=O", "[CX3H0R]=O", "[O!H1]-[CX3H0!R]=O", "[O]-[CX3H0R]=O", "[OH1]-[CX3H0!R]=O",
				"[NH2]-[CX3H0!R]=O", "[NX3H1]-[CX3H0!R]=O", "[NX3H1]-[CX3H0R]=O", "[NX3H0]-[CX3H0!R]=O",
				"[NX3H0]-[CX3H0R]=O", "A-[NH2]", "a-[NH2]", "[NX3H1!R]", "[NX3H1R,nX3H1R]", "[NX3H0!R]",
				"[NX3H0R,nX3H0R]", "[NX3H0!R]-[NX2]=[OX1]", "[NX3H0R,nX3H0R]-[NX2]=[OX1]", "*:[nX2H0]:*", "[NX2H1]",
				"[NX2H0!R]", "[NX2H0R]", "[NX2H0R]-[NX3H1R]", "[NX2H0R]=[CX3H0R]-[NX3H0R]",
				"[NX2H0R]=[CX3H0R]-[NX3H1R]", "[NX2H0R]=[CX3H1R]-[NX3H0R]", "[NX2H0R]=[CX3H1R]-[NX3H1R]",
				"[NX2H0]=[NX2H0]", "[NX2]=C=[SX1]", "[NX2]=[OX1]", "[OX1][NX3][OX1]", "A-C#[NX1]", "a-C#[NX1]", "A-[F]",
				"a-[F]", "[Cl]", "[CX4H2]-[Cl]", "[CX4H1]-[Cl]", "[CX4H0]-[Cl]", "a-[Cl]", "A-[Br]", "a-[Br]", "A-[I]",
				"a-[I]", "[SX2H1]", "a-[SX2H1]", "[SX2H0!R]", "[SX2H0R]", "[SX3][OX1]", "[OX1][SX4][OX1]",
				"[CX3!R]=[SX1]", "[CX3R]=[SX1]", "[SiX4H1]", "[SiX4!R]", "[SiX4H0R]", "[PX3]", "[PX4]=O", "[PX4]=S" });
	}

	private static final ExecutorService exec = Executors.newCachedThreadPool();

	private static DescriptorValue computeDescriptor(IAtomContainer mol, IMolecularDescriptor md) {
		DecriptorComputeCall computer = new DecriptorComputeCall();
		computer.md = md;
		computer.mol = mol;

		Future<DescriptorValue> future = exec.submit(computer);
		DescriptorValue dv = null;
		try {
			dv = future.get(1, TimeUnit.SECONDS);
		} catch (InterruptedException | ExecutionException | TimeoutException e) {

		} finally {
			future.cancel(true);
		}
		return dv;
	}

	private static class DecriptorComputeCall implements Callable<DescriptorValue> {
		public IAtomContainer mol;
		public IMolecularDescriptor md;

		@Override
		public DescriptorValue call() throws Exception {
			return md.calculate(mol);
		}
	}

	private static String atomTypeToString(IAtom at) {
		HashSet<String> mostCommonAtomTypes = new HashSet<String>();
		mostCommonAtomTypes.addAll(Arrays.asList(new String[] { "C.sp3", "C.sp2", "O.sp3", "O.sp2", "N.sp2", "N.amide",
				"N.sp3", "C.sp", "S.3", "O.planar3", "N.nitro", "N.planar3", "S.planar3" }));
		String atomType = at.getAtomTypeName();
		if (mostCommonAtomTypes.contains(atomType.trim())) {
			return atomType.trim();
		} else {
			return atomType.split("\\.")[0];
		}
	}

	private static int[] atomTypeAndCoordinatesToInts(IAtom at) throws CDKException {
		int[] result = new int[3];
		String[] atomTypes = new String[] { "B", "C", "F", "H", "Cl", "I", "S.planar3", "C.sp2", "C.sp3", "N", "O", "P",
				"N.nitro", "S", "X", "S.3", "O.planar3", "N.planar3", "Br", "N.amide", "N.sp2", "Si", "N.sp3", "C.sp",
				"O.sp3", "O.sp2" };
		String atomType = atomTypeToString(at);
		result[0] = -1;
		for (int i = 0; i < atomTypes.length; i++) {
			if (atomTypes[i].equals(atomType)) {
				result[0] = i;
			}
		}
		if (result[0] == -1) {
			throw new CDKException("Unexpected atom type " + at.getAtomTypeName());
		}
		Point2d xy = at.getPoint2d();
		result[1] = Math.round((float) Math.floor(2 * xy.x));
		result[2] = Math.round((float) Math.floor(2 * xy.y));
		return result;
	}

	private static int[] bondTypeAndCoordinatesToInts(IBond b) throws CDKException {
		int[] result = new int[3];
		String order = b.getOrder().name();
		switch (order) {
		case "SINGLE":
			result[0] = 0;
			break;
		case "DOUBLE":
			result[0] = 1;
			break;
		case "TRIPLE":
			result[0] = 2;
			break;
		}
		Point2d xy = b.get2DCenter();
		result[1] = Math.round((float) Math.floor(2 * xy.x));
		result[2] = Math.round((float) Math.floor(2 * xy.y));
		return result;
	}

	private static IAtomContainer calculate2dcoordinates(IAtomContainer mol) throws CDKException {
		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(mol);
		sdg.generateCoordinates();
		IAtomContainer result = sdg.getMolecule();
		Point2d center = GeometryUtil.get2DCenter(result);
		GeometryUtil.translate2D(result, 32.5 - center.x, 32.5 - center.y);
		return result;
	}

	private static boolean isAtomOH(IAtomContainer mol, IAtom atom) {
		int nonHNeigbors = 0;
		int hNeigbors = 0;
		if ((atom.getAtomicNumber() == 8) && (atom.getFormalCharge() == 0)) {
			for (IAtom neighbor : mol.getConnectedAtomsList(atom)) {
				if (neighbor.getAtomicNumber() != 1) {
					nonHNeigbors++;
				} else {
					hNeigbors++;
				}
			}
		}
		return (nonHNeigbors == 1) && ((atom.getImplicitHydrogenCount() == 1) || (hNeigbors == 1));
	}

	private static void replaceOHWithOSiCH33(IAtomContainer mol, IAtom atom) {
		for (IAtom neighbor : mol.getConnectedAtomsList(atom)) {
			if (neighbor.getAtomicNumber() == 1) {
				mol.removeAtom(neighbor);
			}
		}
		atom.setImplicitHydrogenCount(0);
		IAtom at2 = new Atom("Si");
		mol.addAtom(at2);
		mol.addBond(new Bond(atom, at2));
		for (int y = 0; y < 3; y++) {
			IAtom at3 = new Atom("C");
			mol.addAtom(at3);
			mol.addBond(new Bond(at3, at2));
			at3.setImplicitHydrogenCount(3);
		}
	}

}
