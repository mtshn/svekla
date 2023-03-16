package ru.ac.phyche.gcms.svekla;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.lang3.tuple.Pair;
import org.openscience.cdk.exception.CDKException;

/**
 * A molecular descriptors computation manager. This class calls method
 * descriptors from class Chemoinformatics. It allows to cache (store)
 * precomputed decriptors values for SMILES strings, save and load them, find
 * min and max value of a descriptor for a set of SMILES strings.
 *
 */
public class Descriptors {

	private HashMap<String, float[]> precomputed = null;
	private boolean usePrecomputed = true;
	private String[] descriptorsSet = null;
	private float[] min = null;
	private float[] max = null;

	/**
	 * If the descriptor manager is set to compute N descriptors it store arrays min
	 * and max (both float[N]) which contains minimal and maximal values of each
	 * descriptor for some compounds. i-th descriptor can be scaled to [0,1] as
	 * (value-min[i])/(max[i]-min[i])
	 * 
	 * @return pair of min array and max array
	 */
	public Pair<float[], float[]> getMinMaxArray() {
		return Pair.of(this.min, this.max);
	}

	/**
	 * See method get for more detailed explanations
	 * 
	 * @param smiles SMILES string (a compound)
	 * @return descriptors array for the compound, all NaNs are replaced with zeros.
	 * @throws CDKException CDK
	 */
	public float[] getNoNaNs(String smiles) throws CDKException {
		return nansToZero(get(smiles));
	}

	/**
	 * If usePrecomputed is true (this flag can be set while creation of an
	 * instance) it will return precomputed values for this SMILES (which are stored
	 * in the instance). In no precomputed values are available for this SMILES - an
	 * exception will be thrown. If usePrecomputed is false (this flag can be set
	 * while creation of an instance) it will call
	 * Chemoinformatics.descriptors(Chemoinformatics.canonical(smiles, true),
	 * descriptorsSet, min, max) method. See instance method, getMinMaxArray method
	 * for further details. Returned value = (CDK value-min[i])/(max[i]-min[i])
	 * 
	 * @param smiles SMILES string (a compound)
	 * @return descriptors array for the compound.
	 * @throws CDKException CDK error, absence of precomputed value if usage of
	 *                      precomputed descriptors is enabled in this instance.
	 */
	public float[] get(String smiles) throws CDKException {
		if (!usePrecomputed) {
			return Chemoinformatics.descriptors(Chemoinformatics.canonical(smiles, true), descriptorsSet, min, max);
		} else {
			if (precomputed == null) {
				throw new CDKException("No precomputed descriptors for this SMILES string!(1)");
			}
			float[] result = precomputed.get(smiles);
			if (result == null) {
				result = precomputed.get(Chemoinformatics.canonical(smiles, true));
				if (result == null) {
					throw new CDKException("No precomputed descriptors for this SMILES string!(2)");
				}
			}
			return result;
		}
	}

	/**
	 * Precomputing of descriptors. Compute and store descriptors for set of SMILES
	 * strings. This method computes descriptors in parallel, using parallelStream.
	 * usePrecomputed must be true (this flag can be set while creation of an
	 * instance).
	 * 
	 * @param compounds         HashSet of SMILES strings.
	 * @param recalculateMinMax if TRUE - array min will be filled with minimal
	 *                          values of each descriptor for this set of SMILES
	 *                          string and array max will be filled with maximal
	 *                          values of each descriptor. If FALSE - min and max
	 *                          arrays will be unchanged.
	 * @throws CDKException CDK error, thread pool error, usePrecomputed is false.
	 */
	public void precompute(HashSet<String> compounds, boolean recalculateMinMax) throws CDKException {
		if (!usePrecomputed) {
			throw new CDKException("Precomputing is disabled");
		}
		HashSet<String> canonicalCompounds = new HashSet<String>();
		for (String smiles : compounds) {
			canonicalCompounds.add(Chemoinformatics.canonical(smiles, true));
		}
		HashSet<DescriptorsOfOneCompound> descriptors = new HashSet<DescriptorsOfOneCompound>();
		for (String smiles : canonicalCompounds) {
			DescriptorsOfOneCompound d = new DescriptorsOfOneCompound();
			d.smiles = smiles;
			descriptors.add(d);
		}
		try {
			AtomicInteger i = new AtomicInteger(0);
			descriptors.parallelStream().forEach(d -> {
				try {
					i.incrementAndGet();
					if (i.get() % 1000 == 0) {
						System.out.println("Computing descriptors... " + i);
					}
					if (!recalculateMinMax) {
						d.compute(descriptorsSet, min, max);
					} else {
						d.compute(descriptorsSet);
					}
				} catch (CDKException e) {
					throw new RuntimeException(e.getMessage());
				}
			});
		} catch (Throwable e) {
			throw (new CDKException(e.getMessage()));
		}

		if (recalculateMinMax) {
			min = new float[descriptorsSet.length];
			max = new float[descriptorsSet.length];
			for (int i = 0; i < descriptorsSet.length; i++) {
				min[i] = Float.POSITIVE_INFINITY;
				max[i] = Float.NEGATIVE_INFINITY;
			}
			for (DescriptorsOfOneCompound d : descriptors) {
				for (int i = 0; i < descriptorsSet.length; i++) {
					if (d.d[i] < min[i]) {
						min[i] = d.d[i];
					}
					if (d.d[i] > max[i]) {
						max[i] = d.d[i];
					}
				}
			}
			for (int i = 0; i < descriptorsSet.length; i++) {
				if ((min[i] == Float.POSITIVE_INFINITY) || (max[i] == Float.NEGATIVE_INFINITY)) {
					min[i] = Float.NaN;
					max[i] = Float.NaN;
				}
			}
			for (DescriptorsOfOneCompound d : descriptors) {
				this.precomputed.put(d.smiles, Chemoinformatics.scaleMinMax(d.d, min, max));
			}
		} else {
			for (DescriptorsOfOneCompound d : descriptors) {
				this.precomputed.put(d.smiles, d.d);
			}
		}
	}

	/**
	 * Create instance of this class. Always use this method, instead any
	 * constructors!!
	 * 
	 * @param descriptorsSet_ names CDK descriptors names (such as "fragC", "C1SP1",
	 *                        "C2SP1", "C1SP2" etc...).
	 * @param min_            minimal value of each descriptor which is used to
	 *                        scale descriptor values to [0,1]. Can be new
	 *                        float[descriptorsSet_.length], if will be computed
	 *                        later. min_.length must be equal to
	 *                        descriptorsSet_.length.
	 * @param max_            maximal value of each descriptor which is used to
	 *                        scale descriptor values to [0,1]. Can be new
	 *                        float[descriptorsSet_.length], if will be computed
	 *                        later. max_.length must be equal to
	 *                        descriptorsSet_.length.
	 * @param usePrecomputed_ If TRUE - precomputed values will be used. In this
	 *                        case before any call of get(..) or getNoNaNs(...) for
	 *                        any SMILES string must be the call of precompute
	 *                        method for that string. If FALSE - precompute method
	 *                        is disabled (will return exception), get(...) and
	 *                        getNoNaNs(...) can be used directly.
	 * @return an instance of this class
	 */
	public static Descriptors instance(String[] descriptorsSet_, float[] min_, float[] max_, boolean usePrecomputed_) {
		if ((min_.length != descriptorsSet_.length) || (max_.length != descriptorsSet_.length)) {
			throw (new UnsupportedOperationException(
					"min, max arrays have to have the same length as descriptor names array"));
		}
		Descriptors result = new Descriptors();
		result.usePrecomputed = usePrecomputed_;
		result.descriptorsSet = descriptorsSet_;
		result.min = min_;
		result.max = max_;
		if (usePrecomputed_) {
			result.precomputed = new HashMap<String, float[]>();
		}
		return result;
	}

	/**
	 * Empty new instance (usePrecomputed_ = true)
	 * 
	 * @param descriptorsSet_ names CDK descriptors names (such as "fragC", "C1SP1",
	 *                        "C2SP1", "C1SP2" etc...).
	 * @return instance(descriptorsSet_, new float[descriptorsSet_.length], new
	 *         float[descriptorsSet_.length], true);
	 */

	public static Descriptors instance(String[] descriptorsSet_) {
		return instance(descriptorsSet_, new float[descriptorsSet_.length], new float[descriptorsSet_.length], true);
	}

	/**
	 * Save all content of this instance to file. File format. First line: number of
	 * descriptors, space, number of compounds that have precomputed descriptors,
	 * space, space-separated CDK names of descriptors, space separated minimal
	 * values for each descriptor, space separated maximal values for each
	 * descriptor. In the first line 3*(number of descriptors)+2 space separated
	 * values. Then one line per one compound: SMILES string, space, space separated
	 * values of descriptors. The values already are scaled to [0,1] as
	 * (value-min[i])/(max[i]-min[i])
	 * 
	 * 
	 * @param filename file name
	 * @throws IOException IO, only instance with precomputed descriptors can be
	 *                     saved to file, instead - exception will be thrown.
	 */
	public void saveToFile(String filename) throws IOException {
		if (!usePrecomputed) {
			throw (new UnsupportedOperationException(
					"Only instance with precomputed descriptors can be saved to file"));
		}
		FileWriter fw = new FileWriter(filename);
		fw.write(descriptorsSet.length + " " + precomputed.size() + " ");
		for (int i = 0; i < descriptorsSet.length; i++) {
			fw.write(this.descriptorsSet[i].trim() + " ");
		}
		for (int i = 0; i < descriptorsSet.length; i++) {
			fw.write(this.min[i] + " ");
		}
		for (int i = 0; i < descriptorsSet.length; i++) {
			fw.write(this.max[i] + " ");
		}
		fw.write("\n");
		for (Entry<String, float[]> e : precomputed.entrySet()) {
			fw.write(e.getKey() + " ");
			float[] d = e.getValue();
			for (int i = 0; i < d.length; i++) {
				fw.write(d[i] + " ");
			}
			fw.write("\n");
		}
		fw.close();
	}

	/**
	 * Load new instance from file. If no precomputed compounds in file
	 * usePrecomputed_ will be false, else will be true. File format described in
	 * method saveToFile.
	 * 
	 * @param filename file name
	 * @return new instance
	 * @throws IOException IO
	 */
	public static Descriptors readFromFile(String filename) throws IOException {
		BufferedReader inp = new BufferedReader(new InputStreamReader(new FileInputStream(new File(filename))));
		String s = inp.readLine();
		String[] split = s.split("\\s+");
		int nd = Integer.parseInt(split[0]);
		int ncomp = Integer.parseInt(split[1]);
		String[] descriptorsSet_ = new String[nd];
		float[] min_ = new float[nd];
		float[] max_ = new float[nd];

		int j = 2;
		for (int i = 0; i < nd; i++) {
			descriptorsSet_[i] = split[j];
			j++;
		}
		for (int i = 0; i < nd; i++) {
			min_[i] = Float.parseFloat(split[j]);
			j++;
		}
		for (int i = 0; i < nd; i++) {
			max_[i] = Float.parseFloat(split[j]);
			j++;
		}
		Descriptors result = Descriptors.instance(descriptorsSet_, min_, max_, (ncomp != 0));
		for (int i = 0; i < ncomp; i++) {
			s = inp.readLine();
			split = s.split("\\s+");
			String smiles = split[0];
			float[] d = new float[nd];
			for (int k = 0; k < nd; k++) {
				d[k] = Float.parseFloat(split[k + 1]);
			}
			result.precomputed.put(smiles, d);
		}
		inp.close();
		return result;
	}

	/**
	 * 
	 * @param a float[] array
	 * @return how many values are NaN in this array
	 */
	public static int countNaNs(float[] a) {
		int nans = 0;
		for (int i = 0; i < a.length; i++) {
			if (Float.isNaN(a[i])) {
				nans++;
			}
		}
		return nans;
	}

	/**
	 * 
	 * @param a float[] array
	 * @return replace NaNs with zeros (this method returns a copy, doesn't touch
	 *         a).
	 */
	public static float[] nansToZero(float[] a) {
		float[] b = new float[a.length];
		for (int i = 0; i < a.length; i++) {
			if (Float.isNaN(a[i])) {
				b[i] = 0;
			} else {
				b[i] = a[i];
			}
		}
		return b;
	}

	private static class DescriptorsOfOneCompound {
		public String smiles;
		public float[] d;

		public void compute(String[] descriptorsSet, float[] min, float[] max) throws CDKException {
			d = Chemoinformatics.descriptors(smiles, descriptorsSet, min, max);
		}

		public void compute(String[] descriptorsSet) throws CDKException {
			d = Chemoinformatics.descriptors(smiles, descriptorsSet);
		}
	}

	/**
	 * This array stores names of all molecular descriptors supported by CDK but 3D
	 * descriptors (descriptors which require 3D coordinates) and two 2D
	 * descriptors: nAtomLAC and MolIP. This two descriptors are calculated
	 * unreasonable slow for many molecules.
	 */
	public static final String[] descriptors2DBut_nAtomLAC_And_MolIP = new String[] { "WPATH", "WPOL", "MLogP",
			"nSmallRings", "nAromRings", "nRingBlocks", "nAromBlocks", "nRings3", "nRings4", "nRings5", "nRings6",
			"nRings7", "nRings8", "nRings9", "SC-3", "SC-4", "SC-5", "SC-6", "VC-3", "VC-4", "VC-5", "VC-6", "HybRatio",
			"nAtom", "nAromBond", "MDEC-11", "MDEC-12", "MDEC-13", "MDEC-14", "MDEC-22", "MDEC-23", "MDEC-24",
			"MDEC-33", "MDEC-34", "MDEC-44", "MDEO-11", "MDEO-12", "MDEO-22", "MDEN-11", "MDEN-12", "MDEN-13",
			"MDEN-22", "MDEN-23", "MDEN-33", "ALogP", "ALogp2", "AMR", "khs.sLi", "khs.ssBe", "khs.ssssBe", "khs.ssBH",
			"khs.sssB", "khs.ssssB", "khs.sCH3", "khs.dCH2", "khs.ssCH2", "khs.tCH", "khs.dsCH", "khs.aaCH",
			"khs.sssCH", "khs.ddC", "khs.tsC", "khs.dssC", "khs.aasC", "khs.aaaC", "khs.ssssC", "khs.sNH3", "khs.sNH2",
			"khs.ssNH2", "khs.dNH", "khs.ssNH", "khs.aaNH", "khs.tN", "khs.sssNH", "khs.dsN", "khs.aaN", "khs.sssN",
			"khs.ddsN", "khs.aasN", "khs.ssssN", "khs.sOH", "khs.dO", "khs.ssO", "khs.aaO", "khs.sF", "khs.sSiH3",
			"khs.ssSiH2", "khs.sssSiH", "khs.ssssSi", "khs.sPH2", "khs.ssPH", "khs.sssP", "khs.dsssP", "khs.sssssP",
			"khs.sSH", "khs.dS", "khs.ssS", "khs.aaS", "khs.dssS", "khs.ddssS", "khs.sCl", "khs.sGeH3", "khs.ssGeH2",
			"khs.sssGeH", "khs.ssssGe", "khs.sAsH2", "khs.ssAsH", "khs.sssAs", "khs.sssdAs", "khs.sssssAs", "khs.sSeH",
			"khs.dSe", "khs.ssSe", "khs.aaSe", "khs.dssSe", "khs.ddssSe", "khs.sBr", "khs.sSnH3", "khs.ssSnH2",
			"khs.sssSnH", "khs.ssssSn", "khs.sI", "khs.sPbH3", "khs.ssPbH2", "khs.sssPbH", "khs.ssssPb", "nHBDon",
			"BCUTw-1l", "BCUTw-1h", "BCUTc-1l", "BCUTc-1h", "BCUTp-1l", "BCUTp-1h", "nBase", "Fsp3", "PetitjeanNumber",
			"topoShape", "nSpiroAtoms", "bpol", "VABC", "ECCEN", "VAdjMat", "JPLogP", "nAtomLC", "nAcid", "XLogP",
			"apol", "LipinskiFailures", "Zagreb", "WTPT-1", "WTPT-2", "WTPT-3", "WTPT-4", "WTPT-5", "MW", "SP-0",
			"SP-1", "SP-2", "SP-3", "SP-4", "SP-5", "SP-6", "SP-7", "VP-0", "VP-1", "VP-2", "VP-3", "VP-4", "VP-5",
			"VP-6", "VP-7", "nRotB", "FMF", "TopoPSA", "ATSc1", "ATSc2", "ATSc3", "ATSc4", "ATSc5", "SCH-3", "SCH-4",
			"SCH-5", "SCH-6", "SCH-7", "VCH-3", "VCH-4", "VCH-5", "VCH-6", "VCH-7", "C1SP1", "C2SP1", "C1SP2", "C2SP2",
			"C3SP2", "C1SP3", "C2SP3", "C3SP3", "C4SP3", "SPC-4", "SPC-5", "SPC-6", "VPC-4", "VPC-5", "VPC-6",
			"tpsaEfficiency", "nAtomP", "ATSm1", "ATSm2", "ATSm3", "ATSm4", "ATSm5", "fragC", "Kier1", "Kier2", "Kier3",
			"naAromAtom", "ATSp1", "ATSp2", "ATSp3", "ATSp4", "ATSp5", "nHBAcc", "nB", "nA", "nR", "nN", "nD", "nC",
			"nF", "nQ", "nE", "nG", "nH", "nI", "nP", "nL", "nK", "nM", "nS", "nT", "nY", "nV", "nW" };
}
