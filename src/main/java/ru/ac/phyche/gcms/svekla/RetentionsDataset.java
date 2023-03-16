package ru.ac.phyche.gcms.svekla;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Random;

import org.openscience.cdk.exception.CDKException;

/**
 * Retention data set, which can be used for training and testing. It is the
 * core class with routines such as data set splitting and shuffling and so on.
 * Internally it's array of RetentionsEntry
 * 
 */
public class RetentionsDataset {

	private RetentionsEntry[] data;

	/**
	 * 
	 * @return Number of entries (retention index - column type - SMILES) in the
	 *         data set.
	 */
	public int size() {
		return data.length;
	}

	/**
	 * Returns complete data set as array of RetentionsEntry
	 * 
	 * @return RetentionsEntry[]
	 */
	public RetentionsEntry[] getData() {
		return data;
	}

	/**
	 * 
	 * @param data complete data set as array of RetentionsEntry
	 */
	public void setData(RetentionsEntry[] data) {
		this.data = data;
	}

	/**
	 * getData()[i]
	 * 
	 * @param i number of entry to return
	 * @return i-th entry
	 */
	public RetentionsEntry getEntry(int i) {
		return data[i];
	}

	/**
	 * 
	 * @param i number of entry
	 * @return SMILES of i-th entry
	 */
	public String getSmiles(int i) {
		return data[i].getSmiles();
	}

	/**
	 * 
	 * @param i number of entry
	 * @return retention index of i-th entry
	 */
	public float getRetention(int i) {
		return data[i].getRetention();
	}

	/**
	 * 
	 * @param i number of entry
	 * @return column of i-th entry
	 */
	public int getColumn(int i) {
		return data[i].getColumnType();
	}

	/**
	 * data[i] = e;
	 * 
	 * @param e entry
	 * @param i number of entry
	 */
	public void setEntry(RetentionsEntry e, int i) {
		data[i] = e;
	}

	/**
	 * Create new instance from array
	 * 
	 * @param data_ array of RetentionsEntry
	 * @return new instance
	 */
	public static RetentionsDataset create(RetentionsEntry[] data_) {
		RetentionsDataset result = new RetentionsDataset();
		result.data = data_;
		return result;
	}

	/**
	 * Create new instance from ArrayList
	 * 
	 * @param data_ list of entries
	 * @return new instance
	 */
	public static RetentionsDataset create(ArrayList<RetentionsEntry> data_) {
		RetentionsDataset result = new RetentionsDataset();
		result.data = data_.toArray(new RetentionsEntry[data_.size()]);
		return result;
	}

	/**
	 * 
	 * @param t molecular fingerprints type
	 * @param i number of the entry
	 * @return molecular fingerprints for molecule from i-th entry
	 * @throws CDKException Chemoinformatics error
	 */
	public float[] fingerprints(Chemoinformatics.FingerprintsType t, int i) throws CDKException {
		return data[i].fingerprints(t);
	}

	/**
	 * See Chemoinformatics class and Descriptors class
	 * 
	 * @param i                    number of the entry
	 * @param descriptorsGenerator generator of descriptors
	 * @return molecular descriptors
	 * @throws CDKException CDK error, absence of precomputed value if usage of
	 *                      precomputed descriptors is enabled in
	 *                      descriptorsGenerator instance.
	 */
	public float[] descriptors(int i, Descriptors descriptorsGenerator) throws CDKException {
		return data[i].descriptors(descriptorsGenerator);
	}

	/**
	 * See Chemoinformatics class and Descriptors class. NaNs replaced with zeros.
	 * 
	 * @param i                    number of entry
	 * @param descriptorsGenerator generator of descriptors
	 * @return molecular descriptors (NaNs replaced with zeros)
	 * @throws CDKException CDK error, absence of precomputed value if usage of
	 *                      precomputed descriptors is enabled in
	 *                      descriptorsGenerator instance.CDK error, absence of pre
	 */
	public float[] descriptorsNoNaNs(int i, Descriptors descriptorsGenerator) throws CDKException {
		return data[i].descriptorsNoNaNs(descriptorsGenerator);
	}

	/**
	 * See Chemoinformatics class. Abundance of various functional groups in
	 * molecule. Groups enumerated according to work 10.1021/ci600548y, Table 2.
	 * 
	 * @param i number of the entry
	 * @return functional groups count.
	 * @throws CDKException CDK internal error
	 */
	public float[] funcGroups(int i) throws CDKException {
		return data[i].funcGroups();
	}

	/**
	 * See Chemoinformatics class.
	 * 
	 * @param i number of the entry
	 * @return 2D - depiction of structural formula of the molecule
	 * @throws CDKException CDK internal error
	 */
	public float[][] depiction(int i) throws CDKException {
		return data[i].depitction();
	}

	/**
	 * See Chemoinformatics class.
	 * 
	 * @param i number of the entry
	 * @return 2D - representation of molecule for 2D-CNN
	 * @throws CDKException CDK internal error, oversized molecule (typically more
	 *                      than C50)
	 */
	public float[][][] representation2d(int i) throws CDKException {
		return data[i].representation2d();
	}

	/**
	 * 
	 * @param i                    number of the entry
	 * @param descriptorsGenerator descriptors generator
	 * @return concatenated (merged) result of descriptorsNoNaNs(i,
	 *         descriptorsGenerator) and funcGroups(i)
	 * @throws CDKException CDK error, descriptors computation error
	 */
	public float[] extDescriptors(int i, Descriptors descriptorsGenerator) throws CDKException {
		return mergeArrays(this.descriptorsNoNaNs(i, descriptorsGenerator), this.funcGroups(i));
	}

	/**
	 * 
	 * @param i number of the entry
	 * @param t fingerprints type
	 * @return concatenated (merged) result of fingerprints(t, i and funcGroups(i)
	 * @throws CDKException CDK error
	 */
	public float[] extFingerprints(Chemoinformatics.FingerprintsType t, int i) throws CDKException {
		return mergeArrays(this.fingerprints(t, i), this.funcGroups(i));
	}

	/**
	 * Deep copy of the data set (with copying in memory of all records)
	 * 
	 * @return deep copy
	 */
	public RetentionsDataset copy() {
		RetentionsDataset result = new RetentionsDataset();
		result.data = new RetentionsEntry[this.data.length];
		for (int i = 0; i < this.data.length; i++) {
			result.data[i] = this.data[i].deepclone();
		}
		return result;
	}

	/**
	 * Shuffle this data set. Change order of entries.
	 */
	public void shuffle() {
		Random rnd = new Random();
		for (int c = 0; c < 10; c++) {
			for (int i = 0; i < this.data.length; i++) {
				int j = rnd.nextInt(this.data.length);
				RetentionsEntry a = data[j];
				data[j] = data[i];
				data[i] = a;
			}
		}
	}

	/**
	 * This method allows to obtain a set of compounds which are conatined in this
	 * data set. Each compounds is contained only once (SMILES are stored in the set
	 * in "canonical" form). If the data set contains various SMILES strings for
	 * really one compound, it will be contained in the HashSet only once.
	 * 
	 * @param stereochemistry treat cis/trans and optical compounds as different
	 *                        (TRUE) or as identical (FALSE).
	 * @return HashSet of SMILES strings
	 * @throws CDKException CDK error during canonicalization
	 */
	public HashSet<String> compoundsCanonical(boolean stereochemistry) throws CDKException {
		HashSet<String> result = new HashSet<String>();
		for (int i = 0; i < this.data.length; i++) {
			result.add(Chemoinformatics.canonical(data[i].getSmiles(), stereochemistry).trim());
		}
		return result;
	}

	/**
	 * 
	 * @return HashSet of SMILES strings from this data set. No canonicalization is
	 *         used. It can contain actually more than one record for a structure
	 *         due to ambiguity of SMILES representation.
	 */
	public HashSet<String> compounds() {
		HashSet<String> result = new HashSet<String>();
		for (int i = 0; i < this.data.length; i++) {
			result.add(data[i].getSmiles().trim());
		}
		return result;
	}

	/**
	 * 
	 * @return HashSet of InChI strings for compounds from this data set.
	 *         Stereoisomers will be treated as DIFFERENT(!!!) compound if
	 *         stereoisomers are already stored with symbols denoting
	 *         stereochemistry.
	 * @throws CDKException CDK error
	 */
	public HashSet<String> inchiIds() throws CDKException {
		HashSet<String> result = new HashSet<String>();
		for (int i = 0; i < this.data.length; i++) {
			result.add(Chemoinformatics.smilesToInchi(data[i].getSmiles()).trim());
		}
		return result;
	}

	/**
	 * 
	 * @return HashSet of InChI-key strings for compounds from this data set.
	 *         Stereoisomers will be treated as DIFFERENT(!!!) compound if
	 *         stereoisomers are already stored with symbols denoting
	 *         stereochemistry.
	 * @throws CDKException CDK error
	 */
	public HashSet<String> inchiKeys() throws CDKException {
		HashSet<String> result = new HashSet<String>();
		for (int i = 0; i < this.data.length; i++) {
			result.add(Chemoinformatics.smilesToInchiKey(data[i].getSmiles()).trim());
		}
		return result;
	}

	/**
	 * Comparison using InChI strings
	 * 
	 * @param second other data set
	 * @return number of compounds which are contained in both data sets
	 *         simultaneously. Stereoisomers will be treated as DIFFERENT(!!!)
	 *         compound if stereoisomers are already stored with symbols denoting
	 *         stereochemistry. CDK
	 * @throws CDKException CDK
	 */
	public int countIdenticalByInchi(RetentionsDataset second) throws CDKException {
		HashSet<String> a = this.inchiIds();
		HashSet<String> b = second.inchiIds();
		a.retainAll(b);
		return a.size();
	}

	/**
	 * Comparison using InChI-key strings
	 * 
	 * @param second other data set
	 * @return number of compounds which are contained in both data sets
	 *         simultaneously. Stereoisomers will be treated as DIFFERENT(!!!)
	 *         compound if stereoisomers are already stored with symbols denoting
	 *         stereochemistry.
	 * @throws CDKException CDK
	 */
	public int countIdenticalByInchikeys(RetentionsDataset second) throws CDKException {
		HashSet<String> a = this.inchiKeys();
		HashSet<String> b = second.inchiKeys();
		a.retainAll(b);
		return a.size();
	}

	/**
	 * Comparison using "canonical" SMILES strings
	 * 
	 * @param second other data set
	 * @return number of compounds which are contained in both data sets
	 *         simultaneously. Stereoisomers will be treated as DIFFERENT(!!!)
	 *         compound if stereoisomers are already stored with symbols denoting
	 *         stereochemistry.
	 * @throws CDKException CDK error
	 */
	public int countIdenticalByCanonicalSmiles(RetentionsDataset second) throws CDKException {
		HashSet<String> a = this.compoundsCanonical(true);
		HashSet<String> b = second.compoundsCanonical(true);
		a.retainAll(b);
		return a.size();
	}

	/**
	 * This methods calls countIdenticalByInchi method for each data set from array
	 * d1 and each data set from array d2.
	 * 
	 * @param d1 array of data sets
	 * @param d2 array of data sets
	 * @return int[d1.length][d2.length];
	 * @throws CDKException CDK error
	 */
	public static int[][] countIdenticalByInchi(RetentionsDataset[] d1, RetentionsDataset[] d2) throws CDKException {
		int[][] result = new int[d1.length][d2.length];
		for (int i = 0; i < d1.length; i++) {
			for (int j = 0; j < d2.length; j++) {
				result[i][j] = d1[i].countIdenticalByInchi(d2[j]);
			}
		}
		return result;
	}

	private HashSet<String> splitSet(HashSet<String> fullSet, int n) {
		String[] a = fullSet.toArray(new String[fullSet.size()]);
		Random rnd = new Random();
		for (int c = 1; c < 10; c++) {
			for (int i = 0; i < a.length; i++) {
				int j = rnd.nextInt(a.length);
				String b = a[j];
				a[j] = a[i];
				a[i] = b;
			}
		}
		HashSet<String> result = new HashSet<String>();
		for (int i = 0; i < n; i++) {
			result.add(a[i]);
		}
		return result;
	}

	/**
	 * Shuffle and split data sets such way that all entries corresponding to any
	 * compound will be contained in only one of subsets. I.e. will be no compounds
	 * which are contained in both subsets simultaneously. Note, that no
	 * canonicalization of SMILES is used here. If in data set there are cis/trans
	 * or optical isomers they will be treated as different compounds. Use
	 * makeCanonical(false) to remove such occurrences.
	 * 
	 * @param compoundsToSplit number of COMPOUNDS (!!!not data entries!!!) which
	 *                         will be separated to the subset.
	 * @return subset which contains all data entries for compoundsToSplit compounds
	 * @throws CDKException CDK error
	 */
	public RetentionsDataset compoundsBasedSplitAndShuffle(int compoundsToSplit) throws CDKException {
		this.shuffle();
		HashSet<String> split = splitSet(this.compounds(), compoundsToSplit);
		ArrayList<RetentionsEntry> splitData = new ArrayList<RetentionsEntry>();
		ArrayList<RetentionsEntry> retainData = new ArrayList<RetentionsEntry>();
		for (int i = 0; i < this.data.length; i++) {
			if (split.contains(data[i].getSmiles().trim())) {
				splitData.add(this.getEntry(i));
			} else {
				retainData.add(this.getEntry(i));
			}
		}
		RetentionsDataset result = new RetentionsDataset();
		result.data = splitData.toArray(new RetentionsEntry[splitData.size()]);
		this.data = retainData.toArray(new RetentionsEntry[retainData.size()]);
		result.shuffle();
		this.shuffle();
		return result;
	}

	/**
	 * The same as compoundsBasedSplitAndShuffle(int compoundsToSplit) but fraction
	 * of compounds that should be separated is given instead number of compounds.
	 * I.e if all data set contains 1000 compounds and fraction = 0.25 it means that
	 * 250 compounds will be isolated and 750 will be remained.
	 * 
	 * @param fraction number of compounds which will be in split in all compounds
	 * @return subset which contains all data entries for Math.round(fraction *
	 *         ((float) this.compounds().size())) compounds
	 * @throws CDKException CDK error
	 */
	public RetentionsDataset compoundsBasedSplitAndShuffle(float fraction) throws CDKException {
		int n = this.compounds().size();
		int split = Math.round(fraction * ((float) n));
		RetentionsDataset result = this.compoundsBasedSplitAndShuffle(split);
		return result;
	}

	/**
	 * Simple split (based on ENTRIES (records), not compounds. Many records can
	 * correspond to one compound. This method doesn't shuffle the data set. FIRST
	 * toSplit entries will be separated. this.size() - toSplit will be remained
	 * 
	 * @param toSplit how many entries should be separated.
	 * @return data set consisted of first toSplit records. This records will be
	 *         removed from this data set.
	 */
	public RetentionsDataset simpleSplit(int toSplit) {
		ArrayList<RetentionsEntry> splitData = new ArrayList<RetentionsEntry>();
		ArrayList<RetentionsEntry> retainData = new ArrayList<RetentionsEntry>();
		for (int i = 0; i < this.data.length; i++) {
			if (i < toSplit) {
				splitData.add(this.getEntry(i));
			} else {
				retainData.add(this.getEntry(i));
			}
		}
		RetentionsDataset result = new RetentionsDataset();
		result.data = splitData.toArray(new RetentionsEntry[splitData.size()]);
		this.data = retainData.toArray(new RetentionsEntry[retainData.size()]);
		return result;
	}

	/**
	 * Simple split (based on ENTRIES (records), not compounds. Many records can
	 * correspond to one compound. This method doesn't shuffle te data set. FIRST
	 * Math.round(fraction * ((float) this.size())) entries will be separated.
	 * 
	 * @param fraction fraction of the data set to be separated
	 * @return data set with first Math.round(fraction * ((float) this.size()))
	 *         records. This records will be removed from this data set.
	 */
	public RetentionsDataset simpleSplit(float fraction) {
		int n = this.size();
		int split = Math.round(fraction * ((float) n));
		RetentionsDataset result = this.simpleSplit(split);
		return result;
	}

	/**
	 * @param toSplit how many entries will be separated after shuffling
	 * @return set consisted of toSplit entries. This records will be removed from
	 *         this data set.
	 */
	public RetentionsDataset simpleShuffleSplit(int toSplit) {
		this.shuffle();
		RetentionsDataset result = this.simpleSplit(toSplit);
		return result;
	}

	/**
	 * 
	 * @param fraction fraction of the data set to be separated
	 * @return data set with first Math.round(fraction * ((float) this.size()))
	 *         records. This records will be removed from this data set.
	 */
	public RetentionsDataset simpleShuffleSplit(float fraction) {
		this.shuffle();
		RetentionsDataset result = this.simpleSplit(fraction);
		return result;
	}

	/**
	 * Remove from THIS data set all compounds that are contained in second data
	 * set! After this method will be no compounds in both data set simultaneously.
	 * If in data set there are cis/trans or optical isomers, they will be treated
	 * as different compounds. Use makeCanonical(false).
	 * 
	 * @param second another data set.
	 * @throws CDKException CDK errors
	 */
	public void filterIdentical(RetentionsDataset second) throws CDKException {
		ArrayList<RetentionsEntry> retainData = new ArrayList<RetentionsEntry>();
		HashSet<String> b = second.compoundsCanonical(true);
		for (int i = 0; i < this.data.length; i++) {
			if (!b.contains(Chemoinformatics.canonical(data[i].getSmiles(), true).trim())) {
				retainData.add(this.getEntry(i));
			}
		}
		this.data = retainData.toArray(new RetentionsEntry[retainData.size()]);
	}

	/**
	 * Remove from THIS data set all compounds that are contained in second data
	 * set! After this method will be no compounds in both data set simultaneously.
	 * If in data set there are cis/trans or optical isomers, they will be treated
	 * as different compounds. Use makeCanonical(false).
	 * 
	 * @param second another data set.
	 * @throws CDKException CDK errors
	 */
	public void filterIdenticalByInchi(RetentionsDataset second) throws CDKException {
		ArrayList<RetentionsEntry> retainData = new ArrayList<RetentionsEntry>();
		HashSet<String> b = second.inchiIds();
		for (int i = 0; i < this.data.length; i++) {
			if (!b.contains(Chemoinformatics.smilesToInchi(data[i].getSmiles()).trim())) {
				retainData.add(this.getEntry(i));
			}
		}
		this.data = retainData.toArray(new RetentionsEntry[retainData.size()]);
	}

	/**
	 * Save whole data set to file. File format: one line per entry, no empty lines,
	 * no comments. Each line contains SMILES string, retention index (float) and
	 * column type (integer). SMILES, retention index and type of column are
	 * separated by spaces. Example of line: "CCCC 400 0".
	 * 
	 * @param filename name of file (or path) for saving
	 * @throws IOException IO
	 */
	public void saveToFile(String filename) throws IOException {
		FileWriter fw = new FileWriter(filename);
		for (int i = 0; i < this.data.length; i++) {
			fw.write(this.getSmiles(i) + " " + this.getRetention(i) + " " + this.getColumn(i) + "\n");
		}
		fw.close();
	}

	/**
	 * Load data set from file. File format: one line per entry, no empty lines, no
	 * comments. Each line contains SMILES string, retention index (float) and
	 * column type (integer). SMILES, retention index and type of column are
	 * separated by spaces. Example of line: "CCCC 400 0". SMILES strings loaded "as
	 * is", without checking or converting to a canonical form.
	 * 
	 * @param filename file name
	 * @return loaded data set
	 * @throws IOException IO
	 */
	public static RetentionsDataset loadFromFile(String filename) throws IOException {
		BufferedReader inp = new BufferedReader(new InputStreamReader(new FileInputStream(new File(filename))));
		ArrayList<RetentionsEntry> data = new ArrayList<RetentionsEntry>();
		String s = inp.readLine();
		while ((s != null) && (!s.trim().equals(""))) {
			String[] spl = s.split("\\s+");
			data.add(RetentionsEntry.instance(spl[0], Float.parseFloat(spl[1]), Integer.parseInt(spl[2])));
			s = inp.readLine();
		}
		inp.close();
		return RetentionsDataset.create(data);
	}

	/**
	 * Merge multiple data sets into one
	 * 
	 * @param a array of data sets
	 * @return merged data set which contains all entries from all data sets from a
	 */
	public static RetentionsDataset merge(RetentionsDataset[] a) {
		ArrayList<RetentionsEntry> data = new ArrayList<RetentionsEntry>();
		for (int i = 0; i < a.length; i++) {
			data.addAll(Arrays.asList(a[i].getData()));
		}
		return RetentionsDataset.create(data);
	}

	/**
	 * Convert all SMILES strings in this data set to canonical form. See method
	 * canonical from class Chemoinformatics for more information.
	 * 
	 * @param stereochemistry if TRUE - cis/trans isomers and optical isomers will
	 *                        be denoted using special symbols. If false - cis/trans
	 *                        and optical isomers will be identical.
	 * @throws CDKException CDK error
	 */
	public void makeCanoncalAll(boolean stereochemistry) throws CDKException {
		for (int i = 0; i < this.data.length; i++) {
			this.getEntry(i).setSmiles(Chemoinformatics.canonical(this.getSmiles(i), stereochemistry));
		}
	}

	/**
	 * The data set can contain multiple entries for each compound (for each SMILES
	 * string). This method groups all entries for each compound together. For
	 * example if we have a data set with 6 entries "CCC 300 0","CCC 301 1", "CCCC
	 * 400 0","CCC 302 2", "CCCCC 500 0", "CCCC 401 1", this method with return
	 * HashMap with three elements. Key "CCC", ArrayList with three entries ("CCC
	 * 300 0","CCC 301 1", "CCC 302 2"); Key "CCC", ArrayList with two entries
	 * ("CCCC 400 0","CCCC 401 1"); Key "CCCCC", ArrayList with one entry ("CCCCC
	 * 500 0"). SMILES are converted to a canonical form! I.e. (C)CC and CCC will be
	 * grouped together.
	 * 
	 * @param stereochemistry if TRUE - cis/trans isomers and optical isomers will
	 *                        be considered as different compounds. If false -
	 *                        cis/trans and optical isomers will be identical.
	 * @return HashMap. Key - SMILES string, element - all entries with this
	 *         compound.
	 * @throws CDKException error during creation of canonical form.
	 */
	public HashMap<String, ArrayList<RetentionsEntry>> groupByCompounds(boolean stereochemistry) throws CDKException {
		HashMap<String, ArrayList<RetentionsEntry>> result = new HashMap<String, ArrayList<RetentionsEntry>>();
		HashSet<String> compounds = this.compoundsCanonical(stereochemistry);
		for (String smi : compounds) {
			result.put(smi, (new ArrayList<RetentionsEntry>()));
		}
		for (int i = 0; i < this.data.length; i++) {
			ArrayList<RetentionsEntry> c = result.get(Chemoinformatics.canonical(this.getSmiles(i), stereochemistry));
			c.add(this.getEntry(i));
		}
		return result;
	}

	/**
	 * The data set can contain multiple entries for each compound (for each SMILES
	 * string). This method groups all entries for each compound together. For
	 * example if we have a data set with 6 entries "CCC 300 0","CCC 301 1", "CCCC
	 * 400 0","CCC 302 2", "CCCCC 500 0", "CCCC 401 1", this method will return
	 * RetentionsDataset with three entries: "CCCC 301 -1","CCCC 400.5 -1", "CCCCC
	 * 500 -1". Column type everywhere will be -1 (due to averaging of values for
	 * various columns). All values for each compound will be averaged together.
	 * Result will have as many entries how many really different compounds has this
	 * data set.
	 * 
	 * @param stereochemistry if TRUE - cis/trans isomers and optical isomers will
	 *                        be considered as different compounds. If false -
	 *                        cis/trans and optical isomers will be identical.
	 * @return RetentionsDataset that has one entry per compound. All column types
	 *         will be -1.
	 * @throws CDKException error during creation of canonical form.
	 */
	public RetentionsDataset meanByCompounds(boolean stereochemistry) throws CDKException {
		return this.meanOrMedianByCompounds(stereochemistry, false);
	}

	/**
	 * The same as meanByCompounds but median value will be calculated for each
	 * compound instead average.
	 * 
	 * @param stereochemistry if TRUE - cis/trans isomers and optical isomers will
	 *                        be considered as different compounds. If false -
	 *                        cis/trans and optical isomers will be identical.
	 * @return RetentionsDataset that has one entry per compound. All column types
	 *         will be -1.
	 * @throws CDKException error during creation of canonical form.
	 */
	public RetentionsDataset medianByCompounds(boolean stereochemistry) throws CDKException {
		return this.meanOrMedianByCompounds(stereochemistry, true);
	}

	private RetentionsDataset meanOrMedianByCompounds(boolean stereochemistry, boolean median) throws CDKException {
		HashMap<String, ArrayList<RetentionsEntry>> a = this.groupByCompounds(stereochemistry);
		ArrayList<RetentionsEntry> data = new ArrayList<RetentionsEntry>();
		for (Entry<String, ArrayList<RetentionsEntry>> e : a.entrySet()) {
			ArrayList<RetentionsEntry> entries = e.getValue();
			float[] retentions = new float[entries.size()];
			for (int i = 0; i < entries.size(); i++) {
				retentions[i] = entries.get(i).getRetention();
			}
			float meanRet = median ? RetentionsDataset.median(retentions) : RetentionsDataset.mean(retentions);
			data.add(RetentionsEntry.instance(e.getKey(), meanRet));
		}
		RetentionsDataset result = RetentionsDataset.create(data);
		return result;
	}

	/**
	 * 
	 * @param a array of floats
	 * @return average value of elements of a. If a consists of 1.0F, 2.0F, 6.0F
	 *         result will be 3.0F
	 */
	public static float mean(float[] a) {
		float m = 0;
		for (int i = 0; i < a.length; i++) {
			m += a[i];
		}
		return m / ((float) a.length);
	}

	/**
	 * 
	 * @param a array of floats
	 * @return median value of elements of a. If a consists of 1.0F, 2.0F, 6.0F
	 *         result will be 2.0F
	 */
	public static float median(float[] a) {
		if (a.length == 0) {
			return (Float.NaN);
		}
		float[] b = a.clone();
		Arrays.sort(b);
		if (a.length % 2 == 1) {
			int n = a.length / 2;
			return b[n];
		} else {
			int n = a.length / 2;
			return ((b[n] + b[n - 1]) / 2.0F);
		}

	}

	/**
	 * 
	 * @param a array of float
	 * @param b array of float
	 * @return concatenate two arrays. float[a.length+b.length]
	 */
	public static float[] mergeArrays(float a[], float b[]) {
		float[] result = new float[a.length + b.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = (i < a.length) ? a[i] : b[i - a.length];
		}
		return result;
	}

}
