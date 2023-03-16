package ru.ac.phyche.gcms.svekla;

import org.openscience.cdk.exception.CDKException;

/**
 * An entry of a dataset of retentions. It stores information about Kovac
 * retention index, SMILES representation of structure of molecule and type of
 * column for which the index was acquired. For more information about methods
 * related to descriptors, fingerprints, functional groups - see
 * Chemoinformatics class. Here only wrappers.
 *
 */
public class RetentionsEntry {
	private String smiles;
	private float retention;
	private int columnType;

	/**
	 * 
	 * @return SMILES string
	 */
	public String getSmiles() {
		return smiles;
	}

	/**
	 * It is slow. It converts SMILES to InChI internally.
	 * 
	 * @return InChI string (structure of molecule)
	 * @throws CDKException internal CDK error
	 */
	public String getInchi() throws CDKException {
		return Chemoinformatics.smilesToInchi(smiles);
	}

	/**
	 * It is slow. It converts SMILES to InChI-key.
	 * 
	 * @return InChI-key
	 * @throws CDKException internal CDK error
	 */
	public String getInchikey() throws CDKException {
		return Chemoinformatics.smilesToInchiKey(smiles);
	}

	/**
	 * 
	 * @param smiles SMILES string
	 */
	public void setSmiles(String smiles) {
		this.smiles = smiles.trim();
	}

	/**
	 * 
	 * @return Kovac retention index value
	 */
	public float getRetention() {
		return retention;
	}

	/**
	 * 
	 * @param retention Kovac retention index value
	 */
	public void setRetention(float retention) {
		this.retention = retention;
	}

	/**
	 * 0 - "DB-1", 1 - "SE-30", 2 - "OV-101", 3 - "OV-1", 4 - "Methyl_Silicone", 5 -
	 * "CP_Sil_5_CB", 6 - "BP-1", 7 - "HP-1", 8 - "SPB-1", 9 - "RTX-1", 10 -
	 * "Ultra-1", 11 - "Petrocol_DH", 12 - "Polydimethyl_siloxane", 13 -
	 * "OV-1,_SE-30,_Methyl_silicone,_SP-2100,_OV-101,_DB-1,_etc.", 14 - other
	 * standard non polar, 15 - "5_%_Phenyl_methyl_siloxane", 16 - "DB-5", 17 -
	 * "HP-5", 18 - "HP-5MS", 19 - "VF-5MS", 20 - "Squalane", 21 - "HP-5_MS", 22 -
	 * "SE-54", 23 - "DB-5MS", 24 - "Apiezon_L", 25 - "BPX-5", 26 - "SE-52", 27 -
	 * "CP_Sil_8_CB", 28 - "SPB-5", 29 - "RTX-5", 30 - "TR5-MS", 31 - "Ultra-2", 32
	 * - "DB-5_MS", 33 - "ZB-5", 34 - "SLB-5_MS", 35 - other semi-standard nonpolar
	 * 
	 * @return column type
	 */
	public int getColumnType() {
		return columnType;
	}

	/**
	 * 
	 * @param columnType column type
	 */
	public void setColumnType(int columnType) {
		this.columnType = columnType;
	}

	/**
	 * Create new instance with converting SMILES string to canonical form.
	 * 
	 * @param smiles_         SMILES string of the molecule
	 * @param retention_      Kovac retention index
	 * @param columnType_     column type
	 * @param stereochemistry if true - symbols to denote cis/trans and optical
	 *                        isomers will be used
	 * @return newly created entry
	 * @throws CDKException internal CDK error
	 */
	public static RetentionsEntry instanceCanonical(String smiles_, float retention_, int columnType_,
			boolean stereochemistry) throws CDKException {
		RetentionsEntry result = new RetentionsEntry();
		result.columnType = columnType_;
		result.smiles = Chemoinformatics.canonical(smiles_, stereochemistry).trim();
		result.retention = retention_;
		return result;
	}

	/**
	 * 
	 * Create new instance with converting SMILES string to canonical form. Column
	 * type = -1
	 * 
	 * @param smiles_         SMILES string of the molecule
	 * @param retention_      Kovac retention index
	 * @param stereochemistry if true - symbols to denote cis/trans and optical
	 *                        isomers will be used
	 * @return newly created entry
	 * @throws CDKException internal CDK error
	 */
	public static RetentionsEntry instanceCanonical(String smiles_, float retention_, boolean stereochemistry)
			throws CDKException {
		return RetentionsEntry.instanceCanonical(smiles_, retention_, -1, stereochemistry);
	}

	/**
	 * Create new instance without (!) converting SMILES string to canonical form.
	 * 
	 * @param smiles_     SMILES string of the molecule
	 * @param retention_  Kovac retention index
	 * @param columnType_ column type
	 * @return newly created entry
	 */
	public static RetentionsEntry instance(String smiles_, float retention_, int columnType_) {
		RetentionsEntry result = new RetentionsEntry();
		result.columnType = columnType_;
		result.smiles = new String(smiles_.trim());
		result.retention = retention_;
		return result;
	}

	/**
	 * Create new instance without (!) converting SMILES string to canonical form.
	 * Column type =-1
	 * 
	 * @param smiles_    SMILES string of the molecule
	 * @param retention_ Kovac retention index
	 * @return newly created entry
	 */
	public static RetentionsEntry instance(String smiles_, float retention_) {
		return RetentionsEntry.instance(smiles_, retention_, -1);
	}

	@Override
	public RetentionsEntry clone() {
		return RetentionsEntry.instance(smiles, retention, columnType);
	}

	/**
	 * 
	 * @return this.clone();
	 */
	public RetentionsEntry deepclone() {
		return this.clone();
	}

	/**
	 * 
	 * @param e other RetentionsEntry
	 * @return true if this equals e
	 */
	public boolean equals(RetentionsEntry e) {
		return ((e.smiles.trim().equals(this.smiles.trim())) && (this.retention == e.retention)
				&& (this.columnType == e.columnType));
	}

	@Override
	public boolean equals(Object o) {
		if (o.getClass().equals(this.getClass())) {
			RetentionsEntry e = (RetentionsEntry) o;
			return ((e.smiles.trim().equals(this.smiles.trim())) && (this.retention == e.retention)
					&& (this.columnType == e.columnType));
		} else {
			return false;
		}
	}

	/**
	 * See Chemoinformatics class!!!
	 * 
	 * @param t molecular fingerprints type
	 * @return fingerprints
	 * @throws CDKException internal CDK error
	 */
	public float[] fingerprints(Chemoinformatics.FingerprintsType t) throws CDKException {
		return Chemoinformatics.fingerprints(smiles, t);
	}

	/**
	 * See Chemoinformatics class!!!
	 * 
	 * @return numbers - how many functional groups of each type the molecule
	 *         contains.
	 * @throws CDKException internal CDK error
	 */
	public float[] funcGroups() throws CDKException {
		return Chemoinformatics.funcGroups(smiles);
	}

	/**
	 * See Descriptors class!!!
	 * 
	 * @param d Generator of descriptors. See Descriptors class
	 * @return molecular descriptors (appropriate scaled)
	 * @throws CDKException CDK error, absence of precomputed value if usage of
	 *                      precomputed descriptors is enabled in d instance.
	 */
	public float[] descriptors(Descriptors d) throws CDKException {
		return d.get(smiles);
	}

	/**
	 * Returns molecular descriptors. NaNs replaced with zeros.
	 * 
	 * @param d Generator of descriptors. See Descriptors class
	 * @return molecular descriptors (appropriate scaled). NaNs replaced with zeros.
	 * @throws CDKException CDK error, absence of precomputed value if usage of
	 *                      precomputed descriptors is enabled in d instance.
	 */
	public float[] descriptorsNoNaNs(Descriptors d) throws CDKException {
		return d.getNoNaNs(smiles);
	}

	/**
	 * See Chemoinformatics class!!! Create image (depiction of chemical structure)
	 * of molecule. Creates 2-color image with size DEPICTION_SIZE*DEPICTION_SIZE,
	 * each pixel is 0.0F and 1.0F. Blank space (background) 0.0F.
	 * 
	 * @return float[DEPICTION_SIZE][DEPICTION_SIZE] - depiction
	 * @throws CDKException CDK internal errors
	 */
	public float[][] depitction() throws CDKException {
		return Chemoinformatics.smilesToImage(smiles);
	}

	/**
	 * See Chemoinformatics class!!! 2D representation of molecule for 2D-CNN.
	 * 
	 * @return float [29][130][130] input for 2D-CNN
	 * @throws CDKException CDK internal errors, oversized molecules.
	 */
	public float[][][] representation2d() throws CDKException {
		return Chemoinformatics.representation2d(smiles);
	}
}
