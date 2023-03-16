package ru.ac.phyche.gcms.svekla;

import java.util.Arrays;
import java.util.HashSet;

/**
 * 
 * This class works with integer identifiers and names of standard non-polar and
 * semi-standard non-polar columns. 0 -- "DB-1", 1 -- "SE-30", 2 -- "OV-101", 3
 * --"OV-1", 4 -- "Methyl_Silicone", 5 -- "CP_Sil_5_CB", 6 -- "BP-1", 7 --
 * "HP-1", 8 -- "SPB-1", 9 -- "RTX-1", 10 -- "Ultra-1", 11 -- "Petrocol_DH", 12
 * -- "Polydimethyl_siloxane", 13 --
 * "OV-1,_SE-30,_Methyl_silicone,_SP-2100,_OV-101,_DB-1,_etc.", 14 --
 * "Other_non_polar", 15 -- "5_%_Phenyl_methyl_siloxane", 16 -- "DB-5", 17 --
 * "HP-5", 18 -- "HP-5MS", 19 -- "VF-5MS", 20 -- "Squalane", 21 -- "HP-5_MS", 22
 * -- "SE-54", 23 -- "DB-5MS", 24 -- "Apiezon_L", 25 -- "BPX-5", 26 -- "SE-52",
 * 27 --"CP_Sil_8_CB", 28 -- "SPB-5", 29 -- "RTX-5", 30 -- "TR5-MS", 31 --
 * "Ultra-2", 32 -- "DB-5_MS", 33 -- "ZB-5", 34 -- "SLB-5_MS" , 35 --
 * "Other_semi_non_polar". All column names besides "Other_non_polar" and
 * "Other_semi_non_polar" are given according to column names from NIST 17
 * database. Many of them actually describe very similar or equal stationary
 * phases.
 *
 */
public class Columns {
	private static final String[] nonPolarOther = new String[] { "Polydimethyl_siloxanes", "SF-96", "ZB-1",
			"Cross-Linked_Methylsilicone", "PONA", "CP-Sil_5_CB", "HP-101", "CBP-1", "DC-200", "NB-30", "SP-2100",
			"Petrocol_DH-100", "HP-PONA", "Polydimethyl_siloxane:_CP-Sil_5_CB", "E-301", "Normal_alkane_RI", "DB-1MS",
			"RTX-1_PONA", "CP-Sil5_CB_MS", "CP-Sil", "Polymethylsiloxane,_(PMS-20000)", "DB-1HT", "DB-1-MS", "PMS-100",
			"SPB-Sulfur", "Optima-1", "DB-Petro", "CP_Sil_2", "CB-1", "CP-Sil_5", "Ultra-ALLOY-5", "TR-1", "Elite-1",
			"JXR", "RTx-1", "RSL-150", "CP-Sil_PONA_GB", "PB-1", "ZB-1_MS", "DB-1_MS", "LM-1", "EC-1", "CP-Sil_5_Cb",
			"Ultra_1", "SE-33", "DB-Petro_100", "SE-30_MS", "SP_2100", "Optima_1", "PMS-1000", "Se-30", "AT-1",
			"Polidimethyl_siloxane", "Permaphase_DMS", "Solgel-1_(SGE)", "CP_Sil-5_CB", "VF-1_MS", "SSP-1", "SPD-1",
			"Methyl_silicone", "HP_Ultra_1", "TC-1", "DB-1_HT", "HP_Ultra-1", "PE-1HT", "HP-1MS", "BPA-1",
			"Ultra-1_PONA", "GP_SP_2100_DB" };
	private static final String[] nonPolarMain = new String[] { "DB-1", "SE-30", "OV-101", "OV-1", "Methyl_Silicone",
			"CP_Sil_5_CB", "BP-1", "HP-1", "SPB-1", "RTX-1", "Ultra-1", "Petrocol_DH", "Polydimethyl_siloxane",
			"OV-1,_SE-30,_Methyl_silicone,_SP-2100,_OV-101,_DB-1,_etc." };
	private static final String[] semiNonPolarMain = new String[] { "5_%_Phenyl_methyl_siloxane", "DB-5", "HP-5",
			"HP-5MS", "VF-5MS", "Squalane", "HP-5_MS", "SE-54", "DB-5MS", "Apiezon_L", "BPX-5", "SE-52", "CP_Sil_8_CB",
			"SPB-5", "RTX-5", "TR5-MS", "Ultra-2", "DB-5_MS", "ZB-5", "SLB-5_MS" };
	private static final String[] semiNonPolarOther = new String[] { "Apiezon_M", "Vacuum_Grease_Oil_(VM-4)",
			"Polydimethyl_siloxane_with_5_%_Ph_groups", "Rxi-5Sil", "Rtx-5MS", "BP-5", "RSL-200", "Apolane",
			"CP-Sil_8CB-MS", "Rxi-5MS", "VF-5_MS", "RTX-5_MS", "MDN-5", "Rxi-5SilMS", "SLB-5MS", "Mega_5MS",
			"Dexsil_300", "C78,_Branched_paraffin", "LM-5", "OV-3", "PTE-5", "TR-5_MS", "CBP-5", "PE-5", "AT-5", "VF-5",
			"Siloxane,_5_%_Ph", "Polydimethyl_siloxane,_5_%_phenyl", "Apieson_L", "Apieson_L_/_KOH", "Porapack_Q",
			"HG-5", "5_%_Phenyl_polydimethyl_siloxane", "Optima-5", "Nonpolar", "OV-5", "MFE-73", "HT-5", "Apiezon",
			"Elite-5_MS", "Elite-5MS", "Lucopren_G_(silicone_elastomer)", "BPX5", "EC-5",
			"Polydimethyl_siloxane,_unknown_content_of_Ph-groups", "DBP-5", "Apiezon_L_+_KF", "SE-30+Igepal",
			"Apieson_M", "Col-Elite_5MS", "Optima_5", "Rxi-1MS", "RTX-5Sil", "PE-5ht", "NB-54", "CP-Sil_8_CB",
			"methyl_silicone_oil_with_5%_Igepal", "Silicon_High_Vacuum_Grease_(obsolete)", "RTx-5_Sil_MS", "Optima-5MS",
			"SF96+Igepal", "DC-400", "Synachrom", "Rxi-5_MS", "FSOT-RSL-200", "Methyl_phenyl_siloxane_(not_specified)",
			"Mega-5", "PoraPLOT_Q", "Optima-5_MS", "UCW-98", "Ultra_2", "ZB-5_MS", "5_%_Phenyl_methylsiloxane", "DP-5",
			"Methylsiloxane,_5_%_Ph_groups", "Triacontane", "Polydimethyl_siloxane_with_5_%_phenyl_groups", "XTI-5",
			"Apiezon_LH_+_KF", "Equity-5MS", "Equity-5", "PS-255", "VB-5", "OV-73", "ZP-5", "PE-5ms", "MS5",
			"Polydimethyl_siloxane_with_5_%_Ph", "Mega_5_MS", "HP5-MS", "MDN-5S", "Apiezon_N", "NB-5", "Silicone_oil",
			"Methylsiloxane_+_5_%_Ph-groups", "PoraPLOT", "RTX-5Sil_MS", "Equity-5_MS", "RTX-5_Sil_MS",
			"Apiezon_L_+_KOH", "C103H208", "Durabond-5", "CP-Sil8", "SBP-5", "HP-5_MS,_DB-5_MS", "HP_Ultra_2",
			"DB5-30W", "Adamantyl_siloxane", "CP-SIL8", "SLB-5ms", "DB-5HT", "5_%_Phenyl_silicone", "ZEBRON-5",
			"Octacosane", "OV-101_+_Igepal", "DB-5;_CP-Sil_8_CB", "SE-30_(10_%)_+_CW-20M_(1_%)", "MPS-5", "VA-5_MS",
			"n-Dotriacontane", "Dexsil_400", "Chromosorb_101", "RH-5MS", "SE-52/54", "RTV-502", "OV-22",
			"Polydimethyl_siloxane,_5_%_phenyl_groups" };
	private static final HashSet<String> nonPolarOtherSet = new HashSet<String>(Arrays.asList(nonPolarOther));
	private static final HashSet<String> nonPolarMainSet = new HashSet<String>(Arrays.asList(nonPolarMain));
	private static final HashSet<String> semiNonPolarMainSet = new HashSet<String>(Arrays.asList(semiNonPolarMain));
	private static final HashSet<String> semiNonPolarOtherSet = new HashSet<String>(Arrays.asList(semiNonPolarOther));

	/**
	 * 
	 * @param column name of column according to NIST 17
	 * @return corresponding integer identifier
	 */
	public static int columnNum(String column) {
		if (nonPolarMainSet.contains(column)) {
			for (int i = 0; i < nonPolarMain.length; i++) {
				if (nonPolarMain[i].equals(column)) {
					return i;
				}
			}
		}
		if (nonPolarOtherSet.contains(column)) {
			return nonPolarMain.length;
		}
		if (semiNonPolarMainSet.contains(column)) {
			for (int i = 0; i < semiNonPolarMain.length; i++) {
				if (semiNonPolarMain[i].equals(column)) {
					return i + 1 + nonPolarMain.length;
				}
			}
		}
		if (semiNonPolarOtherSet.contains(column)) {
			return 1 + nonPolarMain.length + semiNonPolarMain.length;
		}
		return -1;
	}

	/**
	 * 
	 * @param column column name
	 * @return TRUE if it is standard non-polar column, else FALSE
	 */
	public static boolean isNonPolar(String column) {
		if (column == null) {
			return false;
		}
		if (nonPolarMainSet.contains(column)) {
			return true;
		}
		if (nonPolarOtherSet.contains(column)) {
			return true;
		}
		if (column.equals("Other_non_polar")) {
			return true;
		}
		return false;
	}

	/**
	 * 
	 * @param column column name
	 * @return TRUE if it is semi-standard non-polar column, else FALSE
	 */
	public static boolean isSemiNonPolar(String column) {
		if (column == null) {
			return false;
		}
		if (semiNonPolarMainSet.contains(column)) {
			return true;
		}
		if (semiNonPolarOtherSet.contains(column)) {
			return true;
		}
		if (column.equals("Other_semi_non_polar")) {
			return true;
		}
		return false;
	}

	/**
	 * 
	 * @param column integer id
	 * @return TRUE if it is standard non-polar column, else FALSE
	 */
	public static boolean isNonPolar(int column) {
		if ((column < 1 + nonPolarMain.length) && (column > -1)) {
			return true;
		}
		return false;
	}

	/**
	 * 
	 * @param column integer id
	 * @return TRUE if it is semi-standard non-polar column, else FALSE
	 */
	public static boolean isSemiNonPolar(int column) {
		if ((column > nonPolarMain.length) && (column < 2 + nonPolarMain.length + semiNonPolarMain.length)) {
			return true;
		}
		return false;
	}

	/**
	 * Returns name by id
	 * 
	 * @param n integer id
	 * @return column name
	 */
	public static String column(int n) {
		if ((n < 0) || n > nonPolarMain.length + semiNonPolarMain.length + 1) {
			return "Unknown";
		}
		if (n == nonPolarMain.length) {
			return "Other_non_polar";
		}
		if (n == 1 + nonPolarMain.length + semiNonPolarMain.length) {
			return "Other_semi_non_polar";
		}
		if (n < nonPolarMain.length) {
			return nonPolarMain[n];
		}
		if ((n > nonPolarMain.length) && (n < 1 + nonPolarMain.length + semiNonPolarMain.length)) {
			return semiNonPolarMain[n - nonPolarMain.length - 1];
		}
		return "Unknown";
	}

	/**
	 * Encode column as one-hot float vector
	 * 
	 * @param column integer id
	 * @return float[36] with one 1.0F value and 35 0.0F values
	 */
	public static float[] columnOneHot(int column) {
		float[] result = new float[2 + nonPolarMain.length + semiNonPolarMain.length];
		if ((column > -1) && (column < result.length)) {
			result[column] = 1.0F;
		}
		return result;
	}

	/**
	 * 
	 * @param column integer id
	 * @return new float[]{1.0F, 0.0F} for standard non-polar column, new
	 *         float[]{0.0F, 1.0F} for semi-standard non-polar column
	 */
	public static float[] columnTypeOneHot(int column) {
		float[] result = new float[2];
		result[0] = isNonPolar(column) ? 1.0F : 0.0F;
		result[1] = isSemiNonPolar(column) ? 1.0F : 0.0F;
		return result;
	}

	/**
	 * 
	 * @param column integer id
	 * @return a result of concatenation of results of columnOneHot(column) and
	 *         columnTypeOneHot(column).
	 */
	public static float[] columnAndColumnTypeOneHot(int column) {
		return RetentionsDataset.mergeArrays(columnOneHot(column), columnTypeOneHot(column));
	}
}
