package ru.ac.phyche.gcms.svekla;

import junit.framework.Assert;
import junit.framework.TestCase;

public class ColumnsTest extends TestCase {

	public void testColumnNum() {
		Assert.assertEquals(0, Columns.columnNum("DB-1"));
		Assert.assertEquals(1, Columns.columnNum("SE-30"));
		Assert.assertEquals(2, Columns.columnNum("OV-101"));
		Assert.assertEquals(3, Columns.columnNum("OV-1"));
		Assert.assertEquals(4, Columns.columnNum("Methyl_Silicone"));
		Assert.assertEquals(13, Columns.columnNum("OV-1,_SE-30,_Methyl_silicone,_SP-2100,_OV-101,_DB-1,_etc."));
		Assert.assertEquals(14, Columns.columnNum("Polydimethyl_siloxanes"));
		Assert.assertEquals(14, Columns.columnNum("Cross-Linked_Methylsilicone"));
		Assert.assertEquals(14, Columns.columnNum("SF-96"));
		Assert.assertEquals(14, Columns.columnNum("GP_SP_2100_DB"));
		Assert.assertEquals(14, Columns.columnNum("Ultra-1_PONA"));
		Assert.assertEquals(14, Columns.columnNum("Optima_1"));
		Assert.assertEquals(14, Columns.columnNum("Se-30"));
		Assert.assertEquals(15, Columns.columnNum("5_%_Phenyl_methyl_siloxane"));
		Assert.assertEquals(16, Columns.columnNum("DB-5"));
		Assert.assertEquals(17, Columns.columnNum("HP-5"));
		Assert.assertEquals(18, Columns.columnNum("HP-5MS"));
		Assert.assertEquals(20, Columns.columnNum("Squalane"));
		Assert.assertEquals(33, Columns.columnNum("ZB-5"));
		Assert.assertEquals(34, Columns.columnNum("SLB-5_MS"));
		Assert.assertEquals(35, Columns.columnNum("Apiezon_M"));
		Assert.assertEquals(35, Columns.columnNum("Vacuum_Grease_Oil_(VM-4)"));
		Assert.assertEquals(35, Columns.columnNum("Polydimethyl_siloxane_with_5_%_Ph_groups"));
		Assert.assertEquals(35, Columns.columnNum("Polydimethyl_siloxane,_5_%_phenyl_groups"));
		Assert.assertEquals(35, Columns.columnNum("OV-22"));
		Assert.assertEquals(35, Columns.columnNum("RTV-502"));
		Assert.assertEquals(35, Columns.columnNum("PoraPLOT"));
		Assert.assertEquals(35, Columns.columnNum("Silicone_oil"));
		Assert.assertEquals(35, Columns.columnNum("NB-5"));
		Assert.assertEquals(-1, Columns.columnNum(null));
		Assert.assertEquals(-1, Columns.columnNum(" "));
		Assert.assertEquals(-1, Columns.columnNum(""));
		Assert.assertEquals(-1, Columns.columnNum("aaa"));
	}

	public void testIsNonPolarString() {
		Assert.assertEquals(true, Columns.isNonPolar("DB-1"));
		Assert.assertEquals(true, Columns.isNonPolar("SE-30"));
		Assert.assertEquals(true, Columns.isNonPolar("OV-101"));
		Assert.assertEquals(true, Columns.isNonPolar("OV-1"));
		Assert.assertEquals(true, Columns.isNonPolar("Methyl_Silicone"));
		Assert.assertEquals(true, Columns.isNonPolar("OV-1,_SE-30,_Methyl_silicone,_SP-2100,_OV-101,_DB-1,_etc."));
		Assert.assertEquals(true, Columns.isNonPolar("Polydimethyl_siloxanes"));
		Assert.assertEquals(true, Columns.isNonPolar("Cross-Linked_Methylsilicone"));
		Assert.assertEquals(true, Columns.isNonPolar("SF-96"));
		Assert.assertEquals(true, Columns.isNonPolar("GP_SP_2100_DB"));
		Assert.assertEquals(true, Columns.isNonPolar("Ultra-1_PONA"));
		Assert.assertEquals(true, Columns.isNonPolar("Optima_1"));
		Assert.assertEquals(true, Columns.isNonPolar("Se-30"));
		Assert.assertEquals(false, Columns.isNonPolar("5_%_Phenyl_methyl_siloxane"));
		Assert.assertEquals(false, Columns.isNonPolar("DB-5"));
		Assert.assertEquals(false, Columns.isNonPolar("HP-5"));
		Assert.assertEquals(false, Columns.isNonPolar("HP-5MS"));
		Assert.assertEquals(false, Columns.isNonPolar("Squalane"));
		Assert.assertEquals(false, Columns.isNonPolar("ZB-5"));
		Assert.assertEquals(false, Columns.isNonPolar("SLB-5_MS"));
		Assert.assertEquals(false, Columns.isNonPolar("Apiezon_M"));
		Assert.assertEquals(false, Columns.isNonPolar("Vacuum_Grease_Oil_(VM-4)"));
		Assert.assertEquals(false, Columns.isNonPolar("Polydimethyl_siloxane_with_5_%_Ph_groups"));
		Assert.assertEquals(false, Columns.isNonPolar("Polydimethyl_siloxane,_5_%_phenyl_groups"));
		Assert.assertEquals(false, Columns.isNonPolar("OV-22"));
		Assert.assertEquals(false, Columns.isNonPolar("RTV-502"));
		Assert.assertEquals(false, Columns.isNonPolar("PoraPLOT"));
		Assert.assertEquals(false, Columns.isNonPolar("Silicone_oil"));
		Assert.assertEquals(false, Columns.isNonPolar("NB-5"));
		Assert.assertEquals(false, Columns.isNonPolar(null));
		Assert.assertEquals(false, Columns.isNonPolar(" "));
		Assert.assertEquals(false, Columns.isNonPolar(""));
		Assert.assertEquals(false, Columns.isNonPolar("aaa"));
		Assert.assertEquals(true, Columns.isNonPolar("Other_non_polar"));
		Assert.assertEquals(false, Columns.isNonPolar("Other_semi_non_polar"));
	}

	public void testIsSemiNonPolarString() {
		Assert.assertEquals(false, Columns.isSemiNonPolar("DB-1"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("SE-30"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("OV-101"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("OV-1"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("Methyl_Silicone"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("OV-1,_SE-30,_Methyl_silicone,_SP-2100,_OV-101,_DB-1,_etc."));
		Assert.assertEquals(false, Columns.isSemiNonPolar("Polydimethyl_siloxanes"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("Cross-Linked_Methylsilicone"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("SF-96"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("GP_SP_2100_DB"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("Ultra-1_PONA"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("Optima_1"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("Se-30"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("5_%_Phenyl_methyl_siloxane"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("DB-5"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("HP-5"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("HP-5MS"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("Squalane"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("ZB-5"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("SLB-5_MS"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("Apiezon_M"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("Vacuum_Grease_Oil_(VM-4)"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("Polydimethyl_siloxane_with_5_%_Ph_groups"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("Polydimethyl_siloxane,_5_%_phenyl_groups"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("OV-22"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("RTV-502"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("PoraPLOT"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("Silicone_oil"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("NB-5"));
		Assert.assertEquals(false, Columns.isSemiNonPolar(null));
		Assert.assertEquals(false, Columns.isSemiNonPolar(" "));
		Assert.assertEquals(false, Columns.isSemiNonPolar(""));
		Assert.assertEquals(false, Columns.isSemiNonPolar("aaa"));
		Assert.assertEquals(false, Columns.isSemiNonPolar("Other_non_polar"));
		Assert.assertEquals(true, Columns.isSemiNonPolar("Other_semi_non_polar"));
	}

	public void testIsNonPolarInt() {
		Assert.assertEquals(false, Columns.isNonPolar(-1));
		Assert.assertEquals(false, Columns.isNonPolar(-2));
		for (int i = 0; i < 15; i++) {
			Assert.assertEquals(true, Columns.isNonPolar(i));
		}
		for (int i = 15; i < 150; i++) {
			Assert.assertEquals(false, Columns.isNonPolar(i));
		}
	}

	public void testIsSemiNonPolarInt() {
		for (int i = -20; i < 15; i++) {
			Assert.assertEquals(false, Columns.isSemiNonPolar(i));
		}
		for (int i = 15; i < 36; i++) {
			Assert.assertEquals(true, Columns.isSemiNonPolar(i));
		}
		for (int i = 36; i < 150; i++) {
			Assert.assertEquals(false, Columns.isSemiNonPolar(i));
		}
	}

	public void testColumn() {
		Assert.assertEquals("DB-1", Columns.column(0));
		Assert.assertEquals("SE-30", Columns.column(1));
		Assert.assertEquals("OV-101", Columns.column(2));
		Assert.assertEquals("OV-1", Columns.column(3));
		Assert.assertEquals("Methyl_Silicone", Columns.column(4));
		Assert.assertEquals("OV-1,_SE-30,_Methyl_silicone,_SP-2100,_OV-101,_DB-1,_etc.", Columns.column(13));
		Assert.assertEquals("Other_non_polar", Columns.column(14));
		Assert.assertEquals("5_%_Phenyl_methyl_siloxane", Columns.column(15));
		Assert.assertEquals("DB-5", Columns.column(16));
		Assert.assertEquals("HP-5", Columns.column(17));
		Assert.assertEquals("Squalane", Columns.column(20));
		Assert.assertEquals("ZB-5", Columns.column(33));
		Assert.assertEquals("SLB-5_MS", Columns.column(34));
		Assert.assertEquals("Other_semi_non_polar", Columns.column(35));
		Assert.assertEquals("Unknown", Columns.column(-1));
		Assert.assertEquals("Unknown", Columns.column(-2));
		Assert.assertEquals("Unknown", Columns.column(36));
	}

	public void testColumnOneHot() {
		for (int i = 0; i < 36; i++) {
			float[] r = Columns.columnOneHot(i);
			for (int j = 0; j < 36; j++) {
				Assert.assertEquals(i == j, r[j] == 1F);
				Assert.assertEquals(i != j, r[j] == 0F);
			}
			Assert.assertEquals(r.length, 36);
		}
		float[] r = Columns.columnOneHot(-1);
		for (int j = 0; j < 36; j++) {
			Assert.assertEquals(0.0F, r[j]);
		}
		r = Columns.columnOneHot(-2);
		for (int j = 0; j < 36; j++) {
			Assert.assertEquals(0.0F, r[j]);
		}
	}

	public void testColumnTypeOneHot() {
		for (int i = -100; i < 100; i++) {
			float[] r = Columns.columnTypeOneHot(i);
			Assert.assertEquals(Columns.isNonPolar(i), r[0]==1.0F);
			Assert.assertEquals(Columns.isSemiNonPolar(i), r[1]==1.0F);
			Assert.assertEquals(!Columns.isNonPolar(i), r[0]==0.0F);
			Assert.assertEquals(!Columns.isSemiNonPolar(i), r[1]==0.0F);
			Assert.assertEquals(r.length, 2);
		}
	}

	public void testColumnAndColumnTypeOneHot() {
		for (int i = -100; i < 100; i++) {
			float[] r = Columns.columnAndColumnTypeOneHot(i);
			Assert.assertEquals(Columns.isNonPolar(i), r[36]==1.0F);
			Assert.assertEquals(Columns.isSemiNonPolar(i), r[37]==1.0F);
			Assert.assertEquals(!Columns.isNonPolar(i), r[36]==0.0F);
			Assert.assertEquals(!Columns.isSemiNonPolar(i), r[37]==0.0F);
			Assert.assertEquals(r.length, 38);
		}
		for (int i = 0; i < 36; i++) {
			float[] r = Columns.columnAndColumnTypeOneHot(i);
			for (int j = 0; j < 36; j++) {
				Assert.assertEquals(i == j, r[j] == 1F);
				Assert.assertEquals(i != j, r[j] == 0F);
			}
			Assert.assertEquals(r.length, 38);
		}
		float[] r = Columns.columnAndColumnTypeOneHot(-1);
		for (int j = 0; j < 36; j++) {
			Assert.assertEquals(0.0F, r[j]);
		}
		r = Columns.columnAndColumnTypeOneHot(-2);
		for (int j = 0; j < 36; j++) {
			Assert.assertEquals(0.0F, r[j]);
		}
	}

}
