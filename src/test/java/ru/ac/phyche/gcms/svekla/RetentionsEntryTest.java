package ru.ac.phyche.gcms.svekla;

import junit.framework.TestCase;

import org.openscience.cdk.exception.CDKException;

import junit.framework.Assert;

public class RetentionsEntryTest extends TestCase {

	public void testGetSmiles() {
		RetentionsEntry e1 = RetentionsEntry.instance("CC(C)C", 10, 1);
		RetentionsEntry e2 = RetentionsEntry.instance(" CC(C)C  ", 10);
		Assert.assertEquals("CC(C)C", e1.getSmiles());
		Assert.assertEquals("CC(C)C", e2.getSmiles());
	}

	public void testGetInchi() throws CDKException {
		RetentionsEntry e1 = RetentionsEntry.instance("CC(C)C", 10, 1);
		RetentionsEntry e2 = RetentionsEntry.instance(" C(C)CC  ", 10);
		RetentionsEntry e3 = RetentionsEntry.instanceCanonical("CC(C)C", 10, 1, true);
		RetentionsEntry e4 = RetentionsEntry.instanceCanonical(" CCCC  ", 10, true);
		Assert.assertEquals("InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3", e1.getInchi());
		Assert.assertEquals("InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3", e2.getInchi());
		Assert.assertEquals("InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3", e3.getInchi());
		Assert.assertEquals("InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3", e4.getInchi());
	}

	public void testGetInchikey() throws CDKException {
		RetentionsEntry e1 = RetentionsEntry.instance("CC=CC", 10, 1);
		RetentionsEntry e2 = RetentionsEntry.instanceCanonical("C/C=C/C", 10, 1, true);
		RetentionsEntry e3 = RetentionsEntry.instanceCanonical("C/C=C\\C", 10, 1, true);
		RetentionsEntry e4 = RetentionsEntry.instanceCanonical("C/C=C/C", 10, 1, false);
		RetentionsEntry e5 = RetentionsEntry.instanceCanonical("C/C=C\\C", 10, 1, false);
		Assert.assertEquals("IAQRGUVFOMOMEM-UHFFFAOYSA-N", e1.getInchikey());
		Assert.assertEquals("IAQRGUVFOMOMEM-ONEGZZNKSA-N", e2.getInchikey());
		Assert.assertEquals("IAQRGUVFOMOMEM-ARJAWSKDSA-N", e3.getInchikey());
		Assert.assertEquals("IAQRGUVFOMOMEM-UHFFFAOYSA-N", e4.getInchikey());
		Assert.assertEquals("IAQRGUVFOMOMEM-UHFFFAOYSA-N", e5.getInchikey());
	}

	public void testSetSmiles() {
		RetentionsEntry e1 = RetentionsEntry.instance("CC(C)C", 10, 1);
		RetentionsEntry e2 = RetentionsEntry.instance(" CC(C)C  ", 10);
		e1.setSmiles(" C ");
		e2.setSmiles(" CC ");
		Assert.assertEquals("C", e1.getSmiles());
		Assert.assertEquals("CC", e2.getSmiles());
	}

	public void testGetRetention() throws CDKException {
		RetentionsEntry e1 = RetentionsEntry.instance("CC(C)C", 11, 1);
		RetentionsEntry e2 = RetentionsEntry.instance(" C(C)CC  ", 12);
		RetentionsEntry e3 = RetentionsEntry.instanceCanonical("CC(C)C", 13, 1, true);
		RetentionsEntry e4 = RetentionsEntry.instanceCanonical(" CCCC  ", 14, true);
		Assert.assertEquals(11.0F, e1.getRetention());
		Assert.assertEquals(12.0F, e2.getRetention());
		Assert.assertEquals(13.0F, e3.getRetention());
		Assert.assertEquals(14.0F, e4.getRetention());
	}

	public void testSetRetention() {
		RetentionsEntry e1 = RetentionsEntry.instance("CC(C)C", 10, 1);
		e1.setRetention(20);
		Assert.assertEquals(20.0F, e1.getRetention());
	}

	public void testGetColumnType() {
		RetentionsEntry e1 = RetentionsEntry.instance("CC(C)C", 10, 1);
		RetentionsEntry e2 = RetentionsEntry.instance(" CC(C)C  ", 10);
		Assert.assertEquals(1, e1.getColumnType());
		Assert.assertEquals(-1, e2.getColumnType());
	}

	public void testSetColumnType() {
		RetentionsEntry e2 = RetentionsEntry.instance(" CC(C)C  ", 10);
		e2.setColumnType(1);
		Assert.assertEquals(1, e2.getColumnType());
	}

	public void testInstanceCanonicalStringFloatIntBoolean() throws CDKException {
		RetentionsEntry e1 = RetentionsEntry.instanceCanonical(
				"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O", 11, 1, true);
		RetentionsEntry e2 = RetentionsEntry.instanceCanonical(
				"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O", 12, 1, false);
		RetentionsEntry e3 = RetentionsEntry.instanceCanonical(
				"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O", 13, true);
		RetentionsEntry e4 = RetentionsEntry.instanceCanonical(
				"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O", 14, false);
		Assert.assertEquals("CZMRCDWAGMRECN-UGDNZRGBSA-N", e1.getInchikey());
		Assert.assertEquals("CZMRCDWAGMRECN-UGDNZRGBSA-N", e3.getInchikey());
		Assert.assertEquals("CZMRCDWAGMRECN-UHFFFAOYSA-N", e2.getInchikey());
		Assert.assertEquals("CZMRCDWAGMRECN-UHFFFAOYSA-N", e4.getInchikey());
		Assert.assertEquals(11F, e1.getRetention());
		Assert.assertEquals(12F, e2.getRetention());
		Assert.assertEquals(13F, e3.getRetention());
		Assert.assertEquals(14F, e4.getRetention());
		Assert.assertEquals(1, e1.getColumnType());
		Assert.assertEquals(1, e2.getColumnType());
		Assert.assertEquals(-1, e3.getColumnType());
		Assert.assertEquals(-1, e4.getColumnType());
	}

	public void testClone() {
		RetentionsEntry e1 = RetentionsEntry.instance("CC(C)C", 10, 1);
		RetentionsEntry e2 = e1.clone();
		Assert.assertEquals(false, e1 == e2);
		Assert.assertEquals(false, e1.getSmiles() == e2.getSmiles());
		Assert.assertEquals(e1, e2);
		Assert.assertEquals(true, e1.equals(e2));
		Assert.assertEquals(true, e1.getSmiles().equals(e2.getSmiles()));
		Assert.assertEquals(e1.getSmiles(), e2.getSmiles());
		Assert.assertEquals(1, e2.getColumnType());
		Assert.assertEquals(10F, e2.getRetention());
		Assert.assertEquals(1, e1.getColumnType());
		Assert.assertEquals(10F, e1.getRetention());
	}

	public void testDeepclone() {
		RetentionsEntry e1 = RetentionsEntry.instance("CC(C)C", 10, 1);
		RetentionsEntry e2 = e1.deepclone();
		RetentionsEntry e3 = e2;
		Assert.assertEquals(false, e1 == e2);
		Assert.assertEquals(true, e3 == e2);
		Assert.assertEquals(false, e1.getSmiles() == e2.getSmiles());
		Assert.assertEquals(e1, e2);
		Assert.assertEquals(e3, e2);
	}

	public void testEqualsRetentionsEntry() {
		RetentionsEntry e1 = RetentionsEntry.instance("CC(C)C", 10, 1);
		RetentionsEntry e2 = e1.deepclone();
		Assert.assertEquals(true, e1.equals(e2));
		e2.setRetention(11F);
		Assert.assertEquals(false, e1.equals(e2));
		e2.setRetention(10F);
		Assert.assertEquals(true, e1.equals(e2));
		e2.setColumnType(2);
		Assert.assertEquals(false, e1.equals(e2));
		e2.setColumnType(1);
		Assert.assertEquals(true, e1.equals(e2));
		e2.setSmiles("C");
		Assert.assertEquals(false, e1.equals(e2));
		e2.setSmiles("          CC(C)C   ");
		Assert.assertEquals(true, e1.equals(e2));
	}

	public void testEqualsObject() {
		RetentionsEntry e1 = RetentionsEntry.instance("CC(C)C", 10, 1);
		Object e2 = e1.deepclone();
		Object e3 = new RetentionsDataset();
		Object e4 = "c";
		Assert.assertEquals(true, e1.equals(e2));
		Assert.assertEquals(true, e2.equals(e1));
		Assert.assertEquals(false, e1.equals(e3));
		Assert.assertEquals(false, e3.equals(e1));
		Assert.assertEquals(false, e1.equals(e4));
		Assert.assertEquals(false, e4.equals(e1));
	}

	public void testFingerprints() throws CDKException {
		RetentionsEntry e1 = RetentionsEntry.instance("CC(C)C", 10, 1);
		for (Chemoinformatics.FingerprintsType t : Chemoinformatics.FingerprintsType.values()) {
			float[] fp = e1.fingerprints(t);
			float[] fp2 = Chemoinformatics.fingerprints("CC(C)C", t);
			Assert.assertEquals(fp.length, fp2.length);
			for (int i = 0; i < fp.length; i++) {
				Assert.assertEquals(fp[i], fp2[i]);
			}
		}
	}

	public void testFuncGroups() throws CDKException {
		RetentionsEntry e1 = RetentionsEntry.instance(
				"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O", 10, 1);
		float[] g = e1.funcGroups();
		float[] g2 = Chemoinformatics
				.funcGroups("C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O");
		Assert.assertEquals(g.length, g2.length);
		for (int i = 0; i < g.length; i++) {
			Assert.assertEquals(g[i], g2[i]);
		}
	}

	public void testDepitction() throws CDKException {
		RetentionsEntry e1 = RetentionsEntry.instance(
				"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O", 10, 1);
		float[][] g = e1.depitction();
		float[][] g2 = Chemoinformatics
				.smilesToImage("C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O");
		Assert.assertEquals(g.length, g2.length);
		for (int i = 0; i < g.length; i++) {
			for (int j = 0; j < g[i].length; j++) {
				Assert.assertEquals(g[i][j], g2[i][j]);
			}
		}
	}

}
