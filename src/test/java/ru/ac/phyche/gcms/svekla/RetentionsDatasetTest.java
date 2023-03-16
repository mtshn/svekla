package ru.ac.phyche.gcms.svekla;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import org.openscience.cdk.exception.CDKException;

import junit.framework.Assert;
import junit.framework.TestCase;

public class RetentionsDatasetTest extends TestCase {
	private static final String[] descriptorNames = { "fragC", "C1SP1", "C2SP1", "C1SP2", "C2SP2", "C3SP2", "C1SP3",
			"C2SP3", "C3SP3", "C4SP3", "SCH-3", "SCH-4", "SCH-5", "SCH-6", "SCH-7", "VCH-3", "VCH-4", "VCH-5", "VCH-6",
			"VCH-7", "nAtomLC", "BCUTw-1l", "ATSc2", "ATSc1", "BCUTp-1l", "BCUTc-1h", "ATSc5", "BCUTp-1h", "ATSc4",
			"ATSc3", "BCUTw-1h", "BCUTc-1l", "Kier3" };

	private RetentionsEntry[] a_array() {
		RetentionsEntry a1 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 10, 1);
		RetentionsEntry a2 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 10, 1);
		RetentionsEntry a3 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 12, 1);
		RetentionsEntry a4 = RetentionsEntry.instance("COc1cc(ccc1O)C=O", 14, 1);
		RetentionsEntry a5 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 12, 2);
		RetentionsEntry a6 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 13, 2);
		RetentionsEntry a7 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 10, 3);
		RetentionsEntry a8 = RetentionsEntry.instance("COc1cc(ccc1O)C=O", 15, 3);
		RetentionsEntry a9 = RetentionsEntry.instance("CCCC", 5, 3);
		RetentionsEntry a10 = RetentionsEntry.instance("CCC(C)", 5.5F, 3);
		RetentionsEntry a11 = RetentionsEntry.instance("CCCC", 3.5F, 4);
		RetentionsEntry a12 = RetentionsEntry.instance("C/C=C/C", 7.5F, 1);
		RetentionsEntry a13 = RetentionsEntry.instance("C/C=C/C", 9.5F, 2);
		RetentionsEntry a14 = RetentionsEntry.instance("C/C=C/C", 8.5F, 1);
		RetentionsEntry a15 = RetentionsEntry.instance("CC=CC", 9.5F, 3);
		RetentionsEntry a16 = RetentionsEntry.instance("CC=CC", 8.3F, 1);
		RetentionsEntry a17 = RetentionsEntry.instance("C/C=C\\C", 9.5F, 2);
		RetentionsEntry a18 = RetentionsEntry.instance("C/C=C\\C", 8.3F, 4);
		RetentionsEntry a19 = RetentionsEntry.instance("C/C=C\\C", 8.3F, 3);
		RetentionsEntry a20 = RetentionsEntry.instance("C/C=C\\C", 8.3F, 3);
		RetentionsEntry a21 = RetentionsEntry.instance("C/C=C/CCCCCC(C)CCC", 7.51F, 1);
		RetentionsEntry a22 = RetentionsEntry.instance("C/C=C/CCCCCC(C)CCC", 9.5F, 2);
		RetentionsEntry a23 = RetentionsEntry.instance("C/C=C/CCCCCC(C)CCC", 8.51F, 1);
		RetentionsEntry a24 = RetentionsEntry.instance("CC=CCCCCCC(C)CCC", 9.5F, 3);
		RetentionsEntry a25 = RetentionsEntry.instance("CC=CCCCCCC(C)CCC", 8.31F, 1);
		RetentionsEntry a26 = RetentionsEntry.instance("C/C=C\\CCCCCC(C)CCC", 9.5F, 2);
		RetentionsEntry a27 = RetentionsEntry.instance("C/C=C\\CCCCCC(C)CCC", 8.31F, 4);
		RetentionsEntry a28 = RetentionsEntry.instance("C/C=C\\CCCCCC(C)CCC", 8.31F, 2);
		RetentionsEntry a29 = RetentionsEntry.instance("C/C=C\\CCCCCC(C)CCC", 8.3F, 3);
		RetentionsEntry a30 = RetentionsEntry.instance("C(=N)(N)O", 83F, 5);
		return new RetentionsEntry[] { a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18,
				a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30 };
	}

	private ArrayList<RetentionsEntry> a_list() {
		ArrayList<RetentionsEntry> result = new ArrayList<RetentionsEntry>(Arrays.asList(a_array()));
		return result;
	}

	private RetentionsDataset a() {
		return RetentionsDataset.create(a_list());
	}

	private RetentionsDataset b() {
		RetentionsEntry b1 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 14, 1);
		RetentionsEntry b2 = RetentionsEntry.instance("C(C)C(C)", 5.7F, 3);
		RetentionsEntry b3 = RetentionsEntry.instance("C/C=C/C", 7.5F, 1);
		RetentionsEntry b4 = RetentionsEntry.instance("CC=CCCCCCC(C)CCC", 8.31F, 1);
		RetentionsEntry b5 = RetentionsEntry.instance("CCCCC(=N)(N)O", 83F, 5);
		RetentionsEntry b6 = RetentionsEntry.instance("COc1cc(ccc1O)C=O", 83F, 5);
		return RetentionsDataset.create(new RetentionsEntry[] { b1, b2, b3, b4, b5, b6 });
	}

	private RetentionsDataset c() {
		RetentionsEntry b1 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 14, 1);
		RetentionsEntry b2 = RetentionsEntry.instance("C", 5.7F, 3);
		RetentionsEntry b3 = RetentionsEntry.instance("CC", 7.5F, 1);
		RetentionsEntry b4 = RetentionsEntry.instance("CCC", 8.31F, 1);
		RetentionsEntry b5 = RetentionsEntry.instance("C(=N)(N)O", 83F, 5);
		RetentionsEntry b6 = RetentionsEntry.instance("COc1cc(ccc1O)C=O", 83F, 5);
		return RetentionsDataset.create(new RetentionsEntry[] { b1, b2, b3, b4, b5, b6 });
	}

	private RetentionsDataset d() {
		RetentionsEntry b1 = RetentionsEntry.instance("CO", 5.7F, 3);
		RetentionsEntry b2 = RetentionsEntry.instance("CCO", 7.5F, 1);
		RetentionsEntry b3 = RetentionsEntry.instance("CCCO", 8.31F, 1);
		return RetentionsDataset.create(new RetentionsEntry[] { b1, b2, b3 });
	}

	public void testSize() {
		Assert.assertEquals(30, a().size());
		Assert.assertEquals(30, a().copy().size());
	}

	public void testGetData() {
		RetentionsEntry a[] = a_array();
		RetentionsEntry b[] = a().getData();
		for (int i = 0; i < a.length; i++) {
			Assert.assertEquals(true, a[i].equals(b[i]));
			Assert.assertEquals(true, a[i].getSmiles().equals(b[i].getSmiles()));
			Assert.assertEquals(true, a[i].getRetention() == b[i].getRetention());
		}
	}

	public void testSetData() {
		RetentionsEntry a1 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 10, 1);
		RetentionsEntry a2 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 10, 1);
		RetentionsEntry a3 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 12, 1);
		RetentionsDataset b = RetentionsDataset.create(new RetentionsEntry[] { a1, a2, a3 });
		RetentionsDataset a = a();
		Assert.assertEquals(30, a.size());
		Assert.assertEquals(3, b.size());
		a.setData(b.getData());
		Assert.assertEquals(3, a.size());
		for (int i = 0; i < a.size(); i++) {
			Assert.assertEquals(true, a.getEntry(i) == b.getEntry(i));
		}
	}

	public void testGetEntry() {
		RetentionsDataset a = a();
		for (int i = 0; i < a.size(); i++) {
			Assert.assertEquals(true, a.getEntry(i).equals(a_array()[i]));
			Assert.assertEquals(true, a.getEntry(i).getColumnType() == a_array()[i].getColumnType());
			Assert.assertEquals(true, a.getEntry(i).getRetention() == a_array()[i].getRetention());
			Assert.assertEquals(true, a.getEntry(i).getSmiles().equals(a_array()[i].getSmiles()));
			Assert.assertEquals(false, a.getEntry(i) == (a_array()[i]));
		}
		Assert.assertEquals(true, a.getEntry(0).equals(RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 10, 1)));
	}

	public void testGetSmiles() {
		RetentionsDataset a = a();
		for (int i = 0; i < a.size(); i++) {
			Assert.assertEquals(true, a.getSmiles(i).equals(a_array()[i].getSmiles()));
			Assert.assertEquals(true, a.getSmiles(i).equals(a.getEntry(i).getSmiles()));
		}
		Assert.assertEquals(true, a.getSmiles(0).equals("COC1=C(C=CC(=C1)C=O)O"));
	}

	public void testGetRetention() {
		RetentionsDataset a = a();
		for (int i = 0; i < a.size(); i++) {
			Assert.assertEquals(true, a.getRetention(i) == a_array()[i].getRetention());
			Assert.assertEquals(true, a.getRetention(i) == a.getEntry(i).getRetention());
		}
		Assert.assertEquals(true, a.getRetention(0) == 10F);
	}

	public void testGetColumn() {
		RetentionsDataset a = a();
		for (int i = 0; i < a.size(); i++) {
			Assert.assertEquals(true, a.getColumn(i) == a_array()[i].getColumnType());
			Assert.assertEquals(true, a.getColumn(i) == a.getEntry(i).getColumnType());
		}
		Assert.assertEquals(true, a.getColumn(0) == 1);
	}

	public void testSetEntry() {
		RetentionsDataset a = a();
		RetentionsEntry b = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 10, 1);
		Assert.assertEquals(true, a.getEntry(0).equals(b));
		Assert.assertEquals(false, a.getEntry(0) == b);
		Assert.assertEquals(true, a.getEntry(1).equals(b));
		Assert.assertEquals(false, a.getEntry(2).equals(b));
		a.setEntry(b, 0);
		Assert.assertEquals(true, a.getEntry(0).equals(b));
		Assert.assertEquals(true, a.getEntry(0) == b);
		Assert.assertEquals(false, a.getEntry(2).equals(b));
		a.setEntry(b, 2);
		Assert.assertEquals(true, a.getEntry(0).equals(b));
		Assert.assertEquals(true, a.getEntry(0) == b);
		Assert.assertEquals(true, a.getEntry(2).equals(b));
	}

	public void testCreateRetentionsEntryArray() {
		RetentionsEntry[] p = a_array();
		RetentionsDataset a1 = RetentionsDataset.create(p);
		Assert.assertEquals(p.length, a1.size());
		Assert.assertEquals(p.length, a1.getData().length);
		for (int i = 0; i < p.length; i++) {
			Assert.assertEquals(true, p[i] == a1.getEntry(i));
		}
	}

	public void testCreateArrayListOfRetentionsEntry() {
		ArrayList<RetentionsEntry> p = a_list();
		RetentionsDataset a1 = RetentionsDataset.create(p);
		Assert.assertEquals(p.size(), a1.size());
		Assert.assertEquals(p.size(), a1.getData().length);
		for (int i = 0; i < p.size(); i++) {
			Assert.assertEquals(true, p.get(i) == a1.getEntry(i));
		}
	}

	public void testFingerprints() throws CDKException {
		RetentionsDataset a = a();
		for (Chemoinformatics.FingerprintsType t : Chemoinformatics.FingerprintsType.values()) {
			for (int i = 0; i < a.size(); i++) {
				float fp[] = a.fingerprints(t, i);
				float fp2[] = Chemoinformatics.fingerprints(a_array()[i].getSmiles(), t);
				Assert.assertEquals(fp.length, fp2.length);
				for (int j = 0; j < fp.length; j++) {
					Assert.assertEquals(fp[j], fp2[j]);
				}
			}
		}
	}

	public void testDescriptors() throws CDKException {
		Descriptors d = Descriptors.instance(descriptorNames);
		RetentionsDataset a = a();
		d.precompute(a.compoundsCanonical(true), true);
		for (int i = 0; i < a.size(); i++) {
			float desc[] = a.descriptors(i, d);
			float desc2[] = Chemoinformatics.descriptors(Chemoinformatics.canonical(a_array()[i].getSmiles(), true),
					descriptorNames, d.getMinMaxArray().getLeft(), d.getMinMaxArray().getRight());
			Assert.assertEquals(desc.length, desc2.length);
			for (int j = 0; j < desc.length; j++) {
				Assert.assertEquals(desc[j], desc2[j]);
			}
		}
	}

	public void testDescriptorsNoNaNs() throws CDKException {
		Descriptors d = Descriptors.instance(descriptorNames);
		RetentionsDataset a = a();
		d.precompute(a.compoundsCanonical(true), true);
		for (int i = 0; i < a.size(); i++) {
			float desc[] = a.descriptorsNoNaNs(i, d);
			float desc2[] = d.getNoNaNs(Chemoinformatics.canonical(a.getSmiles(i), true));
			Assert.assertEquals(desc.length, desc2.length);
			for (int j = 0; j < desc.length; j++) {
				Assert.assertEquals(desc[j], desc2[j]);
			}
		}
	}

	public void testFuncGroups() throws CDKException {
		RetentionsDataset a = a();
		for (int i = 0; i < a.size(); i++) {
			float f[] = a.funcGroups(i);
			float f2[] = Chemoinformatics.funcGroups(a_array()[i].getSmiles());
			Assert.assertEquals(f.length, f2.length);
			for (int j = 0; j < f.length; j++) {
				Assert.assertEquals(f[j], f2[j]);
			}
		}
	}

	public void testDepiction() throws CDKException {
		RetentionsDataset a = a();
		for (int i = 0; i < a.size(); i++) {
			float f[][] = a.depiction(i);
			float f2[][] = Chemoinformatics.smilesToImage(a_array()[i].getSmiles());
			Assert.assertEquals(f.length, f2.length);
			for (int j = 0; j < f.length; j++) {
				for (int k = 0; k < f[j].length; k++) {
					Assert.assertEquals(f[j][k], f2[j][k]);
				}
			}
		}
	}

	public void testExtDescriptors() throws CDKException {
		Descriptors d = Descriptors.instance(descriptorNames);
		RetentionsDataset a = a();
		d.precompute(a.compoundsCanonical(true), true);
		for (int i = 0; i < a.size(); i++) {
			float desc[] = a.extDescriptors(i, d);
			float desc2[] = d.getNoNaNs(Chemoinformatics.canonical(a.getSmiles(i), true));
			float[] func = Chemoinformatics.funcGroups(a_array()[i].getSmiles());
			Assert.assertEquals(desc.length, desc2.length + func.length);
			for (int j = 0; j < desc2.length; j++) {
				Assert.assertEquals(desc[j], desc2[j]);
			}
			for (int j = 0; j < func.length; j++) {
				Assert.assertEquals(desc[j + desc2.length], func[j]);
			}
		}
	}

	public void testExtFingerprints() throws CDKException {
		RetentionsDataset a = a();
		for (int i = 0; i < a.size(); i++) {
			float fp[] = a.extFingerprints(Chemoinformatics.FingerprintsType.ADDITIVE_CIRCULAR_4_1024_NO_SCALE, i);
			float fp2[] = Chemoinformatics.fingerprints(a_array()[i].getSmiles(),
					Chemoinformatics.FingerprintsType.ADDITIVE_CIRCULAR_4_1024_NO_SCALE);
			float[] func = Chemoinformatics.funcGroups(a_array()[i].getSmiles());
			Assert.assertEquals(fp.length, fp2.length + func.length);
			for (int j = 0; j < fp2.length; j++) {
				Assert.assertEquals(fp[j], fp2[j]);
			}
			for (int j = 0; j < func.length; j++) {
				Assert.assertEquals(fp[j + fp2.length], func[j]);
			}
		}
	}

	public void testCopy() {
		RetentionsDataset a = a();
		RetentionsDataset b = a.copy();
		RetentionsDataset c = a;
		Assert.assertEquals(false, a == b);
		Assert.assertEquals(true, a == c);
		Assert.assertEquals(true, a.size() == b.size());
		for (int i = 0; i < a.size(); i++) {
			Assert.assertEquals(false, a.getEntry(i) == b.getEntry(i));
			Assert.assertEquals(true, a.getEntry(i).equals(b.getEntry(i)));
			Assert.assertEquals(true, a.getSmiles(i).equals(b.getSmiles(i)));
			Assert.assertEquals(true, a.getRetention(i) == b.getRetention(i));
			Assert.assertEquals(true, a.getColumn(i) == b.getColumn(i));
		}
	}

	public void testShuffle() throws CDKException {
		RetentionsDataset a = a();
		RetentionsDataset b = a.copy();
		a.shuffle();
		Assert.assertEquals(a.size(), b.size());
		HashSet<String> ca = a.compounds();
		HashSet<String> cb = b.compounds();
		Assert.assertEquals(11, ca.size());
		Assert.assertEquals(11, cb.size());
		ca.removeAll(cb);
		Assert.assertEquals(0, ca.size());
		ca = a.inchiIds();
		cb = b.inchiIds();
		Assert.assertEquals(9, ca.size());
		Assert.assertEquals(9, cb.size());
		ca.removeAll(cb);
		Assert.assertEquals(0, ca.size());
		boolean t = true;
		int v = 0;
		for (int i = 0; i < a.size(); i++) {
			if (a.getSmiles(i).equals("C(=N)(N)O")) {
				Assert.assertEquals(83F, a.getRetention(i));
				Assert.assertEquals(5, a.getColumn(i));
			}
			if (!a.getEntry(i).equals(b.getEntry(i))) {
				t = false;
				v++;
			}
		}
		Assert.assertEquals(false, t);
		Assert.assertEquals(true, v > 10);
	}

	public void testCompoundsCanonical() throws CDKException {
		RetentionsDataset a = a();
		HashSet<String> hs = a.compoundsCanonical(true);
		String strs[] = new String[] { "CC=CC", "C/C=C\\CCCCCC(C)CCC", "C/C=C/C", "CC=CCCCCCC(C)CCC", "C/C=C\\C",
				"COc1cc(ccc1O)C=O", "C/C=C/CCCCCC(C)CCC", "C(=N)(N)O", "CCCC" };
		HashSet<String> hs2 = new HashSet<String>();
		for (String s : strs) {
			Assert.assertEquals(true, hs.contains(s));
			hs2.add(s);
		}
		Assert.assertEquals(strs.length, hs.size());
		Assert.assertEquals(strs.length, hs2.size());
		hs2.removeAll(hs);
		Assert.assertEquals(0, hs2.size());

		hs = a.compoundsCanonical(false);
		strs = new String[] { "CC=CC", "CC=CCCCCCC(C)CCC", "COc1cc(ccc1O)C=O", "C(=N)(N)O", "CCCC" };
		hs2 = new HashSet<String>();
		for (String s : strs) {
			Assert.assertEquals(true, hs.contains(s));
			hs2.add(s);
		}
		Assert.assertEquals(strs.length, hs.size());
		Assert.assertEquals(strs.length, hs2.size());
		hs2.removeAll(hs);
		Assert.assertEquals(0, hs2.size());
	}

	public void testCompounds() {
		RetentionsDataset a = a();
		HashSet<String> hs = a.compounds();
		String strs[] = new String[] { "CC=CC", "C/C=C\\CCCCCC(C)CCC", "C/C=C/C", "CC=CCCCCCC(C)CCC", "C/C=C\\C",
				"COc1cc(ccc1O)C=O", "C/C=C/CCCCCC(C)CCC", "C(=N)(N)O", "CCCC", "CCC(C)", "COC1=C(C=CC(=C1)C=O)O" };
		HashSet<String> hs2 = new HashSet<String>();
		for (String s : strs) {
			Assert.assertEquals(true, hs.contains(s));
			hs2.add(s);
		}
		Assert.assertEquals(strs.length, hs.size());
		Assert.assertEquals(strs.length, hs2.size());
		hs2.removeAll(hs);
		Assert.assertEquals(0, hs2.size());
	}

	public void testInchiIds() throws CDKException {
		RetentionsDataset a = a();
		HashSet<String> hs = a.inchiIds();
		String strs[] = new String[] { "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3",
				"InChI=1S/C13H26/c1-4-6-7-8-9-10-12-13(3)11-5-2/h4,6,13H,5,7-12H2,1-3H3",
				"InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3", "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+",
				"InChI=1S/C13H26/c1-4-6-7-8-9-10-12-13(3)11-5-2/h4,6,13H,5,7-12H2,1-3H3/b6-4+",
				"InChI=1S/C13H26/c1-4-6-7-8-9-10-12-13(3)11-5-2/h4,6,13H,5,7-12H2,1-3H3/b6-4-",
				"InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)", "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3-",
				"InChI=1S/C8H8O3/c1-11-8-4-6(5-9)2-3-7(8)10/h2-5,10H,1H3" };
		HashSet<String> hs2 = new HashSet<String>();
		for (String s : strs) {
			Assert.assertEquals(true, hs.contains(s));
			hs2.add(s);
		}
		Assert.assertEquals(strs.length, hs.size());
		Assert.assertEquals(strs.length, hs2.size());
		hs2.removeAll(hs);
		Assert.assertEquals(0, hs2.size());
	}

	public void testInchiKeys() throws CDKException {
		RetentionsDataset a = a();
		HashSet<String> hs = a.inchiKeys();
		String strs[] = new String[] { "IAQRGUVFOMOMEM-UHFFFAOYSA-N", "IJDNQMDRQITEOD-UHFFFAOYSA-N",
				"SWONXRASPPQSKN-UHFFFAOYSA-N", "SWONXRASPPQSKN-XQRVVYSFSA-N", "XSQUKJJJFZCRTK-UHFFFAOYSA-N",
				"SWONXRASPPQSKN-GQCTYLIASA-N", "IAQRGUVFOMOMEM-ARJAWSKDSA-N", "MWOOGOJBHIARFG-UHFFFAOYSA-N",
				"IAQRGUVFOMOMEM-ONEGZZNKSA-N" };
		HashSet<String> hs2 = new HashSet<String>();
		for (String s : strs) {
			Assert.assertEquals(true, hs.contains(s));
			hs2.add(s);
		}
		Assert.assertEquals(strs.length, hs.size());
		Assert.assertEquals(strs.length, hs2.size());
		hs2.removeAll(hs);
		Assert.assertEquals(0, hs2.size());
	}

	public void testCountIdenticalByInchiRetentionsDataset() throws CDKException {
		Assert.assertEquals(4, a().countIdenticalByInchi(b()));
		Assert.assertEquals(4, b().countIdenticalByInchi(a()));
		Assert.assertEquals(2, a().countIdenticalByInchi(c()));
		Assert.assertEquals(2, c().countIdenticalByInchi(a()));
		Assert.assertEquals(0, a().countIdenticalByInchi(d()));
		Assert.assertEquals(0, d().countIdenticalByInchi(a()));
		Assert.assertEquals(1, b().countIdenticalByInchi(c()));
		Assert.assertEquals(1, c().countIdenticalByInchi(b()));
		Assert.assertEquals(0, b().countIdenticalByInchi(d()));
		Assert.assertEquals(0, d().countIdenticalByInchi(b()));
	}

	public void testCountIdenticalByInchikeys() throws CDKException {
		Assert.assertEquals(4, a().countIdenticalByInchikeys(b()));
		Assert.assertEquals(4, b().countIdenticalByInchikeys(a()));
		Assert.assertEquals(2, a().countIdenticalByInchikeys(c()));
		Assert.assertEquals(2, c().countIdenticalByInchikeys(a()));
		Assert.assertEquals(0, a().countIdenticalByInchikeys(d()));
		Assert.assertEquals(0, d().countIdenticalByInchikeys(a()));
		Assert.assertEquals(1, b().countIdenticalByInchikeys(c()));
		Assert.assertEquals(1, c().countIdenticalByInchikeys(b()));
		Assert.assertEquals(0, b().countIdenticalByInchikeys(d()));
		Assert.assertEquals(0, d().countIdenticalByInchikeys(b()));
	}

	public void testCountIdenticalByCanonicalSmiles() throws CDKException {
		Assert.assertEquals(4, a().countIdenticalByCanonicalSmiles(b()));
		Assert.assertEquals(4, b().countIdenticalByCanonicalSmiles(a()));
		Assert.assertEquals(2, a().countIdenticalByCanonicalSmiles(c()));
		Assert.assertEquals(2, c().countIdenticalByCanonicalSmiles(a()));
		Assert.assertEquals(0, a().countIdenticalByCanonicalSmiles(d()));
		Assert.assertEquals(0, d().countIdenticalByCanonicalSmiles(a()));
		Assert.assertEquals(1, b().countIdenticalByCanonicalSmiles(c()));
		Assert.assertEquals(1, c().countIdenticalByCanonicalSmiles(b()));
		Assert.assertEquals(0, b().countIdenticalByCanonicalSmiles(d()));
		Assert.assertEquals(0, d().countIdenticalByCanonicalSmiles(b()));
	}

	public void testCountIdenticalByInchiRetentionsDatasetArrayRetentionsDatasetArray() throws CDKException {
		int[][] q = RetentionsDataset.countIdenticalByInchi(new RetentionsDataset[] { a(), b(), c() },
				new RetentionsDataset[] { a(), b(), c(), d() });
		Assert.assertEquals(9, q[0][0]);
		Assert.assertEquals(4, q[0][1]);
		Assert.assertEquals(2, q[0][2]);
		Assert.assertEquals(0, q[0][3]);
		Assert.assertEquals(4, q[1][0]);
		Assert.assertEquals(5, q[1][1]);
		Assert.assertEquals(1, q[1][2]);
		Assert.assertEquals(0, q[1][3]);
		Assert.assertEquals(2, q[2][0]);
		Assert.assertEquals(1, q[2][1]);
		Assert.assertEquals(5, q[2][2]);
		Assert.assertEquals(0, q[2][3]);
		Assert.assertEquals(3, q.length);
		Assert.assertEquals(4, q[0].length);
	}

	public void testCompoundsBasedSplitAndShuffleInt() throws CDKException {
		RetentionsDataset a = a();
		RetentionsDataset aa = a();
		RetentionsDataset b = a.compoundsBasedSplitAndShuffle(3);
		Assert.assertEquals(3, b.compounds().size());
		Assert.assertEquals(8, a.compounds().size());
		Assert.assertEquals(aa.size(), a.size(), b.size());
		for (int i = 0; i < aa.size(); i++) {
			boolean found = false;
			int foundInDataset = 0;
			for (int j = 0; j < a.size(); j++) {
				if (aa.getEntry(i).equals(a.getEntry(j))) {
					found = true;
				}
			}
			if (found) {
				foundInDataset++;
			}
			found = false;
			for (int k = 0; k < b.size(); k++) {
				if (aa.getEntry(i).equals(b.getEntry(k))) {
					found = true;
				}
			}
			if (found) {
				foundInDataset++;
			}
			Assert.assertEquals(1, foundInDataset);
		}

		a = a();
		a.makeCanoncalAll(true);
		aa = a();
		aa.makeCanoncalAll(true);
		b = a.compoundsBasedSplitAndShuffle(3);
		Assert.assertEquals(3, b.compounds().size());
		Assert.assertEquals(6, a.compounds().size());
		Assert.assertEquals(3, b.compoundsCanonical(true).size());
		Assert.assertEquals(6, a.compoundsCanonical(true).size());
		Assert.assertEquals(0, a.countIdenticalByInchikeys(b));
		Assert.assertEquals(0, b.countIdenticalByInchikeys(a));
		Assert.assertEquals(0, a.countIdenticalByInchi(b));
		Assert.assertEquals(0, b.countIdenticalByInchi(a));
		Assert.assertEquals(aa.size(), a.size(), b.size());

		a = a();
		a.makeCanoncalAll(false);
		aa = a();
		aa.makeCanoncalAll(false);
		b = a.compoundsBasedSplitAndShuffle(3);
		Assert.assertEquals(3, b.compounds().size());
		Assert.assertEquals(2, a.compounds().size());
		Assert.assertEquals(3, b.compoundsCanonical(true).size());
		Assert.assertEquals(2, a.compoundsCanonical(true).size());
		Assert.assertEquals(0, a.countIdenticalByInchikeys(b));
		Assert.assertEquals(0, b.countIdenticalByInchikeys(a));
		Assert.assertEquals(0, a.countIdenticalByInchi(b));
		Assert.assertEquals(0, b.countIdenticalByInchi(a));
		Assert.assertEquals(aa.size(), a.size(), b.size());
	}

	public void testCompoundsBasedSplitAndShuffleFloat() throws CDKException {
		RetentionsDataset a = a();
		RetentionsDataset aa = a();
		RetentionsDataset b = a.compoundsBasedSplitAndShuffle(3F / 11F);
		Assert.assertEquals(3, b.compounds().size());
		Assert.assertEquals(8, a.compounds().size());
		Assert.assertEquals(aa.size(), a.size(), b.size());
		for (int i = 0; i < aa.size(); i++) {
			boolean found = false;
			int foundInDataset = 0;
			for (int j = 0; j < a.size(); j++) {
				if (aa.getEntry(i).equals(a.getEntry(j))) {
					found = true;
				}
			}
			if (found) {
				foundInDataset++;
			}
			found = false;
			for (int k = 0; k < b.size(); k++) {
				if (aa.getEntry(i).equals(b.getEntry(k))) {
					found = true;
				}
			}
			if (found) {
				foundInDataset++;
			}
			Assert.assertEquals(1, foundInDataset);
		}

		a = a();
		a.makeCanoncalAll(true);
		aa = a();
		aa.makeCanoncalAll(true);
		b = a.compoundsBasedSplitAndShuffle(1F / 3F);
		Assert.assertEquals(3, b.compounds().size());
		Assert.assertEquals(6, a.compounds().size());
		Assert.assertEquals(3, b.compoundsCanonical(true).size());
		Assert.assertEquals(6, a.compoundsCanonical(true).size());
		Assert.assertEquals(0, a.countIdenticalByInchikeys(b));
		Assert.assertEquals(0, b.countIdenticalByInchikeys(a));
		Assert.assertEquals(0, a.countIdenticalByInchi(b));
		Assert.assertEquals(0, b.countIdenticalByInchi(a));
		Assert.assertEquals(aa.size(), a.size(), b.size());

		a = a();
		a.makeCanoncalAll(false);
		aa = a();
		aa.makeCanoncalAll(false);
		b = a.compoundsBasedSplitAndShuffle(0.6F);
		Assert.assertEquals(3, b.compounds().size());
		Assert.assertEquals(2, a.compounds().size());
		Assert.assertEquals(3, b.compoundsCanonical(true).size());
		Assert.assertEquals(2, a.compoundsCanonical(true).size());
		Assert.assertEquals(0, a.countIdenticalByInchikeys(b));
		Assert.assertEquals(0, b.countIdenticalByInchikeys(a));
		Assert.assertEquals(0, a.countIdenticalByInchi(b));
		Assert.assertEquals(0, b.countIdenticalByInchi(a));
		Assert.assertEquals(aa.size(), a.size(), b.size());
	}

	public void testSimpleSplitInt() {
		RetentionsEntry a1 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 10, 1);
		RetentionsEntry a2 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 10, 1);
		RetentionsEntry a3 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 12, 1);
		RetentionsEntry a4 = RetentionsEntry.instance("COc1cc(ccc1O)C=O", 14, 1);
		RetentionsEntry a5 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 12, 2);
		RetentionsDataset a = a();
		RetentionsDataset aa = a();
		RetentionsDataset b = a.simpleSplit(5);
		Assert.assertEquals(a1, b.getEntry(0));
		Assert.assertEquals(a2, b.getEntry(1));
		Assert.assertEquals(a3, b.getEntry(2));
		Assert.assertEquals(a4, b.getEntry(3));
		Assert.assertEquals(a5, b.getEntry(4));
		Assert.assertEquals(5, b.size());
		Assert.assertEquals(25, a.size());
		for (int i = 0; i < a.size(); i++) {
			Assert.assertEquals(a.getEntry(i), aa.getEntry(i + 5));
		}
	}

	public void testSimpleSplitFloat() {
		RetentionsEntry a1 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 10, 1);
		RetentionsEntry a2 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 10, 1);
		RetentionsEntry a3 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 12, 1);
		RetentionsEntry a4 = RetentionsEntry.instance("COc1cc(ccc1O)C=O", 14, 1);
		RetentionsEntry a5 = RetentionsEntry.instance("COC1=C(C=CC(=C1)C=O)O", 12, 2);
		RetentionsDataset a = a();
		RetentionsDataset aa = a();
		RetentionsDataset b = a.simpleSplit(1F / 6F);
		Assert.assertEquals(a1, b.getEntry(0));
		Assert.assertEquals(a2, b.getEntry(1));
		Assert.assertEquals(a3, b.getEntry(2));
		Assert.assertEquals(a4, b.getEntry(3));
		Assert.assertEquals(a5, b.getEntry(4));
		Assert.assertEquals(5, b.size());
		Assert.assertEquals(25, a.size());
		for (int i = 0; i < a.size(); i++) {
			Assert.assertEquals(a.getEntry(i), aa.getEntry(i + 5));
		}
	}

	public void testSimpleShuffleSplitInt() {
		RetentionsDataset a = b();
		RetentionsDataset aa = b();
		RetentionsDataset b = a.simpleShuffleSplit(2);
		Assert.assertEquals(aa.size(), a.size(), b.size());
		for (int i = 0; i < aa.size(); i++) {
			boolean found = false;
			int foundInDataset = 0;
			for (int j = 0; j < a.size(); j++) {
				if (aa.getEntry(i).equals(a.getEntry(j))) {
					found = true;
				}
			}
			if (found) {
				foundInDataset++;
			}
			found = false;
			for (int k = 0; k < b.size(); k++) {
				if (aa.getEntry(i).equals(b.getEntry(k))) {
					found = true;
				}
			}
			if (found) {
				foundInDataset++;
			}
			Assert.assertEquals(1, foundInDataset);
		}
	}

	public void testSimpleShuffleSplitFloat() {
		RetentionsDataset a = b();
		RetentionsDataset aa = b();
		RetentionsDataset b = a.simpleShuffleSplit(3F / 8F);
		Assert.assertEquals(aa.size(), a.size(), b.size());
		for (int i = 0; i < aa.size(); i++) {
			boolean found = false;
			int foundInDataset = 0;
			for (int j = 0; j < a.size(); j++) {
				if (aa.getEntry(i).equals(a.getEntry(j))) {
					found = true;
				}
			}
			if (found) {
				foundInDataset++;
			}
			found = false;
			for (int k = 0; k < b.size(); k++) {
				if (aa.getEntry(i).equals(b.getEntry(k))) {
					found = true;
				}
			}
			if (found) {
				foundInDataset++;
			}
			Assert.assertEquals(1, foundInDataset);
		}
	}

	public void testFilterIdentical() throws CDKException {
		RetentionsDataset[] first = new RetentionsDataset[] { a(), b(), c() };
		RetentionsDataset[] second = new RetentionsDataset[] { a(), b(), c(), d() };
		for (int i = 0; i < first.length; i++) {
			for (int j = 0; j < second.length; j++) {
				int a = first[i].countIdenticalByInchi(second[j]);
				RetentionsDataset copy = first[i].copy();
				copy.filterIdentical(second[j]);
				Assert.assertEquals(0, copy.countIdenticalByInchi(second[j]));
				Assert.assertEquals(0, second[j].countIdenticalByInchi(copy));
				Assert.assertEquals(a, first[i].countIdenticalByInchi(second[j]));
				Assert.assertEquals(second[j].size() == 3, 0 == first[i].countIdenticalByInchi(second[j]));
			}
		}
		RetentionsDataset filtered = a();
		filtered.filterIdentical(b());
		Assert.assertEquals(14, filtered.size());
		Assert.assertEquals(5, filtered.compounds().size());
		Assert.assertEquals(5, filtered.compoundsCanonical(true).size());
		Assert.assertEquals(3, filtered.compoundsCanonical(false).size());

		filtered = b();
		filtered.filterIdentical(a());
		Assert.assertEquals(1, filtered.size());
		Assert.assertEquals("CCCCC(=N)(N)O", filtered.getSmiles(0));
		Assert.assertEquals(5, filtered.getColumn(0));
		Assert.assertEquals(83F, filtered.getRetention(0));

		filtered = c();
		filtered.filterIdentical(a());
		Assert.assertEquals(3, filtered.size());
		for (int i = 0; i < filtered.size(); i++) {
			if (filtered.getSmiles(i).equals("CC")) {
				Assert.assertEquals(1, filtered.getColumn(i));
				Assert.assertEquals(7.5F, filtered.getRetention(i));
			}
		}
	}

	public void testFilterIdenticalByInchi() throws CDKException {
		RetentionsDataset[] first = new RetentionsDataset[] { a(), b(), c() };
		RetentionsDataset[] second = new RetentionsDataset[] { a(), b(), c(), d() };
		for (int i = 0; i < first.length; i++) {
			for (int j = 0; j < second.length; j++) {
				int a = first[i].countIdenticalByInchi(second[j]);
				RetentionsDataset copy = first[i].copy();
				copy.filterIdenticalByInchi(second[j]);
				Assert.assertEquals(0, copy.countIdenticalByInchi(second[j]));
				Assert.assertEquals(0, second[j].countIdenticalByInchi(copy));
				Assert.assertEquals(a, first[i].countIdenticalByInchi(second[j]));
				Assert.assertEquals(second[j].size() == 3, 0 == first[i].countIdenticalByInchi(second[j]));
			}
		}
		RetentionsDataset filtered = a();
		filtered.filterIdenticalByInchi(b());
		Assert.assertEquals(14, filtered.size());
		Assert.assertEquals(5, filtered.compounds().size());
		Assert.assertEquals(5, filtered.compoundsCanonical(true).size());
		Assert.assertEquals(3, filtered.compoundsCanonical(false).size());

		filtered = b();
		filtered.filterIdenticalByInchi(a());
		Assert.assertEquals(1, filtered.size());
		Assert.assertEquals("CCCCC(=N)(N)O", filtered.getSmiles(0));
		Assert.assertEquals(5, filtered.getColumn(0));
		Assert.assertEquals(83F, filtered.getRetention(0));

		filtered = c();
		filtered.filterIdenticalByInchi(a());
		Assert.assertEquals(3, filtered.size());
		for (int i = 0; i < filtered.size(); i++) {
			if (filtered.getSmiles(i).equals("CC")) {
				Assert.assertEquals(1, filtered.getColumn(i));
				Assert.assertEquals(7.5F, filtered.getRetention(i));
			}
		}
	}

	public void testSaveToFile() throws IOException {
		RetentionsDataset a = a();
		a.saveToFile("test.txt");
		RetentionsDataset b = RetentionsDataset.loadFromFile("test.txt");
		Assert.assertEquals(false, a == b);
		Assert.assertEquals(true, a.size() == b.size());

		for (int i = 0; i < a.size(); i++) {
			Assert.assertEquals(false, a.getEntry(i) == b.getEntry(i));
			Assert.assertEquals(true, a.getEntry(i).equals(b.getEntry(i)));
			Assert.assertEquals(true, a.getSmiles(i).equals(b.getSmiles(i)));
			Assert.assertEquals(true, a.getRetention(i) == b.getRetention(i));
			Assert.assertEquals(true, a.getColumn(i) == b.getColumn(i));
		}
	}

	public void testLoadFromFileString() throws IOException {
		RetentionsDataset a = RetentionsDataset.create(new RetentionsEntry[] {});
		a.saveToFile("test.txt");
		RetentionsDataset b = RetentionsDataset.loadFromFile("test.txt");
		Assert.assertEquals(true, a.size() == 0);
		Assert.assertEquals(true, b.size() == 0);
	}

	public void testMerge() {
		RetentionsDataset a = a();
		RetentionsDataset aa = a();
		RetentionsDataset b = a.simpleSplit(1F / 6F);
		RetentionsDataset c = a.simpleSplit(1F / 5F);
		Assert.assertEquals(5, b.size());
		Assert.assertEquals(5, c.size());
		Assert.assertEquals(20, a.size());
		RetentionsDataset merged = RetentionsDataset.merge(new RetentionsDataset[] { b, c, a });
		Assert.assertEquals(aa.size(), merged.size());
		for (int i = 0; i < merged.size(); i++) {
			Assert.assertEquals(aa.getEntry(i), merged.getEntry(i));
		}
	}

	public void testMakeCanoncalAll() throws CDKException {
		RetentionsDataset a = a();
		RetentionsDataset b = a();
		a.makeCanoncalAll(true);
		for (int i = 0; i < a.size(); i++) {
			if ((i != 0) && (i != 1) && (i != 2) && (i != 4) && (i != 5) && (i != 6) && (i != 9)) {
				Assert.assertEquals(a.getSmiles(i), b.getSmiles(i));
			}
			if ((i < 8)) {
				Assert.assertEquals("COc1cc(ccc1O)C=O", a.getSmiles(i));
			}
			if (i == 9) {
				Assert.assertEquals("CCCC", a.getSmiles(i));
			}
			Assert.assertEquals(true, a.getRetention(i) == b.getRetention(i));
			Assert.assertEquals(true, a.getColumn(i) == b.getColumn(i));
		}
		a.makeCanoncalAll(false);
		for (int i = 0; i < a.size(); i++) {
			Assert.assertEquals(Chemoinformatics.canonical(b.getSmiles(i), false), a.getSmiles(i));
			Assert.assertEquals(true, a.getRetention(i) == b.getRetention(i));
			Assert.assertEquals(true, a.getColumn(i) == b.getColumn(i));
			if ((i > 10) && (i < 20)) {
				Assert.assertEquals("CC=CC", a.getSmiles(i));
			}
		}

	}

	public void testGroupByCompounds() throws CDKException {
		RetentionsDataset a = a();
		HashMap<String, ArrayList<RetentionsEntry>> groups = a.groupByCompounds(true);
		Assert.assertEquals(9, groups.size());
		Assert.assertEquals(8, groups.get("COc1cc(ccc1O)C=O").size());
		Assert.assertEquals(2, groups.get("CC=CC").size());
		Assert.assertEquals(3, groups.get("C/C=C/C").size());
		Assert.assertEquals(4, groups.get("C/C=C\\C").size());
		Assert.assertEquals(null, groups.get("CC"));
		int i = 0;
		int j = 0;
		for (RetentionsEntry e : groups.get("CC=CC")) {
			if ((e.getRetention() == 9.5F) && (e.getColumnType() == 3)) {
				i += 1;
			}
			if ((e.getRetention() == 8.3F) && (e.getColumnType() == 1)) {
				j += 1;
			}
			Assert.assertEquals(true, ((e.getRetention() == 9.5F) && (e.getColumnType() == 3))
					|| ((e.getRetention() == 8.3F) && (e.getColumnType() == 1)));
		}
		Assert.assertEquals(1, i);
		Assert.assertEquals(1, j);

		groups = a.groupByCompounds(false);
		Assert.assertEquals(5, groups.size());
		Assert.assertEquals(8, groups.get("COc1cc(ccc1O)C=O").size());
		Assert.assertEquals(9, groups.get("CC=CC").size());
		Assert.assertEquals(null, groups.get("C/C=C/C"));
		Assert.assertEquals(null, groups.get("C/C=C\\C"));
		Assert.assertEquals(null, groups.get("CC"));
		i = 0;
		j = 0;
		int k = 0;
		int l = 0;
		int m = 0;
		int n = 0;
		int o = 0;

		for (RetentionsEntry e : groups.get("CC=CC")) {
			if ((e.getRetention() == 7.5F) && (e.getColumnType() == 1)) {
				i += 1;
			}
			if ((e.getRetention() == 9.5F) && (e.getColumnType() == 2)) {
				j += 1;
			}
			if ((e.getRetention() == 8.5F) && (e.getColumnType() == 1)) {
				k += 1;
			}
			if ((e.getRetention() == 9.5F) && (e.getColumnType() == 3)) {
				l += 1;
			}
			if ((e.getRetention() == 8.3F) && (e.getColumnType() == 1)) {
				m += 1;
			}
			if ((e.getRetention() == 8.3F) && (e.getColumnType() == 4)) {
				n += 1;
			}
			if ((e.getRetention() == 8.3F) && (e.getColumnType() == 3)) {
				o += 1;
			}
		}
		Assert.assertEquals(1, i);
		Assert.assertEquals(2, j);
		Assert.assertEquals(1, k);
		Assert.assertEquals(1, l);
		Assert.assertEquals(1, m);
		Assert.assertEquals(1, n);
		Assert.assertEquals(2, o);
	}

	public void testMeanByCompounds() throws CDKException {
		RetentionsDataset a = a();
		RetentionsDataset b = a.meanByCompounds(true);
		Assert.assertEquals(9, b.size());
		Assert.assertEquals(9, b.compounds().size());
		Assert.assertEquals(11, a.compounds().size());
		for (int i = 0; i < b.size(); i++) {
			boolean found = false;
			if (b.getSmiles(i).equals("COc1cc(ccc1O)C=O")) {
				found = true;
				Assert.assertEquals(RetentionsDataset.mean(new float[] { 10, 10, 12, 14, 12, 13, 10, 15 }),
						b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("CCCC")) {
				found = true;
				Assert.assertEquals(14F / 3F, b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("C/C=C/C")) {
				found = true;
				Assert.assertEquals(8.5F, b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("CC=CC")) {
				found = true;
				Assert.assertEquals((8.3F + 9.5F) / 2F, b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("C/C=C\\C")) {
				found = true;
				Assert.assertEquals(RetentionsDataset.mean(new float[] { 9.5F, 8.3F, 8.3F, 8.3F }), b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("C(=N)(N)O")) {
				found = true;
				Assert.assertEquals(RetentionsDataset.mean(new float[] { 83F }), b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (!found) {
				Assert.assertEquals(true, Chemoinformatics.canonical(b.getSmiles(i), false).equals("CC=CCCCCCC(C)CCC"));
			}
		}
		b = a.meanByCompounds(false);
		Assert.assertEquals(5, b.size());
		Assert.assertEquals(5, b.compounds().size());
		for (int i = 0; i < b.size(); i++) {
			boolean found = false;
			if (b.getSmiles(i).equals("COc1cc(ccc1O)C=O")) {
				found = true;
				Assert.assertEquals(RetentionsDataset.mean(new float[] { 10, 10, 12, 14, 12, 13, 10, 15 }),
						b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}

			if (b.getSmiles(i).equals("CC=CC")) {
				found = true;
				Assert.assertEquals(
						RetentionsDataset.mean(new float[] { 7.5F, 9.5F, 8.5F, 9.5F, 8.3F, 9.5F, 8.3F, 8.3F, 8.3F }),
						b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("CCCC")) {
				found = true;
				Assert.assertEquals(RetentionsDataset.mean(new float[] { 5F, 5.5F, 3.5F }), b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("C(=N)(N)O")) {
				found = true;
				Assert.assertEquals(RetentionsDataset.mean(new float[] { 83F }), b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (!found) {
				Assert.assertEquals(true, b.getSmiles(i).equals("CC=CCCCCCC(C)CCC"));
			}
		}
	}

	public void testMedianByCompounds() throws CDKException {
		RetentionsDataset a = a();
		RetentionsDataset b = a.medianByCompounds(true);
		Assert.assertEquals(9, b.size());
		Assert.assertEquals(9, b.compounds().size());
		Assert.assertEquals(11, a.compounds().size());
		for (int i = 0; i < b.size(); i++) {
			boolean found = false;
			if (b.getSmiles(i).equals("COc1cc(ccc1O)C=O")) {
				found = true;
				Assert.assertEquals(RetentionsDataset.median(new float[] { 10, 10, 12, 14, 12, 13, 10, 15 }),
						b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("CCCC")) {
				found = true;
				Assert.assertEquals(5F, b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("C/C=C/C")) {
				found = true;
				Assert.assertEquals(8.5F, b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("CC=CC")) {
				found = true;
				Assert.assertEquals((8.3F + 9.5F) / 2F, b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("C/C=C\\C")) {
				found = true;
				Assert.assertEquals(RetentionsDataset.median(new float[] { 9.5F, 8.3F, 8.3F, 8.3F }),
						b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("C(=N)(N)O")) {
				found = true;
				Assert.assertEquals(RetentionsDataset.median(new float[] { 83F }), b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (!found) {
				Assert.assertEquals(true, Chemoinformatics.canonical(b.getSmiles(i), false).equals("CC=CCCCCCC(C)CCC"));
			}
		}
		b = a.medianByCompounds(false);
		Assert.assertEquals(5, b.size());
		Assert.assertEquals(5, b.compounds().size());
		for (int i = 0; i < b.size(); i++) {
			boolean found = false;
			if (b.getSmiles(i).equals("COc1cc(ccc1O)C=O")) {
				found = true;
				Assert.assertEquals(12F, b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}

			if (b.getSmiles(i).equals("CC=CC")) {
				found = true;
				Assert.assertEquals(
						RetentionsDataset.median(new float[] { 7.5F, 9.5F, 8.5F, 9.5F, 8.3F, 9.5F, 8.3F, 8.3F, 8.3F }),
						b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("CCCC")) {
				found = true;
				Assert.assertEquals(5F, b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (b.getSmiles(i).equals("C(=N)(N)O")) {
				found = true;
				Assert.assertEquals(RetentionsDataset.mean(new float[] { 83F }), b.getRetention(i));
				Assert.assertEquals(-1, b.getColumn(i));
			}
			if (!found) {
				Assert.assertEquals(true, b.getSmiles(i).equals("CC=CCCCCCC(C)CCC"));
			}
		}
	}

	public void testMean() {
		float[] a = {};
		float[] b = { 1.0F };
		float[] c = { 2.0F, 3.0F };
		float[] d = { 3.0F, 4.0F, 5.0F };
		float[] e = { 5.0F, 6.0F, 7.0F, 22.0F };
		float[] f = { 4.0F, 6.0F, 7.0F, 22.0F, 6.0F };
		Assert.assertEquals(Float.NaN, RetentionsDataset.mean(a));
		Assert.assertEquals(1.0F, RetentionsDataset.mean(b));
		Assert.assertEquals(2.5F, RetentionsDataset.mean(c));
		Assert.assertEquals(4.0F, RetentionsDataset.mean(d));
		Assert.assertEquals(10F, RetentionsDataset.mean(e));
		Assert.assertEquals(9F, RetentionsDataset.mean(f));
	}

	public void testMedian() {
		float[] a = {};
		float[] b = { 1.0F };
		float[] c = { 2.0F, 3.0F };
		float[] d = { 3.0F, 4.0F, 5.0F };
		float[] e = { 5.0F, 6.0F, 7.0F, 22.0F };
		float[] f = { 4.0F, 6.0F, 7.0F, 22.0F, 6.0F };
		float[] g = { 4.0F, 6.0F, 7.0F, 22.0F, 36.0F };
		Assert.assertEquals(Float.NaN, RetentionsDataset.median(a));
		Assert.assertEquals(1.0F, RetentionsDataset.median(b));
		Assert.assertEquals(2.5F, RetentionsDataset.median(c));
		Assert.assertEquals(4.0F, RetentionsDataset.median(d));
		Assert.assertEquals(6.5F, RetentionsDataset.median(e));
		Assert.assertEquals(6F, RetentionsDataset.median(f));
		Assert.assertEquals(7F, RetentionsDataset.median(g));
	}

	public void testMergeArrays() {
		float[] a = {};
		float[] b = { 1.0F };
		float[] c = { 2.0F, 3.0F };
		float[] d = { 3.0F, 4.0F, 5.0F };
		float[] e = { 5.0F, 6.0F, 7.0F, 8.0F };
		ArrayList<float[]> f = new ArrayList<float[]>();
		f.add(a);
		f.add(b);
		f.add(c);
		f.add(d);
		f.add(e);
		for (int i = 0; i < f.size(); i++) {
			Assert.assertEquals(f.get(i).length, RetentionsDataset.mergeArrays(a, f.get(i)).length);
			Assert.assertEquals(f.get(i).length, RetentionsDataset.mergeArrays(f.get(i), a).length);
			for (int j = 0; j < f.get(i).length; j++) {
				Assert.assertEquals(f.get(i)[j], RetentionsDataset.mergeArrays(a, f.get(i))[j]);
				Assert.assertEquals(f.get(i)[j], RetentionsDataset.mergeArrays(f.get(i), a)[j]);
			}
		}
		for (int i = 0; i < f.size(); i++) {
			for (int j = 0; j < f.size(); j++) {
				Assert.assertEquals(f.get(i).length + f.get(j).length,
						RetentionsDataset.mergeArrays(f.get(j), f.get(i)).length);
			}
		}
		float[] q = RetentionsDataset.mergeArrays(d, e);
		Assert.assertEquals(q[0], 3.0F);
		Assert.assertEquals(q[1], 4.0F);
		Assert.assertEquals(q[2], 5.0F);
		Assert.assertEquals(q[3], 5.0F);
		Assert.assertEquals(q[4], 6.0F);
		Assert.assertEquals(q[5], 7.0F);
		Assert.assertEquals(q[6], 8.0F);
	}

}
