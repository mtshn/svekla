package ru.ac.phyche.gcms.svekla;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import javax.imageio.ImageIO;
import org.openscience.cdk.exception.CDKException;
import junit.framework.Assert;
import junit.framework.TestCase;

public class ChemoinformaticsTest extends TestCase {

	public void testWeight() throws CDKException {
		Assert.assertEquals(true, Math.abs(Chemoinformatics.weight("c1ccccc1") - 78) < 0.5);
		Assert.assertEquals(true, Math.abs(Chemoinformatics.weight("CCO") - 15 - 14 - 17) < 0.5);

	}

	public void testTokenize() throws CDKException {
		int[] a = Chemoinformatics.tokenize("c1ccccc1CCO");
		Assert.assertEquals(2, a[0]);
		Assert.assertEquals(22, a[1]);
		Assert.assertEquals(2, a[2]);
		Assert.assertEquals(2, a[3]);
		Assert.assertEquals(2, a[4]);
		Assert.assertEquals(2, a[5]);
		Assert.assertEquals(2, a[6]);
		Assert.assertEquals(22, a[7]);
		Assert.assertEquals(1, a[8]);
		Assert.assertEquals(1, a[9]);
		Assert.assertEquals(6, a[10]);

		for (int i = 11; i < 100; i++) {
			Assert.assertEquals(0, a[i]);
		}

		String s = Chemoinformatics.intsToSmiles(a).trim();
		Assert.assertEquals("c1ccccc1CCO", s);
		boolean exception = false;
		try {
			Chemoinformatics.tokenize("c1ccccc1CCmO");
		} catch (CDKException e) {
			exception = true;
		}
		Assert.assertEquals(true, exception);
		exception = false;
		try {
			s = "C";
			for (int i = 0; i < 300; i++) {
				s = s + "C";
			}
			Chemoinformatics.tokenize(s);
		} catch (CDKException e) {
			exception = true;
		}
		Assert.assertEquals(true, exception);
	}

	public void testIntsToSmiles() {
		int[] a = new int[] { 2, 22, 2, 2, 2, 2, 2, 22, 1, 1, 6, 0, 0, 0, 0, 0 };
		String s = Chemoinformatics.intsToSmiles(a).trim();
		Assert.assertEquals("c1ccccc1CCO", s);
	}

	public void testSmilesToInchiKey() throws CDKException {
		Assert.assertEquals("MWOOGOJBHIARFG-UHFFFAOYSA-N",
				Chemoinformatics.smilesToInchiKey("COC1=C(C=CC(=C1)C=O)O").trim());
		Assert.assertEquals("CZMRCDWAGMRECN-UGDNZRGBSA-N",
				Chemoinformatics
						.smilesToInchiKey(
								"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O")
						.trim());
		Assert.assertEquals("CZMRCDWAGMRECN-UGDNZRGBSA-N",
				Chemoinformatics
						.smilesToInchiKey(
								"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O")
						.trim());
		Assert.assertEquals("CZMRCDWAGMRECN-UHFFFAOYSA-N",
				Chemoinformatics.smilesToInchiKey("C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O").trim());
	}

	public void testSmilesToInchi() throws CDKException {
		Assert.assertEquals("InChI=1S/C8H8O3/c1-11-8-4-6(5-9)2-3-7(8)10/h2-5,10H,1H3",
				Chemoinformatics.smilesToInchi("COC1=C(C=CC(=C1)C=O)O").trim());
		Assert.assertEquals(
				"InChI=1S/C12H22O11/c13-1-4-6(16)8(18)9(19)11(21-4)23-12(3-15)10(20)7(17)5(2-14)22-12/h4-11,13-20H,1-3H2/t4-,5-,6-,7-,8+,9-,10+,11-,12+/m1/s1",
				Chemoinformatics
						.smilesToInchi(
								"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O")
						.trim());
		Assert.assertEquals(
				"InChI=1S/C12H22O11/c13-1-4-6(16)8(18)9(19)11(21-4)23-12(3-15)10(20)7(17)5(2-14)22-12/h4-11,13-20H,1-3H2/t4-,5-,6-,7-,8+,9-,10+,11-,12+/m1/s1",
				Chemoinformatics
						.smilesToInchi(
								"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O")
						.trim());
		Assert.assertEquals(
				"InChI=1S/C12H22O11/c13-1-4-6(16)8(18)9(19)11(21-4)23-12(3-15)10(20)7(17)5(2-14)22-12/h4-11,13-20H,1-3H2",
				Chemoinformatics.smilesToInchi("C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O").trim());
	}

	public void testInchiToSmiles() throws CDKException {
		Assert.assertEquals("C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O",
				Chemoinformatics.inchiToSmiles(
						"InChI=1S/C12H22O11/c13-1-4-6(16)8(18)9(19)11(21-4)23-12(3-15)10(20)7(17)5(2-14)22-12/h4-11,13-20H,1-3H2/t4-,5-,6-,7-,8+,9-,10+,11-,12+/m1/s1",
						true));
		Assert.assertEquals(Chemoinformatics.canonical("C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O", true),
				Chemoinformatics.inchiToSmiles(
						"InChI=1S/C12H22O11/c13-1-4-6(16)8(18)9(19)11(21-4)23-12(3-15)10(20)7(17)5(2-14)22-12/h4-11,13-20H,1-3H2/t4-,5-,6-,7-,8+,9-,10+,11-,12+/m1/s1",
						false));
		Assert.assertEquals(Chemoinformatics.inchiToSmiles(
				"InChI=1S/C12H22O11/c13-1-4-6(16)8(18)9(19)11(21-4)23-12(3-15)10(20)7(17)5(2-14)22-12/h4-11,13-20H,1-3H2",
				true),
				Chemoinformatics.inchiToSmiles(
						"InChI=1S/C12H22O11/c13-1-4-6(16)8(18)9(19)11(21-4)23-12(3-15)10(20)7(17)5(2-14)22-12/h4-11,13-20H,1-3H2/t4-,5-,6-,7-,8+,9-,10+,11-,12+/m1/s1",
						false));
	}

	public void testCanonicalizeSmile() throws CDKException {
		Assert.assertEquals("COc1cc(ccc1O)C=O", Chemoinformatics.canonical("COC1=C(C=CC(=C1)C=O)O", true));
		Assert.assertEquals("COc1cc(ccc1O)C=O", Chemoinformatics.canonical("COc1cc(ccc1O)C=O", true));
		Assert.assertEquals("COc1cc(ccc1O)C=O", Chemoinformatics.canonical("COc1cc(ccc1O)C=O", false));
		Assert.assertEquals(Chemoinformatics.canonical("C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O", true),
				Chemoinformatics.canonical(
						"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O",
						false));
	}

	public void testFingerprints() throws CDKException {
		String[] compounds = new String[] { "C", "CC", "COc1cc(ccc1O)C=O",
				"C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O",
				"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O" };
		for (int i = 0; i < 5; i++) {
			Assert.assertEquals(1024, Chemoinformatics.fingerprints(compounds[i],
					Chemoinformatics.FingerprintsType.ADDITIVE_CIRCULAR_4_1024_NO_SCALE).length);
			Assert.assertEquals(1024, Chemoinformatics.fingerprints(compounds[i],
					Chemoinformatics.FingerprintsType.ADDITIVE_CIRCULAR_6_1024_NO_SCALE).length);
			Assert.assertEquals(1024, Chemoinformatics.fingerprints(compounds[i],
					Chemoinformatics.FingerprintsType.CIRCULAR_4_1024).length);
			Assert.assertEquals(1024, Chemoinformatics.fingerprints(compounds[i],
					Chemoinformatics.FingerprintsType.CIRCULAR_6_1024).length);
			Assert.assertEquals(4096, Chemoinformatics.fingerprints(compounds[i],
					Chemoinformatics.FingerprintsType.CIRCULAR_6_4096).length);
			Assert.assertEquals(896,
					Chemoinformatics.fingerprints(compounds[i], Chemoinformatics.FingerprintsType.PUBCHEM).length);
			Assert.assertEquals(0,
					Chemoinformatics.fingerprints(compounds[i], Chemoinformatics.FingerprintsType.NONE).length);
			Assert.assertEquals(192,
					Chemoinformatics.fingerprints(compounds[i], Chemoinformatics.FingerprintsType.MACCS).length);
			Assert.assertEquals(1024,
					Chemoinformatics.fingerprints(compounds[i], Chemoinformatics.FingerprintsType.LINGO).length);
			Assert.assertEquals(4864,
					Chemoinformatics.fingerprints(compounds[i], Chemoinformatics.FingerprintsType.KLEKOTA).length);
		}

		for (Chemoinformatics.FingerprintsType t : Chemoinformatics.FingerprintsType.values()) {
			float prevsum = 0;
			for (int i = 0; i < 5; i++) {
				float sum = 0;
				float[] fp = Chemoinformatics.fingerprints(compounds[i], t);
				for (int j = 0; j < fp.length; j++) {
					if (Math.abs(fp[j] * 2 - 1) != 1.0F) {
						Assert.assertEquals(true,
								((t == Chemoinformatics.FingerprintsType.ADDITIVE_CIRCULAR_6_1024_NO_SCALE)
										|| (t == Chemoinformatics.FingerprintsType.ADDITIVE_CIRCULAR_4_1024_NO_SCALE)));
					}
					sum += fp[j];
				}
				if (t != Chemoinformatics.FingerprintsType.LINGO) {
					if (t != Chemoinformatics.FingerprintsType.NONE) {
						Assert.assertEquals(true, (sum != prevsum) || (i == 4));
						Assert.assertEquals(true, sum != 0);
					}
				}
				prevsum = sum;
			}
		}
	}

	public void testDescriptorsStringStringArrayFloatArrayFloatArray() throws CDKException {
		String[] descriptorNames = { "fragC", "nAtomLC" };
		float[] desc = Chemoinformatics.descriptors("COc1cc(ccc1O)C=O", descriptorNames, new float[] { 0, 0 },
				new float[] { 1000, 70 });
		Assert.assertEquals(true, Math.abs(0.25103 - desc[0]) < 0.000001);
		Assert.assertEquals(true, Math.abs(0.1 - desc[1]) < 0.000001);
	}

	public void testAllDescriptorNames() {
		String[] s = Chemoinformatics.allDescriptorNames();
		HashSet<String> strings = new HashSet<String>();
		for (int i = 0; i < s.length; i++) {
			strings.add(s[i]);
		}
		Assert.assertEquals(310, s.length);
		Assert.assertEquals(310, strings.size());
	}

	public void testDescriptorsStringStringArray() throws CDKException {
		String[] descriptorNames = { "fragC", "C1SP1", "C2SP1", "C1SP2", "C2SP2", "C3SP2", "C1SP3", "C2SP3", "C3SP3",
				"C4SP3", "SCH-3", "SCH-4", "SCH-5", "SCH-6", "SCH-7", "VCH-3", "VCH-4", "VCH-5", "VCH-6", "VCH-7",
				"nAtomLC" };
		String[] compounds = { "COc1cc(ccc1O)C=O", "C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O",
				"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O" };
		float[] answers1 = { 251.03F, 0, 0, 1, 5, 1, 0, 0, 0, 0, 0.0F, 0.0F, 0.0F, 0.06804138F, 0.16426642F, 0.0F, 0.0F,
				0.0F, 0.0240562F, 0.034468F, 7 };
		float[] answers2 = { 1610.11F, 0, 0, 0, 0, 0, 4, 8, 0, 0, 0.0F, 0.0F, 0.06804138F, 0.32578125F, 0.7830142F,
				0.0F, 0.0F, 0.0392837101F, 0.13291862F, 0.2126649075F, 11 };
		float[] desc = Chemoinformatics.descriptors(compounds[0], descriptorNames);
		for (int i = 0; i < desc.length; i++) {
			Assert.assertEquals(true, Math.abs(answers1[i] - desc[i]) < 0.0001);
		}
		Assert.assertEquals(answers1.length, desc.length);
		desc = Chemoinformatics.descriptors(compounds[1], descriptorNames);
		for (int i = 0; i < desc.length; i++) {
			Assert.assertEquals(true, Math.abs(answers2[i] - desc[i]) < 0.0001);
		}
		Assert.assertEquals(answers2.length, desc.length);
		desc = Chemoinformatics.descriptors(compounds[2], descriptorNames);
		for (int i = 0; i < desc.length; i++) {
			Assert.assertEquals(true, Math.abs(answers2[i] - desc[i]) < 0.0001);
		}
		Assert.assertEquals(answers2.length, desc.length);
	}

	public void testSmilesToImage() throws CDKException, IOException {
		String s = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O";

		float[][] d = Chemoinformatics.smilesToImage(s);
		BufferedImage img = Chemoinformatics.smilesToImageAsBufferedImage(s);
		for (int i = 0; i < Chemoinformatics.DEPICTION_SIZE; i++) {
			for (int j = 0; j < Chemoinformatics.DEPICTION_SIZE; j++) {
				if (d[i][j] == 0) {
					Assert.assertEquals(-1, img.getRGB(i, j));
				} else {
					Assert.assertEquals(true, -1 != img.getRGB(i, j));
				}
			}
		}

		Graphics2D graphics = img.createGraphics();
		graphics.setPaint(new Color(0, 0, 0));
		graphics.fillRect(0, 0, img.getWidth(), img.getHeight());
		for (int i = 0; i < Chemoinformatics.DEPICTION_SIZE; i++) {
			for (int j = 0; j < Chemoinformatics.DEPICTION_SIZE; j++) {
				if (d[i][j] != 0) {
					img.setRGB(i, j, -1);
				}
			}
		}
		ImageIO.write(img, "png", new File("test.png"));
	}

	public void testFuncGroups() throws CDKException {
		String s = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O";
		float[] a = Chemoinformatics.funcGroups(s);
		float[] b = new float[] { 0, 3, 0, 0, 8, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 3, 5, 0, 0, 1, 2 };
		for (int i = 0; i < b.length; i++) {
			Assert.assertEquals(b[i], a[i]);
		}
		for (int i = b.length; i < a.length; i++) {
			Assert.assertEquals(0.0F, a[i]);
		}
		s = "O=C(NC)c1cccc(c1)C(c2ccccc2)c3ccccc3";
		a = Chemoinformatics.funcGroups(s);
		b = new float[] { 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 14, 4, 0 };
		for (int i = 0; i < b.length; i++) {
			Assert.assertEquals(b[i], a[i]);
		}
		for (int i = b.length; i < a.length; i++) {
			if ((i != 32) && (i != 26) && (i != 38)) {
				Assert.assertEquals(0.0F, a[i]);
			} else {
				Assert.assertEquals(1.0F, a[i]);
			}
		}
		s = "O=C(NC)c1cccc(c1)C(c2ccccc2)c3cccc(c3)c4ccc5ccccc5(c4)";
		a = Chemoinformatics.funcGroups(s);
		b = new float[] { 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 20, 6, 2 };
		for (int i = 0; i < b.length; i++) {
			Assert.assertEquals(b[i], a[i]);
		}
		for (int i = b.length; i < a.length; i++) {
			if ((i != 32) && (i != 26) && (i != 38)) {
				Assert.assertEquals(0.0F, a[i]);
			} else {
				Assert.assertEquals(1.0F, a[i]);
			}
		}
		s = "c1ccc2ccccc2(c1)";
		a = Chemoinformatics.funcGroups(s);
		b = new float[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 2 };
		for (int i = 0; i < b.length; i++) {
			Assert.assertEquals(b[i], a[i]);
		}

		s = "N1=C(N(C)CCC1)C";
		a = Chemoinformatics.funcGroups(s);
		b = new float[] { 2, 0, 3 };
		for (int i = 0; i < b.length; i++) {
			Assert.assertEquals(b[i], a[i]);
		}
		for (int i = b.length; i < a.length; i++) {
			if ((i != 49) && (i != 11) && (i != 41) && (i != 47)) {
				Assert.assertEquals(0.0F, a[i]);
			} else {
				Assert.assertEquals(1.0F, a[i]);
			}
		}
		s = "O=NN(C)CC1N(N=O)CC(N=O)CC1";
		a = Chemoinformatics.funcGroups(s);
		b = new float[] { 1, 1, 3, 0, 2 };
		for (int i = 0; i < b.length; i++) {
			Assert.assertEquals(b[i], a[i]);
		}
		for (int i = b.length; i < a.length; i++) {
			if ((i != 55) && (i != 41) && (i != 42) && (i != 40) && (i != 43) && (i != 46)) {
				Assert.assertEquals(0.0F, a[i]);
			} else {
				Assert.assertEquals(((i == 46) || (i == 55)) ? 3.0F : 1.0F, a[i]);
			}
		}
		s = "O1NCCCC1";
		a = Chemoinformatics.funcGroups(s);
		b = new float[] { 0, 0, 4 };
		for (int i = 0; i < b.length; i++) {
			Assert.assertEquals(b[i], a[i]);
		}
		for (int i = b.length; i < a.length; i++) {
			if ((i != 39) && (i != 23)) {
				Assert.assertEquals(0.0F, a[i]);
			} else {
				Assert.assertEquals(1.0F, a[i]);
			}
		}
	}

	public void testRepresentation2d() throws CDKException {
		String[] compounds = new String[] { "C", "CC", "CC(C)O", "COc1cc(ccc1O)C=O",
				"C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O",
				"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O",
				"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" };
		int[] atomNums = new int[] { 1, 2, 4, 11, 23, 23, 50 };
		int[] bondNums = new int[] { 0, 1, 3, 11, 24, 24, 49 };

		String[] atomTypes = new String[] { "B", "C", "F", "H", "Cl", "I", "S.planar3", "C.sp2", "C.sp3", "N", "O", "P",
				"N.nitro", "S", "X", "S.3", "O.planar3", "N.planar3", "Br", "N.amide", "N.sp2", "Si", "N.sp3", "C.sp",
				"O.sp3", "O.sp2" };
		for (int s = 0; s < compounds.length; s++) {
			float[][][] a = Chemoinformatics.representation2d(compounds[s]);
			Assert.assertEquals(a.length, atomTypes.length + 3);
			Assert.assertEquals(a[0].length, 130);
			Assert.assertEquals(a[0][0].length, 130);
			int countAtomsBonds = 0;
			int countAtoms = 0;

			for (int i = 0; i < a.length; i++) {
				for (int j = 0; j < a[0].length; j++) {
					for (int k = 0; k < a[0][0].length; k++) {
						if (a[i][j][k] != 0) {
							countAtomsBonds++;
							if (i < atomTypes.length) {
								countAtoms++;
							}
							if (s == 0) {
								Assert.assertEquals(i, 8);
								Assert.assertEquals(j, 65);
								Assert.assertEquals(k, 65);
							}
							if ((s == 0) && (i == 26)) {
								Assert.assertEquals(j, 65);
								Assert.assertEquals(k, 65);
							}
						}
					}
				}
			}
			Assert.assertEquals(atomNums[s], countAtoms);
			Assert.assertEquals(bondNums[s], countAtomsBonds - countAtoms);
		}
		String s = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
		boolean t = false;
		try {
			Chemoinformatics.representation2d(s);
		} catch (CDKException e) {
			t = true;
		}
		Assert.assertEquals(true, t);
	}

	public void testCountOH() throws CDKException {
		String[] compounds = new String[] { "C", "CC", "CC(C)O", "COc1cc(ccc1O)C=O",
				"C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O",
				"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O",
				"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "O=C(O)C1C(O)CCC(=O)C1(OO)",
				"O=C(O)C1C(O)CC(C(=O)C(O)O)C(=O)C1(OO)", "O=P(O)(O)C1C(C)CC1S(=O)(=O)O", "[H][O+](C)CCCOC",
				"[H]O(C)CCCOC", "O" };
		Assert.assertEquals(0, Chemoinformatics.countOH(compounds[0]));
		Assert.assertEquals(0, Chemoinformatics.countOH(compounds[1]));
		Assert.assertEquals(1, Chemoinformatics.countOH(compounds[2]));
		Assert.assertEquals(1, Chemoinformatics.countOH(compounds[3]));
		Assert.assertEquals(8, Chemoinformatics.countOH(compounds[4]));
		Assert.assertEquals(8, Chemoinformatics.countOH(compounds[5]));
		Assert.assertEquals(0, Chemoinformatics.countOH(compounds[6]));
		Assert.assertEquals(3, Chemoinformatics.countOH(compounds[7]));
		Assert.assertEquals(5, Chemoinformatics.countOH(compounds[8]));
		Assert.assertEquals(3, Chemoinformatics.countOH(compounds[9]));
		Assert.assertEquals(0, Chemoinformatics.countOH(compounds[10]));
		Assert.assertEquals(0, Chemoinformatics.countOH(compounds[11]));
		Assert.assertEquals(0, Chemoinformatics.countOH(compounds[12]));
	}

	public void testReplaceOHWithOSiCH33OH() throws CDKException {
		String[] compounds = new String[] { "C", "CC", "CC(C)O", "COc1cc(ccc1O)C=O",
				"C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O",
				"C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@]2(CO)[C@H]([C@@H]([C@@H](CO)O2)O)O)O)O)O)O",
				"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "O=C(O)C1C(O)CCC(=O)C1(OO)",
				"O=C(O)C1C(O)CC(C(=O)C(O)O)C(=O)C1(OO)", "O=P(O)(O)C1C(C)CC1S(=O)(=O)O", "[H][O+](C)CCCOC", "O" };
		for (int i = 0; i < 12; i++) {
			Assert.assertEquals(Chemoinformatics.countOH(compounds[i]) == 0, Chemoinformatics
					.canonical(compounds[i], true).equals(Chemoinformatics.replaceOHWithOSiCH33OH(compounds[i], true)));
		}
		Assert.assertEquals("CC(C)O[Si](C)(C)C", Chemoinformatics.replaceOHWithOSiCH33OH(compounds[2], true));
		Assert.assertEquals("COc1cc(ccc1O[Si](C)(C)C)C=O", Chemoinformatics.replaceOHWithOSiCH33OH(compounds[3], true));
		Assert.assertEquals(
				"C[Si](C)(C)OC1CC(C(=O)C(C1C(=O)O[Si](C)(C)C)OO[Si](C)(C)C)C(=O)C(O[Si](C)(C)C)O[Si](C)(C)C",
				Chemoinformatics.replaceOHWithOSiCH33OH(compounds[8], true));
	}
}
