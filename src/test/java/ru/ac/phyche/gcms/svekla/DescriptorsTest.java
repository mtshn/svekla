package ru.ac.phyche.gcms.svekla;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;

import org.openscience.cdk.exception.CDKException;

import junit.framework.Assert;
import junit.framework.TestCase;

public class DescriptorsTest extends TestCase {

	private String[] smiles = new String[] { "CCCSc1ccc2c(c1)[nH]c(n2)N=C(O)OC", "C1CC1(C(=O)O)N",
			"c1cc(ccc1CC(C(=O)O)N)O", "c1cc(ccc1CC(C(=O)O)N)O", "CC12CCC(CC2CCC3C4CCC(=O)C4(C)CCC31)O",
			"c1nc(c(N)[nH]1)C(=N)O", "c1nc(c(N)[nH]1)C(=N)O", "C(CN)C(=O)O", "C(CN)C(=O)O", "CC(=O)N1c2ccccc2C(=O)C1=O",
			"CC(=O)N1c2ccccc2C(=O)C1=O", "CC(=O)N1c2ccccc2C(=O)C1=O",
			"CC(=CCC[C@@H](C)[C@H]1CC[C@@]2(C)C3=C(CC[C@]12C)[C@@]4(C)CC[C@@H](C(C)(C)[C@@H]4CC3)O)C",
			"CC(C)CCC[C@@H](C)[C@H]1CC[C@@]2(C)C3=C(CC[C@]12C)[C@@]4(C)CC[C@@H](C(C)(C)[C@@H]4CC3)O",
			"C1CCNC(C1)c2cccnc2", "C1CCNC(C1)c2cccnc2", "CC(=Nc1ccccc1)O", "CC(=Nc1ccccc1)O", "CC(=Nc1ccccc1)O",
			"c1ccc(c(c1)N)O", "c1cc(ccc1N)O", "CC(=NC(Cc1c[nH]c2ccccc12)C(=O)O)O", "CC(=NC(Cc1c[nH]c2ccccc12)C(=O)O)O",
			"c1ccc2c(c1)C(=O)c3ccc(c(c3C2=O)O)O", "c1cc(c(cc1O)O)O", "c1ccc(cc1)C(C(=O)c2ccccc2)O",
			"c1ccc2c(c1)c3ccccc3[nH]2", "C(#CCO)CO", "c1ccc(cc1)-c2ccccc2", "c1ccc(c(c1)O)O",
			"CC(CCC(=O)O)C1CCC2C3C(CC(C12C)O)C4(C)CCC(CC4CC3O)O", "C([C@@H](C(=O)O)N)S", "COc1cc(ccc1O)C(C(=O)O)O",
			"c1ccc(cc1)CCN", "CC(C(c1cc(ccc1OC)OC)O)N", "Cc1c(CCO)scn1", "c1c(C(=O)O)nc(nc1O)O",
			"c1ccc(cc1)CC(C(=O)O)O", "c1c(cc(cc1O)O)O", "C(CCC(=O)O)CCC(=O)O", "c1cc(c(c(c1)O)O)O",
			"OP(=O)(O)OP(=O)(O)O", "c1ccc(cc1)CC(=N)O", "CCCCCC/C=C\\CCCCCCCC(=O)OC", "c1cc2ccc(nc2c(c1)O)O",
			"c1cc(c(C(=O)O)nc1)C(=O)O", "CC(C(=O)O)NC", "c1cc(ccc1N(=O)=O)O", "CCCCCCCC=O", "Cc1cc(cc(c1)O)O",
			"CCOP(=O)(OCC)Oc1ccc(cc1)N(=O)=O", "c1ccc(cc1)CC=O", "c1ccc(cc1)CC=O", "c1ccc(cc1)CC=O", "c1ccc(cc1)CC=O",
			"Cc1c(c(c(cn1)CO)C(=O)O)O", "c1cc(C(=O)O)[nH]c1", "c1ccc2ccccc2c1", "c1ccc2c(c1)c(ccn2=O)N(=O)=O",
			"Cc1ccc(cc1)CO", "C(CCC(=O)O)CCO", "C1CC(CCC1C(=O)O)O", "c1(c(c(c(c(c1Cl)Cl)Cl)Cl)Cl)Cl",
			"c1ccc(c(c1)CCC(=O)O)O", "c1cc(ccc1C(C(=O)O)O)O", "CC(C)c1ccc(cc1)C(=O)O", "CC(C)c1ccc(cc1)CO",
			"c1c(nc2c(C(=NC(=N)N2)O)n1)O", "CCCCCCCCCCCC(=O)O", "CC(C)CC(C(=O)O)N", "COc1ccc2c(c1)c(CC(=O)O)c[nH]2",
			"CNCC(c1ccccc1)O", "CNc1ccccc1", "Cc1cccc(c1)CO", "COc1cc2C=CC(=O)Oc2cc1O", "c1ccc(cc1)C(CO)C(=O)O",
			"c1cc(ccc1C(=O)O)C(=O)O", "C(C(C(=O)O)O)(C(=O)O)O", "c1ccc2c(c1)c(CCO)c[nH]2", "Cc1ccc(cc1)S(=O)(=O)O",
			"CC1C2CCC3(C)C=CC(=O)C(=C3C2OC1=O)C", "c1cc(c(c(c1)O)N)C(=O)O", "C1=CN(C2C(C(C(CO)O2)O)O)C(=NC1=N)O",
			"C1CCC(CC1)N=CO", "Cc1cccc(c1)O", "Cc1ccccc1O", "Cc1ccc(cc1)O", "C1CCC(CC1)NS(=O)(=O)O",
			"C1CCC(CC1)NS(=O)(=O)O", "OS(=O)(=O)O", "Cc1cccc(c1O)O", "Cc1ccc(c(c1)O)O", "c1ccc2c(c1)-c3ccccc3C2=O",
			"c1ccc-2c(c1)Cc3ccccc32", "CC(CN)C(=O)O", "CC(CN)C(=O)O", "C(C(COCC(CO)O)O)O", "c1cc(c(c(c1)O)O)C(=O)O",
			"c1cc(c(cc1C(=O)O)O)O", "COc1cc(/C=C/C(=O)O)cc(c1O)OC", "c1cnc(nc1O)O", "C1C(C(CO)OC1n2cnc3c(N)ncnc32)O",
			"CCCCC(CC)COC(=O)c1ccccc1C(=O)OCC(CC)CCCC", "c1cc(c(cc1CC(=O)O)O)O", "c1c(c(nc(n1)O)O)C(=O)O",
			"c1ccc2c(c1)c3ccccc3o2", "C1C(C(COP(=O)(O)O)OC1n2cnc3c(N)ncnc32)O", "COc1cc(ccc1O)/C=C/C(=O)O",
			"CC(=Nc1ccc-2c(Cc3ccccc32)c1)O", "Cc1cc2c(cc1C)nc[nH]2", "c1cc(C(=O)O)oc1", "C(C(=O)O)N=C(C(CS)N)O",
			"C(C(=O)O)N=C(C(CS)N)O", "c1c(cc(c(c1O)O)O)C(=O)O", "CCCCCCCCCCCCCCCCO", "c1ccc(cc1)CCC(=O)O",
			"c1cc(ccc1O)O", "c1ccc(c(c1)CO)O", "c1ccc(c(c1)/C=C/C(=O)O)O", "c1cc(/C=C/C(=O)O)cc(c1)O",
			"COc1cc(ccc1O)CO", "CCCCCCCCCCCCCCCCC(=O)O", "CCCCCCCCCCCCCCCCC(=O)OC", "c1cc(ccc1CC(=O)O)O",
			"c1cc(ccc1CC(C(=O)O)O)O", "c1cc(ccc1CCC(=O)O)O", "c1ccc2c(c1)ccc(n2)O", "C(CO)[C@@H](C(=O)O)N",
			"c1cc(c(cc1O)CC(=O)O)O", "c1cc(CCC(=O)O)cc(c1)O", "c1cc(c(c(c1)O)N)C(=O)O", "C(CO)C(=O)O",
			"COc1cc(ccc1O)C(=O)O", "C(C1C(C(C(n2cnc3c2ncnc3O)O1)O)O)O", "c1ccc(cc1)C(C#N)O", "CC(C(=O)O)C(=O)O",
			"CC1=CC(=O)Oc2cc(ccc12)O", "CNC(Cc1c[nH]c2ccccc12)C(=O)O", "C/C(=C/C(=O)O)/C(=O)O", "CSc1c2c([nH]cn2)ncn1",
			"CSc1c2c([nH]cn2)ncn1", "C(C(=O)O)C(=O)O", "C(C(=O)O)C(=O)O", "C(C(=N)O)C(=N)O", "C(C(=N)O)C(=N)O",
			"C(C(=N)O)C(=N)O", "CC1(CCOC(=O)C1)O", "CCCCCCCCCCCCCCCC", "c1ccc(cc1)CC(=O)O", "C(C(=O)O)OP(=O)(O)O",
			"c1cc(ccc1C(=O)O)O", "C1C(C(C(CC1(C(=O)O)O)O)O)O", "Cc1c(c(CO)c(cn1)CO)O", "c1ccc2c(c1)C(=NS2(=O)=O)O",
			"c1ccc(c(c1)C(=O)O)O", "CCCCCCCCCCCCCCCCCC(=O)O", "C(CC(=O)O)C(=O)O", "CNCC(c1ccc(cc1)O)O",
			"C(C(=O)O)(C(=O)O)O", "CC(C(=O)O)N", "c1cc(cnc1)O", "CC(C)C(C(=O)O)N", "C(=N)(N)O", "C(C(=O)O)N",
			"C(=O)C(C(C(C(CO)O)O)O)O", "C(=O)C(C(C(C(CO)O)O)O)O", "c1cc(c(cc1/C=C/C(=O)O)O)O",
			"CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@H](CC[C@]4(C)[C@H]3CC[C@]12C)O", "C(=C\\C(=O)O)/C=C/C(=O)O",
			"c12c([nH]c(n1)O)nc(nc2O)O", "CC(C)CCCC(C)C1CCC2C3CC(C4(CC(CCC4(C)C3CCC12C)O)O)O",
			"c1cc(c(C(=O)O)nc1)C(=O)O", "CC(C(C(=O)O)N)O", "C(=C/C(=O)O)/C(=O)O", "C(=C\\C(=O)O)/C(=O)O",
			"C(=O)C(C(C(C(CO)O)O)O)O", "C(=O)C(C(C(C(CO)O)O)O)O", "c1cc(cc(c1)O)C(=O)O",
			"C(C(C(C(=O)O)O)O)(C(C(=O)O)O)O", "c1cc(cnc1)C(=N)O", "CCCCCCCCC(=O)O", "CCCCCCCCCCCCCCCC(=O)O",
			"CSCCC(C(=O)O)N", "CCCCC(C(=O)O)N", "c1ccc(cc1)CC(C(=O)O)N", "c1ccc(cc1)CC(C(=O)O)N", "C1CC(C(=O)O)NC1",
			"C(C(C(=O)O)N)O", "c1ccc2c(c1)c(CC(C(=O)O)N)c[nH]2", "c1ccc2c(c1)c(CC(C(=O)O)N)c[nH]2",
			"c1cc(ccc1CC(C(=O)O)N)O", "c1cc(ccc1CC(C(=O)O)N)O", "C(C(C(=O)O)N)C(=O)O", "C(C(C(=O)O)N)C(=O)O",
			"C(CC(=O)O)C(C(=O)O)N", "CCCCCCCCCCCCCCCCCCCCCCCC(=O)O", "c1nc2c([nH]1)ncnc2O", "c1nc2c([nH]1)ncnc2O",
			"C(C(=O)O)O", "C(C(CO)O)O", "C(CN)C(C(=O)O)N", "C=C(CC(=O)O)C(=O)O", "c1ccc(cc1)C(=NCC(=O)O)O",
			"CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O", "C(C(CO)OP(=O)(O)O)O", "CCCCCCCCCC(=O)O", "CCCCCCCC(=O)O",
			"c1cc(c(cc1[C@@H]2[C@H](Cc3c(cc(cc3O2)O)O)O)O)O", "CC(C)CCCC(C)C1CCC2/C(=C/C=C/3\\CC(CCC3=C)O)/CCCC12C",
			"c1ccc(cc1)C2=CC(=O)c3c(cc(cc3O2)O)O", "C1CCC(C1)(C(=O)O)N", "C(CC(C(=O)O)N)CC(=O)O", "CC(CC(=O)O)O",
			"C(=O)C(C(C(C(C(=O)O)O)O)O)O", "C(=O)C(C(C(C(C(=O)O)O)O)O)O", "C(C(C(=O)O)O)OP(=O)(O)O", "C(C(C(=O)O)O)O",
			"c1ccc(cc1)C2=CC(=O)c3ccccc3O2", "c1cc(c(cc1O)C(=O)O)O", "CC(=CCC/C(=C/CO)/C)C",
			"C(C1C(C(C(C(=O)O1)O)O)O)O", "C(C1C(C(C(C(=O)O1)O)O)O)O", "C(C1C(C(C(C(=O)O1)O)O)O)O", "C(CC=O)CC=O",
			"C(CC=O)CC=O", "C(CCC(=O)O)CC(=O)O", "c1ccnc(c1)C(=O)O", "c1ccc(cc1)CN", "CC(CO)(CO)N",
			"CC(CCC(=O)O)C(=O)O", "c1cc(ncc1C(=O)O)O", "CCCCCCB(OC)OC", "CS(Cl)(=O)=O",
			"CC(C)OC(C(C)(C)N(=NOC(C)C)=O)=O", "CC(C)OC(C(C)N(=NOC(C)C)=O)=O", "C(I)(I)I", "S1SSSSSSS1", "CC(S)S",
			"CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC[C@H]4C[C@H](CC[C@]4(C)[C@H]3CC[C@]12C)O", "CC(CC(=O)O)CC(=O)O",
			"c1cc(ccc1C(=O)O)N", "C1=C(CC(C(C1O)O)O)C(=O)O", "c1cc2c(cc1O)cc(C(=O)O)[nH]2",
			"c1cc2c(cc1O)cc(C(=O)O)[nH]2", "CC(C)(CC(=O)O)C(=O)O", "CCCCCCCCCCCCCCCCCCCCCC(=O)O", "C(CC(=O)O)CO",
			"CC(CC(=O)O)(CC(=O)O)O", "c1ccc(cc1)C2=C(C(=O)c3ccccc3O2)O", "CC(C(C)C(=O)O)C(=O)O", "c1ccc(cc1)CO",
			"c1ccc(c(c1)C(=NCC(=O)O)O)O", "CCCCCCCCCCCCCC(=O)O", "CCCCCCCCCCCC(CC(=O)O)O", "C1(C(C(C(C(C1O)O)O)O)O)O",
			"CCCCCCCCCCCCCCCCCCCC(=O)O", "c1cc(c(cc1C(CN)O)O)O", "CNCC(c1ccc(c(c1)O)O)O", "CNCC(c1ccc(c(c1)O)O)O",
			"COc1cc(/C=C/CO)cc(c1O)OC", "COc1cc(/C=C/CO)cc(c1O)OC", "c1ccc2c(c1)ccc3ccccc32", "c1ccc(cc1)C(=N)O",
			"c1cc(cc(c1)O)CC(=O)O", "c1cc(ccc1/C=C/C(=O)O)O", "c1cc2c(cc1O)c(CC(=O)O)c[nH]2",
			"COc1c2C=CC(=O)Oc2cc3c1cco3", "C1CCNC(C1)C(=O)O", "CC(=NCCc1c[nH]c2ccc(cc12)O)O",
			"c1cc2c(cc1O)c(CC(=O)O)c[nH]2", "CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=NCCO)O",
			"CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)OC", "CC/C=C\\C/C=C\\C/C=C\\CCCCCCCC(=O)O", "C(CCN)CN",
			"COc1cc(/C=C/CO)ccc1O", "C1CNC(C1O)C(=O)O", "C1CNC(C1O)C(=O)O", "C(=C\\C(=O)O)/CC(=O)O",
			"C([C@@H]1[C@H]([C@H]([C@H](n2cnc3c(N)ncnc32)O1)O)OP(=O)(O)O)O", "c1cc(c(cc1C(CO)O)O)O", "C(COP(=O)(O)O)N",
			"c1cc(cc(c1)O)O", "c1cc(c(cc1C[C@@H](C(=O)O)N)O)O", "c1cc2c(cc1O)c(CC(C(=O)O)N)c[nH]2",
			"c1ccc2c(c1)c(CC#N)c[nH]2", "c1ccc2c(c1)c(CC(C(=O)O)O)c[nH]2", "c1ccc(c(c1)C(=O)O)C(=O)O",
			"CCCCC(CC)COC(=O)c1ccccc1C(=O)O", "c1ccc(cc1)CSC#N", "COc1ccccc1O", "C(CP(=O)(O)O)N",
			"c1cc(c(cc1C2C(C(=O)c3c(cc(cc3O2)O)O)O)O)O", "c1cc(ccc1CCO)O", "C(CCC(=N)O)CC1CCSS1", "C(CCC(=N)O)CC1CCSS1",
			"c1cc(cnc1)C2CCCN2", "c1cc(c(cc1C(C(=O)O)O)O)O", "c1ccc2c(c1)C(=O)c3cccc(c3C2=O)O", "C(CCCCC(=O)O)CCCCO",
			"CN1C2CCC1CC(C2)OC(=O)C(CO)c3ccccc3", "c1ccc2c(c1)c(CC(=N)O)c[nH]2", "c1ccc2c(c1)c(CC(=N)O)c[nH]2",
			"c1ccc2c(c1)c(CC(=O)O)c[nH]2", "COc1c2c(cco2)cc3C=CC(=O)Oc31", "COc1c2c(cco2)cc3C=CC(=O)Oc31", "C(CO)N",
			"CCCCCCB(OC)OC", "CS(Cl)(=O)=O", "CC(C)OC(C(C)(C)N(=NOC(C)C)=O)=O", "CC(C)OC(C(C)N(=NOC(C)C)=O)=O",
			"C(I)(I)I", "S1SSSSSSS1", "CC(S)S", "CC(C)C1=CC2=CCC3C(C)(CCCC3(C)C(=O)O)C2CC1",
			"CC(C)C1=CC2=CCC3C(C)(CCCC3(C)C(=O)O)C2CC1", "C(C(C(=O)O)N)C(=N)O", "C(C(C(=O)O)N)C(=N)O",
			"C(C(C(=O)O)N)C(=N)O", "COc1cc(ccc1O)CC(=O)O", "CC(=NCCc1c[nH]c2ccc(cc12)OC)O", "c1ccc(cc1)C(C(=O)O)O",
			"COc1ccc(cc1)C2=COc3cc(ccc3C2=O)O", "CNCC(=O)O", "C(CC(=O)O)CN", "C(CC(=O)O)CN", "CC(=O)Oc1ccccc1C(=O)O",
			"c1ccc(cc1)N", "C(CO)C(CO)O", "C(C(CO)O)C(=O)O", "COc1ccc(cc1)C(=O)O", "CCCCCCCC(=O)OC", "CCCCCCCCC(=O)OC",
			"CCCCCCCCCC(=O)OC", "CCCCCCCCCCCC(=O)OC", "CCCCCCCCCCCCCC(=O)OC", "CCCCCCCCCCCCCCCC(=O)OC",
			"CCCCCCCCCCCCCCCCCC(=O)OC", "CCCCCCCCCCCCCCCCCCCC(=O)OC", "CCCCCCCCCCCCCCCCCCCCCC(=O)OC",
			"CCCCCCCCCCCCCCCCCCCCCCCC(=O)OC", "CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)OC", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)OC",
			"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)OC", "C(CCN)CN", "CC[C@H](C)[C@@H](C(=O)O)N",
			"C(C(=O)O)C(CC(=O)O)(C(=O)O)O", "c1cc(ccc1CCN)O", "c1nc2c([nH]1)nc(nc2O)O", "C(=O)(C(=O)O)O", "OP(=O)(O)O",
			"Cc1cnc(nc1O)O", "CCCCCCCCCCCCCCCCCCCCCCCC", "CC(C)c1ccc(C)cc1O",
			"CC(C)CCCC(C)CCCC(C)CCCC1(C)CCc2c(C)c(c(C)c(C)c2O1)OC(=O)C", "CC(C(=O)O)O",
			"C(=C(/CC(=O)O)\\C(=O)O)/C(=O)O", "C(CCN)C[C@@H](C(=O)O)N", "CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=O)O",
			"COP(=O)(O)O", "NO", "C1(C(C(C(C(C1O)O)O)O)O)O", "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)OP(=O)(O)O)O)O)O)O",
			"CCCCCCCCCCCCCCCCCCO", "C(C(=O)O)(C(=O)O)N", "CCCCCCCCCCCCO", "CCCCCCCC/C=C\\CCCCCCCCCC(=O)O",
			"CCCCCCCCCCCCCCCCCC(=O)OCC(CO)O", "CCCCCCCCCCCCCCCC(=O)OCC(CO)O", "CCCCCCCC/C=C/CCCCCCCC(=O)O",
			"CCCCCCCCCCCCCCCC(=O)OC(CO)CO", "CCCCCCCCCCCCCCC(=O)O", "c1ccc(cc1)C(=O)O", "c1cc(cnc1)C(=O)O",
			"C(CCCC(=O)O)CCCC(=O)O", "C(C(=O)O)N", "CCCC(C(=O)O)O", "C(CC(=O)O)CC(=O)O", "C(CS(=O)(=O)O)N",
			"CCCCCCB(OC)OC", "CS(Cl)(=O)=O", "CC(C)OC(C(C)(C)N(=NOC(C)C)=O)=O", "CC(C)OC(C(C)N(=NOC(C)C)=O)=O",
			"C(I)(I)I", "S1SSSSSSS1", "CC(S)S", "COC1=C(C=CC(=C1)C=O)O", "COc1cc(ccc1O)C=O", "CCCC", "C(C)CC" };
	private static final String[] descriptorNames = { "fragC", "C1SP1", "C2SP1", "C1SP2", "C2SP2", "C3SP2", "C1SP3",
			"C2SP3", "C3SP3", "C4SP3", "SCH-3", "SCH-4", "SCH-5", "SCH-6", "SCH-7", "VCH-3", "VCH-4", "VCH-5", "VCH-6",
			"VCH-7", "nAtomLC", "BCUTw-1l", "ATSc2", "ATSc1", "BCUTp-1l", "BCUTc-1h", "ATSc5", "BCUTp-1h", "ATSc4",
			"ATSc3", "BCUTw-1h", "BCUTc-1l", "Kier3" };

	public void testPrecomputeAndGet() throws CDKException {
		float[][] descriptorsNonScaled = new float[smiles.length][];
		float[][] descriptorsScaled = new float[smiles.length][];

		HashSet<String> smilesSet = new HashSet<String>();
		for (int j = 0; j < smiles.length; j++) {
			smilesSet.add(smiles[j]);
		}

		float[] min = new float[descriptorNames.length];
		float[] max = new float[descriptorNames.length];

		for (int i = 0; i < descriptorNames.length; i++) {
			min[i] = Float.POSITIVE_INFINITY;
			max[i] = Float.NEGATIVE_INFINITY;
		}
		for (int j = 0; j < smiles.length; j++) {
			String s = Chemoinformatics.canonical(smiles[j], true);
			float[] desc = Chemoinformatics.descriptors(s, descriptorNames);
			Assert.assertEquals(desc.length, descriptorNames.length);
			for (int i = 0; i < desc.length; i++) {
				if (desc[i] < min[i]) {
					min[i] = desc[i];
				}
				if (desc[i] > max[i]) {
					max[i] = desc[i];
				}
			}
			descriptorsNonScaled[j] = desc;
		}
		for (int i = 0; i < descriptorsNonScaled.length; i++) {
			for (int j = 0; j < descriptorsNonScaled[i].length; j++) {
				Assert.assertEquals(true,
						((descriptorsNonScaled[i][j] <= max[j]) || (Float.isNaN(descriptorsNonScaled[i][j]))));
				Assert.assertEquals(true,
						((descriptorsNonScaled[i][j] >= min[j]) || (Float.isNaN(descriptorsNonScaled[i][j]))));
			}
		}
		for (int j = 0; j < smiles.length; j++) {
			String s = Chemoinformatics.canonical(smiles[j], true);
			float[] desc = Chemoinformatics.descriptors(s, descriptorNames, min, max);
			descriptorsScaled[j] = desc;
		}
		for (int j = 0; j < smiles.length; j++) {
			float[] desc = Chemoinformatics.scaleMinMax(descriptorsNonScaled[j], min, max);
			for (int i = 0; i < desc.length; i++) {
				Assert.assertEquals(descriptorsScaled[j][i], desc[i]);
			}
		}

		Descriptors d = Descriptors.instance(descriptorNames, min, max, true);
		d.precompute(smilesSet, false);
		for (int j = 0; j < smiles.length; j++) {
			String s = smiles[j];
			float[] desc = d.get(s);
			for (int i = 0; i < desc.length; i++) {
				Assert.assertEquals(descriptorsScaled[j][i], desc[i]);
			}
		}
		d = null;
		d = Descriptors.instance(descriptorNames, min, max, false);
		for (int j = 0; j < smiles.length; j++) {
			String s = smiles[j];
			float[] desc = d.get(s);
			for (int i = 0; i < desc.length; i++) {
				Assert.assertEquals(descriptorsScaled[j][i], desc[i]);
			}
		}
		d = null;
		d = Descriptors.instance(descriptorNames, new float[descriptorNames.length], new float[descriptorNames.length],
				true);
		d.precompute(smilesSet, true);
		for (int j = 0; j < smiles.length; j++) {
			String s = smiles[j];
			float[] desc = d.get(s);
			for (int i = 0; i < desc.length; i++) {
				Assert.assertEquals(descriptorsScaled[j][i], desc[i]);
			}
		}
		d = null;
		d = Descriptors.instance(descriptorNames);
		d.precompute(smilesSet, true);
		for (int j = 0; j < smiles.length; j++) {
			String s = smiles[j];
			float[] desc = d.get(s);
			for (int i = 0; i < desc.length; i++) {
				Assert.assertEquals(descriptorsScaled[j][i], desc[i]);
			}
		}
	}

	public void testSaveReadToFromFile() throws CDKException, IOException {
		HashSet<String> smilesSet = new HashSet<String>();
		for (int j = 0; j < smiles.length; j++) {
			smilesSet.add(smiles[j]);
		}

		Descriptors d = Descriptors.instance(descriptorNames);
		d.precompute(smilesSet, true);
		d.saveToFile("test.txt");
		Descriptors d2 = Descriptors.readFromFile("test.txt");
		BufferedReader inp = new BufferedReader(new InputStreamReader(new FileInputStream(new File("test.txt"))));
		String s = inp.readLine().replace("33 324 fragC", "33 0 fragC");
		inp.close();
		FileWriter fw = new FileWriter("test.txt");
		fw.write(s);
		fw.close();
		Descriptors d3 = Descriptors.readFromFile("test.txt");
		for (int i = 0; i < smiles.length; i++) {
			float[] desc1 = d.get(smiles[i]);
			float[] desc2 = d2.get(smiles[i]);
			float[] desc3 = d3.get(smiles[i]);
			for (int j = 0; j < desc1.length; j++) {
				Assert.assertEquals(desc1[j], desc2[j]);
				Assert.assertEquals(desc1[j], desc3[j]);
			}
		}
	}

	public void testNaNs() throws CDKException {
		float[] a = new float[] { Float.NaN, 1, 2, 3, Float.NaN, 3, 4, 6, Float.NaN };
		float[] b = new float[] { 1, 2, 3, Float.NaN, 3, 4, 6, Float.NaN };
		float[] c = new float[] { Float.NaN, 1, 2, 3, 3, 4, 6, Float.NaN };
		float[] d = new float[] { Float.NaN, 1, 2, 3, Float.NaN, 3, 4, 6 };
		float[] e = new float[] { 1, 2, 3, Float.NaN, 3, 4, 6 };
		float[] f = new float[] { 1, 2, 3, 11, 3, 4, 6 };
		float[] g = new float[] { 0, Float.NaN, 1, 2, 3, Float.NaN, 3, 4, 6 };

		Assert.assertEquals(3, Descriptors.countNaNs(a));
		Assert.assertEquals(2, Descriptors.countNaNs(b));
		Assert.assertEquals(2, Descriptors.countNaNs(c));
		Assert.assertEquals(2, Descriptors.countNaNs(d));
		Assert.assertEquals(1, Descriptors.countNaNs(e));
		Assert.assertEquals(0, Descriptors.countNaNs(f));
		Assert.assertEquals(2, Descriptors.countNaNs(g));
		for (int i = 0; i < a.length; i++) {
			Assert.assertEquals(true, ((Descriptors.nansToZero(a)[i] == a[i]) || (Float.isNaN(a[i]))));
		}
		for (int i = 0; i < g.length; i++) {
			Assert.assertEquals(true, ((Descriptors.nansToZero(g)[i] == g[i]) || (Float.isNaN(g[i]))));
		}
		Assert.assertEquals(a.length, Descriptors.nansToZero(a).length);
		Assert.assertEquals(g.length, Descriptors.nansToZero(g).length);
		Assert.assertEquals(0, Descriptors.nansToZero(new float[] {}).length);
		Assert.assertEquals(1, Descriptors.nansToZero(new float[] { Float.NaN }).length);
		Assert.assertEquals(0.0F, Descriptors.nansToZero(new float[] { Float.NaN })[0]);
		HashSet<String> smilesSet = new HashSet<String>();
		for (int j = 0; j < smiles.length; j++) {
			smilesSet.add(smiles[j]);
		}
		Descriptors de = Descriptors.instance(descriptorNames);
		de.precompute(smilesSet, true);
		for (int i = 0; i < smiles.length; i++) {
			float[] desc1 = de.get(smiles[i]);
			float[] desc2 = de.getNoNaNs(smiles[i]);

			for (int j = 0; j < desc1.length; j++) {
				Assert.assertEquals(true, ((desc1[j] == desc2[j]) || (Float.isNaN(desc1[j]) && (desc2[j] == 0.0F))));
			}
		}
	}
}
