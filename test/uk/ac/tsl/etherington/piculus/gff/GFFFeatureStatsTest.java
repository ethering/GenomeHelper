/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.gff;

import java.io.File;
import java.util.HashMap;
import java.util.Set;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author ethering
 */
public class GFFFeatureStatsTest
{
    
    public GFFFeatureStatsTest()
    {
    }
    
    @BeforeClass
    public static void setUpClass()
    {
    }
    
    @AfterClass
    public static void tearDownClass()
    {
    }

    /**
     * Test of getStats method, of class GFFFeatureStats.
     */
    @Test
    public void testGetStats() throws Exception
    {
        System.out.println("getStats");
        String gff = "";
        File refSeq = null;
        String attribute = "";
        GFFFeatureStats instance = new GFFFeatureStats();
        instance.getStats(gff, refSeq, attribute);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getGenomeSizeFromIntArrayHashMap method, of class GFFFeatureStats.
     */
    @Test
    public void testGetGenomeSizeFromIntArrayHashMap()
    {
        System.out.println("getGenomeSizeFromIntArrayHashMap");
        HashMap<String, int[]> genomeMap = null;
        GFFFeatureStats instance = new GFFFeatureStats();
        double expResult = 0.0;
        double result = instance.getGenomeSizeFromIntArrayHashMap(genomeMap);
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFeatureList method, of class GFFFeatureStats.
     */
    @Test
    public void testGetFeatureList() throws Exception
    {
        System.out.println("getFeatureList");
        String gffFile = "";
        GFFFeatureStats instance = new GFFFeatureStats();
        FeatureList expResult = null;
        FeatureList result = instance.getFeatureList(gffFile);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getMeanFeatureLength method, of class GFFFeatureStats.
     */
    @Test
    public void testGetMeanFeatureLength_3args() throws Exception
    {
        System.out.println("getMeanFeatureLength");
        FeatureList fl = null;
        HashMap<String, int[]> genomeMap = null;
        String featureName = "";
        GFFFeatureStats instance = new GFFFeatureStats();
        double expResult = 0.0;
        double result = instance.getMeanFeatureLength(fl, genomeMap, featureName);
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getMeanFeatureLength method, of class GFFFeatureStats.
     */
    @Test
    public void testGetMeanFeatureLength_5args() throws Exception
    {
        System.out.println("getMeanFeatureLength");
        FeatureList fl = null;
        HashMap<String, int[]> genomeMap = null;
        String featureName = "";
        File geneIds = null;
        String attribute = "";
        GFFFeatureStats instance = new GFFFeatureStats();
        double expResult = 0.0;
        double result = instance.getMeanFeatureLength(fl, genomeMap, featureName, geneIds, attribute);
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getMeanFeatureLength method, of class GFFFeatureStats.
     */
    @Test
    public void testGetMeanFeatureLength_FeatureList_String() throws Exception
    {
        System.out.println("getMeanFeatureLength");
        FeatureList fl = null;
        String featureName = "";
        GFFFeatureStats instance = new GFFFeatureStats();
        double expResult = 0.0;
        double result = instance.getMeanFeatureLength(fl, featureName);
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of createNonCodingGenome method, of class GFFFeatureStats.
     */
    @Test
    public void testCreateNonCodingGenome() throws Exception
    {
        System.out.println("createNonCodingGenome");
        FeatureList fl = null;
        File refSeq = null;
        File nonCodingGenome = null;
        GFFFeatureStats instance = new GFFFeatureStats();
        instance.createNonCodingGenome(fl, refSeq, nonCodingGenome);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of createCodingGenome method, of class GFFFeatureStats.
     */
    @Test
    public void testCreateCodingGenome() throws Exception
    {
        System.out.println("createCodingGenome");
        FeatureList fl = null;
        String feature = "";
        File refSeq = null;
        File codingGenome = null;
        GFFFeatureStats instance = new GFFFeatureStats();
        instance.createCodingGenome(fl, feature, refSeq, codingGenome);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getMeanIntronLength method, of class GFFFeatureStats.
     */
    @Test
    public void testGetMeanIntronLength() throws Exception
    {
        System.out.println("getMeanIntronLength");
        FeatureList fl = null;
        String attribute = "";
        double genomeSize = 0.0;
        GFFFeatureStats instance = new GFFFeatureStats();
        double expResult = 0.0;
        double result = instance.getMeanIntronLength(fl, attribute, genomeSize);
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getMeanSecretedIntronLength method, of class GFFFeatureStats.
     */
    @Test
    public void testGetMeanSecretedIntronLength() throws Exception
    {
        System.out.println("getMeanSecretedIntronLength");
        FeatureList fl = null;
        String attribute = "";
        File secretedProteinsFile = null;
        GFFFeatureStats instance = new GFFFeatureStats();
        double expResult = 0.0;
        double result = instance.getMeanSecretedIntronLength(fl, attribute, secretedProteinsFile);
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of calculateCodingRegion method, of class GFFFeatureStats.
     */
    @Test
    public void testCalculateCodingRegion() throws Exception
    {
        System.out.println("calculateCodingRegion");
        FeatureList featList = null;
        HashMap<String, int[]> codingMap = null;
        String attribute = "";
        double genomeLength = 0.0;
        GFFFeatureStats instance = new GFFFeatureStats();
        double expResult = 0.0;
        double result = instance.calculateCodingRegion(featList, codingMap, attribute, genomeLength);
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getBlocks method, of class GFFFeatureStats.
     */
    @Test
    public void testGetBlocks_FeatureList_String()
    {
        System.out.println("getBlocks");
        FeatureList fl = null;
        String featureType = "";
        GFFFeatureStats instance = new GFFFeatureStats();
        HashMap expResult = null;
        HashMap result = instance.getBlocks(fl, featureType);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getBlocks method, of class GFFFeatureStats.
     */
    @Test
    public void testGetBlocks_3args()
    {
        System.out.println("getBlocks");
        FeatureList fl = null;
        String featureType = "";
        String attribute = "";
        GFFFeatureStats instance = new GFFFeatureStats();
        HashMap expResult = null;
        HashMap result = instance.getBlocks(fl, featureType, attribute);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getBlocks method, of class GFFFeatureStats.
     */
    @Test
    public void testGetBlocks_4args()
    {
        System.out.println("getBlocks");
        FeatureList fl = null;
        String featureType = "";
        String attribute = "";
        Set secretedProteins = null;
        GFFFeatureStats instance = new GFFFeatureStats();
        HashMap expResult = null;
        HashMap result = instance.getBlocks(fl, featureType, attribute, secretedProteins);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
}