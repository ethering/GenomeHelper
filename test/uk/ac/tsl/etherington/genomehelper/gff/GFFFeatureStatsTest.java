/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.gff;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.biojava.bio.BioException;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import uk.ac.tsl.etherington.genomehelper.fasta.FastaFeatures;

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
     * Test of getGenomeSizeFromIntArrayHashMap method, of class
     * GFFFeatureStats.
     */
    @Test
    public void testGetGenomeSizeFromIntArrayHashMap() throws FileNotFoundException, BioException, Exception
    {
        System.out.println("getGenomeSizeFromIntArrayHashMap");
        File refSeq = new File("test/test_data_in/piculus_test_refseq.fasta");
        HashMap<String, int[]> genomeMap = new HashMap<>(FastaFeatures.getSequenceAsHashMapIntArray(refSeq));
        GFFFeatureStats instance = new GFFFeatureStats();
        double result = instance.getGenomeSizeFromIntArrayHashMap(genomeMap);
        assertEquals(200, result, 0.0);

    }

    /**
     * Test of getFeatureList method, of class GFFFeatureStats.
     */
    @Test
    public void testGetFeatureList() throws Exception
    {
        System.out.println("getFeatureList");
        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        GFFFeatureStats instance = new GFFFeatureStats();

        FeatureList result = instance.getFeatureList(gffFile);
        assertEquals(10, result.size());

    }

    /**
     * Test of getMeanFeatureLength method, of class GFFFeatureStats.
     */
    @Test
    public void testGetMeanFeatureLength_3args() throws Exception
    {
        System.out.println("getMeanFeatureLength");
        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        File refSeq = new File("test/test_data_in/piculus_test_refseq.fasta");
        String featureName = "exon";
        GFFFeatureStats gffs = new GFFFeatureStats();
        FeatureList fl = gffs.getFeatureList(gffFile);
        HashMap<String, int[]> genomeMap = new HashMap<>(FastaFeatures.getSequenceAsHashMapIntArray(refSeq));
        double result = gffs.getMeanFeatureLength(fl, genomeMap, featureName);
        assertEquals(10, result, 0.0);
    }

    /**
     * Test of getMeanFeatureLength method, of class GFFFeatureStats.
     */
    @Test
    public void testGetMeanFeatureLength_5args() throws Exception
    {
        System.out.println("getMeanFeatureLength");
        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        File refSeq = new File("test/test_data_in/piculus_test_refseq.fasta");
        File geneIds = new File("test/test_data_in/piculus_test_targets.txt");
        String featureName = "exon";
        String attribute = "Parent";

        GFFFeatureStats gffs = new GFFFeatureStats();
        FeatureList fl = gffs.getFeatureList(gffFile);
        HashMap<String, int[]> genomeMap = new HashMap<>(FastaFeatures.getSequenceAsHashMapIntArray(refSeq));
        double result = gffs.getMeanFeatureLength(fl, genomeMap, featureName, geneIds, attribute);
        assertEquals(10, result, 0.0);

    }
    /**
     * Test of getMeanFeatureLength method, of class GFFFeatureStats.
     */
    @Test
    public void testGetMeanFeatureLength_FeatureList_String() throws Exception
    {
        System.out.println("getMeanFeatureLength");
        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        GFFFeatureStats gffs = new GFFFeatureStats();
        FeatureList fl = gffs.getFeatureList(gffFile);
        String featureName = "exon";

        double result = gffs.getMeanFeatureLength(fl, featureName);
        assertEquals(9, result, 0.0);

    }

    /**
     * Test of createNonCodingGenome method, of class GFFFeatureStats.
     */
    @Test
    public void testCreateNonCodingGenome() throws Exception
    {
        System.out.println("createNonCodingGenome");
        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        File refSeq = new File("test/test_data_in/piculus_test_refseq.fasta");
        File nonCodingGenome = new File("test/test_data_out/piculus_non_coding.fasta");
        double refSeqSize = FastaFeatures.getGenomeSize(refSeq);
        System.out.println("Refseq size: " + refSeqSize);

        GFFFeatureStats gffs = new GFFFeatureStats();
        FeatureList fl = gffs.getFeatureList(gffFile);
        gffs.createNonCodingGenome(fl, refSeq, nonCodingGenome);
        double genomeSize = FastaFeatures.getGenomeSize(nonCodingGenome);
        assertEquals(120, genomeSize, 0.0);
    }

    /**
     * Test of createCodingGenome method, of class GFFFeatureStats.
     */
    @Test
    public void testCreateCodingGenome() throws Exception
    {
        System.out.println("createCodingGenome");
        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        File refSeq = new File("test/test_data_in/piculus_test_refseq.fasta");
        File codingGenome = new File("test/test_data_out/piculus_coding.fasta");
        String feature = "mRNA";

        GFFFeatureStats gffs = new GFFFeatureStats();
        FeatureList fl = gffs.getFeatureList(gffFile);
        GFFFeatureStats instance = new GFFFeatureStats();

        instance.createCodingGenome(fl, feature, refSeq, codingGenome);
        double genomeSize = FastaFeatures.getGenomeSize(codingGenome);
        assertEquals(200, genomeSize, 0.0);
    }

    /**
     * Test of getMeanIntronLength method, of class GFFFeatureStats.
     */
    @Test
    public void testGetMeanIntronLength() throws Exception
    {
        System.out.println("getMeanIntronLength");

        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        String attribute = "Parent";
        File refSeq = new File("test/test_data_in/piculus_test_refseq.fasta");

        GFFFeatureStats gffs = new GFFFeatureStats();
        FeatureList fl = gffs.getFeatureList(gffFile);
        double genomeSize = FastaFeatures.getGenomeSize(refSeq);

        double result = gffs.getMeanIntronLength(fl, attribute, genomeSize);
        assertEquals(20, result, 0.0);
    }

    /**
     * Test of getMeanSecretedIntronLength method, of class GFFFeatureStats.
     */
    @Test
    public void testGetMeanSecretedIntronLength() throws Exception
    {
        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        File secretedProteinsFile = new File("test/test_data_in/piculus_test_targets.txt");
        String attribute = "Parent";
        GFFFeatureStats gffs = new GFFFeatureStats();
        FeatureList fl = gffs.getFeatureList(gffFile);

        double result = gffs.getMeanSecretedIntronLength(fl, attribute, secretedProteinsFile);
        assertEquals(20, result, 0.1);
    }

    /**
     * Test of calculateCodingRegion method, of class GFFFeatureStats.
     */
    @Test
    public void testCalculateCodingRegion() throws Exception
    {
        System.out.println("calculateCodingRegion");
        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        String attribute = "Parent";
        File refSeq = new File("test/test_data_in/piculus_test_refseq.fasta");

        GFFFeatureStats gffs = new GFFFeatureStats();
        HashMap<String, int[]> genomeMap = new HashMap<>(FastaFeatures.getSequenceAsHashMapIntArray(refSeq));
        double genomeSize = gffs.getGenomeSizeFromIntArrayHashMap(genomeMap);
        FeatureList fl = gffs.getFeatureList(gffFile);

        double result = gffs.calculateCodingRegion(fl, genomeMap, attribute, genomeSize);
        assertEquals(1, result, 0.0);
    }

    /**
     * Test of getBlocks method, of class GFFFeatureStats.
     */
    @Test
    public void testGetBlocks_FeatureList_String() throws IOException
    {
        System.out.println("getBlocks");
        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        GFFFeatureStats gffs = new GFFFeatureStats();
        FeatureList fl = gffs.getFeatureList(gffFile);

        String featureType = "exon";
        GFFFeatureStats instance = new GFFFeatureStats();
        HashMap result = instance.getBlocks(fl, featureType);
        //result should just contain one key - 'chr1'
        assertEquals(1, result.size());
        assertEquals(true, result.containsKey("chr1"));
        //get the start/stop sites for chr1
        ArrayList<Integer> blocks = (ArrayList<Integer>) result.get("chr1");
        //as there are 8 exons, this means 8 blocks, each with a start and stop site, so 16
        assertEquals(16, blocks.size());
    }

    /**
     * Test of getBlocks method, of class GFFFeatureStats.
     */
    @Test
    public void testGetBlocks_3args() throws IOException
    {
        System.out.println("getBlocks");
        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        GFFFeatureStats gffs = new GFFFeatureStats();
        FeatureList fl = gffs.getFeatureList(gffFile);
        String attribute = "ID";
        String featureType = "mRNA";
        GFFFeatureStats instance = new GFFFeatureStats();
        HashMap<String, ArrayList<Integer>> result = instance.getBlocks(fl, featureType, attribute);
        //result should just contain two keys - '001' and '002'
        assertEquals(2, result.size());
        assertEquals(true, result.containsKey("001"));
        assertEquals(true, result.containsKey("002"));
        //get all keys - there should be 2 in total
        //and count each start/stop site 
        int noFeatures = 0;
        for (Map.Entry<String, ArrayList<Integer>> entry : result.entrySet())
        {
            ArrayList<Integer> al = entry.getValue();
            noFeatures += al.size();
        }
        //as there are 2 mRNAs, this means 2 blocks, each with a start and stop site, so 4
        assertEquals(4, noFeatures);
    }

    /**
     * Test of getBlocks method, of class GFFFeatureStats.
     */
    @Test
    public void testGetBlocks_4args() throws IOException
    {
        System.out.println("getBlocks");
        String gffFile = new File("test/test_data_in/piculus_test_features.gff").toString();
        File secretedProteinsFile = new File("test/test_data_in/piculus_test_targets.txt");
        Set<String> secretedProteins = new HashSet<>();
        BufferedReader br = new BufferedReader(new FileReader(secretedProteinsFile));
        String line = null; //not declared within while loop
        while ((line = br.readLine()) != null)
        {
            secretedProteins.add(line);
            //System.out.println(line);
        }
        GFFFeatureStats gffs = new GFFFeatureStats();
        FeatureList fl = gffs.getFeatureList(gffFile);
        String attribute = "ID";
        String featureType = "mRNA";
        GFFFeatureStats instance = new GFFFeatureStats();
        HashMap<String, ArrayList<Integer>> result = instance.getBlocks(fl, featureType, attribute, secretedProteins);
        //result should just contain one key - '001'
        assertEquals(1, result.size());
        assertEquals(true, result.containsKey("001"));
        assertEquals(false, result.containsKey("002"));
        //get all keys - there should be 1in total
        //and count each start/stop site 
        int noFeatures = 0;
        for (Map.Entry<String, ArrayList<Integer>> entry : result.entrySet())
        {
            ArrayList<Integer> al = entry.getValue();
            noFeatures += al.size();
        }
        //as there is mRNAs, this means 1 blocks, each with a start and stop site, so 2
        assertEquals(2, noFeatures);
    }
}