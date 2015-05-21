/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fasta;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author ethering
 */
public class FastaFeaturesTest
{

    public FastaFeaturesTest()
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
     * Test of getSequenceLengths method, of class FastaFeatures.
     */
    @Test
    public void testGetSequenceLengths() throws Exception
    {
        System.out.println("getSequenceLengths");
        //contains 2000 fasta seqs, all 80 nts long
        File refSeq = new File("test/test_data_in/piculus_test_left.fasta");
        HashMap mp = FastaFeatures.getSequenceLengths(refSeq);
        assertEquals(mp.size(), 2000);
        Iterator it = mp.entrySet().iterator();
        int length = 0;
        while (it.hasNext())
        {
            Map.Entry pairs = (Map.Entry) it.next();
            length += (int) pairs.getValue();
        }
        assertEquals(length, 160000);
    }

    /**
     * Test of getParsedDNASequences method, of class FastaFeatures.
     */
    @Test
    public void testGetParsedDNASequences() throws Exception
    {
        System.out.println("getParsedDNASequences");
        File refSeq = new File("test/test_data_in/piculus_test_left.fasta");
        //check that we get 10 sequences
        LinkedHashMap result = FastaFeatures.getParsedDNASequences(refSeq);
        assertEquals(result.size(), 2000);
    }

    /**
     * Test of getGenomeSize method, of class FastaFeatures.
     */
    @Test
    public void testGetGenomeSize() throws Exception
    {
        System.out.println("getGenomeSize");
        File refSeq = new File("test/test_data_in/piculus_test_left.fasta");
        //the test gemome size should be 160000
        double result = FastaFeatures.getGenomeSize(refSeq);
        assertEquals(result, 160000, 0.0);
    }

    /**
     * Test of getSequenceAsIntArray method, of class FastaFeatures.
     */
    @Test
    public void testGetSequenceAsIntArray() throws Exception
    {
        System.out.println("getSequenceAsIntArray");
        File refSeq = new File("test/test_data_in/piculus_test_left.fasta");
        HashMap result = FastaFeatures.getSequenceAsIntArray(refSeq);
        assertEquals(result.size(), 2000);
    }

    /**
     * Test of getGCContent method, of class FastaFeatures.
     */
    @Test
    public void testGetGCContent() throws Exception
    {
        System.out.println("getGCContent");
        File refSeq = new File("test/test_data_in/piculus_test_gc.fasta");
        FastaFeatures instance = new FastaFeatures();
        double gcCount = instance.getGCContent(refSeq);
        assertEquals(gcCount, 0.5, 0.0);

    }
}