/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fastq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import static uk.ac.tsl.etherington.piculus.fastq.FastqQCTest.countLines;

/**
 *
 * @author ethering
 */
public class FastqMotifFinderTest
{

    public FastqMotifFinderTest()
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
     * Test of findMatches method, of class FastqMotifFinder.
     */
    @Test
    public void testFindMatches() throws Exception
    {
        System.out.println("findMatches");
        File fastqFile = new File("test/test_data_in/piculus_test_left.fastq");
        //motif occurs 67 times (in 65 sequences) on the plus strand and 41 times (in 39 sequences) on the minus (revcom) strand = 108 occurences
        String strPattern = "ATG[TA]GA[TA]";
        File motifCounts = new File("test/test_data_out/piculus_test_motifDNACounts.txt");
        File proteinCounts = new File("test/test_data_out/piculus_test_motifProteinCounts.txt");
        int minCount = 0;
        FastqMotifFinder instance = new FastqMotifFinder();
        instance.findMatches(fastqFile, strPattern, motifCounts, proteinCounts, minCount);
        BufferedReader reader = new BufferedReader(new FileReader(motifCounts));
        String text = null;
        int totalCount = 0;
        while ((text = reader.readLine()) != null)
        {
            String [] array = text.split("\t");
            int count = Integer.parseInt(array[1]);
            totalCount+=count;
        }
        assertEquals(totalCount, 108);
    }

    /**
     * Test of getPEFastqReadsFromMotif method, of class FastqMotifFinder.
     */
    @Test
    public void testGetPEFastqReadsFromMotif() throws Exception
    {
        System.out.println("getPEFastqReadsFromMotif");
        File fastqIn = new File("test/test_data_in/piculus_test_interlaced.fastq");
        String strPattern = "ATG[TA]GA[TA]";
        File fastqOut = new File("test/test_data_out/piculus_test_motifPairs.fastq");
        FastqMotifFinder instance = new FastqMotifFinder();
        instance.getPEFastqReadsFromMotif(fastqIn, strPattern, fastqOut);
        //motif occurs in 65 sequences on the plus strand and in 39 sequences on the minus (revcom) strand = 104 occurences
        //which means 104 pairs of fastq sequences contain the motif = 208 fastq reads * 4 lines per read = 832
        long lineCount = countLines(fastqOut);
        assertEquals(lineCount, 832);
    }
}