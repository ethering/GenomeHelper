/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fasta;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author ethering
 */
public class FastaMotifFinderTest
{
    
    public FastaMotifFinderTest()
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
     * Test of findMatches method, of class FastaMotifFinder.
     */
    @Test
    public void testFindMatches() throws Exception
    {
        System.out.println("findMatches");
        File fastaFile = new File("test/test_data_in/piculus_test_left.fasta");
        //motif occurs 67 times (in 65 sequences) on the plus strand and 41 times (in 39 sequences) on the minus (revcom) strand = 108 occurences
        String strPattern = "ATG[TA]GA[TA]";
        File motifCounts = new File("test/test_data_out/piculus_test_motifDNACounts.txt");
        File proteinCounts = new File("test/test_data_out/piculus_test_motifProteinCounts.txt");
        int minCount = 0;
        FastaMotifFinder instance = new FastaMotifFinder();
        instance.findMatches(fastaFile, strPattern, motifCounts, proteinCounts, minCount);
        BufferedReader reader = new BufferedReader(new FileReader(motifCounts));
        String text = null;
        int totalCount = 0;
        while ((text = reader.readLine()) != null)
        {
            String [] array = text.split("\t");
            int count = Integer.parseInt(array[1]);
            totalCount+=count;
        }
        //the problem here is that once a match is found, matching stops instead of continuing to see if there are more matches
        assertEquals(totalCount, 108);
    }

    /**
     * Test of revcom method, of class FastaMotifFinder.
     */
    @Test
    public void testRevcom() throws Exception
    {
        System.out.println("revcom");
        String seq = "atcg";
        String revcom = FastaMotifFinder.revcom(seq);
        assertEquals("cgat", revcom);
    }

//    /**
//     * Test of revcomFastaFile method, of class FastaMotifFinder.
//     */
//    @Test
//    public void testRevcomFastaFile() throws Exception
//    {
//        System.out.println("revcomFastaFile");
//        File fastaIn = null;
//        File revcomOut = null;
//        FastaMotifFinder instance = new FastaMotifFinder();
//        instance.revcomFastaFile(fastaIn, revcomOut);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
}