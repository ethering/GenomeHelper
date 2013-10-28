/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fasta;

import java.io.File;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import uk.ac.tsl.etherington.piculus.fastq.FastqQCTest;

/**
 *
 * @author ethering
 */
public class RandomFastaTest
{
    
    public RandomFastaTest()
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
     * Test of selectRandomSequences method, of class RandomFasta.
     */
    @Test
    public void testSelectRandomSequences() throws Exception
    {
        System.out.println("selectRandomSequences");
       File fasta = new File("test/test_data_in/piculus_test_left.fasta");
        int numberOfSeqsRequired = 10;
        File outfile = new File("test/test_data_out/piculus_test_random.fasta");
        RandomFasta instance = new RandomFasta();
        instance.selectRandomSequences(fasta, numberOfSeqsRequired, outfile);
        long noLines = FastqQCTest.countLines(outfile);
        assertEquals(noLines, 20);
    }
}