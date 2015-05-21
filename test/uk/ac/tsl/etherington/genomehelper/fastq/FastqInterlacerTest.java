/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fastq;

import java.io.File;
import java.io.IOException;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author ethering
 */
public class FastqInterlacerTest
{

    public FastqInterlacerTest()
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
     * Test of interlace method, of class FastqInterlacer.
     */
    @Test
    public void testInterlace() throws IOException
    {
        System.out.println("interlace");
        File leftReads = new File("test/test_data_in/piculus_test_left_3missing.fastq");
        File rightReads = new File("test/test_data_in/piculus_test_right.fastq");
        File fastqInterlacedFile = new File("test/test_data_out/piculus_test_interlaced_.fastq");
        File fastqSinglesFile = new File("test/test_data_out/piculus_test_interlaced_singles.fastq");
        FastqInterlacer instance = new FastqInterlacer();
        instance.interlace(leftReads, rightReads, fastqInterlacedFile, fastqSinglesFile);
        long interlacedLineCount = FastqQCTest.countLines(fastqInterlacedFile);
        long singlesLineCount = FastqQCTest.countLines(fastqSinglesFile);
        //the interlaced file should contain (2000 - 3)*2*4 = 15976 lines
        assertEquals(interlacedLineCount, 15976);
        //the singles should contain 12 lines
        assertEquals(singlesLineCount, 12);
        
    }
    /**
     * Test of interlaceKnownPairs method, of class FastqInterlacer.
     */
    @Test
    public void testInterlaceKnownPairs() throws Exception
    {
        System.out.println("interlaceKnownPairs");
        File leftReads = new File("test/test_data_in/piculus_test_left.fastq");
        File rightReads = new File("test/test_data_in/piculus_test_right.fastq");
        File fastqInterlacedFile = new File("test/test_data_out/piculus_test_interlaced.fastq");
        FastqInterlacer instance = new FastqInterlacer();
        instance.interlaceKnownPairs(leftReads, rightReads, fastqInterlacedFile);
        long lineCount = FastqQCTest.countLines(fastqInterlacedFile);
        //the interlaced file should contain 2000*2*4=16000 lines
        assertEquals(lineCount, 16000);
    }

    /**
     * Test of deinterlace method, of class FastqInterlacer.
     */
    @Test
    public void testDeinterlace() throws IOException
    {
        System.out.println("deinterlace");
        File fastqFile = new File("test/test_data_in/piculus_test_interlaced_3missing.fastq");
        File leftPairedReads = new File("test/test_data_out/piculus_test_left_fromInterlace.fastq");
        File rightPairedReads = new File("test/test_data_out/piculus_test_right_fromInterlace.fastq");
        File leftSingleReads = new File("test/test_data_out/piculus_test_leftSingles_fromInterlace.fastq");
        File rightSingleReads = new File("test/test_data_out/piculus_test_rightSingles_fromInterlace.fastq");
        FastqInterlacer instance = new FastqInterlacer();
        instance.deinterlace(fastqFile, leftPairedReads, rightPairedReads, leftSingleReads, rightSingleReads);
        long leftLineCount = FastqQCTest.countLines(leftPairedReads);
        long rightLineCount = FastqQCTest.countLines(rightPairedReads);
        long leftSinglesLneCount = FastqQCTest.countLines(leftSingleReads);
        long rightSinglesLneCount = FastqQCTest.countLines(rightSingleReads);
        assertEquals(leftLineCount, 7988);
        assertEquals(rightLineCount, 7988);
        assertEquals(leftSinglesLneCount, 12);
        assertEquals(rightSinglesLneCount, 0);
    }
}