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
import static uk.ac.tsl.etherington.genomehelper.fastq.FastqQCTest.countLines;

/**
 *
 * @author ethering
 */
public class FastqJoinerTest
{
    
    public FastqJoinerTest()
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
     * Test of join method, of class FastqJoiner.
     */
    @Test
    public void testJoin() throws IOException
    {
        System.out.println("join");
        File leftReads = new File("test/test_data_in/piculus_test_left_3missing.fastq");
        File rightReads = new File("test/test_data_in/piculus_test_right.fastq");
        File fastqJoinedFile = new File("test/test_data_out/piculus_test_joined.fastq");
        File fastqSinglesFile = new File("test/test_data_out/piculus_test_joined_singles.fastq");
        FastqJoiner instance = new FastqJoiner();
        instance.join(leftReads, rightReads, fastqJoinedFile, fastqSinglesFile);
        long lineCountJoined = countLines(fastqJoinedFile);
        long lineCountSingles = countLines(fastqSinglesFile);
        assertEquals(lineCountJoined, 7988);
        assertEquals(lineCountSingles, 12);
    }

    /**
     * Test of split method, of class FastqJoiner.
     */
    @Test
    public void testSplit() throws IOException
    {
        System.out.println("split");
        File fastqFile = new File("test/test_data_in/piculus_test_joined.fastq");
        File leftPairdReads = new File("test/test_data_out/piculus_test_split_left.fastq");
        File rightPairedReads = new File("test/test_data_out/piculus_test_split_right.fastq");
        FastqJoiner instance = new FastqJoiner();
        instance.split(fastqFile, leftPairdReads, rightPairedReads);
        long leftCount = countLines(leftPairdReads);
        long rightCount = countLines(rightPairedReads);
        assertEquals(leftCount, 8000);
        assertEquals(rightCount, 8000);
    }
}