/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.bam;

import uk.ac.tsl.etherington.genomehelper.bam.MappedSamRecords;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import uk.ac.tsl.etherington.genomehelper.fastq.FastqQCTest;

/**
 *
 * @author ethering
 */
public class MappedSamRecordsTest
{

    public MappedSamRecordsTest()
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
     * Test of getSingleMappedSamRecords method, of class MappedSamRecords.
     */
    @Test
    public void testGetSingleMappedSamRecords() throws IOException
    {
        System.out.println("getSingleMappedSamRecords");
        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
        File fastqIn = new File("test/test_data_in/piculus_test_left.fastq");
        File fastqOut = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
        boolean isRightHandedReads = false;
        MappedSamRecords instance = new MappedSamRecords();
        instance.getSingleMappedSamRecords(bamFile, fastqIn, fastqOut, isRightHandedReads);
        long lineCount = FastqQCTest.countLines(fastqOut);
        //..so there should be 8000 lines in the output
        assertEquals(lineCount, 0);
    }

    /**
     * Test of getSingleUnmappedSamRecords method, of class MappedSamRecords.
     */
    @Test
    public void testGetSingleUnmappedSamRecords() throws IOException
    {
        System.out.println("getSingleUnmappedSamRecords");
        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
        File fastqIn = new File("test/test_data_in/piculus_test_left.fastq");
        File fastqOut = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
        boolean isRightHandedReads = false;
        MappedSamRecords instance = new MappedSamRecords();
        //none of the 2000 fastq reads should be mapped...
        instance.getSingleUnmappedSamRecords(bamFile, fastqIn, fastqOut, isRightHandedReads);
        long lineCount = FastqQCTest.countLines(fastqOut);
        //..so there should be 8000 lines in the output
        assertEquals(lineCount, 8000);

        isRightHandedReads = true;
        fastqIn = new File("test/test_data_in/piculus_test_right.fastq");
        fastqOut = new File("test/test_data_out/piculus_test_unmapped_right.fastq");
        instance.getSingleUnmappedSamRecords(bamFile, fastqIn, fastqOut, isRightHandedReads);
        lineCount = FastqQCTest.countLines(fastqOut);
        assertEquals(lineCount, 8000);
    }

    /**
     * Test of getOnePairedMappedSamRecords method, of class MappedSamRecords.
     */
    @Test
    public void testGetOnePairedMappedSamRecords() throws IOException
    {
        System.out.println("getOnePairedMappedSamRecords");
        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
        File fastqInLeft = new File("test/test_data_in/piculus_test_left.fastq");
        File fastqInRight = new File("test/test_data_in/piculus_test_right.fastq");
        File fastqOutLeft = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
        File fastqOutRight = new File("test/test_data_out/piculus_test_unmapped_right.fastq");
        MappedSamRecords instance = new MappedSamRecords();
        instance.getOnePairedMappedSamRecords(bamFile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);

        long lineCount = FastqQCTest.countLines(fastqOutLeft);
        assertEquals(lineCount, 0);
        lineCount = FastqQCTest.countLines(fastqOutRight);
        assertEquals(lineCount, 0);

    }

    /**
     * Test of getReads method, of class MappedSamRecords.
     */
    @Test
    public void testGetReads()
    {
        System.out.println("getReads");
        SAMFileReader samReader = new SAMFileReader(new File("test/test_data_in/piculus_test_unmapped.bam"));
        SAMRecordIterator iterator = samReader.iterator();
        boolean getMapped = true;
        MappedSamRecords instance = new MappedSamRecords();
        HashSet result = instance.getReads(samReader, iterator, getMapped);
        assertEquals(result.size(), 0);
        iterator.close();
        iterator = samReader.iterator();
        getMapped = false;
        result = instance.getReads(samReader, iterator, getMapped);
        assertEquals(result.size(), 2000);

    }

    /**
     * Test of getOnePairedUnmappedSamRecords method, of class MappedSamRecords.
     */
    @Test
    public void testGetOnePairedUnmappedSamRecords() throws IOException
    {
        System.out.println("getOnePairedUnmappedSamRecords");
        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
        File fastqInLeft = new File("test/test_data_in/piculus_test_left.fastq");
        File fastqInRight = new File("test/test_data_in/piculus_test_right.fastq");
        File fastqOutLeft = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
        File fastqOutRight = new File("test/test_data_out/piculus_test_unmapped_right.fastq");
        MappedSamRecords instance = new MappedSamRecords();
        instance.getOnePairedUnmappedSamRecords(bamFile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
        long lineCount = FastqQCTest.countLines(fastqOutLeft);
        assertEquals(lineCount, 8000);
        lineCount = FastqQCTest.countLines(fastqOutRight);
        assertEquals(lineCount, 8000);
    }

    /**
     * Test of getBothPairedMappedSamRecords method, of class MappedSamRecords.
     */
    @Test
    public void testGetBothPairedMappedSamRecords() throws IOException
    {
        System.out.println("getBothPairedMappedSamRecords");
        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
        File fastqInLeft = new File("test/test_data_in/piculus_test_left.fastq");
        File fastqInRight = new File("test/test_data_in/piculus_test_right.fastq");
        File fastqOutLeft = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
        File fastqOutRight = new File("test/test_data_out/piculus_test_unmapped_right.fastq");
        MappedSamRecords instance = new MappedSamRecords();
        instance.getBothPairedMappedSamRecords(bamFile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
        long lineCount = FastqQCTest.countLines(fastqOutLeft);
        assertEquals(lineCount, 0);
        lineCount = FastqQCTest.countLines(fastqOutRight);
        assertEquals(lineCount, 0);
    }

    /**
     * Test of getBothPairedUnmappedSamRecords method, of class
     * MappedSamRecords.
     */
    @Test
    public void testGetBothPairedUnmappedSamRecords() throws IOException
    {
        System.out.println("getBothPairedUnmappedSamRecords");
        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
        File fastqInLeft = new File("test/test_data_in/piculus_test_left.fastq");
        File fastqInRight = new File("test/test_data_in/piculus_test_right.fastq");
        File fastqOutLeft = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
        File fastqOutRight = new File("test/test_data_out/piculus_test_unmapped_right.fastq");
        MappedSamRecords instance = new MappedSamRecords();
        instance.getBothPairedUnmappedSamRecords(bamFile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
        long lineCount = FastqQCTest.countLines(fastqOutLeft);
        assertEquals(lineCount, 8000);
        lineCount = FastqQCTest.countLines(fastqOutRight);
        assertEquals(lineCount, 8000);
    }
}