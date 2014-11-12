/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.bam;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Before;

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

//    /**
//     * Test of getSingleMappedSamRecords method, of class MappedSamRecords.
//     */
//    @Test
//    public void testGetSingleMappedSamRecords() throws IOException
//    {
//        System.out.println("getSingleMappedSamRecords");
//        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
//        File fastqIn = new File("test/test_data_in/piculus_test_left.fastq");
//        File fastqOut = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
//        boolean isRightHandedReads = false;
//        MappedSamRecords instance = new MappedSamRecords();
//        instance.getSingleMappedSamRecords(bamFile, fastqIn, fastqOut, isRightHandedReads);
//        long lineCount = FastqQCTest.countLines(fastqOut);
//        //..so there should be 8000 lines in the output
//        assertEquals(lineCount, 0);
//    }
//
//    /**
//     * Test of getSingleUnmappedSamRecords method, of class MappedSamRecords.
//     */
//    @Test
//    public void testGetSingleUnmappedSamRecords() throws IOException
//    {
//        System.out.println("getSingleUnmappedSamRecords");
//        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
//        File fastqIn = new File("test/test_data_in/piculus_test_left.fastq");
//        File fastqOut = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
//        boolean isRightHandedReads = false;
//        MappedSamRecords instance = new MappedSamRecords();
//        //none of the 2000 fastq reads should be mapped...
//        instance.getSingleUnmappedSamRecords(bamFile, fastqIn, fastqOut, isRightHandedReads);
//        long lineCount = FastqQCTest.countLines(fastqOut);
//        //..so there should be 8000 lines in the output
//        assertEquals(lineCount, 8000);
//
//        isRightHandedReads = true;
//        fastqIn = new File("test/test_data_in/piculus_test_right.fastq");
//        fastqOut = new File("test/test_data_out/piculus_test_unmapped_right.fastq");
//        instance.getSingleUnmappedSamRecords(bamFile, fastqIn, fastqOut, isRightHandedReads);
//        lineCount = FastqQCTest.countLines(fastqOut);
//        assertEquals(lineCount, 8000);
//    }
//
//    /**
//     * Test of getOnePairedMappedSamRecords method, of class MappedSamRecords.
//     */
//    @Test
//    public void testGetOnePairedMappedSamRecords() throws IOException
//    {
//        System.out.println("getOnePairedMappedSamRecords");
//        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
//        File fastqInLeft = new File("test/test_data_in/piculus_test_left.fastq");
//        File fastqInRight = new File("test/test_data_in/piculus_test_right.fastq");
//        File fastqOutLeft = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
//        File fastqOutRight = new File("test/test_data_out/piculus_test_unmapped_right.fastq");
//        MappedSamRecords instance = new MappedSamRecords();
//        instance.getOnePairedMappedSamRecords(bamFile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
//
//        long lineCount = FastqQCTest.countLines(fastqOutLeft);
//        assertEquals(lineCount, 0);
//        lineCount = FastqQCTest.countLines(fastqOutRight);
//        assertEquals(lineCount, 0);
//
//    }
//
//    /**
//     * Test of getReads method, of class MappedSamRecords.
//     */
//    @Test
//    public void testGetReads()
//    {
//        System.out.println("getReads");
//        SAMFileReader samReader = new SAMFileReader(new File("test/test_data_in/piculus_test_unmapped.bam"));
//        SAMRecordIterator iterator = samReader.iterator();
//        boolean getMapped = true;
//        MappedSamRecords instance = new MappedSamRecords();
//        HashSet result = instance.getReads(samReader, iterator, getMapped);
//        assertEquals(result.size(), 0);
//        iterator.close();
//        iterator = samReader.iterator();
//        getMapped = false;
//        result = instance.getReads(samReader, iterator, getMapped);
//        assertEquals(result.size(), 2000);
//
//    }
//
//    /**
//     * Test of getOnePairedUnmappedSamRecords method, of class MappedSamRecords.
//     */
//    @Test
//    public void testGetOnePairedUnmappedSamRecords() throws IOException
//    {
//        System.out.println("getOnePairedUnmappedSamRecords");
//        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
//        File fastqInLeft = new File("test/test_data_in/piculus_test_left.fastq");
//        File fastqInRight = new File("test/test_data_in/piculus_test_right.fastq");
//        File fastqOutLeft = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
//        File fastqOutRight = new File("test/test_data_out/piculus_test_unmapped_right.fastq");
//        MappedSamRecords instance = new MappedSamRecords();
//        instance.getOnePairedUnmappedSamRecords(bamFile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
//        long lineCount = FastqQCTest.countLines(fastqOutLeft);
//        assertEquals(lineCount, 8000);
//        lineCount = FastqQCTest.countLines(fastqOutRight);
//        assertEquals(lineCount, 8000);
//    }
//
//    /**
//     * Test of getBothPairedMappedSamRecords method, of class MappedSamRecords.
//     */
//    @Test
//    public void testGetBothPairedMappedSamRecords() throws IOException
//    {
//        System.out.println("getBothPairedMappedSamRecords");
//        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
//        File fastqInLeft = new File("test/test_data_in/piculus_test_left.fastq");
//        File fastqInRight = new File("test/test_data_in/piculus_test_right.fastq");
//        File fastqOutLeft = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
//        File fastqOutRight = new File("test/test_data_out/piculus_test_unmapped_right.fastq");
//        MappedSamRecords instance = new MappedSamRecords();
//        instance.getBothPairedMappedSamRecords(bamFile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
//        long lineCount = FastqQCTest.countLines(fastqOutLeft);
//        assertEquals(lineCount, 0);
//        lineCount = FastqQCTest.countLines(fastqOutRight);
//        assertEquals(lineCount, 0);
//    }
//
//    /**
//     * Test of getBothPairedUnmappedSamRecords method, of class
//     * MappedSamRecords.
//     */
//    @Test
//    public void testGetBothPairedUnmappedSamRecords() throws IOException
//    {
//        System.out.println("getBothPairedUnmappedSamRecords");
//        File bamFile = new File("test/test_data_in/piculus_test_unmapped.bam");
//        File fastqInLeft = new File("test/test_data_in/piculus_test_left.fastq");
//        File fastqInRight = new File("test/test_data_in/piculus_test_right.fastq");
//        File fastqOutLeft = new File("test/test_data_out/piculus_test_unmapped_left.fastq");
//        File fastqOutRight = new File("test/test_data_out/piculus_test_unmapped_right.fastq");
//        MappedSamRecords instance = new MappedSamRecords();
//        instance.getBothPairedUnmappedSamRecords(bamFile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
//        long lineCount = FastqQCTest.countLines(fastqOutLeft);
//        assertEquals(lineCount, 8000);
//        lineCount = FastqQCTest.countLines(fastqOutRight);
//        assertEquals(lineCount, 8000);
//    }

    @Before
    public void setUp() throws Exception
    {
    }

    @After
    public void tearDown() throws Exception
    {
    }

    /**
     * Test of listEitherPairedReadUnmappedFromBam method, of class MappedSamRecords.
     */
    @Test
    public void testListPairedUnmappedReadsFromBam()
    {
        System.out.println("listPairedUnmappedReadsFromBam");
        File bamFile = null;
        MappedSamRecords instance = new MappedSamRecords();
        HashSet expResult = null;
        HashSet result = instance.listEitherPairedReadUnmappedFromBam(bamFile);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of listEitherPairedReadMappedFromBam method, of class MappedSamRecords.
     */
    @Test
    public void testListPairedMappedReadsFromBam()
    {
        System.out.println("listPairedMappedReadsFromBam");
        File bamFile = null;
        MappedSamRecords instance = new MappedSamRecords();
        HashSet expResult = null;
        HashSet result = instance.listEitherPairedReadMappedFromBam(bamFile);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of listSinglePairedReadUnmappedFromBam method, of class MappedSamRecords.
     */
    @Test
    public void testListSingleUnmappedPairedReadsFromBam()
    {
        System.out.println("listSingleUnmappedReadsFromBam");
        File bamFile = null;
        MappedSamRecords instance = new MappedSamRecords();
        HashMap expResult = null;
        HashMap result = instance.listSinglePairedReadUnmappedFromBam(bamFile);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of listMappedSingleEndReadsFromBam method, of class MappedSamRecords.
     */
    @Test
    public void testListSingleMappedReadsFromBam()
    {
        System.out.println("listSingleMappedReadsFromBam");
        File bamFile = null;
        MappedSamRecords instance = new MappedSamRecords();
        HashMap expResult = null;
        HashSet result = instance.listMappedSingleEndReadsFromBam(bamFile);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of writePairedReadsFromHashSet method, of class MappedSamRecords.
     */
    @Test
    public void testWritePairedReadsFromHashSet()
    {
        System.out.println("writePairedReadsFromHashSet");
        HashSet list = null;
        File fastqInLeft = null;
        File fastqInRight = null;
        File fastqOutLeft = null;
        File fastqOutRight = null;
        MappedSamRecords instance = new MappedSamRecords();
        instance.writePairedReadsFromHashSet(list, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of writeSingleReadsFromHashMap method, of class MappedSamRecords.
     */
    @Test
    public void testWriteSingleReadsFromHashMap()
    {
        System.out.println("writeSingleReadsFromHashMap");
        HashMap<String, Boolean> list = null;
        File fastqInLeft = null;
        File fastqInRight = null;
        File fastqSinglesOut = null;
        MappedSamRecords instance = new MappedSamRecords();
        instance.writeSingleReadsFromHashMap(list, fastqInLeft, fastqInRight, fastqSinglesOut);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of printReadsFromBamAndFastq method, of class MappedSamRecords.
     */
    @Test
    public void testPrintReadsFromBamAndFastq()
    {
        System.out.println("printReadsFromBamAndFastq");
        File bamFile = null;
        File fastqInLeft = null;
        File fastqInRight = null;
        MappedSamRecords instance = new MappedSamRecords();
        instance.printReadsFromBamAndFastq(bamFile, fastqInLeft, fastqInRight);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
}