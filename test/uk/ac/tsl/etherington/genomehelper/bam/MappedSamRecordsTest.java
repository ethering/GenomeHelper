/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.bam;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
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
        String fileName = this.getClass().getResource("/test_data_in/piculus_test_unmapped.bam").getFile();
        File bamFile = new File(fileName);
        
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
       String fileName = this.getClass().getResource("/test_data_in/piculus_test_unmapped.bam").getFile();
        File bamFile = new File(fileName);
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
        String fileName = this.getClass().getResource("/test_data_in/piculus_test_unmapped.bam").getFile();
        File bamFile = new File(fileName);
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
        String fileName = this.getClass().getResource("/test_data_in/piculus_test_unmapped.bam").getFile();
        File bamFile = new File(fileName);
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
        String fileName = this.getClass().getResource("/test_data_in/piculus_test_unmapped.bam").getFile();
        File bamFile = new File(fileName);
        File fastqInLeft = null;
        File fastqInRight = null;
        MappedSamRecords instance = new MappedSamRecords();
        instance.printReadsFromBamAndFastq(bamFile, fastqInLeft, fastqInRight);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
}