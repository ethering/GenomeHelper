/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fastq;

import uk.ac.tsl.etherington.genomehelper.fastq.FastqQC;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.util.FastqQualityFormat;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author ethering
 */
public class FastqQCTest
{

    public FastqQCTest()
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

    public static long countLines(java.io.File filename) throws IOException
    {
        InputStream is = new BufferedInputStream(new FileInputStream(filename));
        try
        {
            byte[] c = new byte[1024];
            long count = 0;
            int readChars = 0;
            boolean empty = true;
            while ((readChars = is.read(c)) != -1)
            {
                empty = false;
                for (int i = 0; i < readChars; ++i)
                {
                    if (c[i] == '\n')
                    {
                        ++count;
                    }
                }
            }
            return (count == 0 && !empty) ? 1 : count;
        }
        finally
        {
            is.close();
        }
    }

    /**
     * Test of checkFormat method, of class FastqQC.
     */
    @Test
    public void testCheckFormat()
    {
         System.out.println("Working Directory = " +
              System.getProperty("user.dir"));
        System.out.println("checkFormat");
        System.out.println("checking 'sanger'");
        String format = "sanger";
        FastqQC instance = new FastqQC();
        boolean expResult = true;
        boolean result = instance.checkFormat(format);
        System.out.println(result);
        assertEquals(expResult, result);
        System.out.println("checking 'illumina'");
        format = "illumina";
        expResult = true;
        result = instance.checkFormat(format);
        assertEquals(expResult, result);
        System.out.println("checking 'cornflakes'");

        format = "cornflakes";
        System.out.println(instance.checkFormat(format));
        expResult = false;
        result = instance.checkFormat(format);
        System.out.println(result);
        assertEquals(expResult, result);
    }

    /**
     * Test of checkLengthAndContent method, of class FastqQC.
     */
    @Test
    public void testCheckLengthAndContent()
    {
        System.out.println("qcFastq");
        //a good fastq sequence
        FastqRecord fq1 = new FastqRecord("test_read1", "ATCG", "", "BBB");
        //a bad fastq sequence - contains an N
        FastqRecord fq2 = new FastqRecord("test_read1", "ATCN", "", "BBBB");
        //a bad fastq sequence - too short
        FastqRecord fq3 = new FastqRecord("test_read1", "ATC", "", "BBB");
        //a bad fastq sequence - too long
        FastqRecord fq4 = new FastqRecord("test_read1", "ATCGT", "", "BBBBB");
        System.out.println(fq1.getReadString());
        int singleEndReadLength = 4;
        FastqQC instance = new FastqQC();
        boolean expResult = true;
        boolean result = instance.checkLengthAndContent(fq1, singleEndReadLength);
        assertEquals(expResult, result);
        expResult = false;
        result = instance.checkLengthAndContent(fq2, singleEndReadLength);
        assertEquals(expResult, result);
        result = instance.checkLengthAndContent(fq3, singleEndReadLength);
        assertEquals(expResult, result);
        result = instance.checkLengthAndContent(fq4, singleEndReadLength);
        assertEquals(expResult, result);
    }

    /**
     * Test of qcPairedReads method, of class FastqQC.
     */
    @Test
    public void testQcPairedReads() throws IOException
    {
        System.out.println("qcPairedReads");
        File leftFastqFileIn = new File("test/test_data_in/piculus_test_left.fastq");
        File rightFastqFileIn = new File("test/test_data_in/piculus_test_right.fastq");
        File leftReadsOut = new File("test/test_data_out/piculus_test_left_QC.fastq");
        File rightReadsOut = new File("test/test_data_out/piculus_test_right_QC.fastq");
        File singles = new File("test/test_data_out/piculus_test_singles_QC.fastq");
        int singleEndReadLength = 80;
        String format = "sanger";
        boolean writeBadSeqs = true;
        FastqQC instance = new FastqQC();
        instance.qcPairedReads(leftFastqFileIn, rightFastqFileIn, leftReadsOut, rightReadsOut, singles, singleEndReadLength, format, writeBadSeqs);
        //there are 11 bad reads in the file, so we should have 1989 good reads, or 7956 lines
        long lineCount = countLines(leftReadsOut);
        assertEquals(lineCount, 7956);
    }

    /**
     * Test of qcInterlacedReads method, of class FastqQC.
     */
    @Test
    public void testQcInterlacedReads() throws IOException
    {
        System.out.println("qcInterlacedReads");
        File fastqFileIn = new File("test/test_data_in/piculus_test_interlaced.fastq");
        File fastqFileOut = new File("test/test_data_out/piculus_test_interlaced_QC.fastq");
        File singles = new File("test/test_data_out/piculus_test_singles_QC.fastq");
        int singleEndReadLength = 80;
        String format = "sanger";
        boolean writeBadSeqs = true;
        FastqQC instance = new FastqQC();
        instance.qcInterlacedReads(fastqFileIn, fastqFileOut, singles, singleEndReadLength, format, writeBadSeqs);
        //there are 22 bad reads in the interlaced file, so we should have 3978 good reads, or 15912 lines
        long lineCount = countLines(fastqFileOut);
        assertEquals(lineCount, 15912);
    }

    /**
     * Test of qcInterlacedReadsToPairs method, of class FastqQC.
     */
    @Test
    public void testQcInterlacedReadsToPairs() throws IOException
    {
        System.out.println("qcInterlacedReadsToPairs");
        File fastqFileIn = new File("test/test_data_in/piculus_test_interlaced.fastq");
        File leftReadsOut = new File("test/test_data_out/piculus_test_left_QC.fastq");
        File rightReadsOut = new File("test/test_data_out/piculus_test_right_QC.fastq");
        File singles = new File("test/test_data_out/piculus_test_singles_QC.fastq");
        int singleEndReadLength = 80;
        String format = "sanger";
        boolean writeBadSeqs = true;
        FastqQC instance = new FastqQC();
        instance.qcInterlacedReadsToPairs(fastqFileIn, leftReadsOut, rightReadsOut, singles, singleEndReadLength, format, writeBadSeqs);
        //there are 11 bad reads in the file, so we should have 1989 good reads, or 7956 lines
        long leftLineCount = countLines(leftReadsOut);
        assertEquals(leftLineCount, 7956);
        //there are 11 bad reads in the file, so we should have 1989 good reads, or 7956 lines
        long rightLineCount = countLines(rightReadsOut);
        assertEquals(rightLineCount, 7956);
    }

    /**
     * Test of qcJoinedReads method, of class FastqQC.
     */
    @Test
    public void testQcJoinedReads() throws IOException
    {
        System.out.println("qcJoinedReads");
        File fastqFileIn = new File("test/test_data_in/piculus_test_joined.fastq");
        File leftReadsOut = new File("test/test_data_out/piculus_test_joined_left_QC.fastq");
        File rightReadsOut = new File("test/test_data_out/piculus_test_joined_right_QC.fastq");
        int readLength = 80;
        String format = "sanger";
        boolean writeBadSeqs = true;
        FastqQC instance = new FastqQC();
        instance.qcJoinedReads(fastqFileIn, leftReadsOut, rightReadsOut, readLength, format, writeBadSeqs);
        long leftLineCount = countLines(leftReadsOut);
        assertEquals(leftLineCount, 7956);
        //there are 11 bad reads in the file, so we should have 1989 good reads, or 7956 lines
        long rightLineCount = countLines(rightReadsOut);
        assertEquals(rightLineCount, 7956);
    }

    /**
     * Test of qcSingleEndReads method, of class FastqQC.
     */
    @Test
    public void testQcSingleEndReads() throws IOException
    {
        System.out.println("qcSingleEndReads");
        File fastqFileIn = new File("test/test_data_in/piculus_test_left.fastq");
        File readsOut = new File("test/test_data_out/piculus_test_left_QC.fastq");
        int readLength = 80;
        String format = "sanger";
        boolean writeBadSeqs = true;
        FastqQC instance = new FastqQC();
        instance.qcSingleEndReads(fastqFileIn, readsOut, readLength, format, writeBadSeqs);
        long lineCount = countLines(readsOut);
        assertEquals(lineCount, 7956);
    }

    /**
     * Test of groomRead method, of class FastqQC.
     */
    @Test
    public void testGroomRead()
    {
        System.out.println("groomRead");
        // an illumina quality encoded read
        FastqRecord record = new FastqRecord("@HWI-EAS396_0001:5:1:10468:1298#0/1", "CTTTTAGCAAGATATCTTATCCATTCCATCTTCGATCCACACAATTGAATCATGTAATTCTCCAATGTAACGCAAT",
                "+HWI-EAS396_0001:5:1:10468:1298#0/1", "ddc_cfcccfa[ddab\\_a`cfffdffS_ffc^fYddcWe]`]X^bcbadcffccW^ae[ffffffcdffdfaWcc");

        FastqQC instance = new FastqQC();
        FastqQualityFormat qualFormat = instance.guessFormat(record);
        assertEquals(qualFormat.toString(), "Illumina");

        String format = "illumina";
        FastqRecord groomedRead = instance.groomRead(record, format);
        qualFormat = instance.guessFormat(groomedRead);
        assertEquals(qualFormat.toString(), "Standard");
    }

//    /**
//     * Test of countLengthsAndNs method, of class FastqQC.
//     */
    @Test
    public void testCountLengthsAndNs()
    {
        System.out.println("countLengthsAndNs");
        File fastq = new File("test/test_data_in/piculus_test_left.fastq");
        FastqQC instance = new FastqQC();
        //there should be 11 sequences that are short or contain Ns
        int count = instance.countLengthsAndNs(fastq);
        assertEquals(count, 11);
    }

    /**
     * Test of veryfiyReads method, of class FastqQC.
     */
    @Test
    public void testVeryfiyReads()
    {
        System.out.println("veryfiyReads");
        File fastq = new File("test/test_data_in/piculus_test_left.fastq");
        FastqQC instance = new FastqQC();
        try
        {
            instance.veryfiyReads(fastq);
        }
        catch (Exception e)
        {
            fail("Failed due to exeception: " + e.getLocalizedMessage());
        }
    }

    /**
     * Test of veryfiyPairedReads method, of class FastqQC.
     */
    @Test
    public void testVeryfiyPairedReads()
    {
        System.out.println("veryfiyPairedReads");
        File fastqLeft = new File("test/test_data_in/piculus_test_left.fastq");
        File fastqRight = new File("test/test_data_in/piculus_test_right.fastq");
        FastqQC instance = new FastqQC();
        try
        {
            instance.veryfiyPairedReads(fastqLeft, fastqRight);
        }
        catch (Exception e)
        {
            fail("Failed due to exeception: " + e.getLocalizedMessage());
        }
    }
}