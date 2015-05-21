/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fastq;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedHashMap;
import htsjdk.samtools.fastq.FastqRecord;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.data.sequence.FastaSequence;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author ethering
 */
public class FastqParserTest
{

    public FastqParserTest()
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
     * Test of readNamesToHashSet method, of class FastqParser.
     */
    @Test
    public void testReadNamesToHashSet() throws Exception
    {
        System.out.println("readNamesToHashSet");
        File readNames = new File("test/test_data_in/piculus_test_seqids.txt");
        FastqParser instance = new FastqParser();
        HashSet result = instance.readNamesToHashSet(readNames);
        //there are 7 sequence ids in the file
        assertEquals(result.size(), 7);
    }

    /**
     * Test of getPairedFastqSeqsFromHashSet method, of class FastqParser.
     */
    @Test
    public void testGetPairedFastqSeqsFromHashSet() throws Exception
    {
        System.out.println("getPairedFastqSeqsFromHashSet");
        File readNames = new File("test/test_data_in/piculus_test_seqids.txt");
        File leftFastqFileIn = new File("test/test_data_in/piculus_test_left.fastq");
        File rightFastqFileIn = new File("test/test_data_in/piculus_test_right.fastq");
        File leftReadsOut = new File("test/test_data_out/piculus_test_left_selected.fastq");
        File rightReadsOut = new File("test/test_data_out/piculus_test_right_selected.fastq");
        FastqParser instance = new FastqParser();
        HashSet list = instance.readNamesToHashSet(readNames);
        instance.getPairedFastqSeqsFromHashSet(list, leftFastqFileIn, rightFastqFileIn, leftReadsOut, rightReadsOut);
        //there are 7 readnames in the seqids file, 6 are good, one is bogus
        long leftLineCount = FastqQCTest.countLines(leftReadsOut);
        long rightLineCount = FastqQCTest.countLines(rightReadsOut);
        assertEquals(leftLineCount, 24);
        assertEquals(rightLineCount, 24);
    }

    /**
     * Test of getOneSideFastqSeqsFromList method, of class FastqParser.
     */
    @Test
    public void testGetOneSideFastqSeqsFromList() throws Exception
    {
        System.out.println("getOneSideFastqSeqsFromList");
        File listFile = new File("test/test_data_in/piculus_test_seqids.txt");
        File fastqFileIn = new File("test/test_data_in/piculus_test_left.fastq");
        File fastqFileOut = new File("test/test_data_out/piculus_test_selected.fastq");
        FastqParser instance = new FastqParser();
        instance.getOneSideFastqSeqsFromList(listFile, fastqFileIn, fastqFileOut);
        long linecount = FastqQCTest.countLines(fastqFileOut);
        assertEquals(linecount, 24);
    }

    /**
     * Test of fastqToFastaFile method, of class FastqParser.
     */
    @Test
    public void testFastqToFastaFile() throws Exception
    {
        System.out.println("fastqToFastaFile");
        File fastqIn = new File("test/test_data_in/piculus_test_left.fastq");
        File fastaOut = new File("test/test_data_out/piculus_test_left.fasta");
        FastqParser instance = new FastqParser();
        instance.fastqToFastaFile(fastqIn, fastaOut);
        long fastqCount = FastqQCTest.countLines(fastqIn);
        long fastaCount = FastqQCTest.countLines(fastaOut);
        assertEquals(fastqCount, fastaCount * 2);
        //parse the fasta file to confirm we have valid fasta and count the number of fasta sequences
        LinkedHashMap<String, DNASequence> seqs = FastaReaderHelper.readFastaDNASequence(fastaOut);
        assertEquals(fastaCount / 2, seqs.values().size());
    }

    /**
     * Test of fastqToFastaSeq method, of class FastqParser.
     */
    @Test
    public void testFastqToFastaSeq() throws Exception
    {
        System.out.println("fastqToFastaSeq");
        FastqRecord record = new FastqRecord("test_read1", "ATCG", "", "BBB");
        FastqParser instance = new FastqParser();
        FastaSequence expResult = new FastaSequence("test_read1", "ATCG");
        FastaSequence result = instance.fastqToFastaSeq(record);
        assertEquals(expResult, result);
    }

    /**
     * Test of fastqToFastaSixFrameTranslation method, of class FastqParser.
     */
    @Test
    public void testFastqToFastaSixFrameTranslation() throws Exception
    {
        System.out.println("fastqToFastaSixFrameTranslation");
        File fastqIn = new File("test/test_data_in/piculus_test_left.fastq");
        File fastaOut = new File("test/test_data_out/piculus_test_left_translated.fasta");
        boolean includeDNASeq = false;
        FastqParser instance = new FastqParser();
        instance.fastqToFastaSixFrameTranslation(fastqIn, fastaOut, includeDNASeq);
        long fastqCount = FastqQCTest.countLines(fastqIn);
        long fastaCount = FastqQCTest.countLines(fastaOut);
        //for every fastq seq (4 lines) there should be 6 fasta seqs (12 lines)
        assertEquals(fastqCount, fastaCount / 3);
        LinkedHashMap<String, ProteinSequence> seqs = FastaReaderHelper.readFastaProteinSequence(fastaOut);
        //there should be 6 protein sequences for every fastq
        assertEquals((fastqCount / 4) * 6, seqs.values().size());

        //include the DNA sequence so there should be 7 protein sequences for every fastq
        includeDNASeq = true;
        instance.fastqToFastaSixFrameTranslation(fastqIn, fastaOut, includeDNASeq);
        fastqCount = FastqQCTest.countLines(fastqIn);
        seqs = FastaReaderHelper.readFastaProteinSequence(fastaOut);
        assertEquals((fastqCount / 4) * 7, seqs.values().size());
    }
}