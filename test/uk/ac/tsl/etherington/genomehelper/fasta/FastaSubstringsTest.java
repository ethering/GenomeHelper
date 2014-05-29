/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fasta;

import uk.ac.tsl.etherington.genomehelper.fasta.FastaSubstrings;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import uk.ac.tsl.etherington.genomehelper.fastq.FastqQCTest;

/**
 *
 * @author ethering
 */
public class FastaSubstringsTest
{

    String seqid = "N57638:1:63BW7AAXX:1:1:9536:1046 1:N:0:";
    String seq = "TTTTGAAGAAGGATAGGGATGGAGATATGGGAGCTGAGTAAGCCATGGTTGATGGAAGGAAGGCACGCGAGGAGAGCGGA";
    //subseq is the same as seq, but with the first and last nucleotides missing
    String subseq = "TTTGAAGAAGGATAGGGATGGAGATATGGGAGCTGAGTAAGCCATGGTTGATGGAAGGAAGGCACGCGAGGAGAGCGG";

    public FastaSubstringsTest()
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
     * Test of getSequence method, of class FastaSubstrings.
     */
    @Test
    public void testSeqFromCommandLine_3args() throws Exception
    {
        System.out.println("seqFromCommandLine");
        File fastaIn = new File("test/test_data_in/piculus_test_left.fasta");
        File outfile = new File("test/test_data_out/piculus_test_left.fasta");
        FastaSubstrings instance = new FastaSubstrings();
        //get the required sequence
        instance.getSequence(fastaIn, outfile, seqid);
        LinkedHashMap<String, DNASequence> seqs = FastaReaderHelper.readFastaDNASequence(outfile);
        for (Map.Entry<String, DNASequence> entry : seqs.entrySet())
        {
            String seqName = entry.getKey();
            //test the it has the correct seqid
            assertEquals(seqName, seqid);
            //test it has the correct sequence
            DNASequence dna = entry.getValue();
            assertEquals(dna.getSequenceAsString(), seq);
        }
    }

    /**
     * Test of getSequence method, of class FastaSubstrings.
     */
    @Test
    public void testSeqFromCommandLine_5args() throws Exception
    {
        System.out.println("seqFromCommandLine");
        File fastaIn = new File("test/test_data_in/piculus_test_left.fasta");
        File outfile = new File("test/test_data_out/piculus_test_left.fasta");
        int start = 2;
        int end = 79;
        FastaSubstrings instance = new FastaSubstrings();
        instance.getSubSequence(fastaIn, outfile, seqid, start, end);
        LinkedHashMap<String, DNASequence> seqs = FastaReaderHelper.readFastaDNASequence(outfile);
        for (Map.Entry<String, DNASequence> entry : seqs.entrySet())
        {
            String seqName = entry.getKey();
            //test the it has the correct seqid
            assertEquals(seqName, seqid);
            //test it has the correct sequence
            DNASequence dna = entry.getValue();
            assertEquals(dna.getSequenceAsString(), subseq);
        }
    }

    /**
     * Test of findLongestCommonSequences method, of class FastaSubstrings.
     */
    @Test
    public void testFindLongestCommonSequences() throws Exception
    {
        System.out.println("findLongestCommonSequences");
        File fastaIn = new File("test/test_data_in/piculus_test_common_subseqs.fasta");
        File fastaIn2 = new File("test/test_data_in/piculus_test_common_subseqs2.fasta");
        //there should be 6 substrings: aaa, ttt, tt, aa, t, a
        File outfile = new File("test/test_data_out/piculus_test_subseqs.fasta");
        ArrayList<String> al = new ArrayList<>();
        al.add(fastaIn.toString());
        al.add(fastaIn2.toString());
        FastaSubstrings instance = new FastaSubstrings();
        instance.findLongestCommonSequences(al, outfile);
        //check that there only 6 lines in the outfile
        long noLines = FastqQCTest.countLines(outfile);
        assertEquals(noLines, 6);

    }

    /**
     * Test of getAllSubStrings method, of class FastaSubstrings.
     */
    @Test
    public void testGetAllSubStrings() throws Exception
    {
        System.out.println("getAllSubStrings");
        String fastaIn = new File("test/test_data_in/piculus_test_all_subseqs.fasta").toString();
        //'atcg' has 20 possible substrings (10 on forward, 10 on reverse)
        FastaSubstrings instance = new FastaSubstrings();
        ArrayList result = instance.getAllSubStrings(fastaIn);
        assertEquals(result.size(), 20);
    }
}