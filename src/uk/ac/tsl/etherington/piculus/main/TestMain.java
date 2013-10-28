/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.main;

import java.io.File;
import java.io.IOException;
import uk.ac.tsl.etherington.piculus.fasta.FastaMotifFinder;
import uk.ac.tsl.etherington.piculus.fastq.FastqParser;

/**
 *
 * @author ethering
 */
public class TestMain
{

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, Exception
    {
        File fastqIn = new File("test/test_data_in/piculus_test_left.fastq");
        File fastaOut = new File("test/test_data_in/piculus_test_left.fasta");
        FastqParser fq = new FastqParser();
        fq.fastqToFastaFile(fastqIn, fastaOut);
    }
}
