/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fasta;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.List;
import mathsutils.RandomUtils;
import org.biojava3.data.sequence.FastaSequence;
import org.biojava3.data.sequence.SequenceUtil;

/**
 *
 * @author ethering
 */
public class RandomFasta
{

    /**
     *
     * @param fasta the multi-fasta infile
     * @param numberOfSeqsRequired the number of random sequences required
     * @param outfile the output sequences
     * @throws Exception
     */
    public void selectRandomSequences(File fasta, int numberOfSeqsRequired, File outfile) throws Exception
    {
        //read in the fasta sequences to a List
        List<FastaSequence> seqs = SequenceUtil.readFasta(new FileInputStream(fasta));
        if (seqs.size() < numberOfSeqsRequired)
        {
            System.err.println("The number of random sequences requested is greater than that available.\nThere are " + seqs.size() + " sequences available");
            System.exit(1);
        }
        //get the random ints
        ArrayList<Integer> randomInts = RandomUtils.getRandomArray(numberOfSeqsRequired, seqs.size());
        //a holder
        ArrayList<FastaSequence> dnaSeqs = new ArrayList<>();

        for (Integer i : randomInts)
        {
            dnaSeqs.add(seqs.get(i));
        }
        FileOutputStream fop = new FileOutputStream(outfile);
        SequenceUtil.writeFasta(fop, dnaSeqs);
        fop.close();
    }
}
