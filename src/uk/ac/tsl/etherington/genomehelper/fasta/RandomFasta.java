/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fasta;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import mathsutils.RandomUtils;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava3.data.sequence.FastaSequence;
import org.biojava3.data.sequence.SequenceUtil;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

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

    /**
     * Takes a genome, concatenates it and outputs any number of shuffled versions
     * @param fasta the multi-fasta infile
     * @param numberOfSeqsRequired the number of scrambled genomes required
     * @param filePrefix the output sequences
     * @throws Exception
     */
    public void scrambleGenome(File fasta, int numberOfGenomesRequired, String filePrefix) throws Exception
    {
        //read in the fasta sequences to a List
        BufferedReader br = new BufferedReader(new FileReader(fasta));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");
        //get the reference genome
        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                                                                       alpha.getTokenization("token"), ns);
        String genome = "";
        //concatenate it
        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            String dna = rec.seqString();
            genome = genome.concat(dna);
        }
        
        //create an array list filled with ints from zero to genome legth
        ArrayList<Integer> nos = new ArrayList<>(genome.length());
        for (int i = 0; i < genome.length(); i++)
        { 
            nos.add(i);
        }
        
        //for each scrambled genome needed
        for (int i = 0; i < numberOfGenomesRequired; i++)
        {
            //Fisher–Yates shuffle the array list
            Collections.shuffle(nos);
            //and then go through it from the shuffled start and create a new genome string
            StringBuilder s = new StringBuilder();
            for (int x : nos)
            {
                //System.out.print(alphabet.charAt(x));
                s = s.append(genome.charAt(x));
            }
            //print the new genome to file
            String outfile = filePrefix + "_" + (i + 1) + ".fasta";
            try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(outfile)))
            {
                writer.write(">" + outfile);
                writer.write("\n");
                writer.write(s.toString());
                writer.write("\n");
                writer.close();
            }
        }


    }
}
