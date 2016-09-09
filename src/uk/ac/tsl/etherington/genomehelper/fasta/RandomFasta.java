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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import mathsutils.RandomUtils;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.Symbol;
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
     * Takes a genome, concatenates it and outputs any number of shuffled
     * versions
     *
     * @param fasta the multi-fasta infile
     * @param numberOfGenomesRequired the number of scrambled genomes required
     * @param filePrefix the output sequences
     * @throws Exception
     */
    public void scrambleGenome2(File fasta, int numberOfGenomesRequired, String filePrefix) throws Exception
    {
        //read in the fasta sequences to a List
        BufferedReader br = new BufferedReader(new FileReader(fasta));

        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");
        //get the reference genome
        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br, alpha.getTokenization("token"), ns);

        //calculate the nucleotide composition
        double g = 0;
        double c = 0;
        double t = 0;
        double a = 0;

        while (iterator.hasNext())
        {
            Sequence seq = iterator.nextSequence();

            for (int pos = 1; pos <= seq.length(); pos++)
            {
                Symbol sym = seq.symbolAt(pos);
                if (sym == DNATools.g())
                {
                    g++;
                }

                if (sym == DNATools.c())
                {
                    c++;
                }
                if (sym == DNATools.t())
                {
                    t++;
                }
                if (sym == DNATools.a())
                {
                    a++;
                }
            }
        }

        Double genomeSize = a + t + c + g;
        System.out.println("Genome size = " + genomeSize);
        Double gContent = (g * 100) / genomeSize;
        Double cContent = (c * 100) / genomeSize;
        Double aContent = (a * 100) / genomeSize;
        Double tContent = (t * 100) / genomeSize;

        System.out.println("gContent = " + gContent);
        System.out.println("cContent = " + cContent);
        System.out.println("aContent = " + aContent);
        System.out.println("tContent = " + tContent);

        int gContentI = (int) Math.rint(gContent);
        int cContentI = (int) Math.rint(cContent);
        int aContentI = (int) Math.rint(aContent);
        int tContentI = (int) Math.rint(tContent);
        int genomeContent = gContentI + cContentI + aContentI + tContentI;
        System.out.println("gContentI = " + gContentI);
        System.out.println("cContentI = " + cContentI);
        System.out.println("aContentI = " + aContentI);
        System.out.println("tContentI = " + tContentI);

        System.out.println("Genome array content = " + genomeContent);
        //create a DNA sequence of around 100 nucleotides (totalContent size)

        String[] dnaContent = new String[genomeContent];
        int index = 0;
        for (int i = 0; i < gContentI; i++)
        {
            dnaContent[index] = "g";
            index++;
        }
        for (int i = 0; i < cContentI; i++)
        {
            dnaContent[index] = "c";
            index++;
        }
        for (int i = 0; i < aContentI; i++)
        {
            dnaContent[index] = "a";
            index++;
        }
        for (int i = 0; i < tContentI; i++)
        {
            dnaContent[index] = "t";
            index++;
        }
        //for each scrambled genome needed

        for (int i = 0; i < numberOfGenomesRequired; i++)
        {
            int ntsWritten = 0;
            //create the new genome file
            String prefix = filePrefix + "_" + (i + 1);
            String outfile = prefix + ".fasta";
            try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(outfile)))
            {
                writer.write(">" + prefix);
                writer.write("\n");
                //for the length of the genome 
                for (int x = 0; x < genomeSize; x++)
                {
                    //select a random nucleotide from the dnaContent array and write it to file
                    int idx = new Random().nextInt(dnaContent.length);
                    String random = (dnaContent[idx]);
                    writer.write(random);
                    ntsWritten++;
                }
                writer.write("\n");
                writer.flush();
                writer.close();

                System.out.println("Nts written = " + ntsWritten);

            }
        }

    }

    public void scrambleGenome(File fasta, int numberOfGenomesRequired, String filePrefix) throws Exception
    {
        //read in the fasta sequences to a List
        BufferedReader br = new BufferedReader(new FileReader(fasta));

        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");
        //get the reference genome
        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br, alpha.getTokenization("token"), ns);
        Map<String, Integer> seqLengths = new LinkedHashMap<>();

        //calculate the nucleotide composition
        double g = 0;
        double c = 0;
        double t = 0;
        double a = 0;

        while (iterator.hasNext())
        {
            //Sequence seq = iterator.nextSequence();
            RichSequence seq = iterator.nextRichSequence();
            String seqName = seq.getAccession() + " "+seq.getDescription();
            
            seqLengths.put(seqName, seq.length());

            for (int pos = 1; pos <= seq.length(); pos++)
            {
                Symbol sym = seq.symbolAt(pos);
                if (sym == DNATools.g())
                {
                    g++;
                }

                if (sym == DNATools.c())
                {
                    c++;
                }
                if (sym == DNATools.t())
                {
                    t++;
                }
                if (sym == DNATools.a())
                {
                    a++;
                }
            }
        }

        double genomeSize = a + t + c + g;
        System.out.println("Genome size = " + genomeSize);
        double gContent = (g * 100) / genomeSize;
        double cContent = (c * 100) / genomeSize;
        double aContent = (a * 100) / genomeSize;
        double tContent = (t * 100) / genomeSize;

        System.out.println("gContent = " + gContent);
        System.out.println("cContent = " + cContent);
        System.out.println("aContent = " + aContent);
        System.out.println("tContent = " + tContent);

        int gContentI = (int) Math.rint(gContent);
        int cContentI = (int) Math.rint(cContent);
        int aContentI = (int) Math.rint(aContent);
        int tContentI = (int) Math.rint(tContent);
        int genomeContent = gContentI + cContentI + aContentI + tContentI;
        System.out.println("gContentI = " + gContentI);
        System.out.println("cContentI = " + cContentI);
        System.out.println("aContentI = " + aContentI);
        System.out.println("tContentI = " + tContentI);

        System.out.println("Genome array content = " + genomeContent);
        //create a DNA sequence of around 100 nucleotides (totalContent size)

        String[] dnaContent = new String[genomeContent];
        int index = 0;
        for (int i = 0; i < gContentI; i++)
        {
            dnaContent[index] = "G";
            index++;
        }
        for (int i = 0; i < cContentI; i++)
        {
            dnaContent[index] = "C";
            index++;
        }
        for (int i = 0; i < aContentI; i++)
        {
            dnaContent[index] = "A";
            index++;
        }
        for (int i = 0; i < tContentI; i++)
        {
            dnaContent[index] = "T";
            index++;
        }
        //for each scrambled genome needed

        for (int i = 0; i < numberOfGenomesRequired; i++)
        {
            int ntsWritten = 0;
            //create the new genome file
            String outfile = filePrefix + "_" + (i + 1) + ".fasta";

            try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(outfile)))
            {
                for (Map.Entry<String, Integer> entry : seqLengths.entrySet())
                {
                    String seqName = entry.getKey();
                    Integer seqLength = entry.getValue();
                    writer.write(">" + seqName);
                    writer.write("\n");
                    //for the length of the genome 
                    for (int x = 0; x < seqLength; x++)
                    {
                        //select a random nucleotide from the dnaContent array and write it to file
                        int idx = new Random().nextInt(dnaContent.length);
                        String random = (dnaContent[idx]);
                        writer.write(random);
                        ntsWritten++;
                    }
                    writer.write("\n");
                }
                
                writer.flush();
                writer.close();

                System.out.println("Nts written = " + ntsWritten);

            }

        }
    }
}
