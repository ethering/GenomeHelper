/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fasta;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import mathsutils.RandomUtils;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.DNASequenceCreator;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
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

    public void selectRandomSequences(File fasta, int numberOfSeqsRequired, File outfile) throws Exception
    {
        List<FastaSequence> seqs = SequenceUtil.readFasta(new FileInputStream(fasta));

        ArrayList<Integer> randomInts = RandomUtils.getRandomArray(numberOfSeqsRequired, seqs.size());
        ArrayList<FastaSequence> dnaSeqs = new ArrayList<FastaSequence>();

        for (Integer i : randomInts)
        {
            dnaSeqs.add(seqs.get(i));
        }
        FileOutputStream fop = new FileOutputStream(outfile);
        SequenceUtil.writeFasta(fop, dnaSeqs);
    }
    //if you already know how many seqs you have in the fasta file

    public void selectRandomSequences(File fasta, int numberOfSeqsRequired, File outfile, int noSeqs) throws Exception
    {
        // get a random array of ints
        ArrayList<Integer> randomInts = RandomUtils.getRandomArray(numberOfSeqsRequired, noSeqs);
        //open the fasta file and iterate through it
        FileInputStream inStream = new FileInputStream(fasta);
        FastaReader<DNASequence, NucleotideCompound> fastaReader =
                new FastaReader<DNASequence, NucleotideCompound>(
                inStream,
                new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(),
                new DNASequenceCreator(AmbiguityDNACompoundSet.getDNACompoundSet()));
        LinkedHashMap<String, DNASequence> b = fastaReader.process();
        Integer itCounter = 1;
        BufferedWriter output = new BufferedWriter(new FileWriter(outfile));
        for (Map.Entry<String, DNASequence> entry : b.entrySet())
        {

            if (randomInts.contains(itCounter))
            {
                String id = entry.getKey();
                String seq = entry.getValue().getSequenceAsString();

                output.write(">" + id);
                output.newLine();
                output.write(">" + seq);
                output.newLine();
                randomInts.remove(itCounter);
            }
            itCounter++;
            //System.out.println(itCounter);
        }
        output.close();
    }

    public void selectRandomSequences2(File fasta, int numberOfSeqsRequired, File outfile) throws Exception
    {
        int lines = countLines(fasta.getAbsolutePath());
        int noSeqs = lines / 2;
        System.out.println("There are " + noSeqs + " sequences");
        ArrayList<Integer> randomInts = RandomUtils.getRandomArray(numberOfSeqsRequired, noSeqs);
        BufferedReader br = new BufferedReader(new FileReader(fasta));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                alpha.getTokenization("token"), ns);
        Integer itCounter = 1;
        BufferedWriter output = new BufferedWriter(new FileWriter(outfile));

        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            if (randomInts.contains(itCounter))
            {
                String id = rec.getAccession();
                String seq = rec.seqString();

                output.write(">" + id);
                output.newLine();
                output.write(">" + seq);
                output.newLine();
                randomInts.remove(itCounter);
            }
            itCounter++;
            //System.out.println(itCounter);
        }
        output.close();
    }

    //if you already know how many seqs you have in the fasta file
    public void selectRandomSequences2(File fasta, int numberOfSeqsRequired, File outfile, int noSeqs) throws Exception
    {
        // get a random array of ints
        ArrayList<Integer> randomInts = RandomUtils.getRandomArray(numberOfSeqsRequired, noSeqs);
        //open the fasta file and iterate through it
        BufferedReader br = new BufferedReader(new FileReader(fasta));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                alpha.getTokenization("token"), ns);
        Integer itCounter = 1;
        BufferedWriter output = new BufferedWriter(new FileWriter(outfile));

        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            if (randomInts.contains(itCounter))
            {
                String id = rec.getAccession();
                String seq = rec.seqString();

                output.write(">" + id);
                output.newLine();
                output.write(">" + seq);
                output.newLine();
                randomInts.remove(itCounter);
            }
            itCounter++;
            //System.out.println(itCounter);
        }
        output.close();
    }
    
    public int countLines(String filename) throws IOException
    {
        FileReader fr = new FileReader(filename);
        LineNumberReader lnr = new LineNumberReader(fr);

        int linenumber = 0;

        while (lnr.readLine() != null)
        {
            linenumber++;
        }

        System.out.println("Total number of lines : " + linenumber);

        lnr.close();
        return linenumber;
    }

    public int countLines2(String filename) throws IOException
    {
        InputStream is = new BufferedInputStream(new FileInputStream(filename));
        try
        {
            byte[] c = new byte[1024];
            int count = 0;
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

    
//    public static void main(String args[]) throws Exception
//    {
//        RandomFasta rf = new RandomFasta();
//        File in = new File("/Users/ethering/projects/florian/hamming/baitlibrary.fasta");
//        File out = new File("/Users/ethering/projects/florian/hamming/sample_baitlibrary.fasta");
//        rf.selectRandomSequences2(in, 30, out, 52512);
//    }
}
