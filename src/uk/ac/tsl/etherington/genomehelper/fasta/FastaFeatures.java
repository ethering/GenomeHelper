/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fasta;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.Symbol;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.DNASequenceCreator;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

/**
 * A class of helper methods to obtain information and data structures from
 * fasta files
 *
 * @author ethering
 */
public class FastaFeatures
{

    /**
     *
     * @param refSeq a multi-fasta file of DNA sequences
     * @return a HashMap of sequence names as keys with their sequence-lengths
     * as values
     * @throws FileNotFoundException exception
     * @throws BioException exception
     */
    public static HashMap<String, Integer> getSequenceLengths(File refSeq) throws FileNotFoundException, BioException, Exception
    {

        HashMap<String, Integer> seqLengths = new HashMap<>();
        LinkedHashMap<String, DNASequence> seqs = FastaReaderHelper.readFastaDNASequence(refSeq, true);

        for (Map.Entry<String, DNASequence> entry : seqs.entrySet())
        {
            String seqName = entry.getKey();
            DNASequence dna = entry.getValue();
            seqLengths.put(seqName, dna.getLength());
        }
        return seqLengths;
    }

    /**
     * Parses a multi-fasta file into a LinkedHashMap of sequence names and
     * DNASequences
     *
     * @param refSeq a multi-fasta file of DNA sequences
     * @return a LinkedHashMap of Accession numbers as keys with their DNA
     * sequences as values. All sequence names will have any non-accession
     * information stripped (e.g. "gi|2033910|gb|AA381581.1|AA381581 EST94688
     * Activated T-cells I Homo sapiens" would become
     * "gi|2033910|gb|AA381581.1|AA381581")
     * @throws Exception
     */
    public static LinkedHashMap<String, DNASequence> getParsedDNASequences(File refSeq) throws Exception
    {
        LinkedHashMap<String, DNASequence> tempgenome = FastaReaderHelper.readFastaDNASequence(refSeq);
        LinkedHashMap<String, DNASequence> genome = new LinkedHashMap<>();
        for (Map.Entry<String, DNASequence> entry : tempgenome.entrySet())
        {
            //get the sequence name
            String seqName = entry.getKey();
            //split the name on whitespace to get the sequence ID
            String newSeqName = seqName.split(" ")[0];
            DNASequence dna = entry.getValue();
            genome.put(newSeqName, dna);
        }
        tempgenome.clear();
        return genome;
    }

    /**
     *
     * @param refSeq the reference sequence
     * @return the cumulative length of the provided sequences
     * @throws FileNotFoundException
     * @throws BioException
     */
    public static double getGenomeSize(File refSeq) throws FileNotFoundException, BioException, Exception
    {
        double genomeSize = 0;
        HashMap<String, Integer> seqLengths = new HashMap<>(FastaFeatures.getSequenceLengths(refSeq));
        for (Map.Entry<String, Integer> entry : seqLengths.entrySet())
        {
            int genomeLength = entry.getValue();
            genomeSize += genomeLength;
        }
        return genomeSize;
    }

    /**
     *
     * @param refSeq the reference sequence
     * @return a HashMap where the keys are sequence names and the values are
     * int arrays the length of the corresponding sequence, filed with zeros
     * @throws FileNotFoundException
     * @throws BioException
     */
    public static HashMap getSequenceAsHashMapIntArray(File refSeq) throws FileNotFoundException, BioException, Exception
    {

        HashMap<String, Integer> seqLengths = new HashMap<>(FastaFeatures.getSequenceLengths(refSeq));
        HashMap<String, int[]> codingMap = new HashMap<>();

        for (Map.Entry<String, Integer> entry : seqLengths.entrySet())
        {
            String seqName = entry.getKey();
            int genomeLength = entry.getValue();
            //all elements of an int [] are give value of zero by default
            int[] intArray = new int[genomeLength];
            codingMap.put(seqName, intArray);
        }

        return codingMap;
    }

    /**
     *
     * @param refSeq the reference sequence
     * @return a sorted ArrayList of sequence lengths
     * @throws FileNotFoundException
     * @throws BioException
     */
    public static ArrayList<Integer> getSequenceAsSortedIntArrayList(File refSeq) throws FileNotFoundException, BioException, Exception
    {

        HashMap<String, Integer> seqLengths = new HashMap<>(FastaFeatures.getSequenceLengths(refSeq));
        ArrayList<Integer> lengths = new ArrayList<>();

        for (Map.Entry<String, Integer> entry : seqLengths.entrySet())
        {
            lengths.add(entry.getValue());
        }
        Collections.sort(lengths);
        Collections.reverse(lengths);
        return lengths;
    }

    public static ArrayList<Integer> getSequenceAsSortedIntArrayList(File refSeq, int minContigSize) throws FileNotFoundException, BioException, Exception
    {

        HashMap<String, Integer> seqLengths = new HashMap<>(FastaFeatures.getSequenceLengths(refSeq));
        ArrayList<Integer> lengths = new ArrayList<>();

        for (Map.Entry<String, Integer> entry : seqLengths.entrySet())
        {
            int length = entry.getValue();
            if (length >= minContigSize)
            {
                lengths.add(length);
            }

        }
        Collections.sort(lengths);
        Collections.reverse(lengths);
        return lengths;
    }

    public void getNStats(ArrayList<Integer> sortedSeqLengths)
    {
        double genomeSize = 0;

        for (Integer i : sortedSeqLengths)
        {
            genomeSize += i;
        }

        System.out.printf("Genome size = " + genomeSize);
        System.out.println("\nNo. seqs = " + sortedSeqLengths.size());
        System.out.println("Longest_contig\t" + sortedSeqLengths.get(0));
        System.out.println("N-size\tlength(Mb)");
        double cumulativeSize = 0;
        double mb = 0;
        DecimalFormat decimalFormat = new DecimalFormat("0.00000");

        for (Integer i : sortedSeqLengths)
        {
            cumulativeSize += i;
            mb = (double) i / 1000000.000;

            if ((genomeSize / 100) * 10 < cumulativeSize && (genomeSize / 100) * 10 > (cumulativeSize - i))
            {
                System.out.print("10\t");
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 20 < cumulativeSize && (genomeSize / 100) * 20 > (cumulativeSize - i))
            {
                System.out.print("20\t");
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 30 < cumulativeSize && (genomeSize / 100) * 30 > (cumulativeSize - i))
            {
                System.out.print("30\t");
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 40 < cumulativeSize && (genomeSize / 100) * 40 > (cumulativeSize - i))
            {
                System.out.print("40\t");
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 50 < cumulativeSize && (genomeSize / 100) * 50 > (cumulativeSize - i))
            {
                System.out.print("50\t");
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 60 < cumulativeSize && (genomeSize / 100) * 60 > (cumulativeSize - i))
            {
                System.out.print("60\t");
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 70 < cumulativeSize && (genomeSize / 100) * 70 > (cumulativeSize - i))
            {
                System.out.print("70\t");
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 80 < cumulativeSize && (genomeSize / 100) * 80 > (cumulativeSize - i))
            {
                System.out.print("80\t");
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 90 < cumulativeSize && (genomeSize / 100) * 90 > (cumulativeSize - i))
            {
                System.out.print("90\t");
                System.out.println(decimalFormat.format(mb));
            }
        }

    }

    /*
     Provides:
     N-stats for all seqs
     N-stats for > 1kb
     Size_includeN   
     Size_withoutN   
     Scaffold_Num    
     Genome size
     Scaffold-size distribution (>100, >500, >1kb, >10kb, > 100kb, >1Mb)
     GC content
     */
    public void getGenomeStats(File refSeq) throws Exception
    {
        FileInputStream inStream = new FileInputStream(refSeq);
        FastaReader<DNASequence, NucleotideCompound> fastaReader = new FastaReader<>(inStream,
                                                                                     new GenericFastaHeaderParser<>(),
                                                                                     new DNASequenceCreator(AmbiguityDNACompoundSet.getDNACompoundSet()));
        LinkedHashMap<String, DNASequence> b = fastaReader.process();

        DecimalFormat decimalFormat = new DecimalFormat("0.00000");
        double genomeSize = 0;
        double genomeSizeOver1kb = 0;
        int nseqsOver100b = 0;
        int nseqsOver500b = 0;
        int nseqsOver1kb = 0;
        int nseqsOver10kb = 0;
        int nseqsOver100kb = 0;
        int nseqsOver1mb = 0;

        ArrayList<Integer> seqLengths = new ArrayList<>();
        Map<Character, Integer> numChars = new HashMap<>();

        for (Map.Entry<String, DNASequence> entry : b.entrySet())
        {
            DNASequence dna = entry.getValue();
            int seqLength = dna.getLength();

            seqLengths.add(seqLength);
            genomeSize += seqLength;
            if (seqLength > 100)
            {
                nseqsOver100b++;
            }
            if (seqLength > 500)
            {
                nseqsOver500b++;
            }
            if (seqLength > 1000)
            {
                genomeSizeOver1kb += seqLength;
                nseqsOver1kb++;
            }
            if (seqLength > 10000)
            {

                nseqsOver10kb++;
            }
            if (seqLength > 100000)
            {

                nseqsOver100kb++;
            }
            if (seqLength > 1000000)
            {

                nseqsOver1mb++;
            }
            String seqString = dna.getSequenceAsString().toLowerCase();
            for (int i = 0; i < seqLength; ++i)
            {
                char charAt = seqString.charAt(i);

                if (!numChars.containsKey(charAt))
                {
                    numChars.put(charAt, 1);
                }
                else
                {
                    numChars.put(charAt, numChars.get(charAt) + 1);
                }
            }

        }
        System.out.println("Nucleotide content");
        double n = 0;
        double a = numChars.get('a');
        double t = numChars.get('t');
        double c = numChars.get('c');
        double g = numChars.get('g');
        //just in case we are lucky enough to be N-free!
        if (numChars.containsKey('n'))
        {
             n = numChars.get('n');
        }

        DecimalFormat df = new DecimalFormat("##.##");
        System.out.println("a = " + df.format(a));
        System.out.println("t = " + df.format(t));
        System.out.println("c = " + df.format(c));
        System.out.println("g = " + df.format(g));
        System.out.println("n = " + df.format(n));

        double gc = g + c;
        double all = a + t + c + g;
        float gcContent = (float) ((gc) / all) * 100;

        System.out.println("GC content (G+C)/(A+C+G+T) = " + df.format(gcContent));
        //System.out.println("GC content (G+C)/(A+C+G+T) = " + gcContent);

        Collections.sort(seqLengths);
        Collections.reverse(seqLengths);

        System.out.println("Genome size = " + df.format(genomeSize / 1000000) + " Mb");
        int nCount = 0;
        if (numChars.containsKey('n'))
        {
            nCount = numChars.get('n');
        }

        System.out.println("Genome size without N's = " + (df.format((genomeSize - nCount) / 1000000)) + " Mb");

        System.out.println("Genome size with contigs > 1kb = " + df.format(genomeSizeOver1kb / 1000000) + " Mb");
        System.out.println("\nNo. seqs = " + seqLengths.size());
        System.out.println("\nNo. seqs > 100 bases = " + nseqsOver100b);
        System.out.println("\nNo. seqs > 500 bases = " + nseqsOver500b);
        System.out.println("\nNo. seqs > 1kb = " + nseqsOver1kb);
        System.out.println("\nNo. seqs > 10kb = " + nseqsOver10kb);
        System.out.println("\nNo. seqs > 100kb = " + nseqsOver100kb);
        System.out.println("\nNo. seqs > 1Mb = " + nseqsOver1mb);

        double longestContig = (double) seqLengths.get(0) / 1000000.00;
        System.out.println("Longest_contig\t" + decimalFormat.format(longestContig) + " Mb");
        System.out.println("N-size\tlength(Mb)");
        double cumulativeSize = 0;
        double mb = 0;

        //Get the N-stats
        for (Integer i : seqLengths)
        {
            cumulativeSize += i;

            if ((genomeSize / 100) * 10 < cumulativeSize && (genomeSize / 100) * 10 > (cumulativeSize - i))
            {
                System.out.print("10\t");
                mb = (double) i / 1000000.000;
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 20 < cumulativeSize && (genomeSize / 100) * 20 > (cumulativeSize - i))
            {
                System.out.print("20\t");
                mb = (double) i / 1000000.000;
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 30 < cumulativeSize && (genomeSize / 100) * 30 > (cumulativeSize - i))
            {
                System.out.print("30\t");
                mb = (double) i / 1000000.000;
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 40 < cumulativeSize && (genomeSize / 100) * 40 > (cumulativeSize - i))
            {
                System.out.print("40\t");
                mb = (double) i / 1000000.000;
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 50 < cumulativeSize && (genomeSize / 100) * 50 > (cumulativeSize - i))
            {
                System.out.print("50\t");
                mb = (double) i / 1000000.000;
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 60 < cumulativeSize && (genomeSize / 100) * 60 > (cumulativeSize - i))
            {
                System.out.print("60\t");
                mb = (double) i / 1000000.000;
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 70 < cumulativeSize && (genomeSize / 100) * 70 > (cumulativeSize - i))
            {
                System.out.print("70\t");
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 80 < cumulativeSize && (genomeSize / 100) * 80 > (cumulativeSize - i))
            {
                System.out.print("80\t");
                mb = (double) i / 1000000.000;
                System.out.println(decimalFormat.format(mb));
            }
            if ((genomeSize / 100) * 90 < cumulativeSize && (genomeSize / 100) * 90 > (cumulativeSize - i))
            {
                System.out.print("90\t");
                mb = (double) i / 1000000.000;
                System.out.println(decimalFormat.format(mb));
            }
        }
        System.out.println("N-stats for contigs > 1kb");
        System.out.println("N-size\tlength(Mb)");
        //Get the N-stats > 1kb
        cumulativeSize = 0;
        for (Integer i : seqLengths)
        {
            if (i > 1000)
            {
                cumulativeSize += i;

                if ((genomeSizeOver1kb / 100) * 10 < cumulativeSize && (genomeSizeOver1kb / 100) * 10 > (cumulativeSize - i))
                {
                    System.out.print("10\t");
                    mb = (double) i / 1000000.000;
                    System.out.println(decimalFormat.format(mb));
                }
                if ((genomeSizeOver1kb / 100) * 20 < cumulativeSize && (genomeSizeOver1kb / 100) * 20 > (cumulativeSize - i))
                {
                    System.out.print("20\t");
                    mb = (double) i / 1000000.000;
                    System.out.println(decimalFormat.format(mb));
                }
                if ((genomeSizeOver1kb / 100) * 30 < cumulativeSize && (genomeSizeOver1kb / 100) * 30 > (cumulativeSize - i))
                {
                    System.out.print("30\t");
                    mb = (double) i / 1000000.000;
                    System.out.println(decimalFormat.format(mb));
                }
                if ((genomeSizeOver1kb / 100) * 40 < cumulativeSize && (genomeSizeOver1kb / 100) * 40 > (cumulativeSize - i))
                {
                    System.out.print("40\t");
                    mb = (double) i / 1000000.000;
                    System.out.println(decimalFormat.format(mb));
                }
                if ((genomeSizeOver1kb / 100) * 50 < cumulativeSize && (genomeSizeOver1kb / 100) * 50 > (cumulativeSize - i))
                {
                    System.out.print("50\t");
                    mb = (double) i / 1000000.000;
                    System.out.println(decimalFormat.format(mb));
                }
                if ((genomeSizeOver1kb / 100) * 60 < cumulativeSize && (genomeSizeOver1kb / 100) * 60 > (cumulativeSize - i))
                {
                    System.out.print("60\t");
                    mb = (double) i / 1000000.000;
                    System.out.println(decimalFormat.format(mb));
                }
                if ((genomeSizeOver1kb / 100) * 70 < cumulativeSize && (genomeSizeOver1kb / 100) * 70 > (cumulativeSize - i))
                {
                    System.out.print("70\t");
                    System.out.println(decimalFormat.format(mb));
                }
                if ((genomeSizeOver1kb / 100) * 80 < cumulativeSize && (genomeSizeOver1kb / 100) * 80 > (cumulativeSize - i))
                {
                    System.out.print("80\t");
                    mb = (double) i / 1000000.000;
                    System.out.println(decimalFormat.format(mb));
                }
                if ((genomeSizeOver1kb / 100) * 90 < cumulativeSize && (genomeSizeOver1kb / 100) * 90 > (cumulativeSize - i))
                {
                    System.out.print("90\t");
                    mb = (double) i / 1000000.000;
                    System.out.println(decimalFormat.format(mb));
                }
            }

        }

    }

    public double getGCContent(File fileName) throws FileNotFoundException, NoSuchElementException, BioException
    {
        // Set up sequence iterator
        double noSeqs = 0;
        double combinedGcTally = 0;
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        SimpleNamespace ns = new SimpleNamespace("biojava");

        // You can use any of the convenience methods found in the BioJava 1.6 API
        RichSequenceIterator stream = RichSequence.IOTools.readFastaDNA(br, ns);

        // Iterate over all sequences in the stream
        while (stream.hasNext())
        {
            Sequence seq = stream.nextSequence();
            int gc = 0;
            for (int pos = 1; pos <= seq.length(); ++pos)
            {
                Symbol sym = seq.symbolAt(pos);
                if (sym == DNATools.g() || sym == DNATools.c())
                {
                    ++gc;
                }
            }
            double currentGcCount = (gc * 100.0) / seq.length();
            //System.out.println(seq.getName() + ": "+ currentGcCount+ "%");
            combinedGcTally += currentGcCount;
            noSeqs++;
        }
        double gcCount = (combinedGcTally / noSeqs) / 100;
        System.out.println("Overall gcCount = " + gcCount);
        return gcCount;
    }
}
