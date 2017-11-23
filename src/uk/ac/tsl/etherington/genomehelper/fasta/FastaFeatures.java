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
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.TreeMap;
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

    public Map<String, String[]> compareTwoFastaSequences(File seq1, File seq2) throws Exception
    {

        HashMap<String, String[]> snps = new HashMap<>();
        LinkedHashMap<String, DNASequence> dna1 = FastaReaderHelper.readFastaDNASequence(seq1);
        LinkedHashMap<String, DNASequence> dna2 = FastaReaderHelper.readFastaDNASequence(seq2);

        for (Entry<String, DNASequence> entry : dna1.entrySet())
        {
            DNASequence dna1seq = entry.getValue();
            //System.out.println(entry.getValue().getOriginalHeader() + "=" + entry.getValue().getSequenceAsString());
            if (dna2.containsKey(entry.getValue().getOriginalHeader()))
            {
                DNASequence dna2seq = dna2.get(dna1seq.getOriginalHeader());
                if (dna1seq.getLength() != dna2seq.getLength())
                {
                    System.err.println("Can't compare " + entry.getValue().getOriginalHeader() + " as they are different sizes: " + dna1seq.getLength() + " and  " + dna2seq.getLength());
                }
                else
                {
                    Iterator it1 = dna1seq.iterator();
                    Iterator it2 = dna2seq.iterator();
                    int position = 1;
                    while (it1.hasNext())
                    {
                        NucleotideCompound nt1 = (NucleotideCompound) it1.next();
                        NucleotideCompound nt2 = (NucleotideCompound) it2.next();
                        //System.out.println("Comparing " +nt1.toString()+ " to "+nt2.toString());
                        if (!nt1.equalsIgnoreCase(nt2))
                        {
                            String loci = dna1seq.getOriginalHeader().concat(":").concat(Integer.toString(position));
                            System.out.println("Found " + loci + " " + nt1 + " in File 1 and " + nt2 + " in File 2");
                            String[] snpArray = new String[2];
                            snpArray[0] = nt1.toString();
                            snpArray[1] = nt2.toString();
                            //Arrays.sort(snpArray);
                            snps.put(loci, snpArray);

                        }
                        position++;
                    }
                }
            }
        }
        return snps;
    }

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

        double genomeSizeOver1kb = 0;

        ArrayList<Integer> seqLengths = new ArrayList<>();
        Map<Character, Integer> numChars = new HashMap<>();

        for (Map.Entry<String, DNASequence> entry : b.entrySet())
        {
            DNASequence dna = entry.getValue();
            int seqLength = dna.getLength();
            if (seqLength >= 1000)
            {
                genomeSizeOver1kb += seqLength;
                seqLengths.add(seqLength);
            }
        }

        Collections.sort(seqLengths);
        Collections.reverse(seqLengths);

        double cumulativeSize = 0;
        double mb = 0;
        TreeMap<Integer, Double> nstats = new TreeMap<>();

        //Get the N-stats > 1kb
        cumulativeSize = 0;
        for (Integer i : seqLengths)
        {

            cumulativeSize += i;

            if ((genomeSizeOver1kb / 100) * 10 < cumulativeSize && (genomeSizeOver1kb / 100) * 10 > (cumulativeSize - i))
            {
                //System.out.print("10\t");
                mb = (double) i / 1000000.000;
                //System.out.println(decimalFormat.format(mb));
                nstats.put(10, mb);
            }
            if ((genomeSizeOver1kb / 100) * 20 < cumulativeSize && (genomeSizeOver1kb / 100) * 20 > (cumulativeSize - i))
            {
                //System.out.print("20\t");
                mb = (double) i / 1000000.000;
                //System.out.println(decimalFormat.format(mb));
                nstats.put(20, mb);
            }
            if ((genomeSizeOver1kb / 100) * 30 < cumulativeSize && (genomeSizeOver1kb / 100) * 30 > (cumulativeSize - i))
            {
                //System.out.print("30\t");
                mb = (double) i / 1000000.000;
                //System.out.println(decimalFormat.format(mb));
                nstats.put(30, mb);
            }
            if ((genomeSizeOver1kb / 100) * 40 < cumulativeSize && (genomeSizeOver1kb / 100) * 40 > (cumulativeSize - i))
            {
                //System.out.print("40\t");
                mb = (double) i / 1000000.000;
                //System.out.println(decimalFormat.format(mb));
                nstats.put(40, mb);
            }
            if ((genomeSizeOver1kb / 100) * 50 < cumulativeSize && (genomeSizeOver1kb / 100) * 50 > (cumulativeSize - i))
            {
                //System.out.print("50\t");
                mb = (double) i / 1000000.000;
                //System.out.println(decimalFormat.format(mb));
                nstats.put(50, mb);
            }
            if ((genomeSizeOver1kb / 100) * 60 < cumulativeSize && (genomeSizeOver1kb / 100) * 60 > (cumulativeSize - i))
            {
                //System.out.print("60\t");
                mb = (double) i / 1000000.000;
                //System.out.println(decimalFormat.format(mb));
                nstats.put(60, mb);
            }
            if ((genomeSizeOver1kb / 100) * 70 < cumulativeSize && (genomeSizeOver1kb / 100) * 70 > (cumulativeSize - i))
            {
                //System.out.print("70\t");
                mb = (double) i / 1000000.000;
                //System.out.println(decimalFormat.format(mb));
                nstats.put(70, mb);
            }
            if ((genomeSizeOver1kb / 100) * 80 < cumulativeSize && (genomeSizeOver1kb / 100) * 80 > (cumulativeSize - i))
            {
                //System.out.print("80\t");
                mb = (double) i / 1000000.000;
                //System.out.println(decimalFormat.format(mb));
                nstats.put(80, mb);
            }
            if ((genomeSizeOver1kb / 100) * 90 < cumulativeSize && (genomeSizeOver1kb / 100) * 90 > (cumulativeSize - i))
            {
                //System.out.print("90\t");
                mb = (double) i / 1000000.000;
                //System.out.println(decimalFormat.format(mb));
                nstats.put(90, mb);
            }

        }
        Set set = nstats.entrySet();
        Iterator iterator = set.iterator();
        while (iterator.hasNext())
        {
            Map.Entry mentry = (Map.Entry) iterator.next();
            System.out.print( "\t" + mentry.getKey());
        }
        System.out.println("\n");
        iterator = set.iterator();
        System.out.print(refSeq.getName());
        while (iterator.hasNext())
        {
            Map.Entry mentry = (Map.Entry) iterator.next();
            System.out.print("\t" + mentry.getValue());
        }
        System.out.println("\n");
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

    public void getSeqsOfMatchingLengths(ArrayList<File> files) throws Exception
    {
        //files must be fasta and have matching sequence names to be compared
        //seqLengths contains the names of the sequences and an ArrayList of seqlengths

        HashMap<String, ArrayList<Integer>> seqLengths = new HashMap<>();

        for (File file : files)
        {
            System.out.println("Processing file " + file.toString());
            //HashMap<String, String[]> snps = new HashMap<>();
            LinkedHashMap<String, DNASequence> dna = FastaReaderHelper.readFastaDNASequence(file);

            for (Entry<String, DNASequence> entry : dna.entrySet())
            {
                DNASequence dnaseq = entry.getValue();
                String name = entry.getKey();
                int length = dnaseq.getLength();
                if (seqLengths.containsKey(name))
                {
                    ArrayList al = seqLengths.get(name);
                    al.add(length);
                    seqLengths.put(name, al);
                }
                else
                {
                    ArrayList<Integer> al = new ArrayList<>();
                    al.add(length);
                    seqLengths.put(name, al);
                }

            }
        }
        System.out.println("Found " + seqLengths.size() + " different sequences");
        HashMap<String, Integer> identicalSeqs = new HashMap<>();
        for (Map.Entry<String, ArrayList<Integer>> seqInfo : seqLengths.entrySet())
        {

            //count how many entries are in the lengths array
            int seqCount = seqInfo.getValue().size();
            //if we have the same number of entries as files...
            if (seqCount == files.size())
            {
                //System.out.println(seqInfo.getKey()+" has three values");
                ArrayList<Integer> al = seqInfo.getValue();
                Set<Integer> uniqueLengths = new HashSet<>(al);
                //System.out.println(uniqueLengths.size()+" are unique");
                if (uniqueLengths.size() == 1)
                {
                    String name = seqInfo.getKey();
                    Integer firstElement = al.get(0);
                    identicalSeqs.put(name, firstElement);
                }

            }
        }
        System.out.println("Found " + identicalSeqs.size() + " sequences with identical lengths");
        Map<String, Integer> sortedIdenticalSeqs = sortByValue(identicalSeqs);
        for (Map.Entry<String, Integer> seqInfo : sortedIdenticalSeqs.entrySet())
        {
            System.out.println(seqInfo.getKey() + "\t" + seqInfo.getValue());
        }
    }

    public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map)
    {
        List<Map.Entry<K, V>> list = new LinkedList<>(map.entrySet());
        Collections.sort(list, new Comparator<Map.Entry<K, V>>()
        {
            public int compare(Map.Entry<K, V> o1, Map.Entry<K, V> o2)
            {
                return (o1.getValue()).compareTo(o2.getValue());
            }
        });

        Map<K, V> result = new LinkedHashMap<>();
        for (Map.Entry<K, V> entry : list)
        {
            result.put(entry.getKey(), entry.getValue());
        }
        return result;
    }

}
