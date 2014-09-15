/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fasta;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.Symbol;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

/**
 *A class of helper methods to obtain information and data structures from fasta files 
 * @author ethering
 */
public class FastaFeatures
{
    /**
     * 
     * @param refSeq a multi-fasta file of DNA sequences
     * @return a HashMap of sequence names as keys with their sequence-lengths as values
     * @throws FileNotFoundException exception
     * @throws BioException exception
     */
    public static HashMap getSequenceLengths(File refSeq) throws FileNotFoundException, BioException, Exception
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
     * Parses a multi-fasta file into a LinkedHashMap of sequence names and DNASequences
     * @param refSeq a multi-fasta file of DNA sequences
     * @return a LinkedHashMap of Accession numbers as keys with their DNA sequences as values. All sequence names will have any non-accession information stripped
     * (e.g. "gi|2033910|gb|AA381581.1|AA381581 EST94688 Activated T-cells I Homo sapiens" would become "gi|2033910|gb|AA381581.1|AA381581")
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
     * @param refSeq
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
            genomeSize+=genomeLength;
        }
        return genomeSize;
    }

    /**
     * 
     * @param refSeq
     * @return a HashMap where the keys are sequence names and the values are int arrays the length of the corresponding sequence, filed with zeros
     * @throws FileNotFoundException
     * @throws BioException 
     */
    public static HashMap getSequenceAsIntArray(File refSeq) throws FileNotFoundException, BioException, Exception
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
            combinedGcTally+=currentGcCount;
            noSeqs++;
        }
        double gcCount = (combinedGcTally/noSeqs)/100;
        System.out.println("Overall gcCount = "+gcCount);
        return gcCount;
    }
}
