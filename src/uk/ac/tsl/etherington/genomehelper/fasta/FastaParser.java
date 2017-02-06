/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fasta;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.DNASequenceCreator;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

/**
 *
 * @author ethering
 */
public class FastaParser
{

    /**
     *
     * @param fastaFile
     * @param fastqFile
     * @throws FileNotFoundException
     * @throws BioException
     * @throws Exception
     */
    public void fastaToFastq(File fastaFile, File fastqFile) throws FileNotFoundException, BioException, Exception
    {
        FileWriter fw = new FileWriter(fastqFile.getAbsoluteFile());
        BufferedWriter bw = new BufferedWriter(fw);
        String newLine = System.getProperty("line.separator");

        BufferedReader br = new BufferedReader(new FileReader(fastaFile));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                                                                       alpha.getTokenization("token"), ns);
        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            String fastqSeq = "@";
            String name = rec.getName();
            String dna = rec.seqString();
            fastqSeq = fastqSeq.concat(name + newLine + dna + newLine + "+" + newLine);
            int len = rec.length();
            char[] chars = new char[len];
            Arrays.fill(chars, '#');
            String quality = new String(chars);
            fastqSeq = fastqSeq.concat(quality + newLine);
            bw.write(fastqSeq);

        }

        bw.close();
    }

    /**
     * Concatenates all fasta sequences in multifasta file to a single fasta
     * sequence
     *
     * @param fastaFileIn input fasta file containing multiple fasta sequences
     * @param fastaFileOut output fasta containing single concatenated sequence
     * @throws FileNotFoundException
     * @throws BioException
     * @throws Exception
     */
    public void multiFastaToSingleFasta(File fastaFileIn, File fastaFileOut) throws FileNotFoundException, BioException, Exception
    {
        FileWriter fw = new FileWriter(fastaFileOut.getAbsoluteFile());
        BufferedWriter bw = new BufferedWriter(fw);
        String newLine = System.getProperty("line.separator");

        BufferedReader br = new BufferedReader(new FileReader(fastaFileIn));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                                                                       alpha.getTokenization("token"), ns);
        String name = ">";
        String tempName = fastaFileIn.getName();
        int fileExtension = tempName.lastIndexOf(".");
        String filePrefix = tempName.substring(0, fileExtension);

        name = name.concat(filePrefix + newLine);
        bw.write(name);
        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            String dna = rec.seqString();
            bw.write(dna);

        }
        bw.write(newLine);
        bw.close();
    }

    /**
     * Concatenates all fasta sequences in multifasta file to a single fasta
     * sequence
     *
     * @param fastaFileIn input fasta file containing multiple fasta sequences
     * @param fastaFileOut output fasta containing reformatted sequence
     * @param prefix the prefix for incremental sequence names, e.g. 'Sequence_'
     * will become 'Sequence_1', 'Sequence_2', etc
     * @throws FileNotFoundException
     * @throws BioException
     * @throws Exception
     */
    public void reformatFasta(File fastaFileIn, File fastaFileOut, String prefix) throws FileNotFoundException, BioException, Exception
    {
        FileWriter fw = new FileWriter(fastaFileOut.getAbsoluteFile());
        BufferedWriter bw = new BufferedWriter(fw);
        String newLine = System.getProperty("line.separator");

        BufferedReader br = new BufferedReader(new FileReader(fastaFileIn));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                                                                       alpha.getTokenization("token"), ns);
        int seqId = 1;
        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            String name = ">";
            name = name.concat(prefix + seqId + newLine);
          
            bw.write(name);
            String dna = rec.seqString();
            bw.write(dna);
            bw.write(newLine);
            seqId++;
            System.out.print(name);
            System.out.println(dna);
        }
         bw.close();
    }

    /**
     * Concatenates all fasta sequences in multifasta file to a single fasta
     * sequence
     *
     * @param fastaFileIn input fasta file containing multiple fasta sequences
     * @param fastaFileOut output fasta containing single concatenated sequence
     * @throws FileNotFoundException
     * @throws BioException
     * @throws Exception
     */
    public void multiProtienFastaToSingleFasta(File fastaFileIn, File fastaFileOut) throws FileNotFoundException, BioException, Exception
    {
        FileWriter fw = new FileWriter(fastaFileOut.getAbsoluteFile());
        BufferedWriter bw = new BufferedWriter(fw);
        String newLine = System.getProperty("line.separator");

        BufferedReader br = new BufferedReader(new FileReader(fastaFileIn));
        Alphabet alpha = AlphabetManager.alphabetForName("PROTEIN-TERM");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                                                                       alpha.getTokenization("token"), ns);
        String name = ">";
        String tempName = fastaFileIn.getName();
        int fileExtension = tempName.lastIndexOf(".");
        String filePrefix = tempName.substring(0, fileExtension);

        name = name.concat(filePrefix + newLine);
        bw.write(name);
        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            String dna = rec.seqString();
            bw.write(dna);

        }
        bw.write(newLine);
        bw.close();
    }

    /**
     * Filters a fasta file by length
     *
     * @param fastaIn the fasta file to filter
     * @param fastaOut the filtered fasta file
     * @param minLength the minimum contig length to keep
     * @throws IOException
     * @throws BioException
     */
    public void filterFastaByLength(File fastaIn, File fastaOut, int minLength) throws IOException, BioException
    {
        FileWriter fw = new FileWriter(fastaOut.getAbsoluteFile());
        String newLine = System.getProperty("line.separator");
        FileInputStream inStream = new FileInputStream(fastaIn);
        FastaReader<DNASequence, NucleotideCompound> fastaReader = new FastaReader<>(inStream,
                                                                                     new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(),
                                                                                     new DNASequenceCreator(DNACompoundSet.getDNACompoundSet()));
        LinkedHashMap<String, DNASequence> b = fastaReader.process();
        for (Entry<String, DNASequence> entry : b.entrySet())
        {
            //System.out.println(entry.getValue().getOriginalHeader() + "=" + entry.getValue().getSequenceAsString());
            if (entry.getValue().getSequenceAsString().length() >= minLength)
            {
                fw.write(">" + entry.getValue().getOriginalHeader() + newLine + entry.getValue().getSequenceAsString() + newLine);
            }
        }
        fw.close();
    }

}
