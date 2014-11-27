/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fasta;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojavax.Namespace;
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
        BufferedWriter bw = new BufferedWriter(fw);
        String newLine = System.getProperty("line.separator");
        OutputStream output = new FileOutputStream(fastaOut);

        BufferedReader br = new BufferedReader(new FileReader(fastaIn));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                alpha.getTokenization("token"), ns);
        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();

            if (rec.seqString().length() >= minLength)
            {
                RichSequence.IOTools.writeFasta(output, rec, ns);

            }

        }

    }
}
