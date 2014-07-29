/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fasta;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

/**
 *
 * @author ethering
 */
public class FastaParser
{

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
            
            String fastqSeq = rec.getName();
            String dna = rec.seqString();
            fastqSeq = fastqSeq.concat(newLine + dna + newLine + "+" + newLine);
            int len = rec.length();
            char[] chars = new char[len];
            Arrays.fill(chars, '#');
            String quality = new String(chars);
            fastqSeq = fastqSeq.concat(quality + newLine);
            bw.write(fastqSeq);
            
        }


        bw.close();
    }
}
