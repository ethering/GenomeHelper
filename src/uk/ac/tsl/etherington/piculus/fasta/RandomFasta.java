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
    /**
     * 
     * @param fasta
     * @param numberOfSeqsRequired
     * @param outfile
     * @throws Exception 
     */
    public void selectRandomSequences(File fasta, int numberOfSeqsRequired, File outfile) throws Exception
    {
        //read in the fasta sequences to a List
        List<FastaSequence> seqs = SequenceUtil.readFasta(new FileInputStream(fasta));
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
    }
        
   }
