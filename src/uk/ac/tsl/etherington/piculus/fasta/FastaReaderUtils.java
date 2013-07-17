/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fasta;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.DNASequenceCreator;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;

/**
 *
 * @author ethering
 */
public class FastaReaderUtils
{

    public void FastaOpen(File f) throws IOException
    {
        FileInputStream inStream = new FileInputStream(f);
        FastaReader<DNASequence, NucleotideCompound> fastaReader =
                new FastaReader<DNASequence, NucleotideCompound>(
                inStream,
                new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(),
                new DNASequenceCreator(AmbiguityDNACompoundSet.getDNACompoundSet()));
        LinkedHashMap<String, DNASequence> b = fastaReader.process();
        for (Entry<String, DNASequence> entry : b.entrySet())
        {
            System.out.println(entry.getValue().getOriginalHeader() + "=" + entry.getValue().getSequenceAsString());
        }
    }
}
