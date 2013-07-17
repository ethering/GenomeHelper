/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fastq;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import net.sf.picard.fastq.FastqRecord;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.RNASequence;
import org.biojava3.sequencing.io.fastq.Fastq;
import org.biojava3.sequencing.io.fastq.FastqReader;
import org.biojava3.sequencing.io.fastq.FastqTools;
import org.biojava3.sequencing.io.fastq.SangerFastqReader;

/**
 *
 * @author ethering
 */
public class FastqToPrimer
{

    public void sortReads(File fastqFile)
    {
        net.sf.picard.fastq.FastqReader fq = new net.sf.picard.fastq.FastqReader(fastqFile);
        Iterator it = fq.iterator();
        
        while (it.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it.next();
        }
    }
    
    public void getPrimers(File fastqFile, File seqList) throws IOException
    {
        FastqReader fastqReader = new SangerFastqReader();
        List<DNASequence> sequences = new LinkedList<DNASequence>();
        DNASequence dna;
        RNASequence rna;
        ProteinSequence aa;
        for (Fastq fastq : fastqReader.read(fastqFile))
        {
            dna = FastqTools.createDNASequence(fastq);
            sequences.add(FastqTools.createDNASequence(fastq));
            rna = dna.getRNASequence();
            
            
            // for each frame
            for (int i = 0; i < 3; i++)
            {
                
                String currentRNA = rna.toString().substring(i, rna.toString().length());
                RNASequence newRNA = new RNASequence(currentRNA);
                aa = newRNA.getProteinSequence();
                
                

              
            }

        }
    }
}
