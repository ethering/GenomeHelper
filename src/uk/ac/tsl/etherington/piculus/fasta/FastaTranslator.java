/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fasta;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.Writer;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.RNASequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;

/**
 *
 * @author ethering
 */
public class FastaTranslator
{
    /**
     * Translates a DNA multi-fasta file into ammino acids
     * @param fastaDnaFile the input DNA multi-fasta file
     * @param fastaProteinFile the output protein multi-fasta file
     * @throws Exception 
     */
    public void translateMultiFasta(File fastaDnaFile, File fastaProteinFile) throws Exception
    {
        LinkedHashMap<String, DNASequence> dnaSeqs = FastaReaderHelper.readFastaDNASequence(fastaDnaFile);
        //FastaReaderHelper.readFastaDNASequence for DNA sequences
        Writer out = new BufferedWriter(new FileWriter(fastaProteinFile));
        for (Entry<String, DNASequence> entry : dnaSeqs.entrySet())
        {
            //System.out.println(entry.getValue().getOriginalHeader() + "=" + entry.getValue().getSequenceAsString());
            RNASequence rna = entry.getValue().getRNASequence();
            ProteinSequence aa = rna.getProteinSequence();
            out.write(">" + entry.getKey() + "\n");
            out.write(aa + "\n");
        }
        out.close();
    }
    /**
     * Translate a DNASequence into a ProteinSequence
     * @param dna the DNASequence to be translated
     * @return the translated ProteinSequence
     * @throws Exception 
     */
    public ProteinSequence translateFasta(DNASequence dna) throws Exception
    {
        RNASequence rna = dna.getRNASequence();
        ProteinSequence aa = rna.getProteinSequence();
        System.out.println(dna.getOriginalHeader());
        System.out.println(aa);
        AccessionID id = new AccessionID(dna.getOriginalHeader());
        aa.setAccession(id);
        return aa;
    }
   
}
