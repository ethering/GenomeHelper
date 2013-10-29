/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fasta;

import uk.ac.tsl.etherington.genomehelper.fasta.FastaTranslator;
import java.io.File;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author ethering
 */
public class FastaTranslatorTest
{

    public FastaTranslatorTest()
    {
    }

    @BeforeClass
    public static void setUpClass()
    {
    }

    @AfterClass
    public static void tearDownClass()
    {
    }

      /**
     * Test of translateFasta method, of class FastaTranslator.
     */
    @Test
    public void testTranslateFasta() throws Exception
    {
        System.out.println("translateFasta");
        AccessionID id = new AccessionID("seq1");
        DNASequence dna = new DNASequence("atg");
        dna.setAccession(id);
        FastaTranslator instance = new FastaTranslator();
        ProteinSequence result = instance.translateFasta(dna);
        assertEquals(result.toString(), "M");
    }
}