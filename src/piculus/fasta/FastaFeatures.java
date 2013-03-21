/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package piculus.fasta;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.Writer;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
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
public class FastaFeatures
{

    public static HashMap getSequenceLengths(File filename) throws FileNotFoundException, BioException
    {

        HashMap<String, Integer> seqLengths = new HashMap<String, Integer>();
        BufferedReader br = new BufferedReader(new FileReader(filename));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                alpha.getTokenization("token"), ns);
        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            seqLengths.put(rec.getName(), rec.length());
//            System.out.println(rec.getName());
//            System.out.println(rec.length());
        }
        return seqLengths;
    }
    public static LinkedHashMap<String, DNASequence> getParsedDNASequences(File refSeq) throws Exception
    {
        LinkedHashMap<String, DNASequence> tempgenome = FastaReaderHelper.readFastaDNASequence(refSeq);
        LinkedHashMap<String, DNASequence> genome = new LinkedHashMap<String, DNASequence>();
        for (Map.Entry<String, DNASequence> entry : tempgenome.entrySet())
        {
            String seqName = entry.getKey();
            String newSeqName = new String(seqName.split(" ")[0]);
            DNASequence dna = entry.getValue();
            genome.put(newSeqName, dna);
            System.out.println(newSeqName + " : " + entry.getKey());
        }
        tempgenome.clear();
        return genome;
    }
    
    public static double getGenomeSize(File filename) throws FileNotFoundException, BioException
    {
        double genomeSize = 0;
        HashMap<String, Integer> seqLengths = new HashMap<String, Integer>(FastaFeatures.getSequenceLengths(filename));
        for (Map.Entry<String, Integer> entry : seqLengths.entrySet())
        {
            
            int genomeLength = entry.getValue();
            genomeSize+=genomeLength;
        }
        return genomeSize;
    }

    public static HashMap getSequenceAsIntArray(File filename) throws FileNotFoundException, BioException
    {

        HashMap<String, Integer> seqLengths = new HashMap<String, Integer>(FastaFeatures.getSequenceLengths(filename));
        HashMap<String, int[]> codingMap = new HashMap<String, int[]>();

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

    public void seqFromCommandLine(File fastaIn, File outfile, String seqid) throws FileNotFoundException, Exception
    {
        Writer out = new BufferedWriter(new FileWriter(outfile));

        BufferedReader br = new BufferedReader(new FileReader(fastaIn));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                alpha.getTokenization("token"), ns);
        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            String id = rec.getName();
            if (id.equalsIgnoreCase(seqid))
            {
                String seq = rec.seqString();
                out.write(">" + id + "\n");
                out.write(seq + "\n");
                break;
            }
        }
        out.close();
    }
    public void seqFromCommandLine(File fastaIn, File outfile, String seqid, int start, int end) throws FileNotFoundException, Exception
    {
        int subseqLength = end -start;
        System.out.println("Requested subsequence length is "+subseqLength);
        Writer out = new BufferedWriter(new FileWriter(outfile));

        BufferedReader br = new BufferedReader(new FileReader(fastaIn));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                alpha.getTokenization("token"), ns);
        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            String id = rec.getName();
            if (id.equalsIgnoreCase(seqid))
            {
                String seq = rec.seqString();
                String subseq = seq.substring(start, end);
                subseqLength = subseq.length();
                System.out.println("Provided subsequence length is "+subseqLength);
                
                out.write(">" + id + "\n");
                out.write(subseq + "\n");
                break;
            }
        }
        out.close();
    }

    public static void main(String[] args)
    {
        File gff = new File("/Users/ethering/projects/oomycete_genomes/genome_biology/pinfestans/phytophthora_infestans_t30-4_1_supercontigs-3_edited.fasta");
        FastaFeatures ff = new FastaFeatures();
        try
        {
            HashMap<String, Integer> seqLengths = ff.getSequenceLengths(gff);
        }
        catch (Exception ex)
        {
            Logger.getLogger(FastaFeatures.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
