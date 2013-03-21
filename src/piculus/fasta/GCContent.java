/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package piculus.fasta;

import java.io.*;
import java.util.NoSuchElementException;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.*;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

public class GCContent
{

    public void getGCContent(String fileName) throws FileNotFoundException, NoSuchElementException, BioException
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
    }
//    public static void main (String [] args)
//    {
//        try
//        {
//            String file = args[0];
//            GCContent gcc = new GCContent();
//            gcc.getGCContent(file);
//        }
//        catch (FileNotFoundException ex)
//        {
//            Logger.getLogger(GCContent.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        catch (NoSuchElementException ex)
//        {
//            Logger.getLogger(GCContent.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        catch (BioException ex)
//        {
//            Logger.getLogger(GCContent.class.getName()).log(Level.SEVERE, null, ex);
//        }
//    }
}