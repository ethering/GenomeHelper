package piculus.gff;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author ethering
 */
public class GTF2GFF
{
    //this doesn't work yet!
    public void gtf2gff(File gtfFile, File gffFile) throws FileNotFoundException, IOException
    {
        //<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] 
        BufferedReader in = new BufferedReader(new FileReader(gtfFile));
        BufferedWriter out = new BufferedWriter(new FileWriter(gffFile));
        String line = null;
        while ((line = in.readLine()) != null)
        {
            String[] gtfLine = line.split("\t");
            String gffLine = (gtfLine[0] + "\t" + gtfLine[1] + "\t" + gtfLine[2] + "\t" + gtfLine[3] + "\t" + gtfLine[4] + "\t" + gtfLine[5] + "\t" + gtfLine[6] + "\t" + gtfLine[7] + "\t");
            out.write(gffLine);
            String attString = gtfLine[8].replaceAll("\"", "");
            //System.out.println("attString = " + attString);
            String[] attsArray = attString.split(";");
            for (int i = 0; i < attsArray.length; i++)
            {
                String currentAtt = attsArray[i];
                //System.out.println(currentAtt);
                currentAtt = currentAtt.trim();
                String[] atts = currentAtt.split(" ");

                out.write(atts[0] + "=");
                out.write(atts[1]);
                if (i < attsArray.length-1)
                {
                    out.write("; ");
                }
                

            }
            out.write("\n");
        }
    }

    public static void main(String args[])
    {
        GTF2GFF g2g = new GTF2GFF();
        File gtf = new File("/Users/ethering/temp/crinkler/pi_transcripts_supercontig1.2.gtf");
        File gff = new File("/Users/ethering/temp/crinkler/pi_transcripts_supercontig1.2.gff");
        try
        {
            g2g.gtf2gff(gtf, gff);
        }
        catch (FileNotFoundException ex)
        {
            Logger.getLogger(GTF2GFF.class.getName()).log(Level.SEVERE, null, ex);
        }
        catch (IOException ex)
        {
            Logger.getLogger(GTF2GFF.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
