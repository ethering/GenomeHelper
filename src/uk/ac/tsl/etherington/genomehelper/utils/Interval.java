/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.utils;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import picard.sam.ViewSam;

/**
 *
 * @author ethering
 */
public class Interval
{

    public void gatkToSamInterval(File bamFile, File gatkInterval)
    {

        try
        {

            String outFileString = bamFile.getAbsolutePath();
            int fileDot = outFileString.lastIndexOf(".");
            String startPath = outFileString.substring(0, fileDot);
            File outFile = new File(startPath + ".interval_list");
            try (PrintWriter out = new PrintWriter(outFile))
            {
                ByteArrayOutputStream baos = new ByteArrayOutputStream();
                PrintStream ps = new PrintStream(baos);
                // IMPORTANT: Save the old System.out!
                PrintStream old = System.out;
                // Tell Java to use your special stream
                System.setOut(ps);
                //create a new instance of ViewSam
                ViewSam cmd = new ViewSam();
                
                //invoke the instance
                int rec = cmd.instanceMain(new String[]
                {
                    "INPUT=" + bamFile.getPath(),
                    "HEADER_ONLY=" + true,
                    
                });
                
                // Show what happened
                out.print(baos.toString());
                System.out.flush();
                ps.close();
                baos.close();
                System.setOut(old);
                
                Charset charset = Charset.forName("US-ASCII");
                BufferedReader reader = Files.newBufferedReader(gatkInterval.toPath(), charset);
                int targetCount = 1;
                String line = null;
                while ((line = reader.readLine()) != null)
                {
                    //if the line is not a comment
                    if (!line.startsWith("#"))
                    {
                        String[] lineArray = line.split(":");
                        String chr = lineArray[0];
                        String loc = lineArray[1];
                        String[] locArray = loc.split("-");
                        int start = Integer.parseInt(locArray[0]);
                        
                        out.print(chr + "\t" + start);
                        if (locArray.length > 1)
                        {
                            int end = Integer.parseInt(locArray[1]);
                            out.print("\t" + end);
                        } else
                        {
                            out.print("\t" + start + 1);
                        }
                        
                        out.print("\t+\ttarget_" + targetCount + "\n");
                        targetCount++;
                    }
                }
                
                out.close();
            }

        } catch (IOException | NumberFormatException e)
        {
            System.err.println(e.getMessage());
        }
        

    }

    public static void main(String[] args)
    {
        //File bamFile = new File("/Users/ethering/temp/test_data/genomehelper/SRR1508214_scaffold_unmapped_100_merged_reheader.bam");
        //File outFile = new File("/Users/ethering/temp/test_data/genomehelper/SRR1508214_scaffold_unmapped_100_merged_reheader_out.bam");
        //File gatkInteval = new File("/Users/ethering/temp/test_data/genomehelper/interval.txt");
        //File options = new File("/Users/ethering/temp/test_data/genomehelper/options.txt");
        File bamFile = new File(args[0]);
        File gatkInteval = new File(args[1]);
        Interval i = new Interval();
        i.gatkToSamInterval(bamFile, gatkInteval);
    }

}
