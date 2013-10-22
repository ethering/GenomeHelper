/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fastq;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;

/**
 *
 * @author ethering
 */
public class FastqCompression
{

    public void compressFastq(File fastqFile, File fastqOut) throws FileNotFoundException, IOException, ClassNotFoundException
    {
        FastqReader fqr = new FastqReader(fastqFile);
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodReads = writer.newWriter(fastqOut);

        Iterator<FastqRecord> it = fqr.iterator();

        //loop thru the interlaced file
        while (it.hasNext())
        {
            FastqRecord fastqRecord = it.next();
            if (!fastqRecord.getReadString().equalsIgnoreCase(""))
            {
                goodReads.write(fastqRecord);
            }
        }
        goodReads.close();
    }

    public static void main(String[] args)
    {
        try
        {
            File compFastq = new File("/Users/ethering/temp/solexa/lane1_NoIndex_R1_2000_qc.fastq");
            File compFastqOut = new File("/Users/ethering/temp/solexa/lane1_NoIndex_R1_2000_qc_test.fastq.gz");
            FastqCompression comp = new FastqCompression();
            comp.compressFastq(compFastq, compFastqOut);
        }
        catch (FileNotFoundException ex)
        {
            Logger.getLogger(FastqCompression.class.getName()).log(Level.SEVERE, null, ex);
        }
        catch (IOException | ClassNotFoundException ex)
        {
            Logger.getLogger(FastqCompression.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
}
