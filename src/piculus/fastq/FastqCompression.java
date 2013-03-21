/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package piculus.fastq;

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
        //FastqWriterFactory writer = new FastqWriterFactory();
        FileOutputStream fos = new FileOutputStream(fastqOut);
        GZIPOutputStream gz = new GZIPOutputStream(fos);
        Iterator<FastqRecord> it = fqr.iterator();
        //int peCounter = 0;
        ObjectOutputStream oos = new ObjectOutputStream(gz);
        //loop thru the interlaced file
        while (it.hasNext())
        {
            FastqRecord fastqRecord = it.next();

            oos.writeObject(fastqRecord.getReadHeader());
            oos.writeObject(fastqRecord.getReadString());
            oos.writeObject("+");
            oos.writeObject(fastqRecord.getBaseQualityString());

            System.out.println(fastqRecord.getReadHeader());
            System.out.println(fastqRecord.getReadString());
            System.out.println("+");
            System.out.println(fastqRecord.getBaseQualityString());
        }
        oos.flush();
        oos.close();
        fos.close();
    }

    public void compressFastq2(File fastqFile, File fastqOut) throws FileNotFoundException, IOException, ClassNotFoundException
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
                System.out.println(fastqRecord.getReadHeader());
                System.out.println(fastqRecord.getReadString());
                System.out.println("+");
                System.out.println(fastqRecord.getBaseQualityString());
                goodReads.write(fastqRecord);

            }
        }
    }

    public static void main(String[] args)
    {
        try
        {
            File compFastq = new File("/Users/ethering/temp/solexa/lane1_NoIndex_R1_2000.fastq.gz");
            File compFastqOut = new File("/Users/ethering/temp/solexa/lane1_NoIndex_R1_2000_out.fastq.gz");
            FastqCompression comp = new FastqCompression();
            comp.compressFastq2(compFastq, compFastqOut);
        }
        catch (FileNotFoundException ex)
        {
            Logger.getLogger(FastqCompression.class.getName()).log(Level.SEVERE, null, ex);
        }
        catch (IOException ex)
        {
            Logger.getLogger(FastqCompression.class.getName()).log(Level.SEVERE, null, ex);
        }
        catch (ClassNotFoundException ex)
        {
            Logger.getLogger(FastqCompression.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
}
