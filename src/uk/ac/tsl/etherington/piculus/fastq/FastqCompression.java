/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fastq;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
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
}
