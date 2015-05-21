/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fastq;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;

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
