/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fastq;

import uk.ac.tsl.etherington.genomehelper.fastq.FastqCompression;
import java.io.File;
import java.io.RandomAccessFile;
import java.util.zip.GZIPInputStream;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author ethering
 */
public class FastqCompressionTest
{

    public FastqCompressionTest()
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
     * Tests that a given file is gzipped
     * @param f the file to test
     * @return true if gzipped, otherwise false
     */
    public static boolean isGZipped(File f)
    {
        int magic = 0;
        try
        {
            RandomAccessFile raf = new RandomAccessFile(f, "r");
            magic = raf.read() & 0xff | ((raf.read() << 8) & 0xff00);
            raf.close();
        }
        catch (Throwable e)
        {
            e.printStackTrace(System.err);
        }
        return magic == GZIPInputStream.GZIP_MAGIC;
    }

    /**
     * Test of compressFastq method, of class FastqCompression.
     */
    @Test
    public void testCompressFastq() throws Exception
    {
        System.out.println("compressFastq");
        File fastqFile = new File("test/test_data_in/piculus_test_left.fastq");
        File fastqOut = new File("test/test_data_out/piculus_test_left.fastq.gz");
        FastqCompression instance = new FastqCompression();
        //test that the input isn't compressed
        assertEquals(isGZipped(fastqFile), false);
        instance.compressFastq(fastqFile, fastqOut);
        //test that the output is compressed
        assertEquals(isGZipped(fastqOut), true);
    }
}