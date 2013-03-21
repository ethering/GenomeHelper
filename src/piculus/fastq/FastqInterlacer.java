/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package piculus.fastq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;

/**
 *
 * @author ethering
 */
public class FastqInterlacer
{

    public void interlace(File leftReads, File rightReads, File fastqInterlacedFile, File fastqSinglesFile)
    {
        FastqReader fql = new FastqReader(leftReads);
        FastqReader fqr;
        FastqWriterFactory writer = new FastqWriterFactory();
        //for the matching paired-end reads
        FastqWriter pairedSeqs = writer.newWriter(fastqInterlacedFile);
        //for unpaired singles
        FastqWriter singleSeqs = writer.newWriter(fastqSinglesFile);
        Set<String> pairs = new HashSet<String>();
        int peCounter = 0;
        int singleCounter = 0;

        Iterator lit1 = fql.iterator();
        //loop thru the left hand reads 
        while (lit1.hasNext())
        {
            FastqRecord leftSeqRecord = (FastqRecord) lit1.next();
            //get the read header..
            String leftReadHeader = leftSeqRecord.getReadHeader();
            String[] leftArray = leftReadHeader.split(" ");
            //..and remove the pairing info (e.g. 1:N:0:)
            String leftReadName = leftArray[0];
            fqr = new FastqReader(rightReads);

            Iterator rit1 = fqr.iterator();
            boolean pairFound = false;
            while (rit1.hasNext() && pairFound == false)
            {
                //get the read header and find out the read name
                FastqRecord rightSeqRecord = (FastqRecord) rit1.next();
                String rightReadHeader = rightSeqRecord.getReadHeader();
                String[] rightArray = rightReadHeader.split(" ");
                String rightReadName = rightArray[0];

                if (leftReadName.equalsIgnoreCase(rightReadName))//it's a paired end reads
                {
                    //and write the paired reads to their file
                    pairedSeqs.write(leftSeqRecord);
                    pairedSeqs.write(rightSeqRecord);
                    //remove the current object
                    pairFound = true;
                    pairs.add(leftReadName);
                    peCounter++;
                }
            }
            if (pairFound == false)
            {
                singleSeqs.write(leftSeqRecord);
                singleCounter++;
            }
            fqr.close();
        }
        fql.close();
        //loop through the right hand reads
        fqr = new FastqReader(rightReads);
        Iterator rit2 = fqr.iterator();

        while (rit2.hasNext())
        {
            FastqRecord rightSeqRecord = (FastqRecord) rit2.next();
            String rightReadHeader = rightSeqRecord.getReadHeader();
            String[] rightArray = rightReadHeader.split(" ");
            String rightReadName = rightArray[0];


            if (!pairs.contains(rightReadName))
            {
                singleSeqs.write(rightSeqRecord);
                singleCounter++;
            }
        }

        System.out.println("Completed writing " + peCounter + " paired-reads");
        System.out.println("Completed writing " + singleCounter + " singles");
    }

    //if you know your reads match up, use this one
    public void interlaceKnownPairs(File leftReads, File rightReads, File fastqInterlacedFile) throws IOException
    {

        FastqReader lfq = new FastqReader(leftReads);
        FastqReader rfq = new FastqReader(rightReads);

        FastqWriterFactory writer = new FastqWriterFactory();
        //for the matching paired-end reads
        FastqWriter interlacedSeqs = writer.newWriter(fastqInterlacedFile);
        int peCounter = 0;
        Iterator lit = lfq.iterator();
        Iterator rit = rfq.iterator();
        //loop thru the left hand reads
        while (lit.hasNext())
        {
            FastqRecord leftSeqRecord = (FastqRecord) lit.next();
            if (!rit.hasNext())
            {
                System.err.println("More reads in left file than right");
                System.exit(1);
            }
            else
            {
                FastqRecord rightSeqRecord = (FastqRecord) rit.next();
                //and write the paired reads to their file
                interlacedSeqs.write(leftSeqRecord);
                interlacedSeqs.write(rightSeqRecord);
                peCounter++;
            }
        }
        if (rit.hasNext())
        {
            System.err.println("More reads in right file than left");
            System.exit(1);
        }
        rfq.close();
        lfq.close();
        interlacedSeqs.close();
        System.out.println("Interlaced " + peCounter + " paired-end reads");
    }

    public void printFastq(File fastqFile)
    {
        FastqReader fq = new FastqReader(fastqFile);


        Iterator it = fq.iterator();

        //loop thru the left hand reads 
        while (it.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it.next();
            System.out.println(seqRecord.getReadHeader());

        }

        fq.close();

    }

    public void deinterlace(File fastqFile, File leftPairdReads, File rightPairedReads, File leftSingleReads, File rightSingleReads)
    {
        FastqReader fqr = new FastqReader(fastqFile);
        FastqReader fqr2;
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter leftPairedSeqs = writer.newWriter(leftPairdReads);
        FastqWriter rightPairedSeqs = writer.newWriter(rightPairedReads);
        FastqWriter leftSingleSeqs = writer.newWriter(leftSingleReads);
        FastqWriter rightSingleSeqs = writer.newWriter(rightSingleReads);
        //HashMap<String, PEFastqRead> reads = new HashMap<String, PEFastqRead>();
        Iterator<FastqRecord> it = fqr.iterator();

        int peCounter = 0;
        int lrCounter = 0;
        int rrCounter = 0;
        //loop thru the interlaced file
        while (it.hasNext())
        {

            FastqRecord fastqRecord1 = it.next();
            String readHeader = fastqRecord1.getReadHeader();
            String[] array1 = readHeader.split(" ");
            //get the read name
            String readName1 = array1[0];
            String pair1 = array1[1];
            fqr2 = new FastqReader(fastqFile);
            Iterator<FastqRecord> it2 = fqr2.iterator();
            boolean pairFound = false;
            //System.out.println("\nComparing: "+readName1+" "+pair1);
            while (it2.hasNext() && pairFound == false)
            {
                FastqRecord fastqRecord2 = it2.next();
                String readHeader2 = fastqRecord2.getReadHeader();
                String[] array2 = readHeader2.split(" ");
                //get the read name
                String readName2 = array2[0];
                String pair2 = array2[1];

                if (readName1.equalsIgnoreCase(readName2) && !pair1.equalsIgnoreCase(pair2))//the read names match but they're not the same end
                {

                    if (pair1.startsWith("1:N:0:"))//it's the left hand read
                    {
                        leftPairedSeqs.write(fastqRecord1);
                        rightPairedSeqs.write(fastqRecord2);
                        peCounter++;
                    }
                    pairFound = true;
                }
            }

            if (pairFound == false)
            {
                if (pair1.startsWith("1:N:0:"))//it's a single left hand read
                {
                    leftSingleSeqs.write(fastqRecord1);
                    lrCounter++;
                }
                else //it's a single right read
                {
                    rightSingleSeqs.write(fastqRecord1);
                    rrCounter++;
                }
            }
            fqr2.close();
        }

        System.out.println("Completed writing " + peCounter + " paired-reads");
        System.out.println("Completed writing " + lrCounter + " left-reads");
        System.out.println("Completed writing " + rrCounter + " right-reads");

    }
//    public static void main(String args[])
//    {
//
//        File leftReads = new File("/Users/ethering/temp/fastqapp_data/new_format_paired_left.fastq");
//        File rightReads = new File("/Users/ethering/temp/fastqapp_data/new_format_paired_right.fastq");
//        File zippedFastqInterlacedFile = new File("/Users/ethering/temp/fastqapp_data/interlaced.fq.gz");
//        File log = new File("/Users/ethering/temp/fastqapp_data/log.txt");
//
//
//        Date start = new Date();
//        long l1 = start.getTime();
//        FastqInterlacer fi2 = new FastqInterlacer();
//        try
//        {
//            fi2.interlaceKnownPairs(leftReads, rightReads, zippedFastqInterlacedFile, log);
//        }
//        catch (IOException ex)
//        {
//            Logger.getLogger(FastqInterlacer.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        catch (InterruptedException ex)
//        {
//            Logger.getLogger(FastqInterlacer.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        //fi2.interlaceKnownPairs(leftReads, rightReads, fastqInterlacedFile);
//        Date stop = new Date();
//
//        long l2 = stop.getTime();
//        long diff = l2 - l1;
//        System.out.println("Took " + diff);
//        //fi2.printFastq(fastqInterlacedFile);
//    }
}