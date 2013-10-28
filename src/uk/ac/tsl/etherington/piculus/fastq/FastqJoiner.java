/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fastq;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;

/**
 *A class to join and split fastq files
 * @author ethering
 */
public class FastqJoiner
{

    File leftReads;
    File rightReads;
    File fastqInterlacedFile;
    File fastqSinglesFile;

    /**
     * 
     * @param leftReads the left-handed fastq reads
     * @param rightReads the right-handed fastq reads
     * @param fastqJoinedFile the joined pairs
     * @param fastqSinglesFile single reads that don't have mates
     */
    public void join(File leftReads, File rightReads, File fastqJoinedFile, File fastqSinglesFile)
    {
        FastqReader fql = new FastqReader(leftReads);
        FastqReader fqr;
        FastqWriterFactory writer = new FastqWriterFactory();
        //for the matching paired-end reads
        FastqWriter pairedSeqs = writer.newWriter(fastqJoinedFile);
        //for unpaired singles
        FastqWriter singleSeqs = writer.newWriter(fastqSinglesFile);
        Set<String> pairs = new HashSet<>();
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
                    String joinedSeq = leftSeqRecord.getReadString().concat(rightSeqRecord.getReadString());
                    String joinedQual = leftSeqRecord.getBaseQualityString().concat(rightSeqRecord.getBaseQualityString());
                    FastqRecord joinedRead = new FastqRecord(leftReadName, joinedSeq, "", joinedQual);
                    pairedSeqs.write(joinedRead);

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
        pairedSeqs.close();
        singleSeqs.close();
        System.out.println("Completed writing " + peCounter + " paired-reads");
        System.out.println("Completed writing " + singleCounter + " singles");
    }

    /**
     * 
     * @param fastqFile a fastq file of joined left-handed and right-handed reads
     * @param leftPairdReads the split left-handed reads
     * @param rightPairedReads the split right-handed reads
     */
    public void split(File fastqFile, File leftPairdReads, File rightPairedReads)
    {
        FastqReader fqr = new FastqReader(fastqFile);
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter leftPairedSeqs = writer.newWriter(leftPairdReads);
        FastqWriter rightPairedSeqs = writer.newWriter(rightPairedReads);
        Iterator<FastqRecord> it = fqr.iterator();
        int peCounter = 0;

        //loop thru the interlaced file
        while (it.hasNext())
        {
            FastqRecord fastqRecord = it.next();
            //get the read name
            String readHeader = fastqRecord.getReadHeader();
            //get Sequence
            String seq = fastqRecord.getReadString();
            
            String qualString = fastqRecord.getBaseQualityString();
            int readlength = seq.length();
            int midpoint = readlength/2;

            String leftHeader = readHeader.concat(" 1:N:0:");
            String rightHeader = readHeader.concat(" 2:N:0:");

            String leftRead = seq.substring(0, midpoint);
            String rightRead = seq.substring(midpoint, readlength);

            String leftQual = qualString.substring(0, midpoint);
            String rightQual = qualString.substring(midpoint, readlength);

            FastqRecord leftSeq = new FastqRecord(leftHeader, leftRead, "", leftQual);
            FastqRecord rightSeq = new FastqRecord(rightHeader, rightRead, "", rightQual);
            leftPairedSeqs.write(leftSeq);
            rightPairedSeqs.write(rightSeq);
            peCounter++;

        }
        leftPairedSeqs.close();
        rightPairedSeqs.close();
        System.out.println("Completed spliting " + peCounter + " paired-reads");
    }
    
//    public static void main(String[] args) throws IOException
//    {
//        File leftFastqFileIn = new File("/Users/ethering/NetBeansProjects/Piculus/test/test_data_in/piculus_test_left.fastq");
//        File rightFastqFileIn = new File("/Users/ethering/NetBeansProjects/Piculus/test/test_data_in/piculus_test_right.fastq");
//        File out = new File("/Users/ethering/NetBeansProjects/Piculus/test/test_data_in/piculus_test_joined.fastq");
//        File singles = new File("/Users/ethering/NetBeansProjects/Piculus/test/test_data_out/piculus_test_joined_singles.fastq");
//
//        FastqJoiner fqj = new FastqJoiner();
//        fqj.join(leftFastqFileIn, rightFastqFileIn, out, singles);
//    }
}
