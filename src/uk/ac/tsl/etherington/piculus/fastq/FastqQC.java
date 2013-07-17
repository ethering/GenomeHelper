/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fastq;

import java.io.File;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import net.sf.picard.util.SolexaQualityConverter;

/**
 *
 * @author ethering
 */
public class FastqQC
{

    public void qcPairedReads(File leftFastqFileIn, File rightFastqFileIn, File leftReadsOut, File rightReadsOut, int singleEndReadLength, String format, boolean writeBadSeqs)
    {
        FastqReader fql = new FastqReader(leftFastqFileIn);
        FastqReader fqr = new FastqReader(rightFastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodLeftSeqs = writer.newWriter(leftReadsOut);
        FastqWriter goodRightSeqs = writer.newWriter(rightReadsOut);

        FastqWriter badLeftSeqs = null;
        FastqWriter badRightSeqs = null;

        if (writeBadSeqs == true)
        {
            String badLeftString = "bad_".concat(leftReadsOut.toString());
            String badRightString = "bad_".concat(rightReadsOut.toString());
            File badLeftFile = new File(badLeftString);
            File badRightFile = new File(badRightString);
            badLeftSeqs = writer.newWriter(badLeftFile);
            badRightSeqs = writer.newWriter(badRightFile);
        }


        //create an interator for each file
        Iterator itl = fql.iterator();
        Iterator itr = fqr.iterator();
        int itCounter = 0;
        while (itl.hasNext())
        {
            //get the corresponding reads
            FastqRecord leftSeqRecord = (FastqRecord) itl.next();
            FastqRecord rightSeqRecord = (FastqRecord) itr.next();

            //get the length of them
            int leftSeqLength = leftSeqRecord.getReadString().length();
            int rightSeqLength = rightSeqRecord.getReadString().length();

            //check for N's
            String leftReadString = leftSeqRecord.getReadString();
            String rightReadString = rightSeqRecord.getReadString();
            boolean leftContainsN = leftReadString.contains("N");
            boolean rightContainsN = rightReadString.contains("N");
//            System.out.println("seqLength = "+seqLength);
//            System.out.println("rightSeqLength = "+rightSeqLength);
//            System.out.println("leftContainsN = "+leftContainsN);
//            System.out.println("rightContainsN = "+rightContainsN);

            if (leftSeqLength == singleEndReadLength && rightSeqLength == singleEndReadLength && leftContainsN == false && rightContainsN == false)
            {
                FastqRecord newLeftSeq = groomRead(leftSeqRecord, format);
                FastqRecord newRightSeq = groomRead(rightSeqRecord, format);
                goodLeftSeqs.write(newLeftSeq);
                goodRightSeqs.write(newRightSeq);
                itCounter++;
            }
            else if (writeBadSeqs == true)// && (leftSeqLength != singleEndReadLength || rightSeqLength != singleEndReadLength || leftContainsN == true || rightContainsN == true))
            {
                badLeftSeqs.write(leftSeqRecord);
                badRightSeqs.write(rightSeqRecord);
            }


        }
        System.out.println("Completed writing " + itCounter + " good reads");

    }
    public void qcPairedReadsGalaxy (File leftFastqFileIn, File rightFastqFileIn, File leftReadsOut, File rightReadsOut, int singleEndReadLength)
    {
        FastqReader fql = new FastqReader(leftFastqFileIn);
        FastqReader fqr = new FastqReader(rightFastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodLeftSeqs = writer.newWriter(leftReadsOut);
        FastqWriter goodRightSeqs = writer.newWriter(rightReadsOut);
        //create an interator for each file
        Iterator itl = fql.iterator();
        Iterator itr = fqr.iterator();
        int itCounter = 0;
        while (itl.hasNext())
        {
            //get the corresponding reads
            FastqRecord leftSeqRecord = (FastqRecord) itl.next();
            FastqRecord rightSeqRecord = (FastqRecord) itr.next();

            //get the length of them
            int leftSeqLength = leftSeqRecord.getReadString().length();
            int rightSeqLength = rightSeqRecord.getReadString().length();

            //check for N's
            String leftReadString = leftSeqRecord.getReadString();
            String rightReadString = rightSeqRecord.getReadString();
            boolean leftContainsN = leftReadString.contains("N");
            boolean rightContainsN = rightReadString.contains("N");
//            System.out.println("seqLength = "+seqLength);
//            System.out.println("rightSeqLength = "+rightSeqLength);
//            System.out.println("leftContainsN = "+leftContainsN);
//            System.out.println("rightContainsN = "+rightContainsN);

            if (leftSeqLength == singleEndReadLength && rightSeqLength == singleEndReadLength && leftContainsN == false && rightContainsN == false)
            {
                goodLeftSeqs.write(leftSeqRecord);
                goodRightSeqs.write(rightSeqRecord);
                itCounter++;
            }
        }
        System.out.println("Completed writing " + itCounter + " good reads");

    }

    public void qcInterlacedReads(File fastqFileIn, File leftReadsOut, File rightReadsOut, int readLength, String format, boolean writeBadSeqs)
    {
        FastqReader fq = new FastqReader(fastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodLeftSeqs = writer.newWriter(leftReadsOut);
        FastqWriter goodRightSeqs = writer.newWriter(rightReadsOut);

        FastqWriter badLeftSeqs = null;
        FastqWriter badRightSeqs = null;

        if (writeBadSeqs == true)
        {
            String badLeftString = "bad_".concat(leftReadsOut.toString());
            String badRightString = "bad_".concat(rightReadsOut.toString());
            File badLeftFile = new File(badLeftString);
            File badRightFile = new File(badRightString);
            badLeftSeqs = writer.newWriter(badLeftFile);
            badRightSeqs = writer.newWriter(badRightFile);
        }


        //create an interator for the interlaced file
        Iterator it = fq.iterator();

        int itCounter = 0;
        while (it.hasNext())
        {
            //get the corresponding reads
            FastqRecord leftSeqRecord = (FastqRecord) it.next();
            FastqRecord rightSeqRecord = (FastqRecord) it.next();

            //get the length of them
            int leftSeqLength = leftSeqRecord.getReadString().length();
            int rightSeqLength = rightSeqRecord.getReadString().length();

            //check for N's
            String leftReadString = leftSeqRecord.getReadString();
            String rightReadString = rightSeqRecord.getReadString();
            boolean leftContainsN = leftReadString.contains("N");
            boolean rightContainsN = rightReadString.contains("N");

            if (leftSeqLength == readLength && rightSeqLength == readLength && leftContainsN == false && rightContainsN == false)
            {

                FastqRecord newLeftSeq = groomRead(leftSeqRecord, format);
                FastqRecord newRightSeq = groomRead(rightSeqRecord, format);
                goodLeftSeqs.write(newLeftSeq);
                goodRightSeqs.write(newRightSeq);
                itCounter++;
            }
            if (writeBadSeqs == true && (leftSeqLength != readLength || rightSeqLength != readLength || leftContainsN == true || rightContainsN == true))
            {
                badLeftSeqs.write(leftSeqRecord);
                badRightSeqs.write(rightSeqRecord);
            }


        }
        System.out.println("Completed writing " + itCounter + " good reads");

    }

    public void qcJoinedReads(File fastqFileIn, File leftReadsOut, File rightReadsOut, int readLength, String format, boolean writeBadSeqs)
    {
        FastqReader fq = new FastqReader(fastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodLeftSeqs = writer.newWriter(leftReadsOut);
        FastqWriter goodRightSeqs = writer.newWriter(rightReadsOut);

        FastqWriter badSeqs = null;

        if (writeBadSeqs == true)
        {
            String badString = "bad_".concat(fastqFileIn.toString());
            File badFile = new File(badString);
            badSeqs = writer.newWriter(badFile);
        }


        //create an interator for each file
        Iterator it = fq.iterator();



        int itCounter = 0;
        while (it.hasNext())
        {
            //get the corresponding reads
            FastqRecord seqRecord = (FastqRecord) it.next();
            int seqLength = seqRecord.getReadString().length();
            String readString = seqRecord.getReadString();
            boolean containsN = readString.contains("N");

            if (seqLength / 2 == readLength && containsN == false)
            {
                String qual = seqRecord.getBaseQualityString();
                String leftRead = readString.substring(0, (seqLength / 2));
                String rightRead = readString.substring(seqLength / 2, seqLength);
                String leftQual = qual.substring(0, (seqLength / 2));
                String rightQual = qual.substring(seqLength / 2, seqLength);
                FastqRecord newLeftSeq = new FastqRecord(seqRecord.getReadHeader() + "/1", leftRead, "", leftQual);
                FastqRecord newRightSeq = new FastqRecord(seqRecord.getReadHeader() + "/2", rightRead, "", rightQual);
                goodLeftSeqs.write(newLeftSeq);
                goodRightSeqs.write(newRightSeq);
                itCounter++;
            }


            if (writeBadSeqs == true && (seqLength / 2 != readLength || containsN == true))
            {
                badSeqs.write(seqRecord);

            }


        }
        System.out.println("Completed writing " + itCounter + " good reads");

    }

    public void qcSingleEndReads(File fastqFileIn, File readsOut, int readLength, String format, boolean writeBadSeqs)
    {
        FastqReader fq = new FastqReader(fastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodSeqs = writer.newWriter(readsOut);

        FastqWriter badSeqs = null;


        if (writeBadSeqs == true)
        {
            String badString = "bad_".concat(readsOut.toString());
            File badFile = new File(badString);
            badSeqs = writer.newWriter(badFile);
        }


        //create an interator for the interlaced file
        Iterator it = fq.iterator();

        int itCounter = 0;
        while (it.hasNext())
        {
            //get the corresponding reads
            FastqRecord seqRecord = (FastqRecord) it.next();

            //get the length of them
            int seqLength = seqRecord.getReadString().length();


            //check for N's
            String readString = seqRecord.getReadString();

            boolean containsN = readString.contains("N");


            if (seqLength == readLength && containsN == false)
            {

                FastqRecord newSeq = groomRead(seqRecord, format);
                goodSeqs.write(newSeq);
                itCounter++;
            }
            if (writeBadSeqs == true && (seqLength != readLength || containsN == true))
            {
                badSeqs.write(seqRecord);

            }


        }
        System.out.println("Completed writing " + itCounter + " good reads");

    }

    public FastqRecord groomRead(FastqRecord record, String format)
    {
        String quality = "";


        if (format.equalsIgnoreCase("illumina"))
        {
            byte[] sangerQuals = record.getBaseQualityString().getBytes();
            SolexaQualityConverter.getSingleton().convertSolexa_1_3_QualityCharsToPhredBinary(sangerQuals);
            quality = net.sf.samtools.SAMUtils.phredToFastq(sangerQuals);
        }
        else if (format.equalsIgnoreCase("sanger"))
        {
            quality = record.getBaseQualityString();

        }
        else
        {
            System.err.print("format must be either 'sanger' or 'illumina'");
            System.exit(0);
        }


        FastqRecord newseq = new FastqRecord(record.getReadHeader(), record.getReadString(), "", quality);
        return newseq;

    }

    public void countLengthsAndNs(File fastq)
    {
        FastqReader fq = new FastqReader(fastq);
        int readsWithNs = 0;
        TreeMap<Integer, Integer> lengthDists = new TreeMap<Integer, Integer>();

        //create an interator for each file
        Iterator it = fq.iterator();

        int itCounter = 0;
        while (it.hasNext())
        {
            //get the corresponding reads
            FastqRecord seqRecord = (FastqRecord) it.next();


            //get the length of them
            int seqLength = seqRecord.getReadString().length();
            if (lengthDists.containsKey(seqLength))
            {
                Integer count = lengthDists.get(seqLength);
                count++;
                lengthDists.put(seqLength, count);
            }
            else
            {
                lengthDists.put(seqLength, 1);
            }

            //check for N's
            String leftReadString = seqRecord.getReadString();

            boolean leftContainsN = leftReadString.contains("N");

            if (leftContainsN == true)
            {
                readsWithNs++;
            }
            itCounter++;

        }
        System.out.println("Completed reading " + itCounter + " reads");
        for (Map.Entry<Integer, Integer> entry : lengthDists.entrySet())
        {
            Integer key = entry.getKey();
            Integer value = entry.getValue();

            System.out.println("length: " + key + " count: " + value);
        }
        System.out.println("Found " + readsWithNs + " read with at least one 'N'");

    }
    //just counts the reads and makes sure they look right

    public void veryfiyReads(File fastq)
    {
        System.out.println("I should exit with a read count. If not check the error message");
        FastqReader fq = new FastqReader(fastq);
        Iterator it = fq.iterator();

        int itCounter = 0;
        while (it.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it.next();
            itCounter++;
        }

        System.out.println("Counted " + itCounter + " reads");
    }

    public void veryfiyPairedReads(File fastqLeft, File fastqRight)
    {
        //System.out.println("I should exit with a read count. If not check the error message");
        FastqReader fql = new FastqReader(fastqLeft);
        Iterator itl = fql.iterator();
        FastqReader fqr = new FastqReader(fastqRight);
        Iterator itr = fqr.iterator();

        int itCounter = 0;
        while (itl.hasNext())
        {
            FastqRecord seqRecordLeft = (FastqRecord) itl.next();
            FastqRecord seqRecordRight = (FastqRecord) itr.next();

            String[] leftNameArray = seqRecordLeft.getReadHeader().split(" ");
            String[] rightNameArray = seqRecordRight.getReadHeader().split(" ");
            String leftName = leftNameArray[0];
            String rightName = rightNameArray[0];
            //System.out.println(leftName + " : "+rightName);
            if (!leftName.equalsIgnoreCase(rightName))
            {
                System.out.println("Found an error at read " + itCounter);
                System.err.println(leftName + " is not the same as " + rightName);
                System.exit(0);
            }
            itCounter++;
        }

        System.out.println("Counted " + itCounter + " reads which were all in order");
    }

    public static void main(String[] args)
    {
        File fastqL = new File("/Users/ethering/temp/solexa/lane1_NoIndex_R1_2000.fastq");
        File fastqR = new File("/Users/ethering/temp/solexa/lane1_NoIndex_R2_2000.fastq");
        FastqQC fq = new FastqQC();
        fq.veryfiyPairedReads(fastqL, fastqR);
    }
}
