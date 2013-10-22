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
 * A Quality Control (QC) class that takes fastq sequences in various structure
 * and quality format. The sequences are checked for Ns, long/short reads and
 * format. With paired-end reads, if either of the pairs fails QC then neither
 * of the pair is returned. Reads with the quality format of 'illumina'
 * (pre-Illumina 1.8) are automatically changed to sanger quality format.
 * Currently only accepts Illumina 1.3+
 *
 * @author ethering
 */
public class FastqQC
{

    public boolean checkFormat(String format)
    {
        boolean goodFormat = true;
        if (!format.equalsIgnoreCase("illumina") && !format.equalsIgnoreCase("sanger"))
        {
            goodFormat = false;
            System.err.print("format must be either 'sanger' or 'illumina'");
            System.exit(0);
        }
        return goodFormat;
    }

    /**
     * Removes pairs of reads where at least one of the pair contains short/long
     * reads or a read that contain an 'N', from paired-end fastq files.
     * Optionally writes bad reads to a file.
     *
     * @param leftFastqFileIn the left-handed reads
     * @param rightFastqFileIn the right-handed reads
     * @param leftReadsOut the QCd left-handed reads
     * @param rightReadsOut the QCd right-handed reads
     * @param singleEndReadLength the expected length of a single read
     * @param format the fastq format (can only be 'illumina' or 'sanger')
     * @param writeBadSeqs whether to write the bad reads to a file (bad reads
     * file name will start with 'bad_')
     */
    public void qcPairedReads(File leftFastqFileIn, File rightFastqFileIn, File leftReadsOut, File rightReadsOut, int singleEndReadLength, String format, boolean writeBadSeqs)
    {
        FastqReader fql = new FastqReader(leftFastqFileIn);
        FastqReader fqr = new FastqReader(rightFastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodLeftSeqs = writer.newWriter(leftReadsOut);
        FastqWriter goodRightSeqs = writer.newWriter(rightReadsOut);

        FastqWriter badLeftSeqs = null;
        FastqWriter badRightSeqs = null;
        int itCounter = 0;
        if (checkFormat(format) == true)
        {
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

            while (itl.hasNext())
            {
                //get the corresponding reads
                FastqRecord leftSeqRecord = (FastqRecord) itl.next();
                FastqRecord rightSeqRecord = (FastqRecord) itr.next();
                boolean leftGood = qcFastq(leftSeqRecord, singleEndReadLength);
                boolean rightGood = qcFastq(rightSeqRecord, singleEndReadLength);

                if (leftGood && rightGood)
                {
                    FastqRecord newLeftSeq = groomRead(leftSeqRecord, format);
                    FastqRecord newRightSeq = groomRead(rightSeqRecord, format);
                    goodLeftSeqs.write(newLeftSeq);
                    goodRightSeqs.write(newRightSeq);
                    itCounter++;
                }
                else if (writeBadSeqs == true)
                {
                    badLeftSeqs.write(leftSeqRecord);
                    badRightSeqs.write(rightSeqRecord);
                }
            }
        }
        System.out.println("Completed writing " + itCounter + " good reads");
    }

    public boolean qcFastq(FastqRecord seqRecord, int singleEndReadLength)
    {
        //get the length of them
        int seqLength = seqRecord.getReadString().length();
        String readString = seqRecord.getReadString();
        boolean containsN = readString.contains("N");
        if (seqLength == singleEndReadLength && containsN == false)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    /**
     * Removes pairs of reads where at least one of the pair contains short/long
     * reads or a read that contain an 'N', from an interlaced fastq file
     *
     * @param fastqFileIn an interlaced fastq file
     * @param fastqFileOut the good left-handed reads
     * @param rightReadsOut the good right-handed reads
     * @param readLength the expected length of a single read
     * @param format the fastq format (can only be 'illumina' or 'sanger')
     * @param writeBadSeqs whether to write the bad reads to a file (bad reads
     * file name will start with 'bad_')
     */
    public void qcInterlacedReads(File fastqFileIn, File fastqFileOut, int singleEndReadLength, String format, boolean writeBadSeqs)
    {
        FastqReader fq = new FastqReader(fastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodSeqs = writer.newWriter(fastqFileOut);
        FastqWriter badSeqs = null;
        int itCounter = 0;

        if (checkFormat(format) == true)
        {
            if (writeBadSeqs == true)
            {
                String badString = "bad_".concat(fastqFileOut.toString());
                File badFile = new File(badString);
                badSeqs = writer.newWriter(badFile);
            }

            //create an interator for the interlaced file
            Iterator it = fq.iterator();
            while (it.hasNext())
            {
                //get the corresponding reads
                FastqRecord leftSeqRecord = (FastqRecord) it.next();
                FastqRecord rightSeqRecord = (FastqRecord) it.next();

                boolean leftGood = qcFastq(leftSeqRecord, singleEndReadLength);
                boolean rightGood = qcFastq(rightSeqRecord, singleEndReadLength);

                if (leftGood && rightGood)
                {
                    FastqRecord newLeftSeq = groomRead(leftSeqRecord, format);
                    FastqRecord newRightSeq = groomRead(rightSeqRecord, format);
                    goodSeqs.write(newLeftSeq);
                    goodSeqs.write(newRightSeq);
                    itCounter++;
                }
                else if (writeBadSeqs == true)
                {
                    badSeqs.write(leftSeqRecord);
                    badSeqs.write(rightSeqRecord);
                }
            }
        }

        System.out.println("Completed writing " + itCounter + " good reads");
    }

    /**
     * Removes pairs of reads where at least one of the pair contains short/long
     * reads or a read that contain an 'N', from an interlaced fastq file
     *
     * @param fastqFileIn an interlaced fastq file
     * @param leftReadsOut the good left-handed reads
     * @param rightReadsOut the good right-handed reads
     * @param readLength the expected length of a single read
     * @param format the fastq format (can only be 'illumina' or 'sanger')
     * @param writeBadSeqs whether to write the bad reads to a file (bad reads
     * file name will start with 'bad_')
     */
    public void qcInterlacedReadsToPairs(File fastqFileIn, File leftReadsOut, File rightReadsOut, int singleEndReadLength, String format, boolean writeBadSeqs)
    {
        FastqReader fq = new FastqReader(fastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodLeftSeqs = writer.newWriter(leftReadsOut);
        FastqWriter goodRightSeqs = writer.newWriter(rightReadsOut);

        FastqWriter badLeftSeqs = null;
        FastqWriter badRightSeqs = null;
        int itCounter = 0;

        if (checkFormat(format) == true)
        {
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
            while (it.hasNext())
            {
                //get the corresponding reads
                FastqRecord leftSeqRecord = (FastqRecord) it.next();
                FastqRecord rightSeqRecord = (FastqRecord) it.next();

                 boolean leftGood = qcFastq(leftSeqRecord, singleEndReadLength);
                boolean rightGood = qcFastq(rightSeqRecord, singleEndReadLength);

                if (leftGood && rightGood)
                {
                    FastqRecord newLeftSeq = groomRead(leftSeqRecord, format);
                    FastqRecord newRightSeq = groomRead(rightSeqRecord, format);
                    goodLeftSeqs.write(newLeftSeq);
                    goodRightSeqs.write(newRightSeq);
                    itCounter++;
                }
                else if (writeBadSeqs == true)
                {
                    badLeftSeqs.write(leftSeqRecord);
                    badRightSeqs.write(rightSeqRecord);
                }
            }
        }

        System.out.println("Completed writing " + itCounter + " good reads");
    }

    /**
     * Removes pairs of reads where at least one of the pair contains short/long
     * reads or a read that contain an 'N', from a fastq file where the reads
     * are joined
     *
     * @param fastqFileIn an interlaced fastq file
     * @param leftReadsOut the good left-handed reads
     * @param rightReadsOut the good right-handed reads
     * @param readLength the expected length of a single read (i.e. half the
     * length of one of the joined reads)
     * @param format the fastq format (can only be 'illumina' or 'sanger')
     * @param writeBadSeqs whether to write the bad reads to a file (bad reads
     * file name will start with 'bad_')
     */
    public void qcJoinedReads(File fastqFileIn, File leftReadsOut, File rightReadsOut, int readLength, String format, boolean writeBadSeqs)
    {
        FastqReader fq = new FastqReader(fastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodLeftSeqs = writer.newWriter(leftReadsOut);
        FastqWriter goodRightSeqs = writer.newWriter(rightReadsOut);
        FastqWriter badSeqs = null;
        int itCounter = 0;

        if (checkFormat(format) == true)
        {
            if (writeBadSeqs == true)
            {
                String badString = "bad_".concat(fastqFileIn.toString());
                File badFile = new File(badString);
                badSeqs = writer.newWriter(badFile);
            }
            //create an interator for each file
            Iterator it = fq.iterator();
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
                    FastqRecord leftSeq = new FastqRecord(seqRecord.getReadHeader() + " 1:N:0:", leftRead, "", leftQual);
                    FastqRecord rightSeq = new FastqRecord(seqRecord.getReadHeader() + " 2:N:0", rightRead, "", rightQual);
                    FastqRecord newLeftSeq = groomRead(leftSeq, format);
                    FastqRecord newRightSeq = groomRead(rightSeq, format);
                    goodLeftSeqs.write(newLeftSeq);
                    goodRightSeqs.write(newRightSeq);
                    itCounter++;
                }

                if (writeBadSeqs == true && (seqLength / 2 != readLength || containsN == true))
                {
                    badSeqs.write(seqRecord);
                }
            }
        }
        System.out.println("Completed writing " + itCounter + " good reads");
    }

    /**
     * Removes reads that contain short/long reads or that contain an 'N', from
     * a fastq file
     *
     * @param fastqFileIn the fastq file
     * @param readsOut the good reads
     * @param readLength the expected length of a single read
     * @param format the fastq format (can only be 'illumina' or 'sanger')
     * @param writeBadSeqs whether to write the bad reads to a file (bad reads
     * file name will start with 'bad_')
     */
    public void qcSingleEndReads(File fastqFileIn, File readsOut, int readLength, String format, boolean writeBadSeqs)
    {
        FastqReader fq = new FastqReader(fastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodSeqs = writer.newWriter(readsOut);

        FastqWriter badSeqs = null;
        int itCounter = 0;

        if (checkFormat(format) == true)
        {
            if (writeBadSeqs == true)
            {
                String badString = "bad_".concat(readsOut.toString());
                File badFile = new File(badString);
                badSeqs = writer.newWriter(badFile);
            }

            //create an interator for the interlaced file
            Iterator it = fq.iterator();

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
        }
        System.out.println("Completed writing " + itCounter + " good reads");
    }

    /**
     *
     * @param record the FastqRecord to be groomed
     * @param format the current format. Can only be 'illumina' or 'sanger'
     * @return a FastqRecord in sanger format
     */
    public FastqRecord groomRead(FastqRecord record, String format)
    {
        FastqRecord newseq = new FastqRecord(null, null, null, null);
        if (format.equalsIgnoreCase("illumina"))
        {
            byte[] sangerQuals = record.getBaseQualityString().getBytes();
            SolexaQualityConverter.getSingleton().convertSolexa_1_3_QualityCharsToPhredBinary(sangerQuals);
            String quality = net.sf.samtools.SAMUtils.phredToFastq(sangerQuals);
            newseq = new FastqRecord(record.getReadHeader(), record.getReadString(), "", quality);
            return newseq;
        }
        else if (format.equalsIgnoreCase("sanger"))
        {
            newseq = record;

        }
        else
        {
            System.err.print("format must be either 'sanger' or 'illumina'");
            System.exit(1);
        }

        return newseq;

    }

    /**
     * Takes a fastq file and calculates the number of reads which have a 'N'
     * and the distribution of read lengths
     *
     * @param fastq the fastq file to analyse
     */
    public void countLengthsAndNs(File fastq)
    {
        FastqReader fq = new FastqReader(fastq);
        int readsWithNs = 0;
        TreeMap<Integer, Integer> lengthDists = new TreeMap<>();

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
            String readString = seqRecord.getReadString();

            boolean containsN = readString.contains("N");

            if (containsN == true)
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

    /**
     * Checks a fastq file to verify that all the reads can be parsed into a
     * FastqRecord. Errors should be thrown if any of the reads are not
     * formatted correctly.
     *
     * @param fastq the fastq file to verify
     */
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

    /**
     * Verifies that two files of left/right paired-end reads are in order and
     * contain no fastq format errors
     *
     * @param fastqLeft the left-handed fastq file to verify
     * @param fastqRight the right-handed fastq file to verify
     */
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
}
