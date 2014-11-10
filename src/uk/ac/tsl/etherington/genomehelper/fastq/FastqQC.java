/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.fastq;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import net.sf.picard.util.FastqQualityFormat;
import net.sf.picard.util.QualityEncodingDetector;
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

    /**
     * Checks that a valid format ('illumina' or 'sanger') has been supplied
     *
     * @param format the format (must only be 'illumina' or 'sanger')
     * @return true of the format is valid, otherwise false
     */
    public boolean checkFormat(String format)
    {
        boolean goodFormat = true;
        if (!format.equalsIgnoreCase("illumina") && !format.equalsIgnoreCase("sanger"))
        {
            goodFormat = false;
            System.err.print("format must be either 'sanger' or 'illumina'");

        }
        return goodFormat;
    }

    /**
     * Quality controls a fastq read
     *
     * @param seqRecord a net.sf.picard.fastq.FastqRecord to QC
     * @param singleEndReadLength the expected length of a single end read
     * @return true if the read is the correct length and contains no 'N's,
     * false otherwise
     */
    public boolean checkLengthAndContent(FastqRecord seqRecord, int singleEndReadLength)
    {
        //get the length of them
        int seqLength = seqRecord.getReadString().length();
        String readString = seqRecord.getReadString();
        boolean containsN = readString.toLowerCase().contains("n");

        if (seqLength == singleEndReadLength && containsN == false)
        {
            return true;
        } else
        {
            return false;
        }
    }

    /**
     * Creates a FastqWriter for bad reads that fail Quality Control
     *
     * @param fastqFileOut a FastqWriter to write bad reads to
     * @param writer a FastqWriterFactory. Multiple FastqWriters can be
     * generated from a single FastqWriterFactory
     * @return a FastqWriter to write bad reads to. The bad reads file will
     * start with 'bad_', followed by the outfile name.
     */
    public FastqWriter getBadSeqFastqWriter(File fastqFileOut, FastqWriterFactory writer)
    {
        FastqWriter badSeqWriter = null;
        String parentFilePath = fastqFileOut.getParent();
        String file = fastqFileOut.getName();
        String badStr = "bad_".concat(file);
        File badFile = new File(parentFilePath, badStr);
        badSeqWriter = writer.newWriter(badFile);
        return badSeqWriter;
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
     * @param singles the QC'd reads where the opposite pair has failed QC
     * @param singleEndReadLength the expected length of a single read
     * @param format the fastq format (can only be 'illumina' or 'sanger')
     * @param writeBadSeqs whether to write the bad reads to a file (bad reads
     * file name will start with 'bad_')
     */
    public void qcPairedReads(File leftFastqFileIn, File rightFastqFileIn, File leftReadsOut, File rightReadsOut, File singles, int singleEndReadLength, String format, boolean writeBadSeqs)
    {
        FastqReader fql = new FastqReader(leftFastqFileIn);
        FastqReader fqr = new FastqReader(rightFastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodLeftSeqs = writer.newWriter(leftReadsOut);
        FastqWriter goodRightSeqs = writer.newWriter(rightReadsOut);
        FastqWriter singleSeqs = writer.newWriter(singles);

        FastqWriter badLeftSeqs = null;
        FastqWriter badRightSeqs = null;
        if (writeBadSeqs)
        {
            badLeftSeqs = getBadSeqFastqWriter(leftReadsOut, writer);
            badRightSeqs = getBadSeqFastqWriter(rightReadsOut, writer);
        }

        int itCounter = 0;
        if (checkFormat(format) == true)
        {
            //create an interator for each file
            Iterator itl = fql.iterator();
            Iterator itr = fqr.iterator();

            while (itl.hasNext())
            {
                //get the corresponding reads
                FastqRecord leftSeqRecord = (FastqRecord) itl.next();
                FastqRecord rightSeqRecord = (FastqRecord) itr.next();
                boolean leftGood = checkLengthAndContent(leftSeqRecord, singleEndReadLength);
                boolean rightGood = checkLengthAndContent(rightSeqRecord, singleEndReadLength);

                if (leftGood && rightGood)
                {
                    FastqRecord newLeftSeq = groomRead(leftSeqRecord, format);
                    FastqRecord newRightSeq = groomRead(rightSeqRecord, format);
                    goodLeftSeqs.write(newLeftSeq);
                    goodRightSeqs.write(newRightSeq);
                    itCounter++;
                } else if (leftGood && rightGood == false)
                {
                    FastqRecord seq = groomRead(leftSeqRecord, format);
                    singleSeqs.write(seq);
                } else if (leftGood == false && rightGood)
                {
                    FastqRecord seq = groomRead(rightSeqRecord, format);
                    singleSeqs.write(seq);
                } else if (writeBadSeqs)
                {
                    if (leftGood == false && rightGood == false)
                    {
                        badLeftSeqs.write(leftSeqRecord);
                        badRightSeqs.write(rightSeqRecord);
                    } else if (leftGood && rightGood == false)
                    {
                        badRightSeqs.write(rightSeqRecord);
                    } else if (leftGood == false && rightGood)
                    {
                        badLeftSeqs.write(leftSeqRecord);
                    }

                }
            }
        }
        goodLeftSeqs.close();
        goodRightSeqs.close();
        singleSeqs.close();
        if (writeBadSeqs)
        {
            badLeftSeqs.close();
            badRightSeqs.close();
        }

        System.out.println("Completed writing " + itCounter + " good reads");
    }

    /**
     * Removes pairs of reads where at least one of the pair contains short/long
     * reads or a read that contain an 'N', from an interlaced fastq file
     *
     * @param fastqFileIn an interlaced fastq file
     * @param fastqFileOut the good left-handed reads
     * @param singles the QC'd reads where the opposite pair has failed QC
     * @param singleEndReadLength the expected length of a single read
     * @param format the fastq format (can only be 'illumina' or 'sanger')
     * @param writeBadSeqs whether to write the bad reads to a file (bad reads
     * file name will start with 'bad_')
     */
    public void qcInterlacedReads(File fastqFileIn, File fastqFileOut, File singles, int singleEndReadLength, String format, boolean writeBadSeqs)
    {
        FastqReader fq = new FastqReader(fastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodSeqs = writer.newWriter(fastqFileOut);
        FastqWriter singleSeqs = writer.newWriter(singles);
        FastqWriter badSeqs = null;
        if (writeBadSeqs)
        {
            badSeqs = getBadSeqFastqWriter(fastqFileOut, writer);
        }
        int itCounter = 0;

        if (checkFormat(format) == true)
        {
            //create an interator for the interlaced file
            Iterator it = fq.iterator();
            while (it.hasNext())
            {
                //get the corresponding reads
                FastqRecord leftSeqRecord = (FastqRecord) it.next();
                FastqRecord rightSeqRecord = (FastqRecord) it.next();

                boolean leftGood = checkLengthAndContent(leftSeqRecord, singleEndReadLength);
                boolean rightGood = checkLengthAndContent(rightSeqRecord, singleEndReadLength);

                if (leftGood && rightGood)
                {
                    FastqRecord newLeftSeq = groomRead(leftSeqRecord, format);
                    FastqRecord newRightSeq = groomRead(rightSeqRecord, format);
                    goodSeqs.write(newLeftSeq);
                    goodSeqs.write(newRightSeq);
                    itCounter++;
                } else if (leftGood && rightGood == false)
                {
                    FastqRecord seq = groomRead(leftSeqRecord, format);
                    singleSeqs.write(seq);
                } else if (leftGood == false && rightGood)
                {
                    FastqRecord seq = groomRead(rightSeqRecord, format);
                    singleSeqs.write(seq);
                } 
                
                else if (writeBadSeqs)
                {
                    if (leftGood == false && rightGood == false)
                    {
                        badSeqs.write(leftSeqRecord);
                        badSeqs.write(rightSeqRecord);
                    } else if (leftGood && rightGood == false)
                    {
                        badSeqs.write(rightSeqRecord);
                    } else if (leftGood == false && rightGood)
                    {
                        badSeqs.write(leftSeqRecord);
                    }

                }
            }
        }
        goodSeqs.close();
        singleSeqs.close();
        if (writeBadSeqs)
        {
            badSeqs.close();
        }
        
        System.out.println("Completed writing " + itCounter + " good reads");
    }

    /**
     * Removes pairs of reads where at least one of the pair contains short/long
     * reads or a read that contain an 'N', from an interlaced fastq file and writes
     * the resulting good reads to two separate paired-end files
     *
     * @param fastqFileIn an interlaced fastq file
     * @param leftReadsOut the good left-handed reads
     * @param rightReadsOut the good right-handed reads
     * @param singles any single reads
     * @param singleEndReadLength the expected length of a single read
     * @param format the fastq format (can only be 'illumina' or 'sanger')
     * @param writeBadSeqs whether to write the bad reads to a file (bad reads
     * file name will start with 'bad_')
     */
    public void qcInterlacedReadsToPairs(File fastqFileIn, File leftReadsOut, File rightReadsOut, File singles, int singleEndReadLength, String format, boolean writeBadSeqs)
    {
        FastqReader fq = new FastqReader(fastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodLeftSeqs = writer.newWriter(leftReadsOut);
        FastqWriter goodRightSeqs = writer.newWriter(rightReadsOut);
        FastqWriter singleSeqs = writer.newWriter(singles);
        FastqWriter badLeftSeqs = getBadSeqFastqWriter(leftReadsOut, writer);
        FastqWriter badRightSeqs = getBadSeqFastqWriter(rightReadsOut, writer);
        int itCounter = 0;

        if (checkFormat(format) == true)
        {
            //create an interator for the interlaced file
            Iterator it = fq.iterator();
            while (it.hasNext())
            {
                //get the corresponding reads
                FastqRecord leftSeqRecord = (FastqRecord) it.next();
                FastqRecord rightSeqRecord = (FastqRecord) it.next();

                boolean leftGood = checkLengthAndContent(leftSeqRecord, singleEndReadLength);
                boolean rightGood = checkLengthAndContent(rightSeqRecord, singleEndReadLength);

                if (leftGood && rightGood)
                {
                    FastqRecord newLeftSeq = groomRead(leftSeqRecord, format);
                    FastqRecord newRightSeq = groomRead(rightSeqRecord, format);
                    goodLeftSeqs.write(newLeftSeq);
                    goodRightSeqs.write(newRightSeq);
                    itCounter++;
                } else if (leftGood && rightGood == false)
                {
                    FastqRecord seq = groomRead(leftSeqRecord, format);
                    singleSeqs.write(seq);
                } else if (leftGood == false && rightGood)
                {
                    FastqRecord seq = groomRead(rightSeqRecord, format);
                    singleSeqs.write(seq);
                } else if (writeBadSeqs)
                {
                    if (leftGood == false && rightGood == false)
                    {
                        badLeftSeqs.write(leftSeqRecord);
                        badRightSeqs.write(rightSeqRecord);
                    } else if (leftGood && rightGood == false)
                    {
                        badRightSeqs.write(rightSeqRecord);
                    } else if (leftGood == false && rightGood)
                    {
                        badLeftSeqs.write(leftSeqRecord);
                    }

                }
            }
        }
        goodLeftSeqs.close();
        goodRightSeqs.close();
        singleSeqs.close();
        if (writeBadSeqs)
        {
            badLeftSeqs.close();
            badRightSeqs.close();
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

        FastqWriter badSeqs = getBadSeqFastqWriter(leftReadsOut, writer);

        int itCounter = 0;

        if (checkFormat(format) == true)
        {
            //create an interator for each file
            
            for (FastqRecord seqRecord : fq) {
                int seqLength = seqRecord.getReadString().length();
                String readString = seqRecord.getReadString();
                boolean containsN = readString.toLowerCase().contains("n");

                if (seqLength / 2 == readLength && containsN == false)
                {
                    String qual = seqRecord.getBaseQualityString();
                    String leftRead = readString.substring(0, (seqLength / 2));
                    String rightRead = readString.substring(seqLength / 2, seqLength);
                    String leftQual = qual.substring(0, (seqLength / 2));
                    String rightQual = qual.substring(seqLength / 2, seqLength);
                    FastqRecord leftSeq = new FastqRecord(seqRecord.getReadHeader() + " 1:N:0:", leftRead, "", leftQual);
                    FastqRecord rightSeq = new FastqRecord(seqRecord.getReadHeader() + " 2:N:0:", rightRead, "", rightQual);
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
        goodLeftSeqs.close();
        goodRightSeqs.close();
        badSeqs.close();
        System.out.println("Completed writing " + itCounter + " good reads");
    }

    /**
     * Removes reads that contain short/long reads or that contain an 'N', from
     * a fastq file
     *
     * @param fastqFileIn the fastq file
     * @param readsOut the good reads
     * @param singleEndReadLength the expected length of a single read
     * @param format the fastq format (can only be 'illumina' or 'sanger')
     * @param writeBadSeqs whether to write the bad reads to a file (bad reads
     * file name will start with 'bad_')
     */
    public void qcSingleEndReads(File fastqFileIn, File readsOut, int singleEndReadLength, String format, boolean writeBadSeqs)
    {
        FastqReader fq = new FastqReader(fastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodSeqs = writer.newWriter(readsOut);

        FastqWriter badSeqs = getBadSeqFastqWriter(readsOut, writer);
        int itCounter = 0;

        if (checkFormat(format) == true)
        {
            //create an interator for the interlaced file
            
            for (FastqRecord seqRecord : fq) {
                boolean goodRead = checkLengthAndContent(seqRecord, singleEndReadLength);

                if (goodRead)
                {
                    FastqRecord newSeq = groomRead(seqRecord, format);
                    goodSeqs.write(newSeq);
                    itCounter++;
                } else if (writeBadSeqs == true)
                {
                    badSeqs.write(seqRecord);
                }
            }
        }
        goodSeqs.close();
        badSeqs.close();
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
        } else if (format.equalsIgnoreCase("sanger"))
        {
            newseq = record;

        } else
        {
            System.err.print("format must be either 'sanger' or 'illumina'");
            System.exit(1);
        }
        return newseq;
    }

    /**
     * Guesses the format of a FastqRecord
     *
     * @param record a FastqRecord of unknown format
     * @return the best guess;
     */
    public FastqQualityFormat guessFormat(FastqRecord record)
    {
        QualityEncodingDetector detector = new QualityEncodingDetector();
        detector.add(record);
        FastqQualityFormat qual = detector.generateBestGuess(QualityEncodingDetector.FileContext.FASTQ);
        System.out.println("Qualaty encoding: " + qual.toString());
        return qual;
    }

    /**
     * Counts the number of reads and their cumulative nucleotide content.
     * Useful for when a fastq file has varying read lengths, e.g. after
     * trimming. Prints the number of reads and their cumulative nucleotide
     * content.
     *
     * @param fastq the fastq-formatted sequence file to analyse
     * @return nucleotides the combined number of nucleotides in a fastq file
     */
    public double getNucleotideCount(File fastq)
    {
        FastqReader fq = new FastqReader(fastq);
        double reads = 0;
        double nucleotides = 0;

        for (FastqRecord seqRecord : fq)
        {
            reads++;
            int seqLength = seqRecord.getReadString().length();
            nucleotides += seqLength;
        }
        NumberFormat formatter = new DecimalFormat("###.#####");

        String readsString = formatter.format(reads);
        String ntString = formatter.format(nucleotides);
        System.out.println("No. of reads\tNucleotide count");
        System.out.println(readsString +"\t" + ntString);
        return nucleotides;
    }

    /**
     * Takes a fastq file and calculates the number of reads which have a 'N'
     * and the distribution of read lengths
     *
     * @param fastq the fastq file to analyse
     * @return the number of reads found with 'N' in the read sequence
     */
    public int countLengthsAndNs(File fastq)
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
            } else
            {
                lengthDists.put(seqLength, 1);
            }

            //check for N's
            String readString = seqRecord.getReadString();

            boolean containsN = readString.toLowerCase().contains("n");

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
        return readsWithNs;

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

    /**
     * Removes pairs of reads where at least one of the pair contains short/long
     * reads or a read that contain an 'N', from paired-end fastq files.
     * Optionally writes bad reads to a file.
     *
     * @param leftFastqFileIn the left-handed reads
     * @param rightFastqFileIn the right-handed reads
     * @param leftReadsOut the QCd left-handed reads
     * @param rightReadsOut the QCd right-handed reads
     * @param kmers - an ArrayList of kmers to search for
     */
    public void removePairedReadsWithKmers(File leftFastqFileIn, File rightFastqFileIn, File leftReadsOut, File rightReadsOut, ArrayList<String> kmers)
    {
        FastqReader fql = new FastqReader(leftFastqFileIn);
        FastqReader fqr = new FastqReader(rightFastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodLeftSeqs = writer.newWriter(leftReadsOut);
        FastqWriter goodRightSeqs = writer.newWriter(rightReadsOut);
        System.out.println("kmers: " + kmers.toString());
        int reads = 0;
        int goodReads = 0;

        //create an interator for each file
        Iterator itl = fql.iterator();
        Iterator itr = fqr.iterator();
        FastqRecord leftSeqRecord = null;
        FastqRecord rightSeqRecord = null;

        while (itl.hasNext())
        {
            reads++;
            leftSeqRecord = (FastqRecord) itl.next();
            rightSeqRecord = (FastqRecord) itr.next();
            boolean leftGood = findKmers(kmers, leftSeqRecord);
            boolean rightGood = findKmers(kmers, rightSeqRecord);

            if (leftGood && rightGood)
            {
                goodLeftSeqs.write(leftSeqRecord);
                goodRightSeqs.write(rightSeqRecord);
                goodReads++;
            }
        }

        goodLeftSeqs.close();
        goodRightSeqs.close();
        int removedReads = reads - goodReads;
        System.out.println("Wrote " + goodReads + " good reads from " + reads);
        System.out.println("Removed " + removedReads + " reads");
    }

    public void removeSingleReadsWithKmers(File fastqFileIn, File readsOut, ArrayList<String> kmers)
    {
        FastqReader fql = new FastqReader(fastqFileIn);
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodSeqs = writer.newWriter(readsOut);

        int reads = 0;
        int goodReads = 0;
        FastqRecord seqRecord = null;
        String seq = "";
        //create an interator for each file
        Iterator itl = fql.iterator();
        System.out.println("kmers: " + kmers.toString());
        while (itl.hasNext())
        {
            reads++;
            //get the corresponding reads
            seqRecord = (FastqRecord) itl.next();
            seq = seqRecord.getReadString();
            boolean leftGood = findKmers(kmers, seqRecord);

            if (leftGood)
            {
                goodSeqs.write(seqRecord);
                goodReads++;
            }
        }

        goodSeqs.close();
        int removedReads = reads - goodReads;
        System.out.println("Wrote " + goodReads + " good reads from " + reads);
        System.out.println("Removed " + removedReads + " reads");
    }

    public boolean findKmers(ArrayList<String> kmers, FastqRecord rec)
    {
        boolean goodSeq = true;
        String seq = rec.getReadString();
        for (String kmer : kmers)
        {
            if (seq.toLowerCase().contains(kmer.toLowerCase()))
            {
                goodSeq = false;
                break;

            }
        }
        return goodSeq;
    }

//    public static void main(String[] args)
//    {
//        FastqRecord newseq = new FastqRecord("@HWI-EAS396_0001:5:1:10468:1298#0/1", "CTTTTAGCAAGATATCTTATCCATTCCATCTTCGATCCACACAATTGAATCATGTAATTCTCCAATGTAACGCAAT",
//                "+HWI-EAS396_0001:5:1:10468:1298#0/1", "ddc_cfcccfa[ddab\\_a`cfffdffS_ffc^fYddcWe]`]X^bcbadcffccW^ae[ffffffcdffdfaWcc");
//        FastqQC instance = new FastqQC();
//        instance.guessFormat(newseq);
//
//    }
}
