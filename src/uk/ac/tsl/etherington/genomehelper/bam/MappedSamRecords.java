package uk.ac.tsl.etherington.genomehelper.bam;

/*
 * Takes a SAM File, prints the header and then prints any entries where the
 * read or its mate (or both) are mapped
 */
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import uk.ac.tsl.etherington.genomehelper.fastq.FastqParser;

/**
 *
 * @author ethering
 */
public class MappedSamRecords
{

    /**
     * Returns the fastq sequences that have been mapped to a genome. All fastq
     * reads must either be left-handed or right-handed
     *
     * @param bamFile a sam/bam file
     * @param fastqIn the fastq dataset that was mapped to the sam/bam file
     * @param fastqOut the fastq reads that were actually mapped
     * @param isRightHandedReads whether the reads are right-handed
     */
    public void getSingleMappedSamRecords(File bamFile, File fastqIn, File fastqOut, boolean isRightHandedReads)
    {
        HashSet<String> mappedReads = new HashSet<>();
        final SAMFileReader samReader = new SAMFileReader(bamFile);
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter out = writer.newWriter(fastqOut);

        // Open an iterator for the particular sequence
        SAMRecordIterator iterator = samReader.iterator();
        while (iterator.hasNext())
        {
            SAMRecord samRecord = iterator.next();
            //is the read mapped?
            boolean unMapped = samRecord.getReadUnmappedFlag();

            //..if so...
            if (unMapped == false)
            {
                //get the read name
                String currentRead = samRecord.getReadName();
                boolean isFirstOfPair = samRecord.getFirstOfPairFlag();
                //if we have a right-handed mapped sequence, make sure we have right-handed reads
                if (isRightHandedReads && isFirstOfPair == false)
                {
                    mappedReads.add(currentRead);
                } else if (isRightHandedReads == false && isFirstOfPair)
                {
                    mappedReads.add(currentRead);
                }
            }
        }
        System.out.println("Found " + mappedReads.size() + " mapped reads");
        samReader.close();

        final FastqReader fastqReader = new FastqReader(fastqIn);
        FastqParser.writeRecords(fastqReader, out, mappedReads);
    }

    /**
     * Returns the fastq sequences that have not been mapped to a genome. All
     * fastq reads must either be left-handed or right-handed
     *
     * @param bamFile a sam/bam file
     * @param fastqIn the fastq dataset that was mapped to the sam/bam file
     * @param fastqOut the fastq reads that were not mapped
     * @param isRightHandedReads whether the reads are right-handed
     */
    public void getSingleUnmappedSamRecords(File bamFile, File fastqIn, File fastqOut, boolean isRightHandedReads)
    {
        HashSet<String> unMappedReads = new HashSet<>();
        final SAMFileReader samReader = new SAMFileReader(bamFile);
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter out = writer.newWriter(fastqOut);

        // Open an iterator for the particular sequence
        SAMRecordIterator iterator = samReader.iterator();
        while (iterator.hasNext())
        {
            SAMRecord samRecord = iterator.next();
            //is the read mapped?
            boolean unMapped = samRecord.getReadUnmappedFlag();

            //..if so...
            if (unMapped == true)
            {
                //get the read name
                String currentRead = samRecord.getReadName();
                boolean isFirstOfPair = samRecord.getFirstOfPairFlag();
                //if we have a right-handed mapped sequence, make sure we have right-handed reads
                if (isRightHandedReads && isFirstOfPair == false)
                {
                    unMappedReads.add(currentRead);
                } else if (isRightHandedReads == false && isFirstOfPair)
                {
                    unMappedReads.add(currentRead);
                }

            }
        }
        System.out.println("Found " + unMappedReads.size() + " unmapped reads");
        samReader.close();

        final FastqReader fastqReader = new FastqReader(fastqIn);
        FastqParser.writeRecords(fastqReader, out, unMappedReads);
    }

    /**
     * Returns the pairs of sequences where at least one of the pair has been
     * mapped to a genome.
     *
     * @param bamFile a sam/bam file
     * @param fastqInLeft the left-hand fastq dataset that was mapped to the
     * sam/bam file
     * @param fastqInRight the righ-hand fastq dataset that was mapped to the
     * sam/bam file
     * @param fastqOutLeft the left-hand fastq reads where one of the pair
     * actually mapped
     * @param fastqOutRight the right-hand fastq reads where one of the pair
     * actually mapped
     */
    public void getOnePairedMappedSamRecords(File bamFile, File fastqInLeft, File fastqInRight, File fastqOutLeft, File fastqOutRight)
    {

        final SAMFileReader samReader = new SAMFileReader(bamFile);
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter outLeft = writer.newWriter(fastqOutLeft);
        FastqWriter outRight = writer.newWriter(fastqOutRight);

        // Open an iterator for the particular sequence
        SAMRecordIterator iterator = samReader.iterator();
        HashSet<String> mappedReads = getReads(samReader, iterator, true);

        System.out.println("Found " + mappedReads.size() + " mapped reads");
        samReader.close();

        final FastqReader fastqReaderLeft = new FastqReader(fastqInLeft);
        final FastqReader fastqReaderRight = new FastqReader(fastqInRight);
        while (fastqReaderLeft.hasNext())
        {
            FastqRecord leftRecord = fastqReaderLeft.next();
            FastqRecord rightRecord = fastqReaderRight.next();
            String readName = leftRecord.getReadHeader();
            int hashIndex = readName.indexOf(" ");
            readName = readName.substring(0, hashIndex);
            if (mappedReads.contains(readName))
            {
                outLeft.write(leftRecord);
                outRight.write(rightRecord);
            }
        }
    }

    public HashSet<String> getReads(SAMFileReader samReader, SAMRecordIterator iterator, boolean getMapped)
    {
        HashSet<String> reads = new HashSet<>();
        while (iterator.hasNext())
        {
            SAMRecord samRecord = iterator.next();
            //is the read mapped?
            boolean unMapped = samRecord.getReadUnmappedFlag();

            //..if so...do we want mapped reads..
            if (unMapped == false && getMapped)
            {
                //get the read name
                String currentRead = samRecord.getReadName();
                reads.add(currentRead);
            } else if (unMapped && getMapped == false)
            {
                String currentRead = samRecord.getReadName();
                reads.add(currentRead);
            }
        }
        return reads;
    }

    /**
     * Returns the pairs of sequences where at least one of the pair has not
     * been mapped to a genome.
     *
     * @param bamFile a sam/bam file
     * @param fastqInLeft the left-hand fastq dataset that was mapped to the
     * sam/bam file
     * @param fastqInRight the righ-hand fastq dataset that was mapped to the
     * sam/bam file
     * @param fastqOutLeft the left-hand fastq reads where one of the pair was
     * not mapped
     * @param fastqOutRight the right-hand fastq reads where one of the pair was
     * not mapped
     */
    public void getOnePairedUnmappedSamRecords(File bamFile, File fastqInLeft, File fastqInRight, File fastqOutLeft, File fastqOutRight)
    {

        final SAMFileReader samReader = new SAMFileReader(bamFile);
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter outLeft = writer.newWriter(fastqOutLeft);
        FastqWriter outRight = writer.newWriter(fastqOutRight);

        // Open an iterator for the particular sequence
        SAMRecordIterator iterator = samReader.iterator();
        HashSet<String> unMappedReads = getReads(samReader, iterator, false);

        System.out.println("Found " + unMappedReads.size() + " unmapped reads");
        samReader.close();

        final FastqReader fastqReaderLeft = new FastqReader(fastqInLeft);
        final FastqReader fastqReaderRight = new FastqReader(fastqInRight);
        while (fastqReaderLeft.hasNext())
        {
            FastqRecord leftRecord = fastqReaderLeft.next();
            FastqRecord rightRecord = fastqReaderRight.next();
            String readName = leftRecord.getReadHeader();
            int hashIndex = readName.indexOf(" ");
            readName = readName.substring(0, hashIndex);
            if (unMappedReads.contains(readName))
            {
                outLeft.write(leftRecord);
                outRight.write(rightRecord);
            }
        }
    }

    /**
     * Returns the pairs of sequences where at both of the pair has been mapped
     * to a genome.
     *
     * @param bamFile a sam/bam file
     * @param fastqInLeft the left-hand fastq dataset that was mapped to the
     * sam/bam file
     * @param fastqInRight the righ-hand fastq dataset that was mapped to the
     * sam/bam file
     * @param fastqOutLeft the left-hand fastq reads where both of the pairs
     * actually mapped
     * @param fastqOutRight the right-hand fastq reads where both of the pairs
     * actually mapped
     */
    public void getBothPairedMappedSamRecords(File bamFile, File fastqInLeft, File fastqInRight, File fastqOutLeft, File fastqOutRight)
    {
        HashMap<String, PairedMappedRead> mappedReads = new HashMap<>();
        final SAMFileReader samReader = new SAMFileReader(bamFile);
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter outLeft = writer.newWriter(fastqOutLeft);
        FastqWriter outRight = writer.newWriter(fastqOutRight);

        // Open an iterator for the particular sequence
        SAMRecordIterator iterator = samReader.iterator();
        while (iterator.hasNext())
        {
            SAMRecord samRecord = iterator.next();
            //is the read mapped?
            boolean unMapped = samRecord.getReadUnmappedFlag();

            //..if so...
            if (unMapped == false)
            {
                boolean isLeftRead = samRecord.getFirstOfPairFlag();
                //get the read name
                String currentRead = samRecord.getReadName();
                if (mappedReads.containsKey(currentRead))
                {
                    PairedMappedRead pmr = mappedReads.get(currentRead);

                    if (isLeftRead)
                    {
                        pmr.setLeftMapped(true);
                    } else
                    {
                        pmr.setRightMapped(true);
                    }
                    mappedReads.put(currentRead, pmr);
                } else
                {
                    PairedMappedRead pmr = new PairedMappedRead(false, false);
                    if (isLeftRead)
                    {
                        pmr.setLeftMapped(true);
                    } else
                    {
                        pmr.setRightMapped(true);
                    }
                    mappedReads.put(currentRead, pmr);
                }
            }
        }
        System.out.println("Found " + mappedReads.size() + " mapped reads");
        samReader.close();

        final FastqReader fastqReaderLeft = new FastqReader(fastqInLeft);
        final FastqReader fastqReaderRight = new FastqReader(fastqInRight);
        while (fastqReaderLeft.hasNext())
        {
            FastqRecord leftRecord = fastqReaderLeft.next();
            FastqRecord rightRecord = fastqReaderRight.next();
            String readName = leftRecord.getReadHeader();
            int hashIndex = readName.indexOf(" ");
            readName = readName.substring(0, hashIndex);
            if (mappedReads.containsKey(readName))
            {
                PairedMappedRead pmr = mappedReads.get(readName);
                if (pmr.leftMapped && pmr.rightMapped)
                {
                    outLeft.write(leftRecord);
                    outRight.write(rightRecord);
                }

            }
        }
    }

    /**
     * Returns the pairs of sequences where both of the pair has not been mapped
     * to a genome.
     *
     * @param bamFile a sam/bam file
     * @param fastqInLeft the left-hand fastq dataset that was mapped to the
     * sam/bam file
     * @param fastqInRight the righ-hand fastq dataset that was mapped to the
     * sam/bam file
     * @param fastqOutLeft the left-hand fastq reads where both of the pairs was
     * not mapped
     * @param fastqOutRight the right-hand fastq reads where both of the pairs
     * was not mapped
     */
    public void getBothPairedUnmappedSamRecords(File bamFile, File fastqInLeft, File fastqInRight, File fastqOutLeft, File fastqOutRight)
    {
        HashMap<String, PairedMappedRead> mappedReads = new HashMap<>();
        final SAMFileReader samReader = new SAMFileReader(bamFile);
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter outLeft = writer.newWriter(fastqOutLeft);
        FastqWriter outRight = writer.newWriter(fastqOutRight);

        // Open an iterator for the particular sequence
        SAMRecordIterator iterator = samReader.iterator();
        while (iterator.hasNext())
        {
            SAMRecord samRecord = iterator.next();
            //is the read mapped?
            boolean unMapped = samRecord.getReadUnmappedFlag();

            //..if so...
            if (unMapped == true)
            {
                boolean isLeftRead = samRecord.getFirstOfPairFlag();
                //get the read name
                String currentRead = samRecord.getReadName();
                if (mappedReads.containsKey(currentRead))
                {
                    PairedMappedRead pmr = mappedReads.get(currentRead);

                    if (isLeftRead)
                    {
                        pmr.setLeftMapped(false);
                    } else
                    {
                        pmr.setRightMapped(false);
                    }
                    mappedReads.put(currentRead, pmr);
                } else
                {
                    //we'll be looking for false in both later, so set to true
                    PairedMappedRead pmr = new PairedMappedRead(true, true);
                    if (isLeftRead)
                    {
                        pmr.setLeftMapped(false);
                    } else
                    {
                        pmr.setRightMapped(false);
                    }
                    mappedReads.put(currentRead, pmr);
                }
            }
        }
        System.out.println("Found " + mappedReads.size() + " mapped reads");
        samReader.close();

        final FastqReader fastqReaderLeft = new FastqReader(fastqInLeft);
        final FastqReader fastqReaderRight = new FastqReader(fastqInRight);
        while (fastqReaderLeft.hasNext())
        {
            FastqRecord leftRecord = fastqReaderLeft.next();
            FastqRecord rightRecord = fastqReaderRight.next();
            String readName = leftRecord.getReadHeader();
            int hashIndex = readName.indexOf(" ");
            readName = readName.substring(0, hashIndex);
            if (mappedReads.containsKey(readName))
            {
                PairedMappedRead pmr = mappedReads.get(readName);
                if (pmr.leftMapped == false && pmr.rightMapped == false)
                {
                    outLeft.write(leftRecord);
                    outRight.write(rightRecord);
                }

            }
        }
    }

    /**
     * Returns the pairs of sequences where at both of the pair has been mapped
     * to a genome.
     *
     * @param bamFile a sam/bam file
     * @param fastqInLeft the left-hand fastq dataset that was mapped to the
     * sam/bam file
     * @param fastqInRight the righ-hand fastq dataset that was mapped to the
     * sam/bam file
     * @param fastqOutLeft the left-hand fastq reads where both of the pairs
     * were unmapped
     * @param fastqOutRight the right-hand fastq reads where both of the pairs
     * were unmapped
     * @param singles - single unmapped reads where the mate was mapped
     */
    //public void getUnmappedSamRecords(File bamFile, File fastqInLeft, File fastqInRight, File fastqOutLeft, File fastqOutRight, File singles)
    public void getUnmappedSamRecords(File bamFile)
    {
        Set<String> unmappedPairs = new HashSet<>();
        Set<String> unmappedSingles = new HashSet<>();
        int mappedSingles = 0;
        int mappedPairs = 0;
        try (SAMFileReader samReader = new SAMFileReader(bamFile))
        {
            samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            // Open an iterator for the particular sequence
            SAMRecordIterator iterator = samReader.iterator();
            while (iterator.hasNext())
            {
                SAMRecord samRecord = iterator.next();
                //is the read unmapped?
                boolean unmapped = samRecord.getReadUnmappedFlag();
                //if it is...
                if (unmapped)
                {
                    //is it paired
                    boolean isPaired = samRecord.getReadPairedFlag();
                    //if it is
                    if (isPaired)
                    {
                        //is its mate also unmapped
                        boolean mateUnmapped = samRecord.getMateUnmappedFlag();

                        //if so, add it to the unmappedPairs HashSet
                        if (mateUnmapped)
                        {
                            unmappedPairs.add(samRecord.getReadName());
                        }
                    } //if it's not paired
                    else
                    {
                        //set it to true if it's a leftRead, false if it's not
                        unmappedSingles.add(samRecord.getReadName());
                    }
                } //if mapped
                else
                {
                    boolean isPaired = samRecord.getReadPairedFlag();
                    //is its mate also unmapped
                    if (isPaired)
                    {
                        boolean mateUnmapped = samRecord.getMateUnmappedFlag();

                        //if so, add it to the unmappedPairs HashSet
                        if (mateUnmapped == false)
                        {
                            mappedPairs++;
                        }
                    } //if it's not paired find out if it's a left or right read
                    else
                    {
                        mappedSingles++;
                    }
                }
            }
            System.out.println("Found " + unmappedPairs.size() + " unmapped pairs");
            System.out.println("Found " + unmappedSingles.size() + " unmapped singles");
            System.out.println("Found " + mappedPairs + " mapped pairs");
            System.out.println("Found " + mappedSingles + " mapped singles");

//        FastqWriterFactory writer = new FastqWriterFactory();
//        FastqWriter outLeft = writer.newWriter(fastqOutLeft);
//        FastqWriter outRight = writer.newWriter(fastqOutRight);
//        FastqWriter outSingles = writer.newWriter(singles);
//
//        final FastqReader fastqReaderLeft = new FastqReader(fastqInLeft);
//        final FastqReader fastqReaderRight = new FastqReader(fastqInRight);
//        while (fastqReaderLeft.hasNext())
//        {
//            FastqRecord leftRecord = fastqReaderLeft.next();
//            FastqRecord rightRecord = fastqReaderRight.next();
//            String readName = leftRecord.getReadHeader();
//            int hashIndex = readName.indexOf(" ");
//            readName = readName.substring(0, hashIndex);
//            if (unmappedPairs.contains(readName))
//            {
//                outLeft.write(leftRecord);
//                outRight.write(rightRecord);
//            } else if (unmappedSingles.containsKey(readName))
//            {
//                boolean isLeftRead = unmappedSingles.get(readName);
//                if (isLeftRead)
//                {
//                    outSingles.write(leftRecord);
//                } else
//                {
//                    outSingles.write(rightRecord);
//                }
//            }
//        }
//        outLeft.close();
//        outRight.close();
//        outSingles.close();
        }
    }

    private class PairedMappedRead
    {

        private boolean leftMapped;
        private boolean rightMapped;

        public PairedMappedRead(boolean leftMapped, boolean rightMapped)
        {
            this.leftMapped = leftMapped;
            this.rightMapped = rightMapped;
        }

        public boolean isLeftMapped()
        {
            return leftMapped;
        }

        public boolean isRightMapped()
        {
            return rightMapped;
        }

        public void setLeftMapped(boolean leftMapped)
        {
            this.leftMapped = leftMapped;
        }

        public void setRightMapped(boolean rightMapped)
        {
            this.rightMapped = rightMapped;
        }
    }
}
