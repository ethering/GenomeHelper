package uk.ac.tsl.etherington.genomehelper.bam;

/*
 * Takes a SAM File, prints the header and then prints any entries where the
 * read or its mate (or both) are mapped
 */
import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

/**
 *
 * @author ethering
 */
public class MappedSamRecords
{

    /**
     * Returns a HashSet of unmapped reads where either of the paired reads are
     * unmapped
     *
     * @param bamFile a sam/bam file
     * @return a unique list of unmapped reads
     */
    public HashSet listPairedUnmappedReadsFromBam(File bamFile)
    {
        HashSet<String> unmappedReads = new HashSet<>();

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
                    unmappedReads.add(samRecord.getReadName());
                }
            }
        }
        System.out.println("Found "+unmappedReads.size()+ " unmapped reads");
        return unmappedReads;
    }

    /**
     * Returns a HashSet of unmapped reads where either of the paired reads are
     * unmapped
     *
     * @param bamFile a sam/bam file
     * @return a unique list of mapped reads
     */
    public HashSet listPairedMappedReadsFromBam(File bamFile)
    {
        HashSet<String> mappedReads = new HashSet<>();

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
                if (unmapped == false)
                {
                    mappedReads.add(samRecord.getReadName());
                }
            }
        }
        return mappedReads;
    }

    /**
     * Returns a HashSet of mapped reads where only one of the paired reads are
     * mapped
     *
     * @param bamFile a sam/bam file
     * @return a unique list of mapped reads
     */
    public HashMap listSingleUnmappedReadsFromBam(File bamFile)
    {
        HashMap<String, Boolean> unmappedReads = new HashMap<>();

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
                    //if it is paired
                    if (isPaired)
                    {
                        //is its mate mapped
                        boolean mateMapped = samRecord.getMateUnmappedFlag();
                        //if it's unmapped and its mate is mapped
                        if (mateMapped)
                        {
                            //is it a left or right read
                            boolean isLeft = samRecord.getFirstOfPairFlag();
                            {
                                //add it to the list
                                unmappedReads.put(samRecord.getReadName(), isLeft);
                            }
                        }
                    }
                }
            }
        }
        return unmappedReads;
    }

    /**
     * Returns a HashSet of mapped reads where only one of the paired reads are
     * mapped
     *
     * @param bamFile a sam/bam file
     * @return a unique list of mapped reads
     */
    public HashMap listSingleMappedReadsFromBam(File bamFile)
    {
        HashMap<String, Boolean> mappedReads = new HashMap<>();

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
                if (unmapped == false)
                {
                    //is it paired
                    boolean isPaired = samRecord.getReadPairedFlag();
                    //if it is paired
                    if (isPaired)
                    {
                        //is its mate mapped
                        boolean mateMapped = samRecord.getMateUnmappedFlag();
                        //if it's mapped and its mate is unmapped
                        if (mateMapped == false)
                        {
                            //is it a left or right read
                            boolean isLeft = samRecord.getFirstOfPairFlag();
                            {
                                //add it to the list
                                mappedReads.put(samRecord.getReadName(), isLeft);
                            }
                        }

                    }

                }
            }
        }
        return mappedReads;
    }

    /**
     * 
     * @param list a HashSet of paired reads
     * @param fastqInLeft the left-hand fastq reads which contains the reads in list
     * @param fastqInRight the right-hand fastq reads which contains the reads in list
     * @param fastqOutLeft the output left-hand fastq reads
     * @param fastqOutRight the output right-hand fastq reads
     */
    
    public void writePairedReadsFromHashSet(HashSet list, File fastqInLeft, File fastqInRight, File fastqOutLeft, File fastqOutRight)
    {
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter outLeft = writer.newWriter(fastqOutLeft);
        FastqWriter outRight = writer.newWriter(fastqOutRight);

        final FastqReader fastqReaderLeft = new FastqReader(fastqInLeft);
        final FastqReader fastqReaderRight = new FastqReader(fastqInRight);
        
        int noFound = 0;
        while (fastqReaderLeft.hasNext())
        {
            
            FastqRecord leftRecord = fastqReaderLeft.next();
            FastqRecord rightRecord = fastqReaderRight.next();
            String fullReadName = leftRecord.getReadHeader();
            int hashIndex = fullReadName.indexOf(" ");
            String readName = fullReadName.substring(0, hashIndex);

            if (list.contains(readName))
            {
                noFound++;
                outLeft.write(leftRecord);
                outRight.write(rightRecord);
            }
        }
        outLeft.close();
        outRight.close();
        System.out.println("Found "+noFound + " from "+list.size());
    }
    
        /**
     * 
     * @param list a HashSet of paired reads
     * @param fastqInLeft the left-hand fastq reads which contains the reads in list
     * @param fastqInRight the right-hand fastq reads which contains the reads in list
     * @param fastqSinglesOut the output  fastq reads
     */

    public void writeSingleReadsFromHashMap(HashMap<String, Boolean> list, File fastqInLeft, File fastqInRight, File fastqSinglesOut)
    {
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter out = writer.newWriter(fastqSinglesOut);

        final FastqReader fastqReaderLeft = new FastqReader(fastqInLeft);
        final FastqReader fastqReaderRight = new FastqReader(fastqInRight);

        while (fastqReaderLeft.hasNext())
        {
            FastqRecord leftRecord = fastqReaderLeft.next();
            FastqRecord rightRecord = fastqReaderRight.next();
            String fullReadName = leftRecord.getReadHeader();
            int hashIndex = fullReadName.indexOf(" ");
            String readName = fullReadName.substring(0, hashIndex);

            if (list.containsKey(readName))
            {
                boolean isLeft = list.get(readName);
                if (isLeft)
                {
                    String leftReadName = leftRecord.getReadHeader();
                    String newLeftReadName = leftReadName.concat("1:N:0:");
                    FastqRecord newLeftRecord = new FastqRecord(newLeftReadName, leftRecord.getReadString(), "", leftRecord.getBaseQualityString());
                    out.write(newLeftRecord);
                } else
                {
                    String rightReadName = rightRecord.getReadHeader();
                    String newRightReadName = rightReadName.concat("1:N:0:");
                    FastqRecord newLeftRecord = new FastqRecord(newRightReadName, rightRecord.getReadString(), "", rightRecord.getBaseQualityString());
                    out.write(newLeftRecord);
                }
            }
        }
        out.close();
    }

    /**
     * Returns three files, the two paired-end sequence files where both pairs
     * are unmapped and single unmapped reads where only one of a pair are
     * unmapped
     *
     * @param bamFile a sam/bam file
     * @param fastqInLeft the left-hand fastq dataset that was mapped to the
     * sam/bam file
     * @param fastqInRight the righ-hand fastq dataset that was mapped to the
     * sam/bam file
     */
    public void printReadsFromBamAndFastq(File bamFile, File fastqInLeft, File fastqInRight)
    {
        try (SAMFileReader samReader = new SAMFileReader(bamFile))
        {
            samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            // Open an iterator for the particular sequence
            SAMRecordIterator iterator = samReader.iterator();
            System.out.println("Read Headers");
            while (iterator.hasNext())
            {
                SAMRecord samRecord = iterator.next();
                String readHeader = samRecord.getReadName();
                System.out.println(readHeader);
            }
        }
        System.out.println("\t\t");
        final FastqReader fastqReaderLeft = new FastqReader(fastqInLeft);
        final FastqReader fastqReaderRight = new FastqReader(fastqInRight);
        while (fastqReaderLeft.hasNext())
        {
            FastqRecord leftRecord = fastqReaderLeft.next();
            FastqRecord rightRecord = fastqReaderRight.next();
            String leftRead = leftRecord.getReadHeader();
            System.out.println("Left read " + leftRead);
            String rightRead = rightRecord.getReadHeader();
            System.out.println("Fullreadname " + rightRead);
        }

    }

    private class ScaffoldMappingStats
    {

        private String scaffoldName;
        private int mappedPairs;
        private int mappedSingles;
        private int unmappedPairs;
        private int unmappedSingles;

        public ScaffoldMappingStats(String scaffoldName, int mappedPairs, int mappedSingles, int unmappedPairs, int unmappedSingles)
        {
            this.scaffoldName = scaffoldName;
            this.mappedPairs = mappedPairs;
            this.mappedSingles = mappedSingles;
            this.unmappedPairs = unmappedPairs;
            this.unmappedSingles = unmappedSingles;
        }

        public void setScaffoldName(String scaffoldName)
        {
            this.scaffoldName = scaffoldName;
        }

        public void setMappedPairs(int mappedPairs)
        {
            this.mappedPairs = mappedPairs;
        }

        public void setMappedSingles(int mappedSingles)
        {
            this.mappedSingles = mappedSingles;
        }

        public void setUnmappedPairs(int unmappedPairs)
        {
            this.unmappedPairs = unmappedPairs;
        }

        public void setUnmappedSingles(int unmappedSingles)
        {
            this.unmappedSingles = unmappedSingles;
        }

        public String getScaffoldName()
        {
            return scaffoldName;
        }

        public int getMappedPairs()
        {
            return mappedPairs;
        }

        public int getMappedSingles()
        {
            return mappedSingles;
        }

        public int getUnmappedPairs()
        {
            return unmappedPairs;
        }

        public int getUnmappedSingles()
        {
            return unmappedSingles;
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
