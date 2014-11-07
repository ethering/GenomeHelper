package uk.ac.tsl.etherington.genomehelper.bam;

/*
 * Takes a SAM File, prints the header and then prints any entries where the
 * read or its mate (or both) are mapped
 */
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
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
     * Returns three files, the two paired-end sequence files where both pairs
     * are unmapped and single unmapped reads where only one of a pair are
     * unmapped
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
    public void getUnmappedReadsFromPairedReads(File bamFile, File fastqInLeft, File fastqInRight, File fastqOutLeft, File fastqOutRight, File singles)
    //public void getUnmappedSamRecords(File bamFile, File fastqInLeft, File fastqInRight)
    {
        HashMap<String, ScaffoldMappingStats> mappingStats = new HashMap<>();
        Set<String> unmappedPairs = new HashSet<>();
        HashMap<String, Boolean> unmappedSingles = new HashMap<>();

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter outLeft = writer.newWriter(fastqOutLeft);
        FastqWriter outRight = writer.newWriter(fastqOutRight);
        FastqWriter outSingles = writer.newWriter(singles);
        FastqRecord fq = new FastqRecord(null, null, null, null);
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
                String readHeader = samRecord.getReadName();
                String leftHeader = readHeader.concat(" 1:N:0:");
                String rightHeader = readHeader.concat(" 2:N:0:");
                String scaffold = samRecord.getReferenceName();
                if (!mappingStats.containsKey(scaffold))
                {
                    ScaffoldMappingStats stat = new ScaffoldMappingStats(scaffold, 0, 0, 0, 0);
                    mappingStats.put(scaffold, stat);
                }
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
                        //is its mate also unmapped
                        boolean mateUnmapped = samRecord.getMateUnmappedFlag();
                        //is it a left-hand read (we want to count pairs, so just count the left-hand reads)
                        boolean isLeftRead = samRecord.getFirstOfPairFlag();
                        //if so, add it to the unmappedPairs HashSet
                        if (mateUnmapped && isLeftRead)
                        {
                            fq = new FastqRecord(leftHeader, samRecord.getReadString(), "", samRecord.getBaseQualityString());
                            outLeft.write(fq);
                            unmappedPairs.add(samRecord.getReadName());
                            int unmappedPairsCount = mappingStats.get(scaffold).getUnmappedPairs();
                            unmappedPairsCount++;
                            mappingStats.get(scaffold).setUnmappedPairs(unmappedPairsCount);
                            System.out.println(samRecord.getReadName() + " is unmapped pair");
                        } else if (mateUnmapped && isLeftRead == false)
                        {
                            fq = new FastqRecord(rightHeader, samRecord.getReadString(), "", samRecord.getBaseQualityString());
                            outRight.write(fq);
                        } //if the mate is mapped
                        else if (mateUnmapped == false)
                        {

                            int unmappedSinglesCount = mappingStats.get(scaffold).getUnmappedSingles();
                            unmappedSinglesCount++;
                            mappingStats.get(scaffold).setUnmappedSingles(unmappedSinglesCount);
                            System.out.print(samRecord.getReadName() + " is unmapped single");
                            boolean isLeft = samRecord.getFirstOfPairFlag();

                            if (isLeft)
                            {
                                System.out.println(" left-hand read");
                                fq = new FastqRecord(leftHeader, samRecord.getReadString(), "", samRecord.getBaseQualityString());
                                outSingles.write(fq);
                            } else
                            {
                                fq = new FastqRecord(rightHeader, samRecord.getReadString(), "", samRecord.getBaseQualityString());
                                outSingles.write(fq);
                                System.out.println(" right-hand read");
                            }
                            //add the current read into the unmapped singles
                            unmappedSingles.put(samRecord.getReadName(), isLeft);
                        }
                    } //if it's not paired

                } //if mapped
                else if (unmapped == false)
                {
                    //is it paired
                    boolean isPaired = samRecord.getReadPairedFlag();
                    //if it is paired
                    if (isPaired)
                    {
                        //is its mate also mapped
                        boolean mateUnmapped = samRecord.getMateUnmappedFlag();
                        //is it a left-hand read (we want to count pairs, so just count the left-hand reads)
                        boolean isLeftRead = samRecord.getFirstOfPairFlag();
                        //if the mate is also mapped and it's a left-hand read
                        if (mateUnmapped == false && isLeftRead)
                        {
                            //count it as a mapped pair
                            mappedPairs++;
                            int mappedPairsCount = mappingStats.get(scaffold).getMappedPairs();
                            mappedPairsCount++;
                            mappingStats.get(scaffold).setMappedPairs(mappedPairsCount);
                            //System.out.println(samRecord.getReadName() + " is mapped pair");
                        } //if it's mate is unmapped, add it to the singles
                        else if (mateUnmapped)
                        {
                            mappedSingles++;
                            int mappedSinglesCount = mappingStats.get(scaffold).getMappedSingles();
                            mappedSinglesCount++;
                            mappingStats.get(scaffold).setMappedSingles(mappedSinglesCount);
                            boolean isLeft = samRecord.getFirstOfPairFlag();
                            System.out.print(samRecord.getReadName() + " is mapped single");
                            if (isLeft)
                            {
                                System.out.println(" left-hand read");
                            } else
                            {
                                System.out.println(" right-hand read");
                            }
                        }
                    }

                }
            }
            System.out.println("Found " + mappedPairs + " mapped pairs");
            System.out.println("Found " + mappedSingles + " mapped singles");
            System.out.println("Found " + unmappedPairs.size() + " unmapped pairs");
            System.out.println("Found " + unmappedSingles.size() + " unmapped singles");

            Iterator it = mappingStats.keySet().iterator();

            while (it.hasNext())
            {
                String key = it.next().toString();
                ScaffoldMappingStats stats = mappingStats.get(key);

                System.out.println("Scaffold: " + stats.scaffoldName + "\t" + stats.getMappedPairs() + "\t" + stats.getMappedSingles() + "\t" + stats.getUnmappedPairs() + "\t" + stats.getUnmappedSingles());

            }

            final FastqReader fastqReaderLeft = new FastqReader(fastqInLeft);
            final FastqReader fastqReaderRight = new FastqReader(fastqInRight);
            while (fastqReaderLeft.hasNext())
            {
                FastqRecord leftRecord = fastqReaderLeft.next();
                FastqRecord rightRecord = fastqReaderRight.next();
                String fullReadName = leftRecord.getReadHeader();
                System.out.println("Fullreadname " + fullReadName);
                int hashIndex = fullReadName.indexOf(" ");
                String readName = fullReadName.substring(0, hashIndex);
                System.out.println("Shortreadname " + readName);
                if (unmappedPairs.contains(readName))
                {
                    System.out.println("Paired " + leftRecord.getReadHeader() + "/1");
                    System.out.println("Paired " + rightRecord.getReadHeader() + "/2");
                    //outLeft.write(leftRecord);
                    //outRight.write(rightRecord);
                } else if (unmappedSingles.containsKey(fullReadName))
                {
                    boolean isLeft = unmappedSingles.get(fullReadName);
                    if (isLeft)
                    {
                        System.out.println("Single " + leftRecord.getReadHeader() + "/1");
                    } else
                    {
                        System.out.println("Single " + leftRecord.getReadHeader() + "/2");
                    }
                    //NEED TO FIND OUT WHICH LEFT OF RIGHT TO WRITE
                    //outSingles.write(rightRecord);
                }
            }
            outLeft.close();
            outRight.close();
            outSingles.close();
        }
    }

    /**
     * Returns three files, the two paired-end sequence files where both or
     * either pairs are unmapped and single unmapped reads where only single-end
     * reads are available
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
    public void getAllPairedUnmappedReadsFromPairedReads(File bamFile, File fastqInLeft, File fastqInRight, File fastqOutLeft, File fastqOutRight, File singles)
    //public void getUnmappedSamRecords(File bamFile, File fastqInLeft, File fastqInRight)
    {
        HashMap<String, ScaffoldMappingStats> mappingStats = new HashMap<>();
        Set<String> unmappedPairs = new HashSet<>();
        HashMap<String, Boolean> unmappedSingles = new HashMap<>();

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter outLeft = writer.newWriter(fastqOutLeft);
        FastqWriter outRight = writer.newWriter(fastqOutRight);
        FastqWriter outSingles = writer.newWriter(singles);
        FastqRecord fq = new FastqRecord(null, null, null, null);
        try (SAMFileReader samReader = new SAMFileReader(bamFile))
        {
            samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            // Open an iterator for the particular sequence
            SAMRecordIterator iterator = samReader.iterator();
            while (iterator.hasNext())
            {
                SAMRecord samRecord = iterator.next();
                String readHeader = samRecord.getReadName();
                String leftHeader = readHeader.concat(" 1:N:0:");
                String rightHeader = readHeader.concat(" 2:N:0:");
                String scaffold = samRecord.getReferenceName();
                if (!mappingStats.containsKey(scaffold))
                {
                    ScaffoldMappingStats stat = new ScaffoldMappingStats(scaffold, 0, 0, 0, 0);
                    mappingStats.put(scaffold, stat);
                }
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
                        //is its mate also unmapped
                        //boolean mateUnmapped = samRecord.getMateUnmappedFlag();
                        //is it a left-hand read (we want to count pairs, so just count the left-hand reads)
                        boolean isLeftRead = samRecord.getFirstOfPairFlag();
                        //if so, add it to the unmappedPairs HashSet
                        if (isLeftRead)
                        {
                            fq = new FastqRecord(leftHeader, samRecord.getReadString(), "", samRecord.getBaseQualityString());
                            outLeft.write(fq);
                            unmappedPairs.add(samRecord.getReadName());

                        } else if (isLeftRead == false)
                        {
                            fq = new FastqRecord(rightHeader, samRecord.getReadString(), "", samRecord.getBaseQualityString());
                            outRight.write(fq);
                        }

                    } else//it's not paired
                    {
                        fq = new FastqRecord(samRecord.getReadName(), samRecord.getReadString(), "", samRecord.getBaseQualityString());
                        outSingles.write(fq);
                    }

                }

            }


            final FastqReader fastqReaderLeft = new FastqReader(fastqInLeft);
            final FastqReader fastqReaderRight = new FastqReader(fastqInRight);
            while (fastqReaderLeft.hasNext())
            {
                FastqRecord leftRecord = fastqReaderLeft.next();
                FastqRecord rightRecord = fastqReaderRight.next();
                String fullReadName = leftRecord.getReadHeader();
                System.out.println("Fullreadname " + fullReadName);
                int hashIndex = fullReadName.indexOf(" ");
                String readName = fullReadName.substring(0, hashIndex);
                System.out.println("Shortreadname " + readName);
                if (unmappedPairs.contains(readName))
                {
                    System.out.println("Paired " + leftRecord.getReadHeader() + "/1");
                    System.out.println("Paired " + rightRecord.getReadHeader() + "/2");
                    //outLeft.write(leftRecord);
                    //outRight.write(rightRecord);
                } else if (unmappedSingles.containsKey(fullReadName))
                {
                    boolean isLeft = unmappedSingles.get(fullReadName);
                    if (isLeft)
                    {
                        System.out.println("Single " + leftRecord.getReadHeader() + "/1");
                    } else
                    {
                        System.out.println("Single " + leftRecord.getReadHeader() + "/2");
                    }
                    //NEED TO FIND OUT WHICH LEFT OF RIGHT TO WRITE
                    //outSingles.write(rightRecord);
                }
            }
            outLeft.close();
            outRight.close();
            outSingles.close();
        }
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
     * @param fastqOutLeft the left-hand fastq reads where both of the pairs
     * were unmapped
     * @param fastqOutRight the right-hand fastq reads where both of the pairs
     * were unmapped
     * @param singles - single unmapped reads where the mate was mapped
     */
    public void getMappedReadsFromPairedReads(File bamFile, File fastqInLeft, File fastqInRight, File fastqOutLeft, File fastqOutRight, File singles)
    {
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter outLeft = writer.newWriter(fastqOutLeft);
        FastqWriter outRight = writer.newWriter(fastqOutRight);
        FastqWriter outSingles = writer.newWriter(singles);
        FastqRecord fq = new FastqRecord(null, null, null, null);

        try (SAMFileReader samReader = new SAMFileReader(bamFile))
        {
            samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            // Open an iterator for the particular sequence
            SAMRecordIterator iterator = samReader.iterator();
            while (iterator.hasNext())
            {
                SAMRecord samRecord = iterator.next();
                String readHeader = samRecord.getReadName();
                String leftHeader = readHeader.concat(" 1:N:0:");
                String rightHeader = readHeader.concat(" 2:N:0:");
                
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
                       
                        //is it a left-hand read (we want to count pairs, so just count the left-hand reads)
                        boolean isLeftRead = samRecord.getFirstOfPairFlag();
                        //if so, add it to the unmappedPairs HashSet
                        if (isLeftRead)
                        {
                            fq = new FastqRecord(leftHeader, samRecord.getReadString(), "", samRecord.getBaseQualityString());
                            outLeft.write(fq);
                           
                        } else if (isLeftRead == false)
                        {
                            fq = new FastqRecord(rightHeader, samRecord.getReadString(), "", samRecord.getBaseQualityString());
                            outRight.write(fq);
                        } //if the mate is mapped
                        
                    } else//if it's not paired
                    {
                        fq = new FastqRecord(samRecord.getReadName(), samRecord.getReadString(), "", samRecord.getBaseQualityString());
                        outSingles.write(fq);
                    }

                } 
            }

            outLeft.close();
            outRight.close();
            outSingles.close();
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
