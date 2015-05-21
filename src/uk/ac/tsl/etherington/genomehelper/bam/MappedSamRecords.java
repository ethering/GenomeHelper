package uk.ac.tsl.etherington.genomehelper.bam;

/*
 * Takes a SAM File, prints the header and then prints any entries where the
 * read or its mate (or both) are mapped
 */
import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import htsjdk.samtools.fastq.*;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

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
    public HashSet listEitherPairedReadUnmappedFromBam(File bamFile)
    {
        HashSet<String> unmappedReads = new HashSet<>();
        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        //USE SamReaderFactory
        for (final SAMRecord samRecord : reader)
        {
            //is the read unmapped?
            boolean unmapped = samRecord.getReadUnmappedFlag();
            //if it is...
            if (unmapped)
            {
                unmappedReads.add(samRecord.getReadName());
            }
        }

        System.out.println("Found " + unmappedReads.size() + " unmapped reads");
        return unmappedReads;
    }

    /**
     * Returns a HashSet of unmapped reads where both of the paired reads are
     * unmapped
     *
     * @param bamFile a sam/bam file
     * @return a unique list of unmapped reads where both pairs are unmapped
     */
    public HashSet listBothPairedReadsUnmappedFromBam(File bamFile)
    {
        HashSet<String> unmappedReads = new HashSet<>();

        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        //USE SamReaderFactory
        for (final SAMRecord samRecord : reader)
        {
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
                    boolean mateUnmapped = samRecord.getMateUnmappedFlag();
                    //if both it and it's mate are unmapped
                    if (mateUnmapped)
                    {
                        //add it to the list
                        unmappedReads.add(samRecord.getReadName());
                    }
                }
            }

        }
        System.out.println("Found " + unmappedReads.size() + " unmapped reads");
        return unmappedReads;
    }

    /**
     * Returns a HashSet of reads where both reads are not properly paired
     *
     * @param bamFile a sam/bam file
     * @return a unique list of unmapped reads where both pairs are unmapped
     */
    public HashSet listNonProperlyPairedReadsUnmappedFromBam(File bamFile)
    {
        HashSet<String> unmappedReads = new HashSet<>();

        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        //USE SamReaderFactory
        for (final SAMRecord samRecord : reader)
        {
            //is it paired
            boolean isPaired = samRecord.getReadPairedFlag();
            //if it is paired
            if (isPaired)
            {
                //is the read unmapped?
                boolean properlyPaired = samRecord.getProperPairFlag();
                //if it is...
                if (properlyPaired == false)
                {
                    //add it to the list
                    unmappedReads.add(samRecord.getReadName());
                }
            }
        }

        System.out.println("Found " + unmappedReads.size() + " non-properly paired reads");
        return unmappedReads;
    }

    /**
     * Returns a HashSet of mapped reads where either of the paired reads are
     * mapped
     *
     * @param bamFile a sam/bam file
     * @return a unique list of mapped reads
     */
    public HashSet listEitherPairedReadMappedFromBam(File bamFile)
    {
        HashSet<String> mappedReads = new HashSet<>();
        int unmappedReads = 0;
        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        //USE SamReaderFactory
        for (final SAMRecord samRecord : reader)
        {
            //is the read unmapped?
            boolean unmapped = samRecord.getReadUnmappedFlag();
            //if it is...
            if (unmapped == false)
            {
                mappedReads.add(samRecord.getReadName());
                unmappedReads++;
            }
        }

        System.out.println("Found " + unmappedReads + " mapped reads from " + mappedReads.size() + " different pairs");
        return mappedReads;
    }

    /**
     * Returns a HashSet of mapped reads where both of the paired reads are
     * mapped
     *
     * @param bamFile a sam/bam file
     * @return a unique list of mapped reads where both pairs are mapped
     */
    public HashSet listBothPairedReadsMappedFromBam(File bamFile)
    {
        HashSet<String> mappedReads = new HashSet<>();

        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        //USE SamReaderFactory
        for (final SAMRecord samRecord : reader)
        {
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
                    boolean mateUnmapped = samRecord.getMateUnmappedFlag();
                    //if both it and it's mate are unmapped
                    if (mateUnmapped == false)
                    {
                        //add it to the list
                        mappedReads.add(samRecord.getReadName());
                    }
                }
            }
        }

        System.out.println("Found " + mappedReads.size() + " unmapped reads");
        return mappedReads;
    }

    /**
     * Returns a HashMap of unmapped reads where only one of the paired reads
     * are unmapped
     *
     * @param bamFile a sam/bam file
     * @return a unique list of mapped reads
     */
    public HashMap listSinglePairedReadUnmappedFromBam(File bamFile)
    {
        HashMap<String, Boolean> unmappedReads = new HashMap<>();

        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        //USE SamReaderFactory
        for (final SAMRecord samRecord : reader)
        {
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
                    boolean mateUnmapped = samRecord.getMateUnmappedFlag();
                    //if it's unmapped and its mate is mapped
                    if (mateUnmapped == false)
                    {
                        //is it a left or right read
                        boolean isLeft = samRecord.getFirstOfPairFlag();
                        //add it to the list
                        unmappedReads.put(samRecord.getReadName(), isLeft);
                        //System.out.println(samRecord.getReadName() + " " + isLeft);
                    }
                }
            }
        }



        System.out.println("Found " + unmappedReads.size() + " unmapped reads");
        return unmappedReads;
    }

    /**
     * Returns a HashSet of ummapped unpaired reads
     *
     * @param bamFile a sam/bam file
     * @return a unique list of mapped reads
     */
    public HashSet listUnmappedSingleEndReadsFromBam(File bamFile)
    {
        HashSet<String> unmappedReads = new HashSet<>();

        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        //USE SamReaderFactory
        for (final SAMRecord samRecord : reader)
        {
                //is the read unmapped?
                boolean unmapped = samRecord.getReadUnmappedFlag();
                //if it is...
                if (unmapped)
                {
                    //is it paired
                    boolean isPaired = samRecord.getReadPairedFlag();
                    //if it is paired
                    if (isPaired == false)
                    {
                        //add it to the list
                        unmappedReads.add(samRecord.getReadName());
                    }
                }
            }
        
        System.out.println("Found " + unmappedReads.size() + " mapped reads");
        return unmappedReads;
    }

    /**
     * Returns a HashSet of mapped reads where only one of the paired reads are
     * mapped
     *
     * @param bamFile a sam/bam file
     * @return a unique list of mapped reads
     */
    public HashSet listMappedSingleEndReadsFromBam(File bamFile)
    {
        HashSet<String> mappedReads = new HashSet<>();

       final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        //USE SamReaderFactory
        for (final SAMRecord samRecord : reader)
        {
                //is the read unmapped?
                boolean unmapped = samRecord.getReadUnmappedFlag();
                //if it is...
                if (unmapped == false)
                {
                    //is it paired
                    boolean isPaired = samRecord.getReadPairedFlag();
                    //if it is paired
                    if (isPaired == false)
                    {
                        //add it to the list
                        mappedReads.add(samRecord.getReadName());
                    }
                }
            }
        
        System.out.println("Found " + mappedReads.size() + " mapped reads");
        return mappedReads;
    }

    /**
     *
     * @param list a HashSet of paired reads
     * @param fastqInLeft the left-hand fastq reads which contains the reads in
     * list
     * @param fastqInRight the right-hand fastq reads which contains the reads
     * in list
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
        System.out.println("Found " + noFound + " from " + list.size());
    }

    /**
     *
     * @param list a HashSet of paired reads
     * @param fastqInLeft the left-hand fastq reads which contains the reads in
     * list
     * @param fastqInRight the right-hand fastq reads which contains the reads
     * in list
     * @param fastqSinglesOut the output fastq reads
     */
    public void writeSingleReadsFromHashMap(HashMap<String, Boolean> list, File fastqInLeft, File fastqInRight, File fastqSinglesOut)
    {
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter out = writer.newWriter(fastqSinglesOut);

        final FastqReader fastqReaderLeft = new FastqReader(fastqInLeft);
        final FastqReader fastqReaderRight = new FastqReader(fastqInRight);
        int noFound = 0;
        for (Map.Entry<String, Boolean> entry : list.entrySet())
        {
            String key = entry.getKey();
            Boolean value = entry.getValue();
            System.out.println("Read, " + key + " isLeft " + value);
        }
        while (fastqReaderRight.hasNext())
        {
            FastqRecord leftRecord = fastqReaderLeft.next();
            FastqRecord rightRecord = fastqReaderRight.next();
            String fullReadName = leftRecord.getReadHeader();
            int hashIndex = fullReadName.indexOf(" ");
            String readName = fullReadName.substring(0, hashIndex);

            if (list.containsKey(readName))
            {
                System.out.print("Found " + readName);
                System.out.println(list.get(readName));
                boolean isLeft = list.get(readName);
                System.out.print(" IsLeft? " + isLeft);
                if (isLeft)
                {
                    out.write(leftRecord);
                    System.out.println(" Left");

                } else
                {
                    out.write(rightRecord);
                    System.out.println(" Right");
                }
                noFound++;
            }
        }
        System.out.println("Found " + noFound + " from " + list.size());
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
        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        //USE SamReaderFactory
        for (final SAMRecord samRecord : reader)
        {
                String readHeader = samRecord.getReadName();
                System.out.println(readHeader);
            
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

    public void hashSetToTextFile(HashSet hs, File file) throws FileNotFoundException, UnsupportedEncodingException
    {
        try (PrintWriter writer = new PrintWriter(file, "UTF-8"))
        {
            for (Object hso : hs)
            {
                writer.println(hso.toString());
            }
            writer.close();
        }

    }

    public void hashMapToTextFile(HashMap<String, Boolean> hm, File file) throws FileNotFoundException, UnsupportedEncodingException
    {
        try (PrintWriter writer = new PrintWriter(file, "UTF-8"))
        {
            for (Map.Entry<String, Boolean> entry : hm.entrySet())
            {
                writer.print(entry.getKey());
                writer.print('\t');
                if (entry.setValue(Boolean.TRUE))
                {
                    writer.println("1");
                } else
                {
                    writer.println("2");
                }
            }
            writer.close();
        }
    }

    /**
     * Takes a tab-delimited file in the format <readname><pair-end>, where
     * pair-end is either '1', meaning a left-hand (forward) read or '2',
     * meaning a right-hand (reverse) read.
     *
     * @param file a tab-delimited file containing read-names and pair-end
     * annotation
     * @param leftFastq - the left-hand reads
     * @param rightFastq the right-hand reads
     * @param readsOut the file to write the reads listed in @param file to
     */
    public void getReadsFromList(File file, File leftFastq, File rightFastq, File readsOut)
    {
        HashSet<String> leftReads = new HashSet();
        HashSet<String> rightReads = new HashSet();

        Charset charset = Charset.forName("US-ASCII");
        try (BufferedReader reader = Files.newBufferedReader(file.toPath(), charset))
        {
            String line = null;
            while ((line = reader.readLine()) != null)
            {
                String[] array = line.split("\t");
                String readName = array[0];
                String pair = array[1];
                if (pair.equals("1"))
                {
                    leftReads.add(readName);
                } else
                {
                    rightReads.add(readName);
                }
            }
        } catch (IOException x)
        {
            System.err.format("IOException: %s%n", x);
        }
        int unmappedReads = leftReads.size() + rightReads.size();
        System.out.println("There are " + unmappedReads + " unmapped reads");

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter out = writer.newWriter(readsOut);
        final FastqReader fastqReaderLeft = new FastqReader(leftFastq);
        final FastqReader fastqReaderRight = new FastqReader(rightFastq);

        int noFound = 0;
        while (fastqReaderLeft.hasNext())
        {
            FastqRecord leftRecord = fastqReaderLeft.next();
            FastqRecord rightRecord = fastqReaderRight.next();
            String fullReadName = leftRecord.getReadHeader();
            int hashIndex = fullReadName.indexOf(" ");
            String readName = fullReadName.substring(0, hashIndex);
            if (leftReads.contains(readName))
            {
                out.write(leftRecord);
                noFound++;
            }
            if (rightReads.contains(readName))
            {
                out.write(rightRecord);
                noFound++;
            }

        }
        System.out.println("Wrote " + noFound + " reads");
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
