package piculus.gff;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import org.biojava3.genome.parsers.gff.*;
import piculus.fasta.ContigPosition;


/*
 * To change this template, choose Tools | Templates and open the template in
 * the editor.
 */
/**
 *
 * @author ethering
 */
public class BAMHandler
{

    public File BAM2GFF(String[] args) throws IOException
    {
        String source = "BGI";
        String type = "tRNA";
        Double score = 0.0;
        int frame = 0;

        File outfile = new File(args[1]);
        //HashMap will contain e.g.<Chr1, features>
        HashMap<String, FeatureList> chromFeaturesPlus = new HashMap();
        HashMap<String, FeatureList> chromFeaturesMinus = new HashMap();
        //for each bam file
        for (int i = 2; i < args.length; i++)
        {
            final SAMFileReader reader = new SAMFileReader(new File(args[i]));
            System.out.println("Opening " + args[i]);
            reader.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
            SAMRecordIterator iterator = reader.iterator();
            while (iterator.hasNext())
            {
                //create a gff entry
                SAMRecord samRecord = iterator.next();
                String seqname = samRecord.getReferenceName();
                boolean onMinusStrand = samRecord.getReadNegativeStrandFlag();
                int start = samRecord.getAlignmentStart();
                int end = samRecord.getAlignmentEnd();
                //Location loc = new Location(start, end);
                Location bioloc = Location.fromBio(start, end, '+');

                String readName = samRecord.getReadName();
                //make some attribute
                String attributes = "ID=".concat(readName);
                //make the gff feature
                Feature feat = new Feature(seqname, source, type, bioloc, score, frame, attributes);
                if (onMinusStrand)
                {
                    if (chromFeaturesMinus.containsKey(seqname))
                    {
                        //add the current feature to itPlus                    
                        chromFeaturesMinus.get(seqname).add(feat);
                    }
                    //if itPlus doesn't exist, create a new one
                    else
                    {
                        FeatureList fl = new FeatureList();
                        fl.add(feat);
                        chromFeaturesMinus.put(seqname, fl);
                    }
                }
                else
                {
                    //if the current chromosome or contig already has a FeatureList associated with itPlus..
                    if (chromFeaturesPlus.containsKey(seqname))
                    {
                        //add the current feature to itPlus                    
                        chromFeaturesPlus.get(seqname).add(feat);
                    }
                    //if itPlus doesn't exist, create a new one
                    else
                    {
                        FeatureList fl = new FeatureList();
                        fl.add(feat);
                        chromFeaturesPlus.put(seqname, fl);
                    }
                }
                //System.out.println(feat.toString());

            }
        }
        System.out.println("Finished reading file. Sorting");
        ArrayList<ContigPosition> plusCps = new ArrayList();


        //iterate over all the chromosomes and sort them by start position
        Iterator itPlus = chromFeaturesPlus.entrySet().iterator();
        while (itPlus.hasNext())
        {
            Map.Entry pairs = (Map.Entry) itPlus.next();
            //String ref = (String) pairs.getKey();
            //System.out.println(ref);
            FeatureList bamFeatures = (FeatureList) pairs.getValue();
            //get a new sorted FeatureList, whith all the features in all the bam files sorted
            FeatureList sortedFeatures = bamFeatures.sortByStart();


            //get an iterate for each FeatureList
            Iterator bamFeatIt = sortedFeatures.iterator();
            //get the first feature
            Feature startFi = (Feature) bamFeatIt.next();
            //get the location and the chromosome name
            Location firstLoc = startFi.location();
            int firstStart = firstLoc.bioStart();
            int firstEnd = firstLoc.bioEnd();
            String firstRef = startFi.seqname();
            //make a ContigPosition object from itPlus the first entry
            ContigPosition firstCp = new ContigPosition(firstRef, firstStart, firstEnd);
            //iterate over the remaining features
            while (bamFeatIt.hasNext())
            {
                Feature nextFi = (Feature) bamFeatIt.next();
                Location nextLoc = nextFi.location();
                int nextStart = nextLoc.bioStart();
                int nextEnd = nextLoc.bioEnd();
                String nextRef = nextFi.seqname();
                ContigPosition nextCp = new ContigPosition(nextRef, nextStart, nextEnd);

                if (nextCp.overlaps(firstCp))
                {
                    if (nextCp.getStart() < firstCp.getStart())
                    {
                        firstCp.setStart(nextCp.getStart());
                    }
                    if (nextCp.getEnd() > firstCp.getEnd())
                    {
                        firstCp.setEnd(nextCp.getEnd());
                    }
                }
                else
                {
                    plusCps.add(firstCp);
                    firstCp = nextCp;
                }
            }
            //need to add the last entry
            plusCps.add(firstCp);
        }


        ArrayList<ContigPosition> minusCps = new ArrayList();

        Iterator itMinus = chromFeaturesMinus.entrySet().iterator();
        while (itMinus.hasNext())
        {
            Map.Entry pairs = (Map.Entry) itMinus.next();
            //String ref = (String) pairs.getKey();
            //System.out.println(ref);
            FeatureList bamFeatures = (FeatureList) pairs.getValue();
            //get a new sorted FeatureList, whith all the features in all the bam files sorted
            FeatureList sortedFeatures = bamFeatures.sortByStart();


            //get an iterate for each FeatureList
            Iterator bamFeatIt = sortedFeatures.iterator();
            //get the first feature
            Feature startFi = (Feature) bamFeatIt.next();
            //get the location and the chromosome name
            Location firstLoc = startFi.location();
            int firstStart = firstLoc.bioStart();
            int firstEnd = firstLoc.bioEnd();
            String firstRef = startFi.seqname();
            //make a ContigPosition object from itPlus the first entry
            ContigPosition firstCp = new ContigPosition(firstRef, firstStart, firstEnd);
            //iterate over the remaining features
            while (bamFeatIt.hasNext())
            {
                Feature nextFi = (Feature) bamFeatIt.next();
                Location nextLoc = nextFi.location();
                int nextStart = nextLoc.bioStart();
                int nextEnd = nextLoc.bioEnd();
                String nextRef = nextFi.seqname();
                ContigPosition nextCp = new ContigPosition(nextRef, nextStart, nextEnd);

                if (nextCp.overlaps(firstCp))
                {
                    if (nextCp.getStart() < firstCp.getStart())
                    {
                        firstCp.setStart(nextCp.getStart());
                    }
                    if (nextCp.getEnd() > firstCp.getEnd())
                    {
                        firstCp.setEnd(nextCp.getEnd());
                    }
                }
                else
                {
                    minusCps.add(firstCp);
                    firstCp = nextCp;
                }
            }
            //need to add the last entry
            minusCps.add(firstCp);
        }

        BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
        Iterator it3 = plusCps.iterator();
        int clusterNo = 1;

        while (it3.hasNext())
        {
            ContigPosition cp = (ContigPosition) it3.next();
            Location bioloc = Location.fromBio(cp.getStart(), cp.getEnd(), '+');
            //System.out.println(cp.toString());
            Feature feat = new Feature(cp.getContigid(), source, type, bioloc, score, frame, ("ID=cluster_" + clusterNo));
            GFFObject gff = new GFFObject(feat, '-');
            //System.out.println(cp.toString());
            out.write(gff.toString() + "\n");
            clusterNo++;
        }
        clusterNo = 1;
        Iterator it4 = minusCps.iterator();
        while (it4.hasNext())
        {
            ContigPosition cp = (ContigPosition) it4.next();
            Location bioloc = Location.fromBio(cp.getStart(), cp.getEnd(), '-');
            Feature feat = new Feature(cp.getContigid(), source, type, bioloc, score, frame, ("ID=cluster_" + clusterNo));
            GFFObject gff = new GFFObject(feat, '-');
            //System.out.println(cp.toString());
            out.write(gff.toString() + "\n");
            clusterNo++;
        }
        out.close();
        return outfile;


    }

    //takes a BAM file of mapped reads and compares the mapping co-ordinates to features in a GFF file
    //calculates the number of total mapped reads
    public void createFpkmGffFiles(File bamfile, String gffFile, File fpkmGff) throws IOException, Exception
    {
        //keep account of the total fragments mapped to the bamfile
        double noReadsMappedtoBam = 0;

        String source = "BGI";
        String type = "tRNA";
        Double score = 0.0;
        int frame = 0;
        //the HashMaps follow a <ReadName, No_times_mapped> convention to account for multiple mapped reads
        HashMap<String, Integer> bamReads = new HashMap<String, Integer>();

        //the feature lists will hold the co-ordinates of the mapped reads and the attribute ID will be the read name
        FeatureList fl1 = new FeatureList();

        final SAMFileReader reader1 = new SAMFileReader(bamfile);
        reader1.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
        SAMRecordIterator iterator1 = reader1.iterator();
        String read = null;
        //iterate over the  bamfile
        while (iterator1.hasNext())
        {
            SAMRecord samRecord = iterator1.next();
            boolean unmapped = samRecord.getReadUnmappedFlag();
            //if the read is not unmapped
            if (!unmapped)
            {
                //get the read name (e.g. t00005_250, meaning that the tRNA read t00005 occurs 250 times) 
                read = samRecord.getReadName();

                //split it on the underscore
                String[] readArray = read.split("_");
                Integer integerExpr = Integer.parseInt(readArray[1]);
                //..and get the number of times it occurs
                int expression = integerExpr.intValue();

                //if the current read has already been put into the HashMap, increment it by one
                if (bamReads.containsKey(read))
                {
                    Integer count = bamReads.get(read);
                    bamReads.put(read, count++);
                    noReadsMappedtoBam += expression;
                }
                //else put in a new entry
                else
                {
                    bamReads.put(read, 1);
                    noReadsMappedtoBam += expression;
                }

                //now to make the Feature for the FeatureList
                //get the reference name
                String seqname = samRecord.getReferenceName();
                //find out which strand it's on and assign the correct char for the Location object
                boolean onMinusStrand = samRecord.getReadNegativeStrandFlag();
                char strand;
                if (onMinusStrand)
                {
                    strand = '-';
                }
                else
                {
                    strand = '+';
                }

                int start = samRecord.getAlignmentStart();
                int end = samRecord.getAlignmentEnd();
                Location bioloc = Location.fromBio(start, end, strand);
                //get the read name for the attribute 'ID'
                String readName = samRecord.getReadName();
                String attributes = "ID=".concat(readName);
                //make the new Feature
                Feature feat = new Feature(seqname, source, type, bioloc, score, frame, attributes);
                //add it to the FeatureList
                fl1.add(feat);

            }

        }


        System.out.println(noReadsMappedtoBam + " reads mapped to bam 1");

        //open the new gff files which will have each gff cluster loci with two attributes: ID=original cluster id and FPKM
        BufferedWriter out = new BufferedWriter(new FileWriter(fpkmGff));
        FeatureList gffFeatures = GFF3Reader.read(gffFile);

        for (FeatureI gffFeat : gffFeatures)
        {

            String seqname = gffFeat.seqname();
            Location location = gffFeat.location();
            int length = location.length();
            String clusterNo = gffFeat.getAttribute("ID");

            FeatureList subsetBam1 = fl1.selectOverlapping(seqname, location, false);

            if (subsetBam1.size() != 0)
            {
                double readCount1 = 0;
                for (FeatureI feat : subsetBam1)
                {
                    String id = feat.getAttribute("ID");
                    String[] readArray = id.split("_");
                    Integer integerCount = Integer.parseInt(readArray[1]);
                    //..and get the number of times it occurs
                    int readCount = integerCount.intValue();

                    Integer noTimesMapped = bamReads.get(id);
                    double mappingValue = readCount / noTimesMapped;
                    readCount1 += mappingValue;
                }
                double fpkm = getFPKM(readCount1, noReadsMappedtoBam, length);
                Feature feat = new Feature(seqname, source, type, location, score, frame, ("ID=" + clusterNo + "; FPKM=" + fpkm));
                GFFObject gff = new GFFObject(feat);
                //System.out.println(cp.toString());
                out.write(gff.toString() + "\n");
            }
        }
    }

    //calculate the FPKM
    public double getFPKM(double mappedReadsAtLoci, double totalMappedReads, int lociLength)
    {
        double oneBillion = java.lang.Math.pow(10, 9);
        double fpkm = (oneBillion * mappedReadsAtLoci) / (totalMappedReads * lociLength);
        return fpkm;
    }

    public void createDEGseqFilesFromSam(File bamfile1, File bamfile2, String gffFile, String degSeqFile) throws IOException, Exception
    {

        //the HashMaps follow a <ReadName, No_times_mapped> convention to account for multiple mapped reads
        HashMap<String, Integer> bam1Reads = new HashMap<String, Integer>();
        HashMap<String, Integer> bam2Reads = new HashMap<String, Integer>();

        final SAMFileReader reader1 = new SAMFileReader(bamfile1);
        reader1.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
        SAMRecordIterator iterator1 = reader1.iterator();
        String read = null;
        //iterate over the first bamfile
        while (iterator1.hasNext())
        {
            SAMRecord samRecord = iterator1.next();
            boolean unmapped = samRecord.getReadUnmappedFlag();
            //if the read is not unmapped
            if (!unmapped)
            {
                //get the read name (e.g. t00005_250, meaning that the tRNA read t00005 occurs 250 times) 
                read = samRecord.getReadName();

                //if the current read has already been put into the HashMap, increment it by one
                if (bam1Reads.containsKey(read))
                {
                    Integer count = bam1Reads.get(read);
                    bam1Reads.put(read, count++);
                }
                //else put in a new entry
                else
                {
                    bam1Reads.put(read, 1);
                }
            }
        }
        iterator1.close();

        //do the same for bamfile2
        final SAMFileReader reader2 = new SAMFileReader(bamfile2);
        reader2.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
        SAMRecordIterator iterator2 = reader2.iterator();
        while (iterator2.hasNext())
        {
            SAMRecord samRecord = iterator2.next();
            boolean unmapped = samRecord.getReadUnmappedFlag();
            if (!unmapped)
            {
                read = samRecord.getReadName();

                if (bam2Reads.containsKey(read))
                {
                    Integer count = bam2Reads.get(read);
                    bam2Reads.put(read, count++);
                }
                else
                {
                    bam2Reads.put(read, 1);
                }
            }
        }
        iterator2.close();

        //open the gff file which will have cluster loci with the attribute: ID=cluster id 
        BufferedWriter out = new BufferedWriter(new FileWriter(degSeqFile));
        out.write("\"EnsemblGeneID\"\t\"Chr\"\t\"GeneStart\"\t\"GeneEnd\"\t\"Status\"\t\"ExternalID\"\t\"Col-0-DC3000\"\t\"Scs9-2-DC3000\"\n");
        FeatureList gffFeatures = GFF3Reader.read(gffFile);

        
        for (FeatureI gffFeat : gffFeatures)
        {
            String clusterNo = gffFeat.getAttribute("ID");
            //keep account of the reads mapped to each cluster
            double readCount1 = 0;
            double readCount2 = 0;
            //get the information from the gff file required to query the bam file
            String seqname = gffFeat.seqname();
            Location location = gffFeat.location();
            char strand = location.bioStrand();
            SAMRecordIterator si1 = reader1.query(seqname, location.bioStart(), location.bioEnd(), false);

            while (si1.hasNext())
            {
                SAMRecord sr = si1.next();
                
                boolean unmapped = sr.getReadUnmappedFlag();
                if (!unmapped)
                {
                    //NEED TO ACCOUNT FOR STRANDEDNESS - get strand flag and see if it matches feature strand
                    boolean onMinusStrand = sr.getReadNegativeStrandFlag();
                    if (onMinusStrand && strand == '-' || !onMinusStrand && strand == '+')
                    {
                        //Get the read name..
                        String id = sr.getReadName();
                        //..and parse the read name to get the number of times it occurs
                        String[] readArray = id.split("_");
                        Integer integerCount = Integer.parseInt(readArray[1]);
                        int readCount = integerCount.intValue();
                        //then get its value (No. reads of this sequence/No. times mapped)
                        Integer noTimesMapped = bam1Reads.get(id);
                        double mappingValue = readCount / noTimesMapped;
                        readCount1 += mappingValue;
                    }
                }
            }
            si1.close();

            //same for the second bam file
            SAMRecordIterator si2 = reader2.query(seqname, location.bioStart(), location.bioEnd(), false);
            while (si2.hasNext())
            {
                SAMRecord sr = si2.next();
                boolean unmapped = sr.getReadUnmappedFlag();
                if (!unmapped)
                {
                    boolean onMinusStrand = sr.getReadNegativeStrandFlag();
                    if (onMinusStrand && strand == '-' || !onMinusStrand && strand == '+')
                    {
                        String id = sr.getReadName();
                        String[] readArray = id.split("_");
                        Integer integerCount = Integer.parseInt(readArray[1]);
                        int readCount = integerCount.intValue();
                        Integer noTimesMapped = bam2Reads.get(id);
                        double mappingValue = readCount / noTimesMapped;
                        readCount2 += mappingValue;
                    }
                }
            }
            si2.close();
            out.write(clusterNo + "\t" + seqname + "\t" + location.bioStart() + "\t" + location.bioEnd() + "\tKNOWN\t" + clusterNo + "\t" + readCount1 + "\t" + readCount2 + "\n");
        }
        out.close();
    }
}