package uk.ac.tsl.etherington.genomehelper.gff;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.BioException;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.Strand;
import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.biojava3.genome.parsers.gff.GFF3Reader;
import org.biojava3.genome.parsers.gff.Location;
import uk.ac.tsl.etherington.genomehelper.fasta.FastaFeatures;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;



/*
 *Contains a range of methods that calculate a number of genome statistics for a requested feature type (the 3rd column in a gff file)
 * normally returning the mean length of the feature and printing out lots of stats associated with the feature.
 */
/**
 *
 * @author ethering
 */
public class GFFFeatureStats
{

    /**
     * Combines a number of methods
     *
     * @param gff the GFF annotation file
     * @param refSeq the fasta reference file
     * @param attribute the name of the attribute that will make the genes
     * unique (e.g. "name", gene_id, etc)
     * @throws FileNotFoundException
     * @throws BioException
     * @throws IOException
     */
    public void getStats(String gff, File refSeq, String attribute) throws FileNotFoundException, BioException, IOException, Exception
    {

        HashMap<String, int[]> genomeMap = new HashMap<>(FastaFeatures.getSequenceAsIntArray(refSeq));
        double genomeSize = getGenomeSizeFromIntArrayHashMap(genomeMap);
        double genomeMb = genomeSize / 1048576;
        System.out.println("Genome size = " + genomeSize + "(" + genomeMb + "MB)");
        FeatureList fl = GFF3Reader.read(gff);
        try
        {
            getMeanFeatureLength(fl, "CDS");
            System.out.println("");
            getMeanFeatureLength(fl, "exon");
            System.out.println("");
            getMeanIntronLength(fl, attribute, genomeSize);
            System.out.println("");
            //genomeMap = new HashMap<String, int[]>(fastautils.FastaFeatures.getSequenceAsIntArray(refSeq));
            //getCodingAndNonCodingRegions(fl, genomeMap, attribute, genomeSize);
        }
        catch (Exception ex)
        {
            Logger.getLogger(GFFFeatureStats.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Calculates the length of a genome from collection of int arrays, where
     * each int array is the length of a genome sequence
     *
     * @param genomeMap an integer array of genome sequence lengths
     * @return the size of the genome (in nucleotides)
     */
    public double getGenomeSizeFromIntArrayHashMap(HashMap<String, int[]> genomeMap)
    {
        double genomeSize = 0;
        for (Map.Entry<String, int[]> arraySet : genomeMap.entrySet())
        {
            genomeSize += arraySet.getValue().length;
        }
        return genomeSize;
    }

    /**
     * Gets all the features from a gff file
     *
     * @param gffFile the gff file from which to extract the features
     * @return a org.biojava3.genome.parsers.gff.FeatureList
     * @throws IOException
     */
    public FeatureList getFeatureList(String gffFile) throws IOException
    {
        FeatureList fl = GFF3Reader.read(gffFile);
        return fl;
    }

    /**
     * Prints stats and mean feature length for a given feature
     *
     * @param fl a org.biojava3.genome.parsers.gff.FeatureList
     * @param genomeMap a HashMap where the keys are the reference sequence
     * names and the values a zero-filled int [] the length of the associated
     * reference sequence
     * @param featureName the feature-type to analyse (e.g. mRNA, exon, CDS,
     * etc)
     * @return the mean length of the requested feature name
     * @throws IOException
     * @throws Exception
     */
    public double getMeanFeatureLength(FeatureList fl, HashMap<String, int[]> genomeMap, String featureName) throws IOException, Exception
    {
        FeatureList featTypes = new FeatureList(fl.selectByType(featureName));
        if (featTypes.size() == 0)
        {
            System.out.println("No features to examine:" + featureName);
            return 0;
        }
        for (FeatureI fi : featTypes)
        {
            String seqname = fi.seqname();
            if (genomeMap.containsKey(seqname))
            {
                Location loc = fi.location();
                int start = loc.bioStart();
                int end = loc.bioEnd();
                //get the int array for the current reference sequence
                int[] seqMap = genomeMap.get(seqname);
                //the location is 1-based, but the array is zero-based
                for (int i = start - 1; i < end - 1; i++)
                {
                    //increment the location of the feature in the int [] by one
                    seqMap[i]++;
                }
                genomeMap.put(seqname, seqMap);
            }
        }

        //feats will be an array list of feature lengths so we know the number of features and how long they are
        ArrayList<Integer> feats = new ArrayList<>();
        double genomeLength = 0;
        int altSplicing = 0;
        for (Map.Entry<String, int[]> arraySet : genomeMap.entrySet())
        {
            //for each map, get the int array
            //System.out.println("Chr: "+arraySet.getKey());
            int[] codingPositions = arraySet.getValue();
            genomeLength += codingPositions.length;
            //System.out.println("Length: "+codingPositions.length);
            //iterate over each int array
            for (int i = 0; i < codingPositions.length; i++)
            {
                //if we come across a coding region (1)
                if (codingPositions[i] > 0)
                {
                    //note the start site
                    int featStart = i;
                    //until we don't get a 1, or utnil we get to the end of the contig
                    while (codingPositions[i] > 0 && i < codingPositions.length)
                    {
                        //if the codingPosition has been counted more than once, it must have been alternatively spliced
                        if (codingPositions[i] > 1)
                        {
                            altSplicing += codingPositions[i] - 1;
                        }
                        //advance along the array
                        i++;
                    }
                    //while exited, so we must have reached the end of the feature
                    //note the end site
                    int featEnd = i;
                    //calculate the length of the feature
                    int featLength = featEnd - featStart + 1;
                    //and put it into the feats arraylist
                    feats.add(featLength);
                }

            }
            //System.out.println();
        }
        int codingSpace = 0;
        double[] featArray = new double[feats.size()];
        int i = 0;
        for (Integer featInt : feats)
        {
            codingSpace += featInt;
            featArray[i] = featInt.doubleValue();
            i++;
        }

        System.out.println("Feature calculated: " + featureName);
        double meanLength = codingSpace / feats.size();
        System.out.println("Mean feature length = " + meanLength);
        double stdev = getStandardDeviation(featArray);
        System.out.println("Stddev = " + stdev);
        double stderr = getStandardErrorOfMean(stdev, featArray.length);
        System.out.println("Standard error of mean = "+stderr);
        double codingProportion = codingSpace / genomeLength;
        double genomeMb = genomeLength / 1048576;
        System.out.println("Genome size = " + genomeLength + "(" + genomeMb + "MB)");

        System.out.println("Proportion of genome used for feature = " + codingProportion);
        System.out.println("Amount of genome alternatively spliced = " + altSplicing + " nucleotides");
        return meanLength;
    }

    /**
     * Prints stats and mean feature length for a given feature, but restricts
     * it to a list of geneIds of interest
     *
     * @param fl a org.biojava3.genome.parsers.gff.FeatureList
     * @param genomeMap a HashMap where the keys are the reference sequence
     * names and the values a zero-filled int [] the length of the associated
     * reference sequence
     * @param featureName the feature-type to analyse (e.g. mRNA, exon, CDS,
     * etc)
     * @param geneIds a file of gene IDs, one per line
     * @param attribute the gff attribute that the gene IDs will be listed
     * under, e.g. gene_id, transcript_id, etc
     * @return prints out stats and returns the mean feature length of the
     * featureName
     * @throws IOException
     * @throws Exception
     */
    public double getMeanFeatureLength(FeatureList fl, HashMap<String, int[]> genomeMap, String featureName, File geneIds, String attribute) throws IOException, Exception
    {
        Set<String> secretedProteins = new HashSet<>();
        BufferedReader br = new BufferedReader(new FileReader(geneIds));
        String line = null; //not declared within while loop
        while ((line = br.readLine()) != null)
        {
            secretedProteins.add(line);
        }
        System.out.println("GeneIds contains " + secretedProteins.size() + " names");

        //open the gff file and extract the requested features
        FeatureList featTypes = new FeatureList(fl.selectByType(featureName));

        //iterate over the gff file
        for (FeatureI fi : featTypes)
        {
            if (fi.hasAttribute(attribute))
            {
                String geneId = fi.getAttribute(attribute);
                if (secretedProteins.contains(geneId))
                {
                    String seqname = fi.seqname();
                    if (genomeMap.containsKey(seqname))
                    {
                        Location loc = fi.location();
                        int start = loc.bioStart();
                        int end = loc.bioEnd();
                        int[] seqMap = genomeMap.get(seqname);
                        //the location is 1-based, but the array is zero-based
                        for (int i = start - 1; i < end - 1; i++)
                        {
                            seqMap[i]++;
                        }
                        genomeMap.put(seqname, seqMap);
                    }
                }
            }
        }
        if (genomeMap.isEmpty())
        {
            System.out.println("No features found with attribute " + attribute);
            return 0;
        }
        //feats will be an array list of feature lengths so we know the number of features and how long they are
        ArrayList<Integer> feats = new ArrayList<Integer>();
        double genomeLength = 0;
        int altSplicing = 0;
        for (Map.Entry<String, int[]> arraySet : genomeMap.entrySet())
        {
            //for each map, get the int array
            //System.out.println("Chr: "+arraySet.getKey());
            int[] codingPositions = arraySet.getValue();
            genomeLength += codingPositions.length;
            //System.out.println("Length: "+codingPositions.length);
            //iterate over each int array
            for (int i = 0; i < codingPositions.length; i++)
            {
                //if we come across a coding region (1)
                if (codingPositions[i] == 1)
                {
                    //note the start site
                    int featStart = i;
                    //until we don't get a 1, or utnil we get to the end of the contig
                    while (codingPositions[i] > 0 && i < codingPositions.length)
                    {
                        if (codingPositions[i] > 1)
                        {
                            altSplicing += codingPositions[i] - 1;
                        }
                        //advance along the array
                        i++;
                    }
                    //while exited, so we must have reached the end of the feature
                    //note the end site
                    int featEnd = i;
                    //calculate the length of the feature
                    int featLength = featEnd - featStart + 1;
                    //System.out.println("Feat size: "+featLength);
                    //and put it into the feats arraylist
                    feats.add(featLength);
                }

            }
            //System.out.println();
        }
        int codingSpace = 0;
        double[] featArray = new double[feats.size()];
        int i = 0;
        for (Integer featInt : feats)
        {
            codingSpace += featInt;
            featArray[i] = featInt.doubleValue();
            i++;
        }


        System.out.println("Feature calculated: " + featureName);
        double meanLength = codingSpace / feats.size();
        System.out.println("Mean feature length = " + meanLength);
        double stdev = getStandardDeviation(featArray);
        System.out.println("Stddev = " + stdev);
        double stderr = getStandardErrorOfMean(stdev, featArray.length);
        System.out.println("Standard error of mean = "+stderr);
        double codingProportion = codingSpace / genomeLength;
        double genomeMb = genomeLength / 1048576;
        System.out.println("Genome size = " + genomeLength + "(" + genomeMb + "MB)");

        System.out.println("Proportion of genome used for feature = " + codingProportion);
        System.out.println("Amount of genome alternatively spliced = " + altSplicing + " nucleotides");
        return meanLength;
    }

    /**
     * Calculates the mean feature length of a given feature
     *
     * @param fl a org.biojava3.genome.parsers.gff.FeatureList
     * @param featureName the feature-type to analyse (e.g. mRNA, exon, CDS,
     * etc)
     * @return the mean feature length of the given feature
     * @throws IOException
     */
    public double getMeanFeatureLength(FeatureList fl, String featureName) throws IOException
    {
        FeatureList featTypes = new FeatureList(fl.selectByType(featureName));
        Iterator it = featTypes.iterator();
        double combinedFeatureLength = 0;
        double numberOfFeatures = 0;
        int i = 0;
        double[] featArray = new double[featTypes.size()];
        while (it.hasNext())
        {
            FeatureI fi = (FeatureI) it.next();
            Location loc = fi.location();
            int start = loc.bioStart();
            int end = loc.bioEnd();
            int featureLength = end - start;
            //System.out.print(featureLength + ", ");
            featArray[i] = featureLength;
            //System.out.println("Start = " + start + " End = " + end + " Length = " + featureLength);
            combinedFeatureLength += featureLength;
            numberOfFeatures++;
            i++;

        }
        double meanLength = combinedFeatureLength / numberOfFeatures;
        double stdev = getStandardDeviation(featArray);
        System.out.println("Stddev = " + stdev);
        double stderr = getStandardErrorOfMean(stdev, featArray.length);
        System.out.println("Standard error of mean = "+stderr);
        System.out.println("Feature calculated: " + featureName);
        System.out.println("Mean feature length = " + meanLength);
        return meanLength;
    }

    /**
     * Writes the DNA sequence of the non-coding portion of the genome (i.e.
     * intergenic and intron).
     *
     * @param fl a org.biojava3.genome.parsers.gff.FeatureList
     * @param refSeq the reference sequence (in fasta format)
     * @param nonCodingGenome the non-coding reference geneom (in fasta format)
     * @throws Exception
     */
    public void createNonCodingGenome(FeatureList fl, File refSeq, File nonCodingGenome) throws Exception
    {
        Writer output = new BufferedWriter(new FileWriter(nonCodingGenome));
        LinkedHashMap<String, DNASequence> genome = FastaFeatures.getParsedDNASequences(refSeq);
        //a HashMap for the sequence lengths, e.g. <Chr1><5500>
        HashMap<String, Integer> seqLengths = new HashMap<>(FastaFeatures.getSequenceLengths(refSeq));
        //a HashMap for the non-coding blocks, e.g. <Chr1><[1, 50, 55, 70....]>
        HashMap<String, ArrayList<Integer>> blocks = getBlocks(fl, "exon");



        //iterate over the blocks HashMap and put, as the final array list entry for each sequence
        //the length of the sequence (so as to include the very last bit of intergenic sequence
        for (Map.Entry<String, ArrayList<Integer>> entry : blocks.entrySet())
        {
            String chr = entry.getKey();
            System.out.println("Chr " + chr);
            ArrayList<Integer> al = entry.getValue();
            Collections.sort(al);
            int firstPosition = al.get(0);
            System.out.println("first pos = " + firstPosition);
            int lastPosition = al.get(al.size() - 1);
            System.out.println("last pos = " + lastPosition);

            //check to see if the current contig is included in the seqLengths list (it should be!)

            if (seqLengths.containsKey(chr))
            {
                if (lastPosition != seqLengths.get(chr))
                {
                    al.add(seqLengths.get(chr));
                }
                //if the first block is not '1' add 0 as the first index (this will get incremented to '1' later) so that the first intergenic block starts at the 
                //start of the referene sequence and ends at the start of the first exon 
                if (firstPosition != 1)
                {
                    al.add(0);
                }
                else//the first block starts at '1', so we need to remove it so that the first element is the end of the first exon
                {
                    al.remove(0);
                }

            }
        }

        //iterate over the blocks Hashmap...
        for (Map.Entry<String, ArrayList<Integer>> entry : blocks.entrySet())
        {
            StringBuilder nonCodingContig = new StringBuilder();
            System.out.println(entry.getKey());
            //get each arrayList...
            ArrayList<Integer> al = entry.getValue();
            String contig = entry.getKey();
            nonCodingContig.append(">").append(contig);
            nonCodingContig.append(System.getProperty("line.separator"));

            DNASequence seq = genome.get(contig);
            //sort the array of blocks
            Collections.sort(al);
            Iterator alit = al.iterator();
            //iterate over the blocks
            while (alit.hasNext())
            {
                //get the start of the first block..
                int startBlock = (Integer) alit.next();
                //System.out.println("Start block = "+startBlock);
                //and then the start of the next block
                if (alit.hasNext())
                {
                    int endBlock = (Integer) alit.next();
                    //System.out.println("end block = "+endBlock);
                    System.out.println("Contig = " + contig + " Start = " + (startBlock + 1) + " End = " + (endBlock - 1));
                    //the block start/end coords are both part of the coding sequence,
                    //so to get the intergenic region +1 to the start and -1 to the end
                    //Also, check for consecutive exons (where one exon starts straight after the other)
                    if ((endBlock - (startBlock + 1)) > 0)
                    {
                        String subseq = seq.getSequenceAsString((startBlock + 1), (endBlock - 1), Strand.POSITIVE);
                        nonCodingContig.append(subseq);
                    }
                }
            }
            nonCodingContig.append(System.getProperty("line.separator"));
            output.write(nonCodingContig.toString());
        }
        output.close();
    }

    /**
     *
     * @param fl a org.biojava3.genome.parsers.gff.FeatureList
     * @param feature the feature to use for the coding genome (e.g. exon, cds,
     * mRNA)
     * @param refSeq the reference sequence (in fasta format)
     * @param codingGenome the non-coding reference genome (in fasta format)
     * @throws Exception
     */
    public void createCodingGenome(FeatureList fl, String feature, File refSeq, File codingGenome) throws Exception
    {
        Writer output = new BufferedWriter(new FileWriter(codingGenome));
        LinkedHashMap<String, DNASequence> genome = FastaFeatures.getParsedDNASequences(refSeq);
        //a HashMap for the non-coding blocks, e.g. <Chr1><[1, 50, 55, 70....]>
        HashMap<String, ArrayList<String>> blocks = new HashMap<>();

        //just select the exons
        FeatureList featTypes = new FeatureList(fl.selectByType(feature));
        for (FeatureI currentfi : featTypes)
        {
            String seqName = currentfi.seqname();

            DNASequence dna = genome.get(seqName);
            Location currentLoc = currentfi.location();
            int currentStart = currentLoc.bioStart();
            int currentEnd = currentLoc.bioEnd();
            char bioStrand = currentLoc.bioStrand();
            Strand strand;
            if (bioStrand == '-')
            {
                strand = Strand.NEGATIVE;
            }
            else
            {
                strand = Strand.POSITIVE;
            }
            String subseq = dna.getSequenceAsString(currentStart, currentEnd, strand);
            //System.out.println(seqName +" "+strand.getStringRepresentation());
            //if the current sequence is already in the blocks HashMap, 
            //get the array list and insert the current start and stop positions
            if (blocks.containsKey(seqName))
            {
                ArrayList<String> al = (ArrayList<String>) blocks.get(seqName);
                al.add(subseq);
                blocks.put(seqName, al);
            }
            //if it's a new sequence, get the array list...
            else
            {
                ArrayList<String> al = new ArrayList<>();
                al.add(subseq);
                blocks.put(seqName, al);
            }
        }

        //iterate over the blocks Hashmap...
        for (Map.Entry<String, ArrayList<String>> entry : blocks.entrySet())
        {
            StringBuilder nonCodingContig = new StringBuilder();
            //System.out.println(entry.getKey());
            //get each arrayList...
            ArrayList<String> al = entry.getValue();
            String contig = entry.getKey();
            nonCodingContig.append(">").append(contig);
            nonCodingContig.append(System.getProperty("line.separator"));

            Iterator alit = al.iterator();
            //iterate over the blocks
            while (alit.hasNext())
            {
                String subseq = alit.next().toString();
                nonCodingContig.append(subseq);
            }
            nonCodingContig.append(System.getProperty("line.separator"));
            output.write(nonCodingContig.toString());
        }
        output.close();
    }

    /**
     *
     * @param fl a org.biojava3.genome.parsers.gff.FeatureList
     * @param attribute the name of the attribute that will make your genes
     * unique (e.g. "name", gene_id, etc)
     * @param genomeSize the size of the genome
     * @return the mean intron length
     * @throws IOException
     * @throws FileNotFoundException
     * @throws BioException
     */
    public double getMeanIntronLength(FeatureList fl, String attribute, double genomeSize) throws IOException, FileNotFoundException, BioException
    {
        //String will be the attribute value. The HashMap will contain the strand Character and the ArrayList of exon start/stops
        //this is so we don't calculate introns from overlapping exons on different strands
        
        HashMap<String, ArrayList<Integer>> blocks = getBlocks(fl, "exon", attribute);
        double nGenes = blocks.size();
        double combinedIntronLength = 0;

        ArrayList<Double> intronArray = new ArrayList<>();
        if (blocks.isEmpty())
        {
            System.out.println("No features found with attribute " + attribute);
            return 0;
        }
        for (Map.Entry<String, ArrayList<Integer>> entry : blocks.entrySet())
        {
            //System.out.print(entry.getKey());
            ArrayList<Integer> al = entry.getValue();

            Collections.sort(al);
            Iterator alit = al.iterator();
            //not interested in the first co-ordinate as it's the exon start
            //we need the exon end as our first number, so we'll skip the first one
            alit.next();
            while (alit.hasNext())
            {
                int startBlock = (Integer) alit.next();//ergo the end of the exon
                //System.out.print(startBlock);
                if (alit.hasNext())
                {
                    int endBlock = (Integer) alit.next();//the start of the next exon
                    //System.out.print("\t"+endBlock);
                    //the block start/end coords are both part of the exon sequence,
                    //so to get the intron length calculate end - start -1
                    double blockLength = endBlock - startBlock - 1;
                    intronArray.add(blockLength);
                    //System.out.println("\t" + blockLength);
                    combinedIntronLength += blockLength;

                }
            }
            //System.out.println("");
        }
        double nIntrons = intronArray.size();
        double[] featArray = new double[(int)nIntrons];
        for (int x = 0; x < intronArray.size(); x++)
        {
            featArray[x] = (double) intronArray.get(x);
        }        
        
        double meanLength = combinedIntronLength / nIntrons;
        double proportionIntron = combinedIntronLength / genomeSize;
        double meanIntronGene = nIntrons / nGenes;
        System.out.println("nGenes = " + nGenes);
        System.out.println("nIntrons = " + nIntrons);
        System.out.println("Mean intron length = " + meanLength);
        double stdev = getStandardDeviation(featArray);
        System.out.println("Stddev = " + stdev);
        double stderr = getStandardErrorOfMean(stdev, featArray.length);
        System.out.println("Standard error of mean = "+stderr);
        System.out.println("Proportion of genome intronic = " + proportionIntron);
        System.out.println("Mean No. introns per gene = " + meanIntronGene);
        return meanLength;
    }

    /**
     *
     * @param fl a org.biojava3.genome.parsers.gff.FeatureList
     * @param attribute the name of the attribute that will make the genes
     * unique (e.g. "name", gene_id, etc)
     * @param secretedProteinsFile a list of genes that you want to restrict
     * your analysis to
     * @return the mean intron length of your restricted gene list
     * @throws IOException
     */
    public double getMeanSecretedIntronLength(FeatureList fl, String attribute, File secretedProteinsFile) throws IOException
    {
        Set<String> secretedProteins = new HashSet<>();
        BufferedReader br = new BufferedReader(new FileReader(secretedProteinsFile));
        String line = null; //not declared within while loop
        while ((line = br.readLine()) != null)
        {
            secretedProteins.add(line);
            //System.out.println(line);
        }
        //String will be the attribute value. The HashMap will contain the strand Character and the ArrayList of exon start/stops
        //this is so we don't calculate introns from overlapping exons on different strands
        HashMap<String, ArrayList<Integer>> blocks = getBlocks(fl, "exon", attribute, secretedProteins);
        ArrayList<Double> intronArray = new ArrayList<>();
        if (blocks.isEmpty())
        {
            System.out.println("No features found with attribute " + attribute);
            return 0;
        }
        double combinedIntronLength = 0;


        for (Map.Entry<String, ArrayList<Integer>> entry : blocks.entrySet())
        {
            //System.out.print(entry.getKey());
            ArrayList<Integer> al = entry.getValue();

            Collections.sort(al);
            Iterator alit = al.iterator();

            //not interested in the first co-ordinate as it's the exon start
            //we need the exon end as our first number, so we'll skip the first one
            alit.next();
            while (alit.hasNext())
            {
                int startBlock = (Integer) alit.next();//ergo the end of the exon
                if (alit.hasNext())
                {
                    int endBlock = (Integer) alit.next();//the start of the next exon
                    //the block start/end coords are both part of the exon sequence,
                    //so to get the intron length calculate end - start -1
                    double blockLength = endBlock - startBlock - 1;
                    intronArray.add(blockLength);
                    combinedIntronLength += blockLength;
                }
            }
        }
        double nIntrons = intronArray.size();
        double[] featArray = new double[(int)nIntrons];
        for (int x = 0; x < intronArray.size(); x++)
        {
            featArray[x] = (double) intronArray.get(x);
        }
        double nGenes = blocks.size();

        double meanLength = combinedIntronLength / nIntrons;

        double meanIntronGene = nIntrons / nGenes;
        System.out.println("nGenes = " + nGenes);
        System.out.println("nIntrons = " + nIntrons);

        System.out.println("Mean No. introns per gene = " + meanIntronGene);
        System.out.println("Mean intron length = " + meanLength);
        double stdev = getStandardDeviation(featArray);
        System.out.println("Stddev = " + stdev);
        double stderr = getStandardErrorOfMean(stdev, featArray.length);
        System.out.println("Standard error of mean = "+stderr);
        System.out.println("Stddev = " + stdev);
        return meanLength;
    }

    /**
     * Calculate the proportion of the genome that contains coding regions (from
     * the start of the first exon to the end of the last exon)
     *
     * @param featList a org.biojava3.genome.parsers.gff.FeatureList
     * @param codingMap a HashMap where the keys are the reference sequence
     * names and the values a zero-filled int [] the length of the associated
     * reference sequence
     * @param attribute the name of the attribute that will make the genes
     * unique (e.g. "name", gene_id, etc)
     * @param genomeLength the length of the genome
     * @return the proportion of the genome which contains exons
     * @throws IOException
     * @throws Exception
     */
    public double calculateCodingRegion(FeatureList featList, HashMap<String, int[]> codingMap, String attribute, double genomeLength) throws IOException, Exception
    {
        HashMap<String, HashMap<String, ArrayList<Integer>>> features = new HashMap<>();
        //go across the gff file, note the sequence name, gene id and start/stop positions of each feature
        //store info e.g. <Chr1><gene_5><[5, 10, 30, 40, 100, 120]>
        FeatureList fl = featList.selectByType("exon");
        for (FeatureI fi : fl)
        {
            if (fi.hasAttribute(attribute))
            {
                //get the sequence name and the start and end position of the  exon
                String seqName = fi.seqname();
                //System.out.println(seqName);
                Location loc = fi.location();
                int start = loc.bioStart();
                int end = loc.bioEnd();
                String geneId = fi.getAttribute(attribute);
                //System.out.println(geneId);
                //if the  sequence is already in the blocks HashMap, 
                //get the array list and insert the  start and stop positions
                if (features.containsKey(seqName))
                {
                    HashMap<String, ArrayList<Integer>> hm = features.get(seqName);//e.g. <gene_5><[5, 10, 20, 30]>

                    ArrayList<Integer> al;
                    //if there's an entry for the current geneId, e.g. <gene_5>
                    if (hm.containsKey(geneId))
                    {
                        //get the arrayList for the geneId
                        al = (ArrayList<Integer>) hm.get(geneId);//e.g. <[5, 10, 20, 30]>
                        al.add(start);//e.g. 50
                        al.add(end);//e.g. 62
                        hm.put(geneId, al);//e.g. <gene_5><[5, 10, 20, 30, 50, 62]>

                    }
                    //if it's a new geneId, e.g. <gene_6>
                    else
                    {
                        //create a new ArrayList
                        al = new ArrayList<>();
                        //add the start and stop positions to it
                        al.add(start);//e.g. 10
                        al.add(end);//e.g. 25
                        // and put it into hm
                        hm.put(geneId, al);//e.g. <gene_6><[10, 25]>
                    }
                    //add hm 
                    features.put(seqName, hm);
                }
                //if it's a new sequence, create a new HashMap<String, ArrayList<Integer>> and ArrayList<Integer>
                else
                {
                    HashMap<String, ArrayList<Integer>> hm = new HashMap<>();
                    ArrayList<Integer> al = new ArrayList<>();
                    al.add(start);
                    al.add(end);
                    hm.put(geneId, al);
                    features.put(seqName, hm);
                }
            }
        }
        if (features.isEmpty())
        {
            System.out.println("No features found with attribute " + attribute);
            return 0;
        }
        System.out.println("Finished reading gff file");

        //iterate over the stored features, sort the array of start/stop positions
        // and get the highest and lowest positions for each gene (i.e. start and end of the gene)
        for (Map.Entry<String, HashMap<String, ArrayList<Integer>>> entry : features.entrySet())
        {
            //for each chromosome
            String seqName = entry.getKey();
            //System.out.println("SeqName = "+seqName);
            //get the int [] array from the codingMap
            int[] seqMap = codingMap.get(seqName);
            //get the hashMap of gene ids and their associated arraylist of start/stop features
            HashMap<String, ArrayList<Integer>> hm = features.get(seqName);
            //iterate over that hashMap
            for (Map.Entry<String, ArrayList<Integer>> geneArray : hm.entrySet())
            {
                //System.out.println("Feature id = "+geneArray.getKey());
                //and get the arrayList
                ArrayList<Integer> al = geneArray.getValue();
                //sort the arrayList and get the first and last element.
                //this should be the start and end of each coding region
                Collections.sort(al);
                int lowest = al.get(0);
                int highest = al.get(al.size() - 1);
                //System.out.println("lowest = "+lowest);
                //System.out.println("highest = "+highest);
                //and set the coding region to '1'
                for (int i = lowest - 1; i < highest; i++)
                {
                    seqMap[i] = 1;
                }

            }
            //put the int [] array back into codingMap
            codingMap.put(seqName, seqMap);
        }

        // iterate over the codingMap and add up all the 1's. This is our coding space
        int codingSpace = 0;
        for (Map.Entry<String, int[]> arraySet : codingMap.entrySet())
        {
            //String chr = arraySet.getKey();
            int[] codingPositions = arraySet.getValue();
            // System.out.println(chr);//and set the coding region to '1'
            for (int i = 0; i < codingPositions.length; i++)
            {
                //System.out.print(" " + codingPositions[i]);
                codingSpace += codingPositions[i];
            }
            //System.out.println();
        }

        System.out.println("Coding space = " + codingSpace);
        System.out.println("Genome Length = " + genomeLength);
        //sequence length - coding space is the non-coding length
        double codingProportion = codingSpace / genomeLength;
        System.out.println("Coding proportion = " + codingProportion);
        //coding space / sequence length is the proportion of coding space
        double nonCodingSpace = genomeLength - codingSpace;
        double nonCodingProportion = nonCodingSpace / genomeLength;
        System.out.println("Non-coding proportion = " + nonCodingProportion);
        //non-coding space / sequence length is the proportion of non-coding space
        return codingProportion;
    }

    /**
     * Creates start/end blocks of features
     *
     * @param fl a org.biojava3.genome.parsers.gff.FeatureList
     * @param featureType the type of feature to make blocks from
     * @return a HashMap of blocks where the keys are the names of the reference
     * sequences and the values are arrays of Integers with the consecutive
     * start and then finish co-ordinates of all the requested featureTypes
     */
    public HashMap<String, ArrayList<Integer>> getBlocks(FeatureList fl, String featureType)
    {
        HashMap<String, ArrayList<Integer>> blocks = new HashMap<>();
        FeatureList featTypes = fl.selectByType(featureType);
        for (FeatureI currentfi : featTypes)
        {
            String seqName = currentfi.seqname();
            Location currentLoc = currentfi.location();
            int currentStart = currentLoc.bioStart();
            int currentEnd = currentLoc.bioEnd();
            //System.out.println(geneId + " " + currentStart + " " + currentEnd);
            //see if the attribute name exists
            if (blocks.containsKey(seqName))
            {
                //if it does, get the HashMap and test to see if the current strand exists
                ArrayList<Integer> al = (ArrayList<Integer>) blocks.get(seqName);
                al.add(currentStart);
                al.add(currentEnd);
                blocks.put(seqName, al);
            }
            //if the attribute doesn't exist
            else
            {
                //create a new ArrayList
                ArrayList<Integer> al = new ArrayList<>();
                //add the coords to the arraylist
                al.add(currentStart);
                al.add(currentEnd);
                //add the new geneId to the blocks hashmap
                blocks.put(seqName, al);
            }
        }
        return blocks;
    }

    /**
     * Creates start/end blocks of features
     *
     * @param fl a org.biojava3.genome.parsers.gff.FeatureList
     * @param featureType the type of feature to make blocks from
     * @param attribute the name of the attribute that will make the genes
     * unique (e.g. "name", gene_id, etc)
     * @return a HashMap of blocks where the keys are the names of the reference
     * sequences and the values are arrays of Integers with the consecutive
     * start and then finish co-ordinates of all the requested featureTypes
     */
    public HashMap<String, ArrayList<Integer>> getBlocks(FeatureList fl, String featureType, String attribute)
    {
        HashMap<String, ArrayList<Integer>> blocks = new HashMap<>();
        FeatureList featTypes = fl.selectByType(featureType);
        for (FeatureI currentfi : featTypes)
        {
            if (currentfi.hasAttribute(attribute))
            {
                String geneId = currentfi.getAttribute(attribute);
                Location currentLoc = currentfi.location();
                int currentStart = currentLoc.bioStart();
                int currentEnd = currentLoc.bioEnd();
                //System.out.println(geneId + " " + currentStart + " " + currentEnd);
                //see if the attribute name exists
                if (blocks.containsKey(geneId))
                {
                    //if it does, get the HashMap and test to see if the current strand exists
                    ArrayList<Integer> al = (ArrayList<Integer>) blocks.get(geneId);
                    al.add(currentStart);
                    al.add(currentEnd);
                    blocks.put(geneId, al);
                }
                //if the attribute doesn't exist
                else
                {
                    //create a new ArrayList
                    ArrayList<Integer> al = new ArrayList<>();
                    //add the coords to the arraylist
                    al.add(currentStart);
                    al.add(currentEnd);
                    //add the new geneId to the blocks hashmap
                    blocks.put(geneId, al);
                }
            }
        }
        return blocks;
    }

    /**
     * Creates start/end blocks of features
     *
     * @param fl a org.biojava3.genome.parsers.gff.FeatureList
     * @param featureType the type of feature to make blocks from
     * @param attribute the name of the attribute that will make the genes
     * unique (e.g. "name", gene_id, etc)
     * @param secretedProteins a list of genes that you want to restrict your
     * analysis to
     * @return a HashMap of blocks where the keys are the names of the reference
     * sequences and the values are arrays of Integers with the consecutive
     * start and then finish co-ordinates of all the requested featureTypes
     */
    public HashMap<String, ArrayList<Integer>> getBlocks(FeatureList fl, String featureType, String attribute, Set secretedProteins)
    {
        HashMap<String, ArrayList<Integer>> blocks = new HashMap<>();
        FeatureList featTypes = fl.selectByType(featureType);
        for (FeatureI currentfi : featTypes)
        {
            if (currentfi.hasAttribute(attribute))
            {
                String geneId = currentfi.getAttribute(attribute);
                if (secretedProteins.contains(geneId))
                {
                    Location currentLoc = currentfi.location();
                    int currentStart = currentLoc.bioStart();
                    int currentEnd = currentLoc.bioEnd();
                    //System.out.println(geneId + " " + currentStart + " " + currentEnd);
                    //see if the attribute name exists
                    if (blocks.containsKey(geneId))
                    {
                        //if it does, get the HashMap and test to see if the current strand exists
                        ArrayList<Integer> al = (ArrayList<Integer>) blocks.get(geneId);
                        al.add(currentStart);
                        al.add(currentEnd);
                        blocks.put(geneId, al);
                    }
                    //if the attribute doesn't exist
                    else
                    {
                        //create a new ArrayList
                        ArrayList<Integer> al = new ArrayList<>();
                        //add the coords to the arraylist
                        al.add(currentStart);
                        al.add(currentEnd);
                        //add the new geneId to the blocks hashmap
                        blocks.put(geneId, al);
                    }
                }

            }
        }
        return blocks;
    }
    
    public double getStandardDeviation(double [] featArray)
    {
        StandardDeviation sd = new StandardDeviation();
        double stdev = sd.evaluate(featArray);
        return stdev;
    }
    
    public double getStandardErrorOfMean(double stdev, double sampleSize)
    {
        //Do this by dividing the standard deviation by the square root of the sample size. 
        double stderr = stdev/ Math.sqrt(sampleSize);
        return stderr;
    }
}

