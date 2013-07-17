/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fasta;

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
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.RNASequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;

/**
 *
 * @author ethering
 */
public class FastaTranslator
{

    public LinkedHashMap<String, ProteinSequence> translateMultiFasta(File fastaDnaFile) throws Exception
    {
        LinkedHashMap<String, DNASequence> dnaSeqs = FastaReaderHelper.readFastaDNASequence(fastaDnaFile);
        //FastaReaderHelper.readFastaDNASequence for DNA sequences
        LinkedHashMap<String, ProteinSequence> aaSeqs = new LinkedHashMap<String, ProteinSequence>();
        for (Entry<String, DNASequence> entry : dnaSeqs.entrySet())
        {
            //System.out.println(entry.getValue().getOriginalHeader() + "=" + entry.getValue().getSequenceAsString());
            RNASequence rna = entry.getValue().getRNASequence();
            ProteinSequence aa = rna.getProteinSequence();
            aaSeqs.put(entry.getKey(), aa);
        }
        return aaSeqs;
    }

    public ProteinSequence translateFasta(DNASequence dna) throws Exception
    {

        RNASequence rna = dna.getRNASequence();
        ProteinSequence aa = rna.getProteinSequence();
        System.out.println(dna.getOriginalHeader());
        System.out.println(aa);
        AccessionID id = new AccessionID(dna.getOriginalHeader());
        aa.setAccession(id);
        return aa;
    }

    public void translateCrinklerDNACountFile(File dnaCounts, File proteinCounts) throws FileNotFoundException, IOException
    {

        HashMap<String, double[]> aaMotifs = new HashMap<String, double[]>();
        ArrayList<String> aaSequences = new ArrayList<String>();
        //open the file
        BufferedReader input = new BufferedReader(new FileReader(dnaCounts));

        String line = null;
        //get the first line. It contains the dataset info
        String header = line = input.readLine();
        //split it and then count the number of array indexes to calculate the number of datasets
        String[] numDatasetsArray = line.split("\t");
        int numDatasets = numDatasetsArray.length - 1;
        /*
         * readLine is a bit quirky :
         * it returns the content of a line MINUS the newline.
         * it returns null only for the END of the stream.
         * it returns an empty String if two newlines appear in a row.
         */
        while ((line = input.readLine()) != null)
        {
            //split each line into an array. The first index will be the DNA sequence, the rest will be the counts
            //System.out.println(line);
            String[] array = line.split("\t");
            DNASequence dna = new DNASequence(array[0]);
            double[] doubleArray = new double[array.length - 1];

            //put the counts into an int array
            for (int i = 1; i < array.length; i++)
            {
                int count = Integer.parseInt(array[i]);
                doubleArray[i - 1] = count;
            }
            //make a protein sequence from the DNA sequence
            RNASequence rna = dna.getRNASequence();
            ProteinSequence aa = rna.getProteinSequence();
            //if it ends in an 'X', remove the 'X'
            if (aa.toString().endsWith("X"))
            {
                String aaString = aa.getSequenceAsString();
                String aaSubStr = aaString.substring(0, aaString.length() - 1);
                aa = new ProteinSequence(aaSubStr);
            }
            String finalAAString = aa.getSequenceAsString();
            //if the current protein seq is in the aaMotifs collection
            if (aaMotifs.containsKey(finalAAString))
            {
                //get the counts from aaMotifs and add the current counts to them
                double[] counts = aaMotifs.get(finalAAString);
                for (int i = 0; i < counts.length; i++)
                {
                    double count = counts[i];
                    count += doubleArray[i];
                    counts[i] = count;
                }
                aaMotifs.put(finalAAString, counts);
            }
            //if it's not in, add it
            else
            {
                aaMotifs.put(finalAAString, doubleArray);
                //keep track of the unique sequences
                aaSequences.add(finalAAString);
            }
        }
        input.close();
        //print out the new counts
        Iterator it = aaMotifs.entrySet().iterator();
        System.out.println("\n\nPrinting aaMotifs");
        while (it.hasNext())
        {
            //get the protein name and the counts
            Map.Entry pairs = (Map.Entry) it.next();
            String motif = (String) pairs.getKey();
            double[] doubleArray = (double[]) pairs.getValue();
            System.out.print(motif);
            for (int i = 0; i < doubleArray.length; i++)
            {
                System.out.print("\t" + doubleArray[i]);
            }
            System.out.println();
        }

        System.out.println("\n\n");
        //sort the arraylist by sequence length
        Collections.sort(aaSequences, new StringLengthComparator());
        //sort the aaMotifs map
        LinkedHashMap<String, double[]> sortedAaMotifs = new LinkedHashMap<String, double[]>();
        for (String aaseq : aaSequences)
        {
            sortedAaMotifs.put(aaseq, aaMotifs.get(aaseq));
        }
 //       Iterator it2 = sortedAaMotifs.entrySet().iterator();
//        System.out.println("\n\nprinting sortedAaMotifs");
//        while (it2.hasNext())
//        {
//            //get the protein name and the counts
//            Map.Entry pairs = (Map.Entry) it2.next();
//            String motif = (String) pairs.getKey();
//            double[] doubleArray = (double[]) pairs.getValue();
//            System.out.print(motif);
//            for (int i = 0; i < doubleArray.length; i++)
//            {
//                System.out.print("\t" + doubleArray[i]);
//            }
//            System.out.println();
//        }
//
//        System.out.println("\n\n");

        //clone the unique seqs
        Collections.sort(aaSequences, new StringLengthComparator());
        ArrayList<String> aaSeqClone = new ArrayList<String>();
        for (String str : aaSequences)
        {
            aaSeqClone.add(str);
        }
        ArrayList<String> aaSeqClone2 = new ArrayList<String>();
        for (String str : aaSequences)
        {
            aaSeqClone2.add(str);
        }
        //Collections.sort(aaSeqClone, new StringLengthComparator());
        //System.out.println("aaSequences.size = " + aaSequences.size());
        //go thru the unique seqs and get each sequence (as a charSeq)
        for (int i = 0; i < aaSequences.size(); i++)
        {
            //System.out.println(i);
            String str1 = aaSequences.get(i);
            int strLength = str1.length();
            //CharSequence cs = str1.subSequence(0, str1.length());
            //System.out.println("Looking at " + str1);
            //go through the clone
            boolean found = false;
            for (int j = i+1; j < aaSeqClone2.size(); j++)
            {
                String str2 = aaSeqClone2.get(j).substring(0, strLength);
                //System.out.println("Comparing to " + str2);

                //test to see if cs is a substring of str2
                if (str1.equalsIgnoreCase(str2))
                {
                    found = true;
                    //System.out.println(str1 + " equals "+str2);
                }
            }
            if (found == true)
            {
                //if it is, remove it from str1
                aaSeqClone.remove(str1);
                //System.out.println("Removing " + str1);

            }

        }
        //aaSequences should be a list of the longest unique sequences
        //Collections.sort(aaSequences, new StringLengthComparator());
        //System.out.println("Printing sorted aaSeqs. There should be no subseqs");
        //printArray(aaSeqClone);
        //System.out.println("\n\n");
        //finalAaMotifs will contain the final list of motif counts for the longest unique sequences
        LinkedHashMap<String, double[]> finalAaMotifs = new LinkedHashMap<String, double[]>();
        //put the 
        for (String aaseq:aaSeqClone)
        {
            finalAaMotifs.put(aaseq, aaMotifs.get(aaseq));
        }
        //iterate over the protein seqs again
        Iterator itm = sortedAaMotifs.entrySet().iterator();
        while (itm.hasNext())
        {
            //get the protein name and the counts
            Map.Entry pairs = (Map.Entry) itm.next();
            String motif = (String) pairs.getKey();
            double[] doubleArray = (double[]) pairs.getValue();

            CharSequence cs = motif.subSequence(0, motif.length());

            //count how many times the current motif occurs as a substring in the longest unique seq list
            int noTimesFound = 0;
            for (int j = 0; j < aaSeqClone.size(); j++)
            {
                String longestSeq = aaSeqClone.get(j);
                if (!longestSeq.equalsIgnoreCase(cs.toString()) &&  longestSeq.contains(cs))
                {
                    noTimesFound++;
                }
            }


            //go thru the longest motifs again
            for (int j = 0; j < aaSeqClone.size(); j++)
            {
                String longestSeq = aaSeqClone.get(j);

                //if the current protein occurs as a substring in longestSeq
                if (!longestSeq.equalsIgnoreCase(cs.toString()) &&  longestSeq.contains(cs))
                {
                    //check to see if it exists in the final list
                    if (finalAaMotifs.containsKey(longestSeq))
                    {
                        //if it does, get the counts from the final list, then take the counts from the current protein
                        //divide each count by the number of times it's found, round the number and add it to the final list

                        double[] counts = finalAaMotifs.get(longestSeq);
                        for (int i = 0; i < counts.length; i++)
                        {
                            double count = counts[i];
                            double countToDist = doubleArray[i];
                            double share = countToDist / noTimesFound;
                            count += share;
                            counts[i] = count;
                        }
                        finalAaMotifs.put(longestSeq, counts);
                    }
                    else
                    {
                        System.out.println("Ooops this should exist "+longestSeq);
                    }
                }
            }
        }
        Writer output = new BufferedWriter(new FileWriter(proteinCounts));

        System.out.println(header);
        output.write(header + "\n");
        Iterator itm2 = finalAaMotifs.entrySet().iterator();
        int[][] tally = new int[finalAaMotifs.entrySet().size()][numDatasets];
        int datasetNo = 0;
        while (itm2.hasNext())
        {
            Map.Entry pairs = (Map.Entry) itm2.next();
            String motif = (String) pairs.getKey();
            double[] doubleArray = (double[]) pairs.getValue();
            System.out.print(motif);
            output.write(motif);
            for (int i = 0; i < doubleArray.length; i++)
            {
                System.out.print("\t" + Math.round(doubleArray[i]));
                output.write("\t" + Math.round(doubleArray[i]));
                if (Math.round(doubleArray[i]) > 0)
                {
                    tally[datasetNo][i] = 1;
                }
                else
                {
                    tally[datasetNo][i] = 0;
                }
            }
            System.out.println();
            output.write("\n");
            datasetNo++;
        }
        output.close();

        if (numDatasets == 2)
        {
            int inBoth = 0;
            int in1not2 = 0;
            int in2not1 = 0;
            for (int i = 0; i < tally.length; i++)
            {
                if (tally[i][0] == 1 && tally[i][1] == 1)
                {
                    inBoth++;
                }
                else if (tally[i][0] == 1 && tally[i][1] == 0)
                {
                    in1not2++;
                }
                else if (tally[i][0] == 0 && tally[i][1] == 1)
                {
                    in2not1++;
                }
                else
                {
                    System.out.println("Oppps something has gone awry!");
                }
            }
            System.out.println("In both = " + inBoth + ". In 1 only " + in1not2 + ". In 2 only " + in2not1);
        }
        if (numDatasets == 3)
        {
            int inAll = 0;
            int in1only = 0;
            int in2only = 0;
            int in3only = 0;
            int in1and2not3 = 0;
            int in1and3not2 = 0;
            int in2and3not1 = 0;
            for (int i = 0; i < tally.length; i++)
            {
                if (tally[i][0] == 1 && tally[i][1] == 1 && tally[i][2] == 1)
                {
                    inAll++;
                }
                else if (tally[i][0] == 1 && tally[i][1] == 1 && tally[i][2] == 0)
                {
                    in1and2not3++;
                }
                else if (tally[i][0] == 1 && tally[i][1] == 0 && tally[i][2] == 1)
                {
                    in1and3not2++;
                }
                else if (tally[i][0] == 0 && tally[i][1] == 1 && tally[i][2] == 1)
                {
                    in2and3not1++;
                }
                else if (tally[i][0] == 1 && tally[i][1] == 0 && tally[i][2] == 0)
                {
                    in1only++;
                }
                else if (tally[i][0] == 0 && tally[i][1] == 1 && tally[i][2] == 0)
                {
                    in2only++;
                }
                else if (tally[i][0] == 0 && tally[i][1] == 0 && tally[i][2] == 1)
                {
                    in3only++;
                }
                else
                {
                    System.out.println("Oppps something has gone awry!");
                    System.out.println(tally[i][0] + " " + tally[i][1] + " " + tally[i][2]);
                }
            }
            System.out.println("In all = " + inAll + ". in1and2not3 " + in1and2not3 + ". In in1and3not2 " + in1and3not2 + " in2and3not1 " + in2and3not1);
            System.out.println("in1only " + in1only + " in2only " + in2only + " in3only " + in3only);
        }

    }

    private void printArray(ArrayList<String> array)
    {
        for (String str : array)
        {
            System.out.println(str);
        }
    }

    public class StringLengthComparator implements Comparator<String>
    {

        @Override
        public int compare(String o1, String o2)
        {
            if (o1.length() > o2.length())
            {
                return 1;
            }
            else if (o1.length() < o2.length())
            {
                return -1;
            }
            return o1.compareTo(o2);
        }
    }

    public static void main(String[] args)
    {
        File f = new File("/Users/ethering/temp/crinkler/ec_GTGCTGGTGG[TCA][ATCG][TG][TC][ACTG]CC_dna_counts_distributed.txt");
        File out = new File("/Users/ethering/temp/crinkler/ec_GTGCTGGTGG[TCA][ATCG][TG][TC][ACTG]CC_aa_counts_distributed.txt");
        //File f = new File("/Users/ethering/temp/crinkler/usuknl_GTGCTGGTGG[TC]G[TG]T[TG]CC_dna_counts_tab.txt");
        //File out = new File("/Users/ethering/temp/crinkler/usuknl_GTGCTGGTGG[TC]G[TG]T[TG]CC_aa_counts_tab.txt");
        FastaTranslator ft = new FastaTranslator();
        try
        {
            ft.translateCrinklerDNACountFile(f, out);


        }
        catch (Exception ex)
        {
            Logger.getLogger(FastaTranslator.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
    }
}
