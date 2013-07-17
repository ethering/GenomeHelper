/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fasta;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Map.Entry;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import uk.ac.tsl.etherington.piculus.fastq.FastqMotifFinder;

/**
 *
 * @author ethering
 */
public class FastaSubstrings
{
    //al will contain an list of paths to fasta seqs

    public void findLongestCommonSequences(ArrayList<String> al, File outfile) throws FileNotFoundException, Exception
    {
        Writer out = new BufferedWriter(new FileWriter(outfile));
        int numberOfDataSets = al.size();
        //this will contain the String fileName and LinkedHashSet of the ordered-by-length subseqs
        HashMap<String, LinkedHashSet> fastaSeqs = new HashMap<String, LinkedHashSet>();
        //for each fasta file
        for (String fastaIn : al)
        {
            System.out.println(fastaIn);
            ArrayList<String> substr = new ArrayList<String>();
            //open the file
            LinkedHashMap<String, DNASequence> a = FastaReaderHelper.readFastaDNASequence(new File(fastaIn));
            //read in the dna seqs
            for (Entry<String, DNASequence> entry : a.entrySet())
            {
                //System.out.println(entry.getValue().getOriginalHeader() + "=" + entry.getValue().getSequenceAsString());
                String fastaSeq = entry.getValue().getSequenceAsString();
                //starting with the first nt, create a subseq of the length of the seq, then reduce the lenght of the seq by one 
                //from the end of the seq. When that's finished, start at the second nt and do the same until we have every possible subseq
                for (int i = 0; i < fastaSeq.length(); i++)
                {
                    int start = i;
                    for (int j = fastaSeq.length(); j > i; j--)
                    {
                        int end = j;
                        String str = fastaSeq.substring(start, end);
                        //System.out.println(str);
                        //add the string to an arrayList
                        substr.add(str);

                    }
                }
                //get the reverse complement and do the same
                String revcom = FastqMotifFinder.revcom(fastaSeq);
                for (int i = 0; i < revcom.length(); i++)
                {
                    int start = i;
                    for (int j = revcom.length(); j > i; j--)
                    {
                        int end = j;
                        String str = revcom.substring(start, end);
                        //System.out.println(str);
                        substr.add(str);

                    }
                }
            }
            //sort the arrayList by length
            Collections.sort(substr, new StringLengthComparator());
            //and create a unique set of subseqs
            LinkedHashSet<String> lhs = new LinkedHashSet<String>(substr);
            //add it to the top level 'store' - fastaSeqs
            fastaSeqs.put(fastaIn, lhs);


        }

        //go through the fastaSeqs store
        Iterator it = fastaSeqs.entrySet().iterator();
        while (it.hasNext())
        {
            Map.Entry currentPairs = (Map.Entry) it.next();
            //get the name of the current file and the subseq contents
            String currentFile = (String) currentPairs.getKey();
            LinkedHashSet<String> currentFastaSubs = (LinkedHashSet<String>) currentPairs.getValue();
            System.out.println("Current File: " + currentFile);
            //for each subsequence
            for (String fastaSub : currentFastaSubs)
            {
                int subSeqsFound = 0;
                //System.out.println(fastaSub);
                //loop thru the other files
                Iterator it2 = fastaSeqs.entrySet().iterator();
                while (it2.hasNext())
                {
                    Map.Entry pairsInLoop = (Map.Entry) it2.next();
                    //get the name of the current file in the loop
                    String fileInLoop = (String) pairsInLoop.getKey();
                    //if it's the same as the file we're looking at, go on to the next one
                    if (!fileInLoop.equalsIgnoreCase(currentFile))
                    {
                        LinkedHashSet<String> fastaSubsInLoop = (LinkedHashSet<String>) pairsInLoop.getValue();
                        if (fastaSubsInLoop.contains(fastaSub))
                        {
                            subSeqsFound++;
                        }
                    }
                    //if it's not, get the collection and check to see if the current subseq is contained in the collection in
                    
                }
                //if the subseq is found in all the other datasets
                if (subSeqsFound == (numberOfDataSets - 1))
                {
                    out.write(fastaSub+"\n");
                }

            }
        }
    }

    public class StringLengthComparator implements Comparator<String>
    {

        public int compare(String o1, String o2)
        {
            if (o1.length() > o2.length())
            {
                return -1;
            }
            else if (o1.length() < o2.length())
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }
    }

    
}
