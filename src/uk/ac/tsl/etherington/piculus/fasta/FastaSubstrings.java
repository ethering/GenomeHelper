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
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.jtr.transliterate.CharacterParseException;
import uk.ac.tsl.etherington.piculus.fastq.FastqMotifFinder;

/**
 *
 * @author ethering
 */
public class FastaSubstrings
{
    /**
     * Extract a single fasta sequence from a multi-fasta file
     * @param fastaIn the input multi-fasta file
     * @param outfile the output file for the requested fasta sequence
     * @param seqid the name of the sequence to be extracted from the multi-fasta file
     * @throws FileNotFoundException
     * @throws Exception 
     */
    public void seqFromCommandLine(File fastaIn, File outfile, String seqid) throws FileNotFoundException, Exception
    {
        boolean seqFound = false;
        Writer out = new BufferedWriter(new FileWriter(outfile));

        BufferedReader br = new BufferedReader(new FileReader(fastaIn));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                alpha.getTokenization("token"), ns);
        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            String id = rec.getName();
            if (id.equalsIgnoreCase(seqid))
            {
                String seq = rec.seqString();
                out.write(">" + id + "\n");
                out.write(seq + "\n");
                seqFound = true;
                break;
            }
        }
        out.close();
        if (seqFound == false)
        {
            System.out.println("Couldn't find a sequence with that name");
        }
    }
    /**
     * Extract a sub-sequence of a single fasta sequence from a multi-fasta file
     * @param fastaIn the input multi-fasta file
     * @param outfile the output file for the requested fasta sequence
     * @param seqid the name of the sequence to be extracted from the multi-fasta file
     * @param start the 1-based inclusive start position of the requested sub-sequence
     * @param end the 1-based inclusive end position of the requested sub-sequence
     * @throws FileNotFoundException
     * @throws Exception 
     */
    public void seqFromCommandLine(File fastaIn, File outfile, String seqid, int start, int end) throws FileNotFoundException, Exception
    {
        boolean seqFound = false;
        int subseqLength = end - start+1;
        //the String.substring method is zero-based and is inclusive for the start and exclusive for the end.
        //Changing it to 1-based means that we need to subtract 1 from the start index, but leave the end index alone.
        start--;
        System.out.println("Requested subsequence length is " + subseqLength);
        Writer out = new BufferedWriter(new FileWriter(outfile));

        BufferedReader br = new BufferedReader(new FileReader(fastaIn));
        Alphabet alpha = AlphabetManager.alphabetForName("DNA");
        SimpleNamespace ns = new SimpleNamespace("biojava");

        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
                alpha.getTokenization("token"), ns);
        while (iterator.hasNext())
        {
            RichSequence rec = iterator.nextRichSequence();
            String id = rec.getName();
            if (id.equalsIgnoreCase(seqid))
            {
                String seq = rec.seqString();
                String subseq = seq.substring(start, end);
                subseqLength = subseq.length();
                System.out.println("Provided subsequence length is " + subseqLength);

                out.write(">" + id + "\n");
                out.write(subseq + "\n");
                seqFound = true;
                break;
            }
        }
        out.close();
        if (seqFound == false)
        {
            System.out.println("Couldn't find a sequence with that name");
        }
    }

    /**
     *Takes a list of fasta files and finds the longest subsequences that are common to all files
     * @param al an ArrayList of fasta file paths
     * @param outfile the list of longest subsequences found in all files
     * @throws FileNotFoundException
     * @throws Exception
     */
    public void findLongestCommonSequences(ArrayList<String> al, File outfile) throws FileNotFoundException, Exception
    {
        Writer out = new BufferedWriter(new FileWriter(outfile));
        int numberOfDataSets = al.size();
        //this will contain the String fileName and LinkedHashSet of the ordered-by-length subseqs
        HashMap<String, LinkedHashSet> fastaSeqs = new HashMap<>();
        //for each fasta file
        for (String fastaIn : al)
        {
            System.out.println(fastaIn);
            ArrayList<String> substr = getAllSubStrings(fastaIn);
            //sort the arrayList by length
            Collections.sort(substr, new StringLengthComparator());
            //and create a unique set of subseqs
            LinkedHashSet<String> uniqueSubSeqs = new LinkedHashSet<>(substr);
            //add it to the top level 'store' - fastaSeqs
            fastaSeqs.put(fastaIn, uniqueSubSeqs);
        }

        //go through the fastaSeqs store
        Iterator it = fastaSeqs.entrySet().iterator();
        while (it.hasNext())
        {
            Map.Entry currentPairs = (Map.Entry) it.next();
            //get the name of the current file and the subseq contents
            String currentFile = (String) currentPairs.getKey();
            //get the subsequences from the current file
            LinkedHashSet<String> currentFastaSubs = (LinkedHashSet<String>) currentPairs.getValue();
            System.out.println("Current File: " + currentFile);
            //for each subsequence
            for (String fastaSub : currentFastaSubs)
            {
                int subSeqsFound = 0;
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
                        //if it's not, get the collection and check to see if the current subseq is contained in the collection in
                        LinkedHashSet<String> fastaSubsInLoop = (LinkedHashSet<String>) pairsInLoop.getValue();
                        //if the subseq is found in another file increment the count
                        if (fastaSubsInLoop.contains(fastaSub))
                        {
                            subSeqsFound++;
                        }
                    }
                }
                //if the subseq is found in all the other datasets, write it to file
                if (subSeqsFound == (numberOfDataSets - 1))
                {
                    out.write(fastaSub + "\n");
                }
            }
        }
    }
    /**
     * 
     * @param fastaIn the path to a fasta file
     * @return all possible substrings (including reverse complements) found within the fasta file
     * @throws CharacterParseException
     * @throws Exception 
     */
    public ArrayList<String> getAllSubStrings(String fastaIn) throws CharacterParseException, Exception
    {
        ArrayList<String> substr = new ArrayList<>();
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
            String revcom = FastaMotifFinder.revcom(fastaSeq);
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
        return substr;
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
