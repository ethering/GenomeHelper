/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package piculus.fasta;

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
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import org.biojava3.core.sequence.*;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.transcription.Frame;
import org.biojava3.data.sequence.FastaSequence;
import org.jtr.transliterate.CharacterParseException;
import org.jtr.transliterate.CharacterReplacer;
import org.jtr.transliterate.Perl5Parser;

/**
 *
 * @author ethering
 */
public class FastaMotifFinder
{

    public void findMatches(File fastaFile, String strPattern, File motifCounts, File proteinCounts, int minCount) throws IOException, CharacterParseException, Exception
    {
        DNASequence dna;
        RNASequence rna;
        ProteinSequence aa;
        //hash maps to count the dna and aa motifs
        HashMap<String, Integer> dnaMotifs = new HashMap<String, Integer>();
        HashMap<String, Integer> aaMotifs = new HashMap<String, Integer>();
        ValueComparator dbvc = new ValueComparator(dnaMotifs);
        ValueComparator abvc = new ValueComparator(aaMotifs);
        TreeMap<String, Integer> sorted_dna_map = new TreeMap<String, Integer>(dbvc);
        TreeMap<String, Integer> sorted_aa_map = new TreeMap<String, Integer>(abvc);
        //a compiler for out
        Pattern pattern = Pattern.compile(strPattern);
        //Writers for the dna motif counts and protein motif counts

        Writer dnaMotifWriter = new BufferedWriter(new FileWriter(motifCounts));
        Writer aaMotifWriter = new BufferedWriter(new FileWriter(proteinCounts));
        //open the fastq file and loop thru it
        LinkedHashMap<String, DNASequence> dnaSeqs = FastaReaderHelper.readFastaDNASequence(fastaFile);
        // Get a set of the entries
        Set seqs = dnaSeqs.entrySet();
        // Get an iterator
        Iterator i = seqs.iterator();
        int matchCounter = 0;
        int matchesOnPlus = 0;
        int matchesOnMinus = 0;
        // Display elements
        while (i.hasNext())
        {
            Map.Entry me = (Map.Entry) i.next();
            dna = (DNASequence) me.getValue();
            String seq = dna.getSequenceAsString();   
            //store the sequence in a matcher
            Matcher matcher = pattern.matcher(seq);

            //if a match is found
            if (matcher.find())
            {
                matchesOnPlus++;
                //get the occurance of the first match and create a subsequence from there
                String subs = seq.substring(matcher.start(), matcher.end());
                matchCounter++;
                //if the dna sequence alread exists, increment it by one
                if (dnaMotifs.containsKey(subs))
                {
                    Integer count = dnaMotifs.get(subs);
                    count = count + 1;
                    dnaMotifs.put(subs, count);
                }
                //else put in a new entry
                else
                {
                    dnaMotifs.put(subs, 1);
                }
            }
            else
            {
                //get the reverse complement
                String revSeq = revcom(seq);
                Matcher revMatcher = pattern.matcher(revSeq);
                //if a reverse complement match is found is found
                if (revMatcher.find())
                {
                    //get the occurance of the first match and create a subsequence from there
                    String subs = revSeq.substring(revMatcher.start(), revMatcher.end());
                    matchCounter++;
                    matchesOnMinus++;
                    //if the dna sequence alread exists, increment it by one
                    if (dnaMotifs.containsKey(subs))
                    {
                        Integer count = dnaMotifs.get(subs);
                        count = count + 1;
                        dnaMotifs.put(subs, count);
                    }
                    //else put in a new entry
                    else
                    {
                        dnaMotifs.put(subs, 1);
                    }
                }
            }
        }
        System.out.println("Found " + matchCounter + " matches");
        System.out.println("Found " + matchesOnPlus + " matches on plus");
        System.out.println("Found " + matchesOnMinus + " matches on minus");
        sorted_dna_map.putAll(dnaMotifs);
        Iterator itm = sorted_dna_map.entrySet().iterator();
        while (itm.hasNext())
        {
            Map.Entry pairs = (Map.Entry) itm.next();
            String motif = (String) pairs.getKey();
            Integer count = (Integer) pairs.getValue();
            if (count >= minCount)
            {
                dnaMotifWriter.write(motif + "\t");
                dnaMotifWriter.write(count + "\n");
                dna = new DNASequence(motif);
                rna = dna.getRNASequence();
                aa = rna.getProteinSequence();
                String aaString = aa.toString();
                if (aaMotifs.containsKey(aaString))
                {
                    Integer aacount = aaMotifs.get(aaString);
                    aacount = aacount + 1;
                    aaMotifs.put(aaString, aacount);
                }
                else
                {
                    aaMotifs.put(aaString, 1);
                }
            }
        }
        dnaMotifWriter.close();

        sorted_aa_map.putAll(aaMotifs);
        Iterator ita = sorted_aa_map.entrySet().iterator();
        while (ita.hasNext())
        {
            Map.Entry pairs = (Map.Entry) ita.next();
            String motif = (String) pairs.getKey();
            Integer count = (Integer) pairs.getValue();
            if (count >= minCount)
            {
                aaMotifWriter.write(motif + "\t");
                aaMotifWriter.write(count + "\n");
            }
        }
        aaMotifWriter.close();
    }

   
    public static void main(String[] args)
    {
            FastaMotifFinder fmf = new FastaMotifFinder();
            File fqfile = new File("/Users/ethering/temp/solexa/uk_us/US_from_UK_US_shared.fastq");
            File fqOut = new File("/Users/ethering/temp/solexa/uk_us/US_from_UK_US_shared_found.fastq");

            String pattern = "GTGCTGGTGG[TC]G[TG]T[TG]CC";

            
    }

    class ValueComparator implements Comparator<String>
    {

        Map<String, Integer> base;

        public ValueComparator(Map<String, Integer> base)
        {
            this.base = base;
        }

        // Note: this comparator imposes orderings that are inconsistent with equals.    
        public int compare(String a, String b)
        {
            if (base.get(a) >= base.get(b))
            {
                return -1;
            }
            else
            {
                return 1;
            } // returning 0 would merge keys
        }
    }

    public static String revcom(String seq) throws CharacterParseException
    {
        //System.out.println("seqIn: "+seq);
        CharacterReplacer cReplacer = Perl5Parser.makeReplacer("tr/ATCG/TAGC/");
        seq = reverseString(seq);
        String newSeq = cReplacer.doReplacement(seq);
        //System.out.println("seqOut: "+newSeq);
        return newSeq;
    }

    public static String reverseString(String s)
    {
        return new StringBuffer(s).reverse().toString();


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