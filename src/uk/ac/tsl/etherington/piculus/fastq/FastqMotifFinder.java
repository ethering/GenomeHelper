/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fastq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import org.biojava3.core.sequence.*;
import org.jtr.transliterate.CharacterParseException;
import uk.ac.tsl.etherington.piculus.fasta.FastaMotifFinder;

/**
 * A class to identify motifs in FASTQ formatted reads
 *
 * @author ethering
 */
public class FastqMotifFinder
{

    /**
     * Looks through a DNA fastq file for a given motif. The motif can contain
     * regular expressions, such as ATG[CG]G[AT], which will search for ATGCGA,
     * ATGCGT, ATGGGA and ATGGGT. Refer to java.util.regex.Pattern for more
     * details. The output is two files, one a tab-delimited file of DNA motif
     * counts, the other a count of the corresponding ammino acid translations.
     * The reverse complement is also examined for the motif.
     *
     * @param fastqFile the input file of DNA sequences in fastq format
     * @param strPattern the pattern to be searched for
     * @param motifCounts a tab-delimited file of the occurrence of each pattern
     * @param proteinCounts a tab-delimited file of the occurrence of each
     * pattern when translated into ammino acids
     * @param minCount the minimum number of times a motif must be found to be
     * included in the results
     * @throws IOException
     * @throws CharacterParseException
     */
    public void findMatches(File fastqFile, String strPattern, File motifCounts, File proteinCounts, int minCount) throws IOException, CharacterParseException
    {
        DNASequence dna;
        RNASequence rna;
        ProteinSequence aa;
        //hash maps to count the dna and aa motifs
        HashMap<String, Integer> dnaMotifs = new HashMap<>();
        HashMap<String, Integer> aaMotifs = new HashMap<>();
        ValueComparator dbvc = new ValueComparator(dnaMotifs);
        ValueComparator abvc = new ValueComparator(aaMotifs);
        TreeMap<String, Integer> sorted_dna_map = new TreeMap<>(dbvc);
        TreeMap<String, Integer> sorted_aa_map = new TreeMap<>(abvc);
        //a compiler for out
        Pattern pattern = Pattern.compile(strPattern.toLowerCase());
        //Writers for the dna motif counts and protein motif counts

        Writer dnaMotifWriter = new BufferedWriter(new FileWriter(motifCounts));
        Writer aaMotifWriter = new BufferedWriter(new FileWriter(proteinCounts));
        //open the fastq file and loop thru it
        FastqReader fq = new FastqReader(fastqFile);
        Iterator it = fq.iterator();
        int matchCounter = 0;
        int matchesOnPlus = 0;
        int matchesOnMinus = 0;
        while (it.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it.next();
            String seq = seqRecord.getReadString().toLowerCase();
            //store the sequence in a matcher
            Matcher matcher = pattern.matcher(seq);

            //if a match is found
            while (matcher.find())
            {
                matchesOnPlus++;
                //get the occurance of the first match and create a subsequence from there
                String subs = seq.substring(matcher.start(), seq.length());
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

            //get the reverse complement
            String revSeq = FastaMotifFinder.revcom(seq).toLowerCase();
            Matcher revMatcher = pattern.matcher(revSeq);
            //if a reverse complement match is found is found
            while (revMatcher.find())
            {
                //get the occurance of the first match and create a subsequence from there
                String subs = revSeq.substring(revMatcher.start(), revSeq.length());
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

    /**
     *
     * @param fastqIn single interlaced or concatenated paired-end fastq
     * sequence file
     * @param strPattern the pattern for which you are searching
     * @param fastqOut the file of paired-end sequences which have a match in
     * one of the pairs
     * @throws IOException
     * @throws CharacterParseException
     */
    public void getPEFastqReadsFromMotif(File fastqIn, String strPattern, File fastqOut) throws IOException, CharacterParseException
    {
        //hash maps to count the dna and aa motifs
        HashSet<String> matches = new HashSet<>();

        //a compiler for out
        Pattern pattern = Pattern.compile(strPattern.toLowerCase());
        //Writers for matching reads

        //open the fastq file and loop thru it
        FastqReader fq = new FastqReader(fastqIn);
        Iterator it = fq.iterator();
        int matchCounter = 0;
        int matchesOnPlus = 0;
        int matchesOnMinus = 0;
        while (it.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it.next();
            String seq = seqRecord.getReadString().toLowerCase();
            //store the sequence in a matcher
            Matcher matcher = pattern.matcher(seq);

            //if a match is found
            while (matcher.find())
            {
                String read = seqRecord.getReadHeader();
                int spaceIndex = read.indexOf(" ");
                String readName = read.substring(0, spaceIndex);
                matches.add(readName);
                matchCounter++;
                matchesOnPlus++;
            }

            //get the reverse complement
            String revSeq = FastaMotifFinder.revcom(seq).toLowerCase();
            Matcher revMatcher = pattern.matcher(revSeq);
            //if a reverse complement match is found is found
            while (revMatcher.find())
            {
                String read = seqRecord.getReadHeader();
                //System.out.println(read);
                int spaceIndex = read.indexOf(" ");
                String readName = read.substring(0, spaceIndex);
                matches.add(readName);
                matchCounter++;
                matchesOnMinus++;
            }

        }
        fq.close();


        FastqWriterFactory writer = new FastqWriterFactory();
        //for the matching paired-end reads
        FastqWriter fqWriter = writer.newWriter(fastqOut);
        //fq = new FastqReader(fastqIn);
        FastqReader fq2 = new FastqReader(fastqIn);
        Iterator it2 = fq2.iterator();
        while (it2.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it2.next();
            String read = seqRecord.getReadHeader();
            int spaceIndex = read.indexOf(" ");
            String readName = read.substring(0, spaceIndex);
            if (matches.contains(readName))
            {
                fqWriter.write(seqRecord);
            }
        }
        System.out.println("Found " + matchCounter + " matches");
        System.out.println("Found " + matchesOnPlus + " matches on plus");
        System.out.println("Found " + matchesOnMinus + " matches on minus");

        fqWriter.close();
    }

    class ValueComparator implements Comparator<String>
    {

        Map<String, Integer> base;

        public ValueComparator(Map<String, Integer> base)
        {
            this.base = base;
        }

        // Note: this comparator imposes orderings that are inconsistent with equals.    
        @Override
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

    public class StringLengthComparator implements Comparator<String>
    {

        @Override
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