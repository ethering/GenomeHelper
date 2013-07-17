/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.fastq;

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
import java.util.Map;
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
import org.biojava3.core.sequence.transcription.Frame;
import org.biojava3.data.sequence.FastaSequence;
import org.jtr.transliterate.CharacterParseException;
import org.jtr.transliterate.CharacterReplacer;
import org.jtr.transliterate.Perl5Parser;

/**
 *
 * @author ethering
 */
public class FastqMotifFinder
{

    public void findMatches(File fastqFile, String strPattern, File motifCounts, File proteinCounts, int minCount) throws IOException, CharacterParseException
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
        FastqReader fq = new FastqReader(fastqFile);
        Iterator it = fq.iterator();
        int matchCounter = 0;
        int matchesOnPlus = 0;
        int matchesOnMinus = 0;
        while (it.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it.next();
            String seq = seqRecord.getReadString();
            //store the sequence in a matcher
            Matcher matcher = pattern.matcher(seq);

            //if a match is found
            if (matcher.find())
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
            else
            {
                //get the reverse complement
                String revSeq = revcom(seq);
                Matcher revMatcher = pattern.matcher(revSeq);
                //if a reverse complement match is found is found
                if (revMatcher.find())
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

    public void findMatchesDNAOnly(File fastqFile, String strPattern, File motifCounts, int minCount) throws IOException, CharacterParseException
    {
        //hash maps to count the dna and aa motifs
        HashMap<String, Integer> dnaMotifs = new HashMap<String, Integer>();

        ValueComparator dbvc = new ValueComparator(dnaMotifs);

        TreeMap<String, Integer> sorted_dna_map = new TreeMap<String, Integer>(dbvc);

        //a compiler for out
        Pattern pattern = Pattern.compile(strPattern);
        //Writers for the dna motif counts and protein motif counts

        Writer dnaMotifWriter = new BufferedWriter(new FileWriter(motifCounts));

        //open the fastq file and loop thru it
        FastqReader fq = new FastqReader(fastqFile);
        Iterator it = fq.iterator();
        int matchCounter = 0;
        int matchesOnPlus = 0;
        int matchesOnMinus = 0;
        while (it.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it.next();
            String seq = seqRecord.getReadString();
            //store the sequence in a matcher
            Matcher matcher = pattern.matcher(seq);

            //if a match is found
            if (matcher.find())
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
            else
            {
                //get the reverse complement
                String revSeq = revcom(seq);
                Matcher revMatcher = pattern.matcher(revSeq);
                //if a reverse complement match is found is found
                if (revMatcher.find())
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
            }
        }
        dnaMotifWriter.close();
    }
    
    

    public void getPEFastqReadsFromMotif(File fastqIn, String strPattern, File fastqOut) throws IOException, CharacterParseException
    {

        //hash maps to count the dna and aa motifs
        HashSet<String> matches = new HashSet<String>();

        //a compiler for out
        Pattern pattern = Pattern.compile(strPattern);
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
            String seq = seqRecord.getReadString();
            //store the sequence in a matcher
            Matcher matcher = pattern.matcher(seq);

            //if a match is found
            if (matcher.find())
            {
                String read = seqRecord.getReadHeader();
                int hashIndex = read.indexOf("#");
                String readName = read.substring(0, hashIndex);
                matches.add(readName);
                matchCounter++;
                matchesOnPlus++;
            }
            else
            {
                //get the reverse complement
                String revSeq = revcom(seq);
                Matcher revMatcher = pattern.matcher(revSeq);
                //if a reverse complement match is found is found
                if (revMatcher.find())
                {
                    String read = seqRecord.getReadHeader();
                    //System.out.println(read);
                    int hashIndex = read.indexOf("#");
                    String readName = read.substring(0, hashIndex);
                    matches.add(readName);
                    matchCounter++;
                    matchesOnMinus++;
                }
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
            int hashIndex = read.indexOf("#");
            String readName = read.substring(0, hashIndex);
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

    public void getReadNamesOfProteinMatches(File fastqFile, String strPattern, File motifsToFlag, File genotypeMotifs, File readNames) throws IOException, CharacterParseException
    {
        //include the individual aa_adjusted_count.txt files for the geneotype
        DNASequence dna;
        RNASequence rna;
        ProteinSequence aa;
        ArrayList<String> superMotifs = new ArrayList<String>();
        ArrayList<String> subMotifs = new ArrayList<String>();
        ArrayList<String> readNamesArray = new ArrayList<String>();
        BufferedReader motifReader = new BufferedReader(new FileReader(motifsToFlag));
        String motif;
        //read the supersequence motifs and put them into the motifs AL
        while ((motif = motifReader.readLine()) != null)
        {
            superMotifs.add(motif);
        }
        System.out.println("Super Motifs: " + superMotifs);
        //read in the individual geneotype motifs and find which motif matches, or is the longest subsequence for each supersequence
        BufferedReader genoMotifReader = new BufferedReader(new FileReader(genotypeMotifs));
        String motifLine;
        while ((motifLine = genoMotifReader.readLine()) != null)
        {
            String[] mofitString = motifLine.split("\\t");
            subMotifs.add(mofitString[0]);
        }

        Collections.sort(subMotifs, new StringLengthComparator());

        //if it's a subsequence, replace the supersequence with it
        ArrayList<String> newSuperMotifs = new ArrayList<String>();
        for (String superStr : superMotifs)
        {
            for (String subStr : subMotifs)
            {
                CharSequence cs = subStr.subSequence(0, subStr.length());
                if (superStr.contains(cs))
                {
                    newSuperMotifs.add(subStr);
                    break;
                }
            }
        }

        System.out.println("New Super Motifs: " + newSuperMotifs);
        //a compiler for out
        Pattern pattern = Pattern.compile(strPattern);
        //Writers for the dna motif counts and protein motif counts
        Writer readNameWriter = new BufferedWriter(new FileWriter(readNames));
        //open the fastq file and loop thru it
        FastqReader fq = new FastqReader(fastqFile);
        Iterator it = fq.iterator();
        int matchCounter = 0;

        while (it.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it.next();
            String seq = seqRecord.getReadString();
            //store the sequence in a matcher
            Matcher matcher = pattern.matcher(seq);

            //if a match is found
            if (matcher.find())
            {
                //get the occurance of the first match and create a subsequence from there
                String subs = seq.substring(matcher.start(), seq.length());
                //translate the sequence
                dna = new DNASequence(subs);
                rna = dna.getRNASequence();
                aa = rna.getProteinSequence();
                String aaString = aa.toString();

                //if the dna sequence alread exists, increment it by one
                if (newSuperMotifs.contains(aaString))
                {
                    readNamesArray.add(seqRecord.getReadHeader());
                    matchCounter++;
                }
            }

            //get the reverse complement
            String revSeq = revcom(seq);
            Matcher revMatcher = pattern.matcher(revSeq);
            //if a reverse complement match is found is found
            if (revMatcher.find())
            {
                //get the occurance of the first match and create a subsequence from there
                String subs = revSeq.substring(revMatcher.start(), revSeq.length());
                //translate the sequence
                dna = new DNASequence(subs);
                rna = dna.getRNASequence();
                aa = rna.getProteinSequence();
                String aaString = aa.toString();

                if (newSuperMotifs.contains(aaString))
                {
                    readNamesArray.add(seqRecord.getReadHeader());
                    matchCounter++;
                }
            }
        }

        System.out.println("Found " + matchCounter + " matches");
        Iterator itm = readNamesArray.iterator();

        while (itm.hasNext())
        {
            String read = (String) itm.next();
            readNameWriter.write(read + "\n");
        }

        readNameWriter.close();
    }

    public void getReadNamesOfDNAMatches(File fastqFile, String strPattern, File motifsToFlag, File genotypeMotifs, File readNames) throws IOException, CharacterParseException
    {
        //include the individual aa_adjusted_count.txt files for the geneotype
        ArrayList<String> superMotifs = new ArrayList<String>();
        ArrayList<String> subMotifs = new ArrayList<String>();
        ArrayList<String> readNamesArray = new ArrayList<String>();
        BufferedReader motifReader = new BufferedReader(new FileReader(motifsToFlag));
        String motif;
        //read the supersequence motifs and put them into the motifs AL
        while ((motif = motifReader.readLine()) != null)
        {
            superMotifs.add(motif);
        }
        System.out.println("Super Motifs: " + superMotifs);
        //read in the individual geneotype motifs and find which motif matches, or is the longest subsequence for each supersequence
        BufferedReader genoMotifReader = new BufferedReader(new FileReader(genotypeMotifs));
        String motifLine;
        while ((motifLine = genoMotifReader.readLine()) != null)
        {
            String[] mofitString = motifLine.split("\\t");
            subMotifs.add(mofitString[0]);
        }

        Collections.sort(subMotifs, new StringLengthComparator());

        //if it's a subsequence, replace the supersequence with it
        ArrayList<String> newSuperMotifs = new ArrayList<String>();
        for (String superStr : superMotifs)
        {
            for (String subStr : subMotifs)
            {
                CharSequence cs = subStr.subSequence(0, subStr.length());
                if (superStr.contains(cs))
                {
                    newSuperMotifs.add(subStr);
                    break;
                }
            }
        }

        System.out.println("New Super Motifs: " + newSuperMotifs);
        //a compiler for out
        Pattern pattern = Pattern.compile(strPattern);
        //Writers for the dna motif counts and protein motif counts
        Writer readNameWriter = new BufferedWriter(new FileWriter(readNames));
        //open the fastq file and loop thru it
        FastqReader fq = new FastqReader(fastqFile);
        Iterator it = fq.iterator();
        int matchCounter = 0;

        while (it.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it.next();
            String seq = seqRecord.getReadString();
            //store the sequence in a matcher
            Matcher matcher = pattern.matcher(seq);

            //if a match is found
            if (matcher.find())
            {
                //get the occurance of the first match and create a subsequence from there
                String subs = seq.substring(matcher.start(), seq.length());
                //translate the sequence


                //if the dna sequence alread exists, increment it by one
                if (newSuperMotifs.contains(subs))
                {
                    readNamesArray.add(seqRecord.getReadHeader());
                    matchCounter++;
                }
            }

            //get the reverse complement
            String revSeq = revcom(seq);
            Matcher revMatcher = pattern.matcher(revSeq);
            //if a reverse complement match is found is found
            if (revMatcher.find())
            {
                //get the occurance of the first match and create a subsequence from there
                String subs = revSeq.substring(revMatcher.start(), revSeq.length());

                if (newSuperMotifs.contains(subs))
                {
                    readNamesArray.add(seqRecord.getReadHeader());
                    matchCounter++;
                }
            }
        }

        System.out.println("Found " + matchCounter + " matches");
        Iterator itm = readNamesArray.iterator();

        while (itm.hasNext())
        {
            String read = (String) itm.next();
            readNameWriter.write(read + "\n");
        }

        readNameWriter.close();
    }

    public void getGoodAAPrimers(File fastqIn, File readsFile, File primersOut, String strPattern) throws FileNotFoundException, IOException
    {
        //strPattern will be e.g.VLV[VA][VL]P
        HashMap<String, Integer> readPairs = new HashMap<String, Integer>();
        BufferedReader reader = new BufferedReader(new FileReader(readsFile));
        Writer out = new BufferedWriter(new FileWriter(primersOut));
        String read;
        int readsReadIn = 0;
        int matchesFound = 0;
        //read the required sequences and their sidedness into a hashmap
        while ((read = reader.readLine()) != null)
        {
            int hashIndex = read.indexOf("#");
            String readName = read.substring(0, hashIndex);
            String pair = read.substring(read.length() - 1, read.length());
            int end = Integer.parseInt(pair);
            readPairs.put(readName, end);
            readsReadIn++;
        }
        reader.close();
        System.out.println("Reads read in = " + readsReadIn);
        DNASequence dna;
        RNASequence rna;
        ProteinSequence aa;
        Pattern pattern = Pattern.compile(strPattern);

        HashMap<String, PEFastqRead> reads = new HashMap<String, PEFastqRead>();
        final FastqReader fastqReader1 = new FastqReader(fastqIn);

        while (fastqReader1.hasNext())
        {
            FastqRecord record = fastqReader1.next();
            String pairedReadName = record.getReadHeader();

            int hashIndex = pairedReadName.indexOf("#");
            String readName = pairedReadName.substring(0, hashIndex);
            String pair = pairedReadName.substring(pairedReadName.length() - 1, pairedReadName.length());
            int end = Integer.parseInt(pair);

            if (reads.containsKey(readName))
            {
                PEFastqRead pefr = reads.get(readName);

                if (end == 1)
                {
                    pefr.setLeftRead(record);
                }
                else if (end == 2)
                {
                    pefr.setRightRead(record);
                }
                else
                {
                    System.err.println("Couldn't determine read pair");
                    System.exit(0);
                }
                reads.put(readName, pefr);
            }
            else
            {
                PEFastqRead pefr = new PEFastqRead(readName, null, null);
                if (end == 1)
                {
                    pefr.setLeftRead(record);
                }
                else if (end == 2)
                {
                    pefr.setRightRead(record);
                }
                else
                {
                    System.err.println("Couldn't determine read pair");
                    System.exit(0);
                }
                reads.put(readName, pefr);
            }
        }

        final FastqReader fastqReader = new FastqReader(fastqIn);
        //loop through the fastq seqs
        while (fastqReader.hasNext())
        {
            FastqRecord record = fastqReader.next();
            String fqReadName = record.getReadHeader();
            FastqParser fp = new FastqParser();
            FastaSequence fasta = fp.fastqToFastaSeq(record);
            dna = new DNASequence(fasta.getSequence());

            //translate the dna into the six reading frames
            Frame[] frames = Frame.getAllFrames();
            int matchesFoundPerRecord = 0;
            for (Frame frame : frames)
            {
                rna = dna.getRNASequence(frame);
                aa = rna.getProteinSequence();
                //look for the VLVVP/VLVALP motif
                Matcher matcher = pattern.matcher(aa.toString());

                //if a match is found
                if (matcher.find())
                {
                    matchesFound++;
                    matchesFoundPerRecord++;
                    //get the readname and sidedness of the current read
                    int hashIndex = fqReadName.indexOf("#");
                    String readName = fqReadName.substring(0, hashIndex);
                    String pair = fqReadName.substring(fqReadName.length() - 1, fqReadName.length());
                    int currentEnd = Integer.parseInt(pair);

                    //get the sidedness from the hash map
                    Integer readSideInt = readPairs.get(readName);
                    int readSide = readSideInt.intValue();
                    //and see if they're the same. If they are...
                    if (currentEnd == readSide)
                    {
                        //get the frame number off the current frame
                        int frameNumber = frame.ordinal();
                        //if we have a left-hand read, then we want to make sure the reading from is on the positive strand (0-2)
                        if ((readSide == 1 && frameNumber <= 2) || (readSide == 2 && frameNumber >= 3))
                        {
                            PEFastqRead pfr = reads.get(readName);
                            String leftRead = pfr.getLeftRead().getReadString();
                            String rightRead = pfr.getRightRead().getReadString();
                            String leftSeq = (">primer_" + aa.toString() + "_" + readName + "#0/1" + "\n" + leftRead + "\n");
                            String rightSeq = (">primer_" + aa.toString() + "_" + readName + "#0/2" + "\n" + rightRead + "\n");
                            out.write(leftSeq);
                            out.write(rightSeq);
                        }
                    }
                    if (matchesFoundPerRecord > 1)
                    {
                        System.out.println("Found more than one match in " + fqReadName);
                    }
                }
            }


        }
        out.write("\n");
        out.close();
        System.out.println("Found " + matchesFound + " matches");
    }

    public void getGoodDNAPrimers(File fastqIn, File readsFile, File primersOut, String strPattern) throws FileNotFoundException, IOException, CharacterParseException
    {
        //strPattern will be e.g.GTGCTGGTGG[TC]G[TG]T[TG]CC
        HashMap<String, Integer> readPairs = new HashMap<String, Integer>();
        BufferedReader reader = new BufferedReader(new FileReader(readsFile));
        Writer out = new BufferedWriter(new FileWriter(primersOut));
        String read;
        int readsReadIn = 0;
        int matchesFound = 0;
        //read the required sequences and their sidedness into a hashmap
        while ((read = reader.readLine()) != null)
        {
            int hashIndex = read.indexOf("#");
            String readName = read.substring(0, hashIndex);
            String pair = read.substring(read.length() - 1, read.length());
            int end = Integer.parseInt(pair);
            readPairs.put(readName, end);
            readsReadIn++;
        }
        reader.close();
        System.out.println("Reads read in = " + readsReadIn);
        DNASequence dna;

        Pattern pattern = Pattern.compile(strPattern);

        HashMap<String, PEFastqRead> reads = new HashMap<String, PEFastqRead>();
        final FastqReader fastqReader1 = new FastqReader(fastqIn);

        while (fastqReader1.hasNext())
        {
            FastqRecord record = fastqReader1.next();
            String pairedReadName = record.getReadHeader();

            int hashIndex = pairedReadName.indexOf("#");
            String readName = pairedReadName.substring(0, hashIndex);
            String pair = pairedReadName.substring(pairedReadName.length() - 1, pairedReadName.length());
            int end = Integer.parseInt(pair);

            if (reads.containsKey(readName))
            {
                PEFastqRead pefr = reads.get(readName);

                if (end == 1)
                {
                    pefr.setLeftRead(record);
                }
                else if (end == 2)
                {
                    pefr.setRightRead(record);
                }
                else
                {
                    System.err.println("Couldn't determine read pair");
                    System.exit(0);
                }
                reads.put(readName, pefr);
            }
            else
            {
                PEFastqRead pefr = new PEFastqRead(readName, null, null);
                if (end == 1)
                {
                    pefr.setLeftRead(record);
                }
                else if (end == 2)
                {
                    pefr.setRightRead(record);
                }
                else
                {
                    System.err.println("Couldn't determine read pair");
                    System.exit(0);
                }
                reads.put(readName, pefr);
            }
        }

        final FastqReader fastqReader = new FastqReader(fastqIn);
        //loop through the fastq seqs
        while (fastqReader.hasNext())
        {
            FastqRecord record = fastqReader.next();
            String fqReadName = record.getReadHeader();
            String seq = record.getReadString();

            //look for the  motif
            Matcher matcher = pattern.matcher(seq);

            //if a match is found on the forward strand
            if (matcher.find())
            {
                matchesFound++;
                //get the readname and sidedness of the current read
                int hashIndex = fqReadName.indexOf("#");
                String readName = fqReadName.substring(0, hashIndex);
                String pair = fqReadName.substring(fqReadName.length() - 1, fqReadName.length());
                int currentEnd = Integer.parseInt(pair);

                //get the sidedness from the hash map
                Integer readSideInt = readPairs.get(readName);
                int readSide = readSideInt.intValue();
                //and see if they're the same. If they are...
                if (currentEnd == 1 && readSide == 1)
                {

                    PEFastqRead pfr = reads.get(readName);
                    String leftRead = pfr.getLeftRead().getReadString();
                    String rightRead = pfr.getRightRead().getReadString();
                    String leftSeq = (">primer_left_" + readName + "\n" + leftRead + "\n");
                    String rightSeq = (">primer_right_" + readName + "\n" + rightRead + "\n");
                    out.write(leftSeq);
                    out.write(rightSeq);

                }
            }
            else
            {
                String revcom = revcom(seq);
                Matcher revmatcher = pattern.matcher(revcom);

                //if a match is found on the forward strand
                if (revmatcher.find())
                {
                    matchesFound++;
                    //get the readname and sidedness of the current read
                    int hashIndex = fqReadName.indexOf("#");
                    String readName = fqReadName.substring(0, hashIndex);
                    String pair = fqReadName.substring(fqReadName.length() - 1, fqReadName.length());
                    int currentEnd = Integer.parseInt(pair);

                    //get the sidedness from the hash map
                    Integer readSideInt = readPairs.get(readName);
                    int readSide = readSideInt.intValue();
                    //and see if they're the same. If they are...
                    if (currentEnd == 2 && readSide == 2)
                    {

                        PEFastqRead pfr = reads.get(readName);
                        String leftRead = pfr.getLeftRead().getReadString();
                        String rightRead = pfr.getRightRead().getReadString();
                        String leftSeq = (">primer_left_" + readName + "\n" + leftRead + "\n");
                        String rightSeq = (">primer_right_" + readName + "\n" + rightRead + "\n");
                        out.write(leftSeq);
                        out.write(rightSeq);

                    }
                }
            }
        }
        out.write("\n");
        out.close();
        System.out.println("Found " + matchesFound + " matches");
    }

    public static void main(String[] args)
    {
        try
        {
            FastqMotifFinder fmf = new FastqMotifFinder();
            File fqfile = new File("/Users/ethering/temp/solexa/uk_us/US_from_UK_US_shared.fastq");
            File fqOut = new File("/Users/ethering/temp/solexa/uk_us/US_from_UK_US_shared_found.fastq");

            String pattern = "GTGCTGGTGG[TC]G[TG]T[TG]CC";

            fmf.getPEFastqReadsFromMotif(fqfile, pattern, fqOut);
        }
        catch (IOException ex)
        {
            Logger.getLogger(FastqMotifFinder.class.getName()).log(Level.SEVERE, null, ex);
        }
        catch (CharacterParseException ex)
        {
            Logger.getLogger(FastqMotifFinder.class.getName()).log(Level.SEVERE, null, ex);
        }

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