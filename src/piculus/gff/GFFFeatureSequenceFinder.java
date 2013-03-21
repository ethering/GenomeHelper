package piculus.gff;


import piculus.fasta.ContigPosition;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.RNASequence;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.template.SequenceView;
import org.biojava3.core.sequence.transcription.Frame;
import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.biojava3.genome.parsers.gff.GFF3Reader;
import org.biojava3.genome.parsers.gff.GeneMarkGTFReader;
import org.biojava3.genome.parsers.gff.Location;
import piculus.fastq.FastqMotifFinder;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author ethering
 */
public class GFFFeatureSequenceFinder
{

    public void findDNASequences(String gffFile, File fastaFile, String featureType, String strPattern) throws Exception
    {
        int noHitsFound = 0;
        int noHitsNotFound = 0;
        Pattern pattern = Pattern.compile(strPattern);
        Matcher forwardMatcher;
        Matcher reverseMatcher;
        LinkedHashMap<String, DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(fastaFile);
        LinkedHashMap<String, DNASequence> newSequenceMap = new LinkedHashMap<String, DNASequence>();
        for (Map.Entry<String, DNASequence> entry : sequenceMap.entrySet())
        {
            String seqName = entry.getKey();
            String[] seqNamesArray = seqName.split(" ");
            String newSeqName = seqNamesArray[0];
            newSequenceMap.put(newSeqName, entry.getValue());
        }
        FeatureList fl = GFF3Reader.read(gffFile);
        FeatureList featTypes = new FeatureList(fl.selectByType(featureType));
        //ArrayList<String> al = new ArrayList<String>(featTypes.attributeValues("ID"));
        //FeatureList atts = gffFile.selectByAttribute("ID");
        Iterator it = featTypes.iterator();
        while (it.hasNext())
        {
            FeatureI fi = (FeatureI) it.next();

            String chr = fi.seqname();
            Location loc = fi.location();

            int start = loc.bioStart();
            int end = loc.bioEnd();
            if (newSequenceMap.containsKey(chr))
            {
                String fastaSeq = newSequenceMap.get(chr).getSequenceAsString();
                String forwardSeq = fastaSeq.substring(start, end);
                String reverseSeq = FastqMotifFinder.revcom(forwardSeq);

                forwardMatcher = pattern.matcher(forwardSeq);
                reverseMatcher = pattern.matcher(reverseSeq);

                //if a match is found
                if (forwardMatcher.find() || reverseMatcher.find())
                {
                    noHitsFound++;
                }
                else
                {
                    noHitsNotFound++;
                }

            }
            else
            {
                System.out.println("Can't find " + chr + " in fasta file");
            }
        }
        System.out.println("Sequences with motifs: " + noHitsFound);
        System.out.println("Sequences without motifs: " + noHitsNotFound);

    }

    public void findAllDNAHitsInGenome(String crnGffFile, String transcriptGffFile, File fastaFile, String strPattern) throws Exception
    {
        int noHitsInCrns = 0;
        int noHitsInTx = 0;
        int noHitsNotFound = 0;
        ArrayList<ContigPosition> unFoundCps = new ArrayList<ContigPosition>();

        ArrayList<ContigPosition> hits = new ArrayList<ContigPosition>();
        Pattern pattern = Pattern.compile(strPattern);
        System.out.println("Pattern: " + pattern.toString());
        Matcher matcher;

        DNASequence dna;

        //open up the fasta files
        LinkedHashMap<String, DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(fastaFile);

        //loop through the fasta files
        for (Map.Entry<String, DNASequence> entry : sequenceMap.entrySet())
        {
            String seqName = entry.getKey();

            // strip off any text after the sequence id
            String[] seqNamesArray = seqName.split(" ");
            String newSeqName = seqNamesArray[0];
            //get the dna sequence
            dna = entry.getValue();

            //check for a match
            matcher = pattern.matcher(dna.toString().toUpperCase());

            //if a match is found
            while (matcher.find())
            {
                //get the location of the start and end and create a ContigPostion object.
                int matchStart = matcher.start();
                int matchEnd = matcher.end();
                ContigPosition cp = new ContigPosition(newSeqName, matchStart, matchEnd, '+');
                //add it to hits
                hits.add(cp);
            }
            SequenceView<NucleotideCompound> revcom = dna.getReverseComplement();
            String dnaRevSeq = revcom.getSequenceAsString();
            //String dnaRevSeq = fastqapp.FastqMotifFinder.revcom(dna.toString().toUpperCase());
            int revSeqLength = dnaRevSeq.length();
            //System.out.println("Rev seq lenth "+revSeqLength);
            matcher = pattern.matcher(dnaRevSeq);
            while (matcher.find())
            {
                //get the location of the start and end and create a ContigPostion object.
                int matchStart = matcher.start();
                int matchEnd = matcher.end();
                int newMatchStart = revSeqLength - matchEnd + 1;
                int newMatchEnd = revSeqLength - matchStart + 1;
                ContigPosition cp = new ContigPosition(newSeqName, newMatchStart, newMatchEnd, '-');
                //add it to hits
                hits.add(cp);
            }
        }

        for (ContigPosition cp : hits)
        {
            //System.out.println(cp.toString());
            FeatureList crnfl = GFF3Reader.read(crnGffFile);
            //iterate over the crn annotation and see if the match is in a crn
            Iterator crnit = crnfl.iterator();
            boolean foundInCrn = false;

            while (crnit.hasNext() && foundInCrn == false)
            {
                FeatureI fi = (FeatureI) crnit.next();
                foundInCrn = featureOverlapsContigPosition(cp, fi);
                //if the cp is in a crn, count it
                if (foundInCrn == true)
                {
                    noHitsInCrns++;
                }
            }
            //if the feature isn't in a crn
            if (foundInCrn == false)
            {
                boolean foundInTx = false;
                //open up the genome transcript annotation file and see if it's found in any other transcripts
                FeatureList txfl = GeneMarkGTFReader.read(transcriptGffFile);
                Iterator txit = txfl.iterator();
                while (txit.hasNext() && foundInTx == false)
                {
                    FeatureI txfi = (FeatureI) txit.next();
                    foundInTx = featureOverlapsContigPosition(cp, txfi);
                    if (foundInTx == true)
                    {
                        noHitsInTx++;
                    }
                }
                if (foundInTx == false)
                {
                    unFoundCps.add(cp);
                    noHitsNotFound++;
                }
            }
        }
        if (unFoundCps.size() > 0)
        {
            System.out.println("Couldn't find the following hits:");
            for (ContigPosition cp : unFoundCps)
            {
                System.out.println(cp.toString());
            }
        }

        System.out.println("No hits found in geneome: " + hits.size());
        System.out.println("No hits found in crn= " + noHitsInCrns);
        System.out.println("No hits found in tx= " + noHitsInTx);
        System.out.println("No hits not found in any annotation= " + noHitsNotFound);
    }

    public void findAllDNAHitsInGenome(String crnGffFile, String transcriptGffFile, File fastaFile, File strPatternFastaFile) throws Exception
    {
        int noHitsInCrns = 0;
        int noHitsInTx = 0;
        int noHitsNotFound = 0;
        LinkedHashMap<String, DNASequence> patternSeqs = FastaReaderHelper.readFastaDNASequence(strPatternFastaFile);

        ArrayList<ContigPosition> hits = new ArrayList<ContigPosition>();
        Matcher matcher;
        DNASequence dna;
        DNASequence searchSeq;
        //open up the fasta files
        LinkedHashMap<String, DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(fastaFile);

        //loop through the fasta files
        for (Map.Entry<String, DNASequence> entry : sequenceMap.entrySet())
        {
            String seqName = entry.getKey();
            // strip off any text after the sequence id
            String[] seqNamesArray = seqName.split(" ");
            String newSeqName = seqNamesArray[0];
            //get the dna sequence
            dna = entry.getValue();

            //loop through the pattern file check for a match
            for (Map.Entry<String, DNASequence> searchSeqs : patternSeqs.entrySet())
            {
                searchSeq = searchSeqs.getValue();
                Pattern pattern = Pattern.compile(searchSeq.toString());
                matcher = pattern.matcher(dna.toString());

                //if a match is found
                while (matcher.find())
                {
                    //get the location of the start and end and create a ContigPostion object.
                    int matchStart = matcher.start();
                    int matchEnd = matcher.end();
                    ContigPosition cp = new ContigPosition(newSeqName, matchStart, matchEnd);
                    //add it to hits
                    hits.add(cp);
                }
            }
            //do the same for the reverse strand
            String dnaRevSeq = FastqMotifFinder.revcom(dna.toString());
            for (Map.Entry<String, DNASequence> searchSeqs : patternSeqs.entrySet())
            {
                searchSeq = searchSeqs.getValue();
                Pattern pattern = Pattern.compile(searchSeq.toString());
                matcher = pattern.matcher(dnaRevSeq);

                //if a match is found
                while (matcher.find())
                {
                    //get the location of the start and end and create a ContigPostion object.
                    int matchStart = matcher.start();
                    int matchEnd = matcher.end();
                    ContigPosition cp = new ContigPosition(newSeqName, matchStart, matchEnd);
                    //add it to hits
                    hits.add(cp);
                }
            }
        }

        for (ContigPosition cp : hits)
        {
            //System.out.println(cp.toString());
            FeatureList crnfl = GFF3Reader.read(crnGffFile);
            //iterate over the crn annotation and see if the match is in a crn
            Iterator crnit = crnfl.iterator();
            boolean foundInCrn = false;

            while (crnit.hasNext() && foundInCrn == false)
            {
                FeatureI fi = (FeatureI) crnit.next();
                foundInCrn = featureOverlapsContigPosition(cp, fi);
                //if the cp is in a crn, count it
                if (foundInCrn == true)
                {
                    noHitsInCrns++;
                }
            }
            //if the feature isn't in a crn
            if (foundInCrn == false)
            {
                boolean foundInTx = false;
                //open up the genome transcript annotation file and see if it's found in any other transcripts
                FeatureList txfl = GeneMarkGTFReader.read(transcriptGffFile);
                Iterator txit = txfl.iterator();
                while (txit.hasNext() && foundInTx == false)
                {
                    FeatureI txfi = (FeatureI) txit.next();
                    foundInTx = featureOverlapsContigPosition(cp, txfi);
                    if (foundInTx == true)
                    {
                        noHitsInTx++;
                    }
                }
                if (foundInTx == false)
                {
                    System.out.println("Couldn't find " + cp.toString());
                    noHitsNotFound++;
                }
            }
        }
        System.out.println("No hits found in geneome: " + hits.size());
        System.out.println("No hits found in crn= " + noHitsInCrns);
        System.out.println("No hits found in tx= " + noHitsInTx);
        System.out.println("No hits not found in any annotation= " + noHitsNotFound);
    }

    public void findAllDNAHitsInGenome2(String crnGffFile, String transcriptGffFile, File fastaFile, String strPattern) throws Exception
    {
        int noHitsInCrns = 0;
        int noHitsInTx = 0;
        int noHitsNotFound = 0;
        ArrayList<ContigPosition> unFoundCps = new ArrayList<ContigPosition>();

        ArrayList<ContigPosition> hits = new ArrayList<ContigPosition>();
        Pattern pattern = Pattern.compile(strPattern);
        System.out.println("Pattern: " + pattern.toString());
        Matcher matcher;

        DNASequence dna;

        //open up the fasta files
        LinkedHashMap<String, DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(fastaFile);

        //loop through the fasta files
        for (Map.Entry<String, DNASequence> entry : sequenceMap.entrySet())
        {
            String seqName = entry.getKey();

            // strip off any text after the sequence id
            String[] seqNamesArray = seqName.split(" ");
            String newSeqName = seqNamesArray[0];
            //get the dna sequence
            dna = entry.getValue();

            //check for a match
            matcher = pattern.matcher(dna.toString().toUpperCase());

            //if a match is found
            while (matcher.find())
            {
                //get the location of the start and end and create a ContigPostion object.
                int matchStart = matcher.start();
                int matchEnd = matcher.end();
                ContigPosition cp = new ContigPosition(newSeqName, matchStart, matchEnd, '+');
                //add it to hits
                hits.add(cp);
            }
            SequenceView<NucleotideCompound> revcom = dna.getReverseComplement();
            String dnaRevSeq = revcom.getSequenceAsString();
            //String dnaRevSeq = fastqapp.FastqMotifFinder.revcom(dna.toString().toUpperCase());
            int revSeqLength = dnaRevSeq.length();
            //System.out.println("Rev seq lenth "+revSeqLength);
            matcher = pattern.matcher(dnaRevSeq);
            while (matcher.find())
            {
                //get the location of the start and end and create a ContigPostion object.
                int matchStart = matcher.start();
                int matchEnd = matcher.end();
                int newMatchStart = revSeqLength - matchEnd + 1;
                int newMatchEnd = revSeqLength - matchStart + 1;
                ContigPosition cp = new ContigPosition(newSeqName, newMatchStart, newMatchEnd, '-');
                //add it to hits
                hits.add(cp);
            }
        }

        FeatureList crnfl = GFF3Reader.read(crnGffFile);
        FeatureList txfl = GeneMarkGTFReader.read(transcriptGffFile);
        for (ContigPosition cp : hits)
        {
            //System.out.println(cp.toString());

            //iterate over the crn annotation and see if the match is in a crn
            boolean foundInCrn = false;
            for (FeatureI fi : crnfl)
            {
                if (foundInCrn == false)
                {
                    foundInCrn = featureOverlapsContigPosition(cp, fi);
                    //if the cp is in a crn, count it
                    if (foundInCrn == true)
                    {
                        noHitsInCrns++;
                    }
                }
            }

            //if the feature isn't in a crn
            if (foundInCrn == false)
            {
                boolean foundInTx = false;
                //open up the genome transcript annotation file and see if it's found in any other transcripts
                for (FeatureI txfi : txfl)
                {
                    if (foundInTx == false)
                    {
                        foundInTx = featureOverlapsContigPosition(cp, txfi);
                        if (foundInTx == true)
                        {
                            noHitsInTx++;
                        }
                    }
                }
                
                if (foundInTx == false)
                {
                    unFoundCps.add(cp);
                    noHitsNotFound++;
                }
            }
        }
        if (unFoundCps.size() > 0)
        {
            System.out.println("Couldn't find the following hits:");
            for (ContigPosition cp : unFoundCps)
            {
                System.out.println(cp.toString());
            }
        }

        System.out.println("No hits found in geneome: " + hits.size());
        System.out.println("No hits found in crn= " + noHitsInCrns);
        System.out.println("No hits found in tx= " + noHitsInTx);
        System.out.println("No hits not found in any annotation= " + noHitsNotFound);
    }

    public void findProteinSequences(String gffFile, File fastaFile, String featureType, String strPattern) throws Exception
    {
        int noHitsFound = 0;
        int noHitsNotFound = 0;
        Pattern pattern = Pattern.compile(strPattern);
        System.out.println("Pattern: " + pattern.toString());
        Matcher matcher;
        DNASequence dna;
        RNASequence rna;
        ProteinSequence aa;

        System.out.println("The reference sequences should be DNA sequences");
        LinkedHashMap<String, DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(fastaFile);
        LinkedHashMap<String, DNASequence> newSequenceMap = new LinkedHashMap<String, DNASequence>();
        for (Map.Entry<String, DNASequence> entry : sequenceMap.entrySet())
        {
            String seqName = entry.getKey();

            String[] seqNamesArray = seqName.split(" ");
            String newSeqName = seqNamesArray[0];

            newSequenceMap.put(newSeqName, entry.getValue());
        }

        FeatureList fl;
        if (gffFile.endsWith(".gtf"))
        {
            fl = GeneMarkGTFReader.read(gffFile);
        }
        else
        {
            fl = GFF3Reader.read(gffFile);
        }

        FeatureList featTypes = new FeatureList(fl.selectByType(featureType));

        Iterator it = featTypes.iterator();
        while (it.hasNext())
        {
            FeatureI fi = (FeatureI) it.next();
            String chr = fi.seqname();
            //System.out.println(chr);
            Location loc = fi.location();

            char strandChar = loc.bioStrand();
            Strand strand = Strand.POSITIVE;
            if (strandChar == '-')
            {
                strand = Strand.NEGATIVE;
            }

            int start = loc.bioStart();
            int end = loc.bioEnd();
            if (newSequenceMap.containsKey(chr))
            {
                DNASequence mapChr = newSequenceMap.get(chr);
                String dnaSubSeq = mapChr.getSequenceAsString(start, end, strand);
                if (dnaSubSeq.length() > 2)
                {
                    dna = new DNASequence(dnaSubSeq);
                    //System.out.println(dnaSubSeq.toString());
                    rna = dna.getRNASequence();
                    aa = rna.getProteinSequence();
                    String forwardAAString = aa.toString();
                    matcher = pattern.matcher(forwardAAString);

                    //if a match is found
                    if (matcher.find())
                    {
                        noHitsFound++;
                    }
                    else
                    {
                        noHitsNotFound++;
                    }
                }
            }
            else
            {
                System.out.println("Can't find " + chr + " in fasta file");
            }


        }

        System.out.println("Sequences with motifs: " + noHitsFound);
        System.out.println("Sequences without motifs: " + noHitsNotFound);
    }

    public void findProteinSequencesFrom6Frame(String gffFile, File fastaFile, String featureType, String strPattern) throws Exception
    {
        int noHitsFound = 0;
        int noHitsNotFound = 0;
        Pattern pattern = Pattern.compile(strPattern);
        System.out.println("Pattern: " + pattern.toString());
        Matcher matcher;

        DNASequence dna;
        RNASequence rna;
        ProteinSequence aa;

        LinkedHashMap<String, DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(fastaFile);
        LinkedHashMap<String, DNASequence> newSequenceMap = new LinkedHashMap<String, DNASequence>();
        for (Map.Entry<String, DNASequence> entry : sequenceMap.entrySet())
        {
            String seqName = entry.getKey();
            String[] seqNamesArray = seqName.split(" ");
            String newSeqName = seqNamesArray[0];
            newSequenceMap.put(newSeqName, entry.getValue());
        }

        FeatureList fl;
        if (gffFile.endsWith(".gtf"))
        {
            fl = GeneMarkGTFReader.read(gffFile);
        }
        else
        {
            fl = GFF3Reader.read(gffFile);
        }
        FeatureList featTypes = new FeatureList(fl.selectByType(featureType));
        //ArrayList<String> al = new ArrayList<String>(featTypes.attributeValues("ID"));
        //FeatureList atts = gffFile.selectByAttribute("ID");
        Iterator it = featTypes.iterator();
        while (it.hasNext())
        {
            FeatureI fi = (FeatureI) it.next();

            String chr = fi.seqname();
            //System.out.println(chr);
            Location loc = fi.location();

            int start = loc.bioStart();
            int end = loc.bioEnd();
            if (newSequenceMap.containsKey(chr))
            {
                DNASequence mapChr = newSequenceMap.get(chr);
                String dnaStr = mapChr.getSubSequence(start, end).getSequenceAsString();
                if (dnaStr.length() > 2)
                {
                    dna = new DNASequence(dnaStr);
                    Frame[] frames = Frame.getAllFrames();
                    boolean matchFound = false;
                    for (Frame frame : frames)
                    {
                        rna = dna.getRNASequence(frame);
                        aa = rna.getProteinSequence();
                        matcher = pattern.matcher(aa.toString());

                        //if a match is found
                        if (matcher.find())
                        {
                            matchFound = true;
                        }
                    }
                    if (matchFound == true)
                    {
                        noHitsFound++;
                    }
                    else
                    {
                        noHitsNotFound++;
                    }
                }
            }
            else
            {
                System.out.println("Can't find " + chr + " in fasta file");
            }


        }
        System.out.println(
                "Sequences with motifs: " + noHitsFound);
        System.out.println(
                "Sequences without motifs: " + noHitsNotFound);
    }

    public boolean featureOverlapsContigPosition(ContigPosition cp, FeatureI feat)
    {
        boolean found = false;

        String cpId = cp.getContigid();
        String featId = feat.seqname();
        if (featId.equalsIgnoreCase(cpId))
        {
            Location featLoc = feat.location();

            if (featLoc.bioStrand() == cp.getStrand())
            {
                Location cpLoc = new Location(cp.getStart(), cp.getEnd());
                if (cp.getStrand() == '-')
                {
                    cpLoc = cpLoc.minus();
                }

                if (cpLoc.overlaps(featLoc))
                {
                    found = true;
                }
            }
        }
        return found;
    }

    //a temp method to look for VLVVVP in all genes
    public void findProteinSequencesInGenome(String gffFile, File fastaFile, String featureType, String strPattern) throws Exception
    {
        HashSet<String> genesWithPattern = new HashSet<String>();
        HashSet<String> genesWithoutPattern = new HashSet<String>();

        Pattern pattern = Pattern.compile(strPattern);
        System.out.println("Pattern: " + pattern.toString());
        Matcher matcher;
        DNASequence dna;
        RNASequence rna;
        ProteinSequence aa;

        System.out.println("The reference sequences should be DNA sequences");
        LinkedHashMap<String, DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(fastaFile);
        LinkedHashMap<String, DNASequence> newSequenceMap = new LinkedHashMap<String, DNASequence>();
        for (Map.Entry<String, DNASequence> entry : sequenceMap.entrySet())
        {
            String seqName = entry.getKey();

            String[] seqNamesArray = seqName.split(" ");
            String newSeqName = seqNamesArray[0];

            newSequenceMap.put(newSeqName, entry.getValue());
        }

        FeatureList fl;
        if (gffFile.endsWith(".gtf"))
        {
            fl = GeneMarkGTFReader.read(gffFile);
        }
        else
        {
            fl = GFF3Reader.read(gffFile);
        }

        FeatureList featTypes = new FeatureList(fl.selectByType(featureType));

        Iterator it = featTypes.iterator();
        while (it.hasNext())
        {
            FeatureI fi = (FeatureI) it.next();
            String chr = fi.seqname();
            //System.out.println(chr);
            Location loc = fi.location();

            char strandChar = loc.bioStrand();
            Strand strand = Strand.POSITIVE;
            if (strandChar == '-')
            {
                strand = Strand.NEGATIVE;
            }

            int start = loc.bioStart();
            int end = loc.bioEnd();
            if (newSequenceMap.containsKey(chr))
            {
                DNASequence mapChr = newSequenceMap.get(chr);
                String dnaSubSeq = mapChr.getSequenceAsString(start, end, strand);
                if (dnaSubSeq.length() > 2)
                {
                    dna = new DNASequence(dnaSubSeq);
                    //System.out.println(dnaSubSeq.toString());
                    rna = dna.getRNASequence();
                    aa = rna.getProteinSequence();
                    String forwardAAString = aa.toString();
                    matcher = pattern.matcher(forwardAAString);
                    String geneId = fi.getAttribute("gene_id");
                    //if a match is found
                    if (matcher.find())
                    {

                        genesWithPattern.add(geneId);
                    }
                    else
                    {
                        genesWithoutPattern.add(geneId);
                    }
                }
            }
            else
            {
                System.out.println("Can't find " + chr + " in fasta file");
            }


        }

        System.out.println("Sequences with motifs: " + genesWithPattern.size());
        System.out.println("Sequences without motifs: " + genesWithoutPattern.size());
    }

    public static void main(String[] args)
    {
        String crnGffFile = "/Users/ethering/temp/crinkler/crn_1.2_annoation_genes.gff";
        String transcriptGffFile = "/Users/ethering/temp/crinkler/pi_transcripts_supercontig1.2.gtf";
        File fastaFile = new File("/Users/ethering/temp/crinkler/pi_supercontig_1.2.fasta");
        //String strPattern = "GTGCTGGTGG";
        String strPattern = "GTGCTGGTGG[TC]G[TG]T[TG]CC";
        //File fastaSearchFile = new File("/Users/ethering/temp/crinkler/New_cds_unique.fas");


        GFFFeatureSequenceFinder fsf = new GFFFeatureSequenceFinder();
        try
        {
            fsf.findAllDNAHitsInGenome2(crnGffFile, transcriptGffFile, fastaFile, strPattern);
        }
        catch (Exception ex)
        {
            Logger.getLogger(GFFFeatureSequenceFinder.class
                    .getName()).log(Level.SEVERE, null, ex);
        }


    }
}
