package piculus.gff;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.ArrayList;
import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.biojava3.genome.parsers.gff.GFF3Reader;
import org.biojava3.genome.parsers.gff.Location;
import piculus.fasta.ContigPosition;



/*
 * To change this template, choose Tools | Templates and open the template in
 * the editor.
 */
/**
 *
 * @author ethering
 */
public class CompareGtfWithCuffdiff
{

    public ArrayList<ContigPosition> getContigPositionsFromCuffdiff(String cuffdiffFile) throws FileNotFoundException, IOException
    {
        int totalTranscriptLength = 0;
        ArrayList<ContigPosition> cuffTranscripts = new ArrayList<ContigPosition>();
        
        BufferedReader input = new BufferedReader(new FileReader(cuffdiffFile));
        try
        {
            String line = input.readLine();//skip the first header line //not declared within while loop

            while ((line = input.readLine()) != null)
            {
                String[] array = line.split("\t");
                String locus = array[3];
                //System.out.println("Locus = " + locus);
                String chr = locus.substring(0, locus.indexOf(":"));
                String coords = locus.substring(locus.indexOf(":") + 1, locus.length());
                String[] cordsArray = coords.split("-");
                int start = Integer.parseInt(cordsArray[0]);
                int end = Integer.parseInt(cordsArray[1]);
                ContigPosition cp = new ContigPosition(chr, start, end);
                if (!cp.exists(cuffTranscripts))
                {
                    cuffTranscripts.add(cp);
                    int transcriptLength = end - start + 1;
                    totalTranscriptLength += transcriptLength;
                }

                //System.out.println("Chr = " + chr + " start = " + start + " end = " + end);
            }
        }
        finally
        {
            input.close();
        }

        System.out.println("No. of cuffdiff genes = " + cuffTranscripts.size());
        double averageTranscriptLength = totalTranscriptLength / cuffTranscripts.size();
        System.out.println("Transcript space = " + totalTranscriptLength);
        System.out.println("Average transcript length = " + averageTranscriptLength);
        return cuffTranscripts;
    }

    public HashMap<String, ContigPosition> getContigPositionsFromCuffcompare(String cuffCompareFile) throws IOException
    {
        int totalTranscriptLength = 0;
        HashMap<String, ContigPosition> gtfTranscripts = new HashMap<String, ContigPosition>();
        FeatureList gffFile = GFF3Reader.read(cuffCompareFile);

        Iterator it = gffFile.iterator();
        while (it.hasNext())
        {
            FeatureI fi = (FeatureI) it.next();
            String chr = fi.seqname();//the scaffold id
            Location loc = fi.location();
            int start = loc.bioStart();
            int end = loc.bioEnd();
            char strand = loc.bioStrand();
            String xloc = fi.getAttribute("gene_id");//the xloc number

            //find the lowest start and highest end position for each xloc
            if (gtfTranscripts.containsKey(xloc))
            {
                ContigPosition cpc = gtfTranscripts.get(xloc);
                if (start < cpc.getStart())
                {
                    cpc.setStart(start);
                }
                if (end > cpc.getEnd())
                {
                    cpc.setEnd(end);
                }
                gtfTranscripts.put(xloc, cpc);
            }
            else
            {
                ContigPosition cp = new ContigPosition(chr, start, end, strand);
                gtfTranscripts.put(xloc, cp);
            }
        }
        
        //remove duplicate transcripts
        HashMap<String, ContigPosition> uniqueGtfTranscripts = new HashMap<String, ContigPosition>();
        Iterator gtfIt = gtfTranscripts.entrySet().iterator();
        while (gtfIt.hasNext())
        {
            Map.Entry pairs = (Map.Entry) gtfIt.next();
            String gtfId = (String) pairs.getKey();
            ContigPosition gtfCp = (ContigPosition) pairs.getValue();
            boolean contained = false;

            Iterator novelIt = uniqueGtfTranscripts.entrySet().iterator();
            while (novelIt.hasNext())
            {
                Map.Entry novelPairs = (Map.Entry) novelIt.next();
                //String id = (String) gtfPairs.getKey();
                ContigPosition novelCp = (ContigPosition) novelPairs.getValue();
                if (gtfCp.equals(novelCp))
                {
                    contained = true;
                    //System.out.println("Found identical genes "+gtfId+" "+gtfCp.toString());
                }
            }
            if(contained == false)
            {
                uniqueGtfTranscripts.put(gtfId, gtfCp);
            }


        }
        
        
        Iterator novelGtfIt = uniqueGtfTranscripts.entrySet().iterator();
        while (novelGtfIt.hasNext())
        {
            Map.Entry pairs = (Map.Entry) novelGtfIt.next();
            String id = (String) pairs.getKey();
            ContigPosition gtfCp = (ContigPosition) pairs.getValue();
            String cpId = gtfCp.getContigid();
            int start = gtfCp.getStart();
            int end = gtfCp.getEnd();
            int transcriptLength = end - start + 1;
            totalTranscriptLength += transcriptLength;
            //System.out.println(id+"\t"+cpId+":"+start+"-"+end);
            
        }
        
        System.out.println("No of cuffcompare transcripts = " + uniqueGtfTranscripts.size());

        double averageTranscriptLength = totalTranscriptLength / uniqueGtfTranscripts.size();
        System.out.println("Cuffcompare transcript space = " + totalTranscriptLength);
        System.out.println("Average cuffcompare transcript length = " + averageTranscriptLength);

        return uniqueGtfTranscripts;
    }

    public HashMap<String, ContigPosition> getContigPositionsFromGtf(String gtfFile) throws IOException
    {
        int totalTranscriptLength = 0;
        HashMap<String, ContigPosition> gtfTranscripts = new HashMap<String, ContigPosition>();
        FeatureList gffFile = GFF3Reader.read(gtfFile);

        Iterator it = gffFile.iterator();
        while (it.hasNext())
        {
            FeatureI fi = (FeatureI) it.next();
            String chr = fi.seqname();
            Location loc = fi.location();
            int start = loc.bioStart();
            int end = loc.bioEnd();
            char strand = loc.bioStrand();
            String name = fi.getAttribute("name");
            //System.out.println("Name = " + name + " Att = " + att);

            if (gtfTranscripts.containsKey(name))
            {
                ContigPosition cpc = gtfTranscripts.get(name);
                if (start < cpc.getStart())
                {
                    cpc.setStart(start);
                }
                if (end > cpc.getEnd())
                {
                    cpc.setEnd(end);
                }
                gtfTranscripts.put(name, cpc);
            }
            else
            {
                ContigPosition cp = new ContigPosition(chr, start, end, strand);
                gtfTranscripts.put(name, cp);
            }
        }

        Iterator gtfIt = gtfTranscripts.entrySet().iterator();
        while (gtfIt.hasNext())
        {
            Map.Entry pairs = (Map.Entry) gtfIt.next();
            String id = (String) pairs.getKey();
            ContigPosition gtfCp = (ContigPosition) pairs.getValue();
            int start = gtfCp.getStart();
            int end = gtfCp.getEnd();
            int transcriptLength = end - start + 1;
            totalTranscriptLength += transcriptLength;
        }
        System.out.println("No of gtf transcripts = " + gtfTranscripts.size());

        double averageTranscriptLength = totalTranscriptLength / gtfTranscripts.size();
        System.out.println("Transcript space = " + totalTranscriptLength);
        System.out.println("Average transcript length = " + averageTranscriptLength);



        return gtfTranscripts;
    }

    public void compareContigPositions(HashMap<String, ContigPosition> cuffCps, HashMap<String, ContigPosition> gtfCps)
    {
        HashMap<ContigPosition, ArrayList<String>> uniqueGtfIds = new HashMap<ContigPosition, ArrayList<String>>();
        int identialTranscripts = 0;
        int overlappingTranscripts = 0;
        int combinedLengthOfUniqueCuffTranscripts = 0;
        double totalPercentOverlaps = 0;
        int novelCuffTranscripts = 0;
        int cuffSize = 0;

        //iterate over the cuffcompare contig positions (which  now represent the longest transcripts)
        Iterator cuffIt = cuffCps.entrySet().iterator();
        while (cuffIt.hasNext())
        {
            Map.Entry pairs = (Map.Entry) cuffIt.next();
            //get the xloc number (cuffcompare gene id)
            //String xloc = (String) pairs.getKey();

            ContigPosition cuffCp = (ContigPosition) pairs.getValue();
            //get the cuffcompare scaffold number and start and end position
            //String cuffName = cuffCp.getContigid();
            int cuffstart = cuffCp.getStart();
            int cuffend = cuffCp.getEnd();
            //calculate the length of the current scaffold transcript
            cuffSize = cuffend - cuffstart + 1;
            boolean overlaps = false;

            //iterate over the ehux transcriptome annotation contig positions
            Iterator gtfIt = gtfCps.entrySet().iterator();
            while (gtfIt.hasNext())
            {
                Map.Entry gtfPairs = (Map.Entry) gtfIt.next();
                //get the transcript id and the scaffold id and start and end position of where it is found
                String id = (String) gtfPairs.getKey();//transcript id
                ContigPosition gtfCp = (ContigPosition) gtfPairs.getValue();
                int gtfstart = gtfCp.getStart();
                int gtfend = gtfCp.getEnd();

                if (gtfCp.equals(cuffCp))
                {
                    System.out.println(id + " " + cuffCp.toString() + " are identical transcript");
                    identialTranscripts++;
                }
                if (gtfCp.overlaps(cuffCp))
                {
                    //calculate the overlap between the two contig positions
                    int highestEnd;//the highest end position of both sequences
                    int lowestStart;//the lowest start positon of both sequences
                    int lowestEnd; //the lowest end position (i.e. the end position of the sequence that doesn't have the highest end)
                    int highestStart;//the higest start position (i.e. the start position of the sequence that doesn't have the lowest start)

//                    Seq1                  highestStart _______________________________________________ highestEnd
//                    Seq2   lowestStart _________________________________________ lowestEnd
//
//                                                       xxxxxxxxxxxxxxxxxxxxxxxxx overlap (lowestEnd - highestStart)
//                    
//                                       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx length (highestEnd - lowestStart)


                    if (cuffend >= gtfend)
                    {
                        highestEnd = cuffend;
                        lowestEnd = gtfend;
                    }
                    else
                    {
                        highestEnd = gtfend;
                        lowestEnd = cuffend;

                    }

                    if (cuffstart <= gtfstart)
                    {
                        lowestStart = cuffstart;
                        highestStart = gtfstart;
                    }
                    else
                    {
                        lowestStart = gtfstart;
                        highestStart = cuffstart;
                    }

                    double length = highestEnd - lowestStart;
                    double overlap = lowestEnd - highestStart;
                    //System.out.println("Length = "+length + " overlap = "+overlap);
                    double percentOverlap = (overlap / length) * 100;
                    //System.out.println(id + " " + gtfCp.toString() + " overlaps with cuffCp " + cuffCp.toString()+ " percent overlap =  = "+percentOverlap);
                    totalPercentOverlaps += percentOverlap;

                    overlaps = true;
                    if (uniqueGtfIds.containsKey(cuffCp))
                    {
                        ArrayList al = uniqueGtfIds.get(cuffCp);
                        al.add(id);
                        uniqueGtfIds.put(cuffCp, al);
                    }
                    else
                    {
                        ArrayList<String> al = new ArrayList<String>();
                        al.add(id);
                        uniqueGtfIds.put(cuffCp, al);

                    }
                }
            }
            if (overlaps == false)
            {
                //System.out.println(id + " doesn't overlap");
                novelCuffTranscripts++;
                combinedLengthOfUniqueCuffTranscripts += cuffSize;
            }
            else
            {
                overlappingTranscripts++;

            }
            //System.out.println("Name = " + id + " Location: chr " + gtfCp.getContigid() + " start = " + gtfCp.getStart() + " end = " + gtfCp.getEnd());
        }
        double averageOverlap = totalPercentOverlaps / overlappingTranscripts;
        System.out.println("Overlapping transcripts = " + overlappingTranscripts);
        System.out.println("Identical identialTranscripts = " + identialTranscripts);
        System.out.println("Novel cuff transcripts = " + novelCuffTranscripts);
        System.out.println("Combined length of Novel cuff transcripts = " + combinedLengthOfUniqueCuffTranscripts);
        System.out.println("Average overlap = " + averageOverlap);
        System.out.println("Transcripts with overlaps:");


        Iterator it = uniqueGtfIds.entrySet().iterator();
        while (it.hasNext())
        {
            Map.Entry entry = (Map.Entry) it.next();
            ContigPosition mapCuffCp = (ContigPosition) entry.getKey();
            String cuffId = mapCuffCp.getContigid();
            int mapCuffStart = mapCuffCp.getStart();
            int mapCuffEnd = mapCuffCp.getEnd();
            System.out.print(cuffId + ":" + mapCuffStart + "-" + mapCuffEnd);
            ArrayList<String> gtfIdAl = (ArrayList<String>) entry.getValue();
            for (String gtfId : gtfIdAl)
            {
                System.out.print("\t" + gtfId);
            }
            System.out.println("");

        }
    }
    public void compareContigPositions(ArrayList<ContigPosition> cuffCps, HashMap<String, ContigPosition> gtfCps)
    {
        HashMap<ContigPosition, ArrayList<String>> uniqueGtfIds = new HashMap<ContigPosition, ArrayList<String>>();
        int identialTranscripts = 0;
        int overlappingTranscripts = 0;
        int combinedLengthOfUniqueCuffTranscripts = 0;
        double totalPercentOverlaps = 0;
        int novelCuffTranscripts = 0;
        int cuffSize = 0;

        //iterate over the cuffcompare contig positions (which  now represent the longest transcripts)
        Iterator cuffIt = cuffCps.iterator();
        while (cuffIt.hasNext())
        {

            ContigPosition cuffCp = (ContigPosition) cuffIt.next();
            //get the cuffcompare scaffold number and start and end position
            //String cuffName = cuffCp.getContigid();
            int cuffstart = cuffCp.getStart();
            int cuffend = cuffCp.getEnd();
            //calculate the length of the current scaffold transcript
            cuffSize = cuffend - cuffstart + 1;
            boolean overlaps = false;

            //iterate over the ehux transcriptome annotation contig positions
            Iterator gtfIt = gtfCps.entrySet().iterator();
            while (gtfIt.hasNext())
            {
                Map.Entry gtfPairs = (Map.Entry) gtfIt.next();
                //get the transcript id and the scaffold id and start and end position of where it is found
                String id = (String) gtfPairs.getKey();//transcript id
                ContigPosition gtfCp = (ContigPosition) gtfPairs.getValue();
                int gtfstart = gtfCp.getStart();
                int gtfend = gtfCp.getEnd();

                if (gtfCp.equals(cuffCp))
                {
                    System.out.println(id + " " + cuffCp.toString() + " are identical transcript");
                    identialTranscripts++;
                }
                if (gtfCp.overlaps(cuffCp))
                {
                    //calculate the overlap between the two contig positions
                    int highestEnd;//the highest end position of both sequences
                    int lowestStart;//the lowest start positon of both sequences
                    int lowestEnd; //the lowest end position (i.e. the end position of the sequence that doesn't have the highest end)
                    int highestStart;//the higest start position (i.e. the start position of the sequence that doesn't have the lowest start)

//                    Seq1                  highestStart _______________________________________________ highestEnd
//                    Seq2   lowestStart _________________________________________ lowestEnd
//
//                                                       xxxxxxxxxxxxxxxxxxxxxxxxx overlap (lowestEnd - highestStart)
//                    
//                                       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx length (highestEnd - lowestStart)


                    if (cuffend >= gtfend)
                    {
                        highestEnd = cuffend;
                        lowestEnd = gtfend;
                    }
                    else
                    {
                        highestEnd = gtfend;
                        lowestEnd = cuffend;

                    }

                    if (cuffstart <= gtfstart)
                    {
                        lowestStart = cuffstart;
                        highestStart = gtfstart;
                    }
                    else
                    {
                        lowestStart = gtfstart;
                        highestStart = cuffstart;
                    }

                    double length = highestEnd - lowestStart;
                    double overlap = lowestEnd - highestStart;
                    //System.out.println("Length = "+length + " overlap = "+overlap);
                    double percentOverlap = (overlap / length) * 100;
                    //System.out.println(id + " " + gtfCp.toString() + " overlaps with cuffCp " + cuffCp.toString()+ " percent overlap =  = "+percentOverlap);
                    totalPercentOverlaps += percentOverlap;

                    overlaps = true;
                    if (uniqueGtfIds.containsKey(cuffCp))
                    {
                        ArrayList al = uniqueGtfIds.get(cuffCp);
                        al.add(id);
                        uniqueGtfIds.put(cuffCp, al);
                    }
                    else
                    {
                        ArrayList<String> al = new ArrayList<String>();
                        al.add(id);
                        uniqueGtfIds.put(cuffCp, al);

                    }
                }
            }
            if (overlaps == false)
            {
                //System.out.println(id + " doesn't overlap");
                novelCuffTranscripts++;
                combinedLengthOfUniqueCuffTranscripts += cuffSize;
            }
            else
            {
                overlappingTranscripts++;

            }
            //System.out.println("Name = " + id + " Location: chr " + gtfCp.getContigid() + " start = " + gtfCp.getStart() + " end = " + gtfCp.getEnd());
        }
        double averageOverlap = totalPercentOverlaps / overlappingTranscripts;
        System.out.println("Overlapping transcripts = " + overlappingTranscripts);
        System.out.println("Identical identialTranscripts = " + identialTranscripts);
        System.out.println("Novel cuff transcripts = " + novelCuffTranscripts);
        System.out.println("Combined length of Novel cuff transcripts = " + combinedLengthOfUniqueCuffTranscripts);
        System.out.println("Average overlap = " + averageOverlap);
        System.out.println("Transcripts with overlaps:");


        Iterator it = uniqueGtfIds.entrySet().iterator();
        while (it.hasNext())
        {
            Map.Entry entry = (Map.Entry) it.next();
            ContigPosition mapCuffCp = (ContigPosition) entry.getKey();
            String cuffId = mapCuffCp.getContigid();
            int mapCuffStart = mapCuffCp.getStart();
            int mapCuffEnd = mapCuffCp.getEnd();
            System.out.print(cuffId + ":" + mapCuffStart + "-" + mapCuffEnd);
            ArrayList<String> gtfIdAl = (ArrayList<String>) entry.getValue();
            for (String gtfId : gtfIdAl)
            {
                System.out.print("\t" + gtfId);
            }
            System.out.println("");

        }
    }
}
