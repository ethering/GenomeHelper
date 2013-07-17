package uk.ac.tsl.etherington.piculus.gff;


import java.io.*;
import java.util.*;
import java.util.Map.Entry;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.biojava3.genome.parsers.gff.GFF3Reader;
import org.biojava3.genome.parsers.gff.Location;
import uk.ac.tsl.etherington.piculus.fasta.ContigPosition;


/*
 * To change this template, choose Tools | Templates and open the template in
 * the editor.
 */
/**
 *
 * @author ethering
 */
public class GffFileHandler
{

    public void printAttributes(String file) throws IOException
    {


        FeatureList fl = GFF3Reader.read(file);
        FeatureList featTypes = new FeatureList(fl.selectByType("mRNA"));
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
            String type = fi.type();
            String group = fi.group();
            String[] ids = group.split("=");
            String id = ids[1];


            System.out.println(chr + "\t" + start + "\t" + end + "\t" + type + "\t" + id);
        }
    }

    public void makeBins(File fastaFile, File gffFile, int binSize) throws IOException, Exception
    {
        BufferedWriter out = new BufferedWriter(new FileWriter(gffFile));
        LinkedHashMap<String, DNASequence> a = FastaReaderHelper.readFastaDNASequence(fastaFile);
        //FastaReaderHelper.readFastaDNASequence for DNA sequences
        out.write("##gff-version 3\n");
        for (Entry<String, DNASequence> entry : a.entrySet())
        {
            String seqId = entry.getValue().getOriginalHeader();
            int seqLength = entry.getValue().getSequenceAsString().length();
            //System.out.println("SeqLength = "+seqId + " "+seqLength);
            for (int i = 1; i < seqLength + 1; i = i + binSize)
            {
                int seqEnd = i + binSize - 1;
                if (seqEnd > seqLength)
                {
                    seqEnd = seqLength;
                }

                //Feature feat = new Feature(seqId, "tsl_bioinformatics", "genomic_location", l, 0.0, 0, "ID="+seqId+"_"+i);
                out.write(seqId + "\ttsl_bioinformatics" + "\tDNA\t" + i + "\t" + seqEnd + "\t.\t+\t.\t" + "ID=" + seqId + "_" + i);
                out.write("\n");
            }

        }
        out.close();
    }

    public void findNonOverlappingFeatures(String gffFileString1, String gffFileString2) throws IOException, Exception
    {
        HashMap<String, ContigPosition> cps = new HashMap<String, ContigPosition>();
        File f = new File(gffFileString1);
        System.out.println("File exists = " + f.exists());
        FeatureList gffFile = GFF3Reader.read(gffFileString1);


        Iterator it = gffFile.iterator();
        while (it.hasNext())
        {
            FeatureI fi = (FeatureI) it.next();
            String chr = fi.seqname();
            Location loc = fi.location();
            int start = loc.bioStart();
            int end = loc.bioEnd();
            char strand = loc.bioStrand();
            //String group = fi.group();
            //String[] ids = group.split("=");
            //String geneid = ids[1];
            //System.out.println("Group = "+group);
            String name = fi.getAttribute("name");
            String att = fi.getAttribute("transcriptId");
            //System.out.println("Name = " + name + " Att = " + att);

            if (cps.containsKey(name))
            {
                ContigPosition cpc = cps.get(name);
                if (start < cpc.getStart())
                {
                    cpc.setStart(start);
                }
                if (end > cpc.getEnd())
                {
                    cpc.setEnd(end);
                }
                cps.put(name, cpc);
            }
            else
            {
                ContigPosition cp = new ContigPosition(chr, start, end, strand);
                cps.put(name, cp);
            }

        }
        HashMap<String, ContigPosition> cps2 = new HashMap<String, ContigPosition>();
        File f2 = new File(gffFileString2);
        System.out.println("File exists = " + f2.exists());
        FeatureList gffFile2 = GFF3Reader.read(gffFileString2);


        Iterator it2 = gffFile2.iterator();
        while (it2.hasNext())
        {
            FeatureI fi = (FeatureI) it2.next();
            String chr = fi.seqname();
            Location loc = fi.location();
            int start = loc.bioStart();
            int end = loc.bioEnd();
            char strand = loc.bioStrand();
            //String group = fi.group();
            //String[] ids = group.split("=");
            //String geneid = ids[1];
            //System.out.println("Group = "+group);
            String name = fi.getAttribute("gene_id");
            String att = fi.getAttribute("transcriptId");
            //System.out.println("Name = " + name + " Att = " + att);

            if (cps2.containsKey(name))
            {
                ContigPosition cpc = cps2.get(name);
                if (start < cpc.getStart())
                {
                    cpc.setStart(start);
                }
                if (end > cpc.getEnd())
                {
                    cpc.setEnd(end);
                }
                cps2.put(name, cpc);
            }
            else
            {
                ContigPosition cp = new ContigPosition(chr, start, end, strand);
                cps2.put(name, cp);
            }

        }
        Iterator it3 = cps.entrySet().iterator();
        while (it3.hasNext())
        {
            Map.Entry pairs = (Map.Entry) it3.next();
            String id = (String) pairs.getKey();
            ContigPosition cpm = (ContigPosition) pairs.getValue();
            System.out.println("Name = " + id + " Location: chr " + cpm.getContigid() + " start = " + cpm.getStart() + " end = " + cpm.getEnd());

        }
    }

    public void idNewTranscripts(String gtf, File cuffdiff, File outfile) throws IOException
    {

        BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
        Set<String> loci = new HashSet<String>();
        HashMap<String, ContigPosition> gtfCps = new HashMap<String, ContigPosition>();
        FeatureList gffFile = GFF3Reader.read(gtf);
        //find the start and end of each transcript
        Iterator it = gffFile.iterator();
        while (it.hasNext())
        {
            FeatureI fi = (FeatureI) it.next();
            //System.out.println("feature"+ fi.toString());
            String chr = fi.seqname();
            //System.out.println("chr = "+chr);
            Location loc = fi.location();
            int start = loc.bioStart();
            //System.out.println("start = "+start);
            int end = loc.bioEnd();
            char strand = loc.bioStrand();

            String name = fi.getAttribute("name");
            //System.out.println("Name = "+name);

            if (gtfCps.containsKey(name))
            {
                ContigPosition cpc = gtfCps.get(name);
                if (start < cpc.getStart())
                {
                    cpc.setStart(start);
                }
                if (end > cpc.getEnd())
                {
                    cpc.setEnd(end);
                }
                gtfCps.put(name, cpc);
            }
            else
            {
                ContigPosition cp = new ContigPosition(chr, start, end, strand);
                gtfCps.put(name, cp);
            }
        }
        System.out.println("Identified " + gtfCps.size() + " contigposs");

        HashMap<String, ContigPosition> cuffCps = new HashMap<String, ContigPosition>();
        FileInputStream fstream = new FileInputStream(cuffdiff);
        // Get the object of DataInputStream
        DataInputStream in = new DataInputStream(fstream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));

        String strLine = br.readLine();//skip the first line
        //Read File Line By Line
        while ((strLine = br.readLine()) != null)
        {
            //System.out.println(strLine);
            String[] data = strLine.split("\\t");
            String geneId = data[0];
            String locus = data[3];
            //System.out.println("Locus "+locus);

            String[] locusArray = locus.split(":");
            String scaffoldID = locusArray[0];
            //System.out.println("scafoldid "+scaffoldID);

            String[] coords = locusArray[1].split("-");
            int start = Integer.parseInt(coords[0]);
            int end = Integer.parseInt(coords[1]);
            //System.out.println("Start "+start + " end "+end);
            if (cuffCps.containsKey(geneId))
            {
                ContigPosition cpc = cuffCps.get(geneId);
                if (start < cpc.getStart())
                {
                    cpc.setStart(start);
                }
                if (end > cpc.getEnd())
                {
                    cpc.setEnd(end);
                }
                cuffCps.put(scaffoldID, cpc);
            }
            else
            {
                ContigPosition cp = new ContigPosition(scaffoldID, start, end);
                cuffCps.put(geneId, cp);
            }
        }
        //Close the input stream
        in.close();
        System.out.println("Identified " + cuffCps.size() + " cuff poss");

        Iterator it3 = cuffCps.entrySet().iterator();

        while (it3.hasNext())
        {
            boolean overlaps = false;
            Map.Entry cuffpairs = (Map.Entry) it3.next();
            //System.out.println(cuffpairs.getKey() + " = " + cuffpairs.getValue());
            String xloc = (String) cuffpairs.getKey();
            ContigPosition cuffcp = (ContigPosition) cuffpairs.getValue();
            String cuffId = cuffcp.getContigid();
            Iterator it2 = gtfCps.entrySet().iterator();
            while (it2.hasNext())
            {
                Map.Entry gtfpairs = (Map.Entry) it2.next();
                //System.out.println(gtfpairs.getKey() + " = " + gtfpairs.getValue());
                //String gtfId = (String) gtfpairs.getKey();
                ContigPosition cp = (ContigPosition) gtfpairs.getValue();
                String gtfId = cp.getContigid();

                if (gtfId.equalsIgnoreCase(cuffId))
                {


                    if (cp.overlaps(cuffcp))
                    {

                        //System.out.println(cuffId + " = " + cuffcp.getStart() + ":" + cuffcp.getEnd());
                        //System.out.println("Overlaps with " + gtfId + " = " + cp.getStart() + ":" + cp.getEnd() + "\n");
                        overlaps = true;
                    }
                }
            }
            if (overlaps == false)
            {
                //System.out.println(cuffId + ":" + cuffcp.getStart() + "-" + cuffcp.getEnd() + " doesn't overlap");
                //int length = cuffcp.getEnd() - cuffcp.getStart() +1;
                out.write(cuffId + "\tcufflinks\tgene\t" + cuffcp.getStart() + "\t" + cuffcp.getEnd() + "\t.\t+\t.\tname \"" + xloc + "\"\n");
                loci.add(xloc);
            }

        }
        out.close();
        System.out.println("Found " + loci.size() + " loci");

    }

    //given a list of IDs (contained in attribute 'ID='), identify the co-ordinated of the ID on the Arabidopsis genome
    //and test to see if they overlap a coding region
    public void identifyExpressedCodingRegions(File idFile, String clusterGffFile, String arabidopsisGeneFile) throws FileNotFoundException, IOException, Exception
    {
        HashSet<String> genesHit = new HashSet<String>();
        //something to put the clusterids in
        HashSet<String> clusters = new HashSet<String>();
        //open up the infile of cluster ids
        FileInputStream fstream = new FileInputStream(idFile);
        // Get the object of DataInputStream
        DataInputStream in = new DataInputStream(fstream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));
        String strLine = br.readLine();//skip the first line
        //Read File Line By Line
        while ((strLine = br.readLine()) != null)
        {
            //System.out.println(strLine);
            String[] data = strLine.split("\\t");
            String clusterId = data[0];
            //add the cluserID into the clusters HashSet
            clusters.add(clusterId);
        }
        FeatureList arabidopsisFl = GFF3Reader.read(arabidopsisGeneFile).selectByType("gene");
        FeatureList clusterFl = GFF3Reader.read(clusterGffFile);
        //loop thru the clusters gff file
        for (FeatureI clustFeat : clusterFl)
        {
            //get the clusterID
            String clusterNo = clustFeat.getAttribute("ID");
            //check to see if the current clusterID is in the clusters HashSet
            if (clusters.contains(clusterNo))
            {
                //if it is, get the seqname and location objects
                String seqname = clustFeat.seqname();
                Location location = clustFeat.location();
                //get any Arabidopsis genes overlapping this location and strand
                FeatureList arabidopsisFlSubset = arabidopsisFl.selectOverlapping(seqname, location, false);
                if (!arabidopsisFlSubset.isEmpty())
                {
                    for (FeatureI atFeat : arabidopsisFlSubset)
                    {
                        //System.out.println(atFeat.toString());
                        String id = atFeat.getAttribute("ID");
                        genesHit.add(id);
                    }
                }
            }
        }
        System.out.println("No genes hit = " + genesHit.size());
    }

    

//    public static void main(String[] args)
//    {
//        GffFileHandler hand = new GffFileHandler();
//        String loc = "/Users/ethering/temp/nblrr/nblrr.gff";
//        try
//        {
//            hand.gffToInterval(loc);
//        }
//        catch (IOException ex)
//        {
//            Logger.getLogger(GffFileHandler.class.getName()).log(Level.SEVERE, null, ex);
//        }
//    }
}
