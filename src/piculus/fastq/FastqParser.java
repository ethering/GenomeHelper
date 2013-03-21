package piculus.fastq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.RNASequence;
import org.biojava3.core.sequence.transcription.Frame;
import org.biojava3.data.sequence.FastaSequence;


/*
 * To change this template, choose Tools | Templates and open the template in
 * the editor.
 */
/**
 *
 * @author ethering
 */
public class FastqParser
{

    /**
     * prints the sequence id's from a fastq file
     *
     * @param fastqfile
     *
     */
    public void getFastqReadNames(File fastqfile)
    {
        final FastqReader reader = new FastqReader(fastqfile);
        while (reader.hasNext())
        {
            FastqRecord record = reader.next();
            System.out.println(record.getReadHeader());

        }

    }

    /**
     * takes a Hashset of fastq sequence ids (with no paired end info), goes
     * through the pairs of the two fastq files, chops off any 'handendness'
     * info and writes any sequences found in the HashSet to a new file
     *
     * @param list
     * @param fastqFileInLeftString
     * @param fastqFileInRightString
     * @param fastqFileOutLeftString
     * @param fastqFileOutRightString
     *
     */
    public HashSet readNamesToHashSet(File readNames) throws FileNotFoundException, IOException
    {
        HashSet reads = new HashSet();
        BufferedReader reader = new BufferedReader(new FileReader(readNames));
        String motif;
        while ((motif = reader.readLine()) != null)
        {
            reads.add(motif);
        }
        return reads;
    }

    public void getPairedFastqSeqsFromHashSet(HashSet list, String fastqFileInLeftString, String fastqFileInRightString,
            String fastqFileOutLeftString, String fastqFileOutRightString) throws FileNotFoundException, IOException
    {


        File fastqLeft = new File(fastqFileInLeftString);
        File fastqRight = new File(fastqFileInRightString);

        File fastqLeftOut = new File(fastqFileOutLeftString);
        File fastqRightOut = new File(fastqFileOutRightString);

        final FastqReader readerLeft = new FastqReader(fastqLeft);
        final FastqReader readerRight = new FastqReader(fastqRight);


        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter outLeft = writer.newWriter(fastqLeftOut);
        FastqWriter outRight = writer.newWriter(fastqRightOut);
        int recordsfound = 0;

        while (readerLeft.hasNext())
        {
            FastqRecord recordLeft = readerLeft.next();
            FastqRecord recordRight = readerRight.next();
            String leftRead = recordLeft.getReadHeader();

            int hashIndex = leftRead.indexOf("#");
            leftRead = leftRead.substring(0, hashIndex);
            if (list.contains(leftRead))
            {
                recordsfound++;
                outLeft.write(recordLeft);
                outRight.write(recordRight);
            }
        }
        System.out.println("Processed " + list.size() + " , found " + recordsfound);

    }

    /**
     * takes a File of fastq sequence ids (with no paired end info), goes
     * through a single fastq file, chops off any 'handendness' info and writes
     * any sequences found in the HashSet to a new file
     *
     * @param listFile
     * @param fastqFileIn
     * @param fastqFileOut
     *
     */
    public void getOneSideFastqSeqsFromList(File listFile, File fastqFileIn, File fastqFileOut) throws FileNotFoundException, IOException
    {
        HashSet<String> list = new HashSet<String>();
        BufferedReader input = new BufferedReader(new FileReader(listFile));

        String line = null;
        //System.out.println("List:");
        while ((line = input.readLine()) != null)
        {
            int hashIndex = line.indexOf("#");
            line = line.substring(0, hashIndex);
            list.add(line);
            //System.out.println(line);
        }
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter out = writer.newWriter(fastqFileOut);

        final FastqReader reader = new FastqReader(fastqFileIn);

        int recordsfound = 0;
        //System.out.println("Reads:");
        while (reader.hasNext())
        {
            FastqRecord record = reader.next();

            String readName = record.getReadHeader();
            int hashIndex = readName.indexOf("#");
            readName = readName.substring(0, hashIndex);
            //System.out.println(readName);
            if (list.contains(readName))
            {
                recordsfound++;
                out.write(record);
            }
        }
        System.out.println("Processed " + list.size() + " , found " + recordsfound);
        //return fastqFileOut;
    }

    public void sortInterleavedFastq(File fastqFileIn, File fastqFileOut)
    {
        HashMap<String, PEFastqRead> reads = new HashMap<String, PEFastqRead>();
        final FastqReader reader = new FastqReader(fastqFileIn);

        while (reader.hasNext())
        {
            FastqRecord record = reader.next();

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
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter out = writer.newWriter(fastqFileOut);

        Iterator itm = reads.entrySet().iterator();
        while (itm.hasNext())
        {
            Map.Entry pairs = (Map.Entry) itm.next();
            PEFastqRead fqpair = (PEFastqRead) pairs.getValue();
            if (fqpair.getLeftRead() != null && fqpair.getRightRead() != null)
            {
                out.write(fqpair.getLeftRead());
                out.write(fqpair.getRightRead());
            }
            else if (fqpair.getRightRead() == null)
            {
                System.out.println("Couldn't find left read for " + pairs.getKey().toString());
            }
            else
            {
                System.out.println("Couldn't find right read for " + pairs.getKey().toString());
            }
        }
    }

    public void fastqToFastaFile(File fastqIn, File fastaOut) throws IOException
    {
        final FastqReader reader = new FastqReader(fastqIn);
        Writer out = new BufferedWriter(new FileWriter(fastaOut));
        while (reader.hasNext())
        {
            FastqRecord record = reader.next();
            out.write(fastqToFastaSeq(record).toString());
        }
        out.close();
    }

    public FastaSequence fastqToFastaSeq(FastqRecord record) throws IOException
    {
        String readName = record.getReadHeader();
        String seq = record.getReadString();
        FastaSequence fasta = new FastaSequence(readName, seq);

        return fasta;
    }

    public void fastqToFastaSixFrameTranslation(File fastqIn, File fastaOut, boolean includeDNASeq) throws IOException
    {
        DNASequence dna;
        RNASequence rna;
        ProteinSequence aa;
        String strPattern = "VLV[VA][VL]P";
        Pattern pattern = Pattern.compile(strPattern);
        final FastqReader reader = new FastqReader(fastqIn);
        Writer out = new BufferedWriter(new FileWriter(fastaOut));
        while (reader.hasNext())
        {
            FastqRecord record = reader.next();
            FastaSequence fasta = fastqToFastaSeq(record);
            dna = new DNASequence(fasta.getSequence());
            Frame[] frames = Frame.getAllFrames();
            if (includeDNASeq == true)
            {
                out.write(fasta.toString());
            }
            for (Frame frame : frames)
            {
                rna = dna.getRNASequence(frame);
                aa = rna.getProteinSequence();
                Matcher matcher = pattern.matcher(aa.toString());

                //if a match is found
                if (matcher.find())
                {
                    String aaSeq = (">" + record.getReadHeader() + "_" + (frame.ordinal() + 1) + "\n" + aa.toString() + "\n");
                    out.write(aaSeq);
                }
            }
            out.write("\n");
        }
        out.close();
    }

    public void catFastqFiles(File outfile, File[] files)
    {
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter out = writer.newWriter(outfile);
        for (File file : files)
        {
            System.out.println(file);
            FastqReader reader = new FastqReader(file);
            while (reader.hasNext())
            {
                FastqRecord record = reader.next();
                out.write(record);

            }
            reader.close();
        }
    }

//    public static void main(String[] args)
//    {
//        FastqParser fp = new FastqParser();
//        File fq1 = new File("/Users/ethering/temp/solexa/lane1_NoIndex_R1_2000.fastq");
//        File fq2 = new File("/Users/ethering/temp/lane1_NoIndex_R2_2000.fastq");
//        File [] fileArray = new File [args.length];
//        for (int i = 0; i < args.length; i++)
//        {
//            fileArray[i]=new File(args[i]);
//        }
//        
//            fp.catFastqFiles(fileArray);
//       
//    }
}
