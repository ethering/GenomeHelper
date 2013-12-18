package uk.ac.tsl.etherington.genomehelper.fastq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.HashSet;
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
    /**
     * Returns a HashSet of fastq read names.
     *
     * @param readNames a text file of read names, one per line
     * @return HashSet of read names
     * @throws FileNotFoundException
     * @throws IOException
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
    
       /**
     * 
     * @param fastqReader
     * @param out
     * @param mappedReads 
     */
    public static void writeRecords(FastqReader fastqReader, FastqWriter out, HashSet<String> mappedReads)
    {
        while (fastqReader.hasNext())
        {
            FastqRecord record = fastqReader.next();
            String readName = record.getReadHeader();
            int hashIndex = readName.indexOf(" ");
            readName = readName.substring(0, hashIndex);
            if (mappedReads.contains(readName))
            {
                out.write(record);
            }
        }
        out.close();
    }


    /**
     * Extract paired-end fastq sequences from a list of sequence names
     *
     * @param list a HashSet of read names
     * @param fastqLeft all left-handed reads
     * @param fastqRight all right-handed reads
     * @param fastqLeftOut the left-handed reads in the list
     * @param fastqRightOut the right-handed reads in the list
     * @throws FileNotFoundException
     * @throws IOException
     */
    public void getPairedFastqSeqsFromHashSet(HashSet list, File fastqLeft, File fastqRight,
            File fastqLeftOut, File fastqRightOut) throws FileNotFoundException, IOException
    {
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

            int hashIndex = leftRead.indexOf(" ");
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
     * Extract fastq sequences in a list from a file of single-end fastq
     * sequences
     *
     * @param listFile a file of read names, one per line
     * @param fastqFileIn the fastq file to search for the read names
     * @param fastqFileOut the fastq file of found read names
     *
     */
    public void getOneSideFastqSeqsFromList(File listFile, File fastqFileIn, File fastqFileOut) throws FileNotFoundException, IOException
    {
        HashSet<String> list = new HashSet<>();
        BufferedReader input = new BufferedReader(new FileReader(listFile));

        String line = null;
        //System.out.println("List:");
        while ((line = input.readLine()) != null)
        {
            String[] array = line.split(" ");
            //if there's a space in the seq id then it will contain the 'handedness' of the read, e.g. 1:N:0:
            if (array.length > 1)
            {
                //split it up and get the sequence id
                int hashIndex = line.indexOf(" ");
                line = line.substring(0, hashIndex);
                if (line.startsWith("@"))
                {
                    //get rid of the @ sign
                    line = line.substring(1, line.length());
                }
                list.add(line);
            }
            else
            {
                list.add(line);
            }

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
            int hashIndex = readName.indexOf(" ");
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

    /**
     * Turns a fastq file into fasta format
     *
     * @param fastqIn the fastq intput file
     * @param fastaOut the fasta output file
     * @throws IOException
     */
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

    /**
     * Takes a net.sf.picard.fastq.FastqRecord and changes it into a
     * org.biojava3.data.sequence.FastqSequence
     *
     * @param record the FastqRecord to change
     * @return a FastaSequence object
     * @throws IOException
     */
    public FastaSequence fastqToFastaSeq(FastqRecord record) throws IOException
    {
        String readName = record.getReadHeader();
        String seq = record.getReadString();
        FastaSequence fasta = new FastaSequence(readName, seq);
        return fasta;
    }

    /**
     * Creates a file of six-frame translations (in fasta format) from a fastq
     * file
     *
     * @param fastqIn the fastq file to translate
     * @param fastaOut the translation file
     * @param includeDNASeq set to true if the original DNA sequence is required
     * in the outfile
     * @throws IOException
     */
    public void fastqToFastaSixFrameTranslation(File fastqIn, File fastaOut, boolean includeDNASeq) throws IOException
    {
        DNASequence dna;
        RNASequence rna;
        ProteinSequence aa;
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
                String aaSeq = (">" + record.getReadHeader() + "_" + (frame.ordinal() + 1) + "\n" + aa.toString() + "\n");
                out.write(aaSeq);
            }
        }
        out.close();
    }
}
