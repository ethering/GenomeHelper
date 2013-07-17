package uk.ac.tsl.etherington.piculus.fastq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import org.jtr.transliterate.CharacterParseException;
import org.jtr.transliterate.CharacterReplacer;
import org.jtr.transliterate.Perl5Parser;

/*
 * To change this template, choose Tools | Templates and open the template in
 * the editor.
 */
/**
 *
 * @author ethering
 */
public class FastqReformater
{

    public void fastqSplitReadsInHalf(File fastqFileIn, File leftHalfOut, File rightHalfOut)
    {
        FastqReader fq = new FastqReader(fastqFileIn);

        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter goodLeftSeqs = writer.newWriter(leftHalfOut);
        FastqWriter goodRightSeqs = writer.newWriter(rightHalfOut);


        //create an interator for each file
        Iterator it = fq.iterator();

        int itCounter = 0;
        while (it.hasNext())
        {
            //get the corresponding reads
            FastqRecord seqRecord = (FastqRecord) it.next();
            int seqLength = seqRecord.getReadString().length();
            String readString = seqRecord.getReadString();
            String qual = seqRecord.getBaseQualityString();
            String leftRead = readString.substring(0, (seqLength / 2));
            String rightRead = readString.substring(seqLength / 2, seqLength);
            String leftQual = qual.substring(0, (seqLength / 2));
            String rightQual = qual.substring(seqLength / 2, seqLength);

            int index = seqRecord.getReadHeader().indexOf("#");
            String newName = seqRecord.getReadHeader().substring(0, index).concat("_split#0");

            FastqRecord newLeftSeq = new FastqRecord(newName.concat("/1"), leftRead, "", leftQual);
            FastqRecord newRightSeq = new FastqRecord(newName.concat("/2"), rightRead, "", rightQual);

            goodLeftSeqs.write(newLeftSeq);
            goodRightSeqs.write(newRightSeq);

            itCounter++;

        }
        System.out.println("Completed writing " + itCounter + " good reads");

    }

    public void reformatPEinfo(File fastqFileIn, File fastqFileOut)
    {

        FastqReader fqr = new FastqReader(fastqFileIn);
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter newSeqs = writer.newWriter(fastqFileOut);

        Iterator it = fqr.iterator();
        int itCounter = 0;
        while (it.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it.next();
            String id = seqRecord.getReadHeader();

            StringBuilder b = new StringBuilder(id);
            b.replace(id.lastIndexOf("/1_"), id.lastIndexOf("/1_") + 3, "/");
            id = b.toString();


            FastqRecord newSeq = new FastqRecord(id, seqRecord.getReadString(), "", seqRecord.getBaseQualityString());

            newSeqs.write(newSeq);


            itCounter++;
        }
        System.out.println("Completed writing " + itCounter + " reads");

    }

    public void fastq2Fasta(File fastqIn, File fastaOut)
    {
        FileWriter fw = null;
        int itCounter = 0;
        try
        {
            FastqReader fqr = new FastqReader(fastqIn);
            fw = new FileWriter(fastaOut.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);
            Iterator it = fqr.iterator();

            while (it.hasNext())
            {
                FastqRecord seqRecord = (FastqRecord) it.next();
                String id = seqRecord.getReadHeader();

                String newId = id.replaceAll(" ", "_");
                String seq = seqRecord.getReadString();

                bw.write(">");
                bw.write(newId + "\n");
                bw.write(seq + "\n");
                itCounter++;
            }
            bw.close();
        }
        catch (IOException ex)
        {
            Logger.getLogger(FastqReformater.class.getName()).log(Level.SEVERE, null, ex);
        }
        finally
        {
            try
            {
                fw.close();
            }
            catch (IOException ex)
            {
                Logger.getLogger(FastqReformater.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        System.out.println("Processed " + itCounter + " reads.");
    }

    public void swapBases(File fastqIn, File fastqOut) throws CharacterParseException
    {
        CharacterReplacer cReplacer = Perl5Parser.makeReplacer("tr/ATCG/CGAT/");
        int itCounter = 0;

        FastqReader fqr = new FastqReader(fastqIn);
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter out = writer.newWriter(fastqOut);


        Iterator it = fqr.iterator();

        while (it.hasNext())
        {
            FastqRecord seqRecord = (FastqRecord) it.next();
            
            String seq = seqRecord.getReadString();
            String newSeq = cReplacer.doReplacement(seq);
            
            FastqRecord newSeqRecord = new FastqRecord(seqRecord.getReadHeader(), newSeq, "", seqRecord.getBaseQualityString());
            out.write(newSeqRecord);
            itCounter++;
        }
        out.close();
        System.out.println("Processed " + itCounter + " reads.");
    }

    
}
