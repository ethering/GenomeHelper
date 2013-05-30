/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package piculus.fastq;

import java.io.File;
import java.util.HashSet;
import java.util.Iterator;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

/**
 *
 * @author ethering
 */
public class GetUnmappedSeqs
{

    //get a HashSet of all the sequences
    public HashSet getListOfSequences(File inLeft, File inRight)
    {
        HashSet<String> seqList = new HashSet<String>();

        FastqReader fql = new FastqReader(inLeft);
        FastqReader fqr = new FastqReader(inRight);

        Iterator itl = fql.iterator();
        Iterator itr = fqr.iterator();

        while (itl.hasNext())
        {
            FastqRecord leftSeqRecord = (FastqRecord) itl.next();
            FastqRecord rightSeqRecord = (FastqRecord) itr.next();
            seqList.add(leftSeqRecord.getReadHeader());
            seqList.add(rightSeqRecord.getReadHeader());

        }
        fql.close();
        fqr.close();
        System.out.println("There are " + seqList.size() + " sequences");
        return seqList;
    }

    //takes a bamfile and a hashset of reads in fastq file and removes the names of any reads that are mapped
    public HashSet getListOfUnmappedSeqs(HashSet readList, File bamFile)
    {

        final SAMFileReader reader = new SAMFileReader(bamFile);
        reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        // Open an iterator for the particular sequence
        SAMRecordIterator iterator = reader.iterator();

        while (iterator.hasNext())
        {
            SAMRecord samRecord = iterator.next();
            boolean unMapped = samRecord.getReadUnmappedFlag();
            if (!unMapped)
            {
                String readName = samRecord.getReadName();
                boolean firstMate = samRecord.getFirstOfPairFlag();
                if (firstMate)
                {
                    readName = readName.concat("/1");
                }
                else
                {
                    readName = readName.concat("/2");
                }
                //System.out.println(readName);
                if (readList.contains(readName))
                {
                    readList.remove(readName);
                }
            }

        }
        //should be a list of unmapped sequences remaining
        System.out.println("There are " + readList.size() + " unmapped sequences");
        return readList;
    }

    public void getUnmappedSeqs(HashSet unmappedSeqs, File inLeft, File inRight, File outLeft, File outRight)
    {
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter leftSeqs = writer.newWriter(outLeft);
        FastqWriter rightSeqs = writer.newWriter(outRight);
        FastqReader fql = new FastqReader(inLeft);
        FastqReader fqr = new FastqReader(inRight);
        Iterator itl = fql.iterator();
        Iterator itr = fqr.iterator();
        int writtenSeqs = 0;

        while (itl.hasNext())
        {
            FastqRecord leftSeqRecord = (FastqRecord) itl.next();
            FastqRecord rightSeqRecord = (FastqRecord) itr.next();
            String leftReadName = leftSeqRecord.getReadHeader();
            String rightReadName = rightSeqRecord.getReadHeader();
            if (unmappedSeqs.contains(leftReadName))
            {
                leftSeqs.write(leftSeqRecord);
                writtenSeqs++;
            }
            if (unmappedSeqs.contains(rightReadName))
            {
                rightSeqs.write(rightSeqRecord);
                writtenSeqs++;
            }
        }
        System.out.println("Written " + writtenSeqs + " to files");
    }

    public void getUnmappedPairedSeqs(HashSet unmappedSeqs, File inLeft, File inRight, File outLeft, File outRight)
    {
        FastqWriterFactory writer = new FastqWriterFactory();
        FastqWriter leftSeqs = writer.newWriter(outLeft);
        FastqWriter rightSeqs = writer.newWriter(outRight);
        FastqReader fql = new FastqReader(inLeft);
        FastqReader fqr = new FastqReader(inRight);
        FastqRecord leftSeqRecord;
        FastqRecord rightSeqRecord;
        String readName;
        String[] readNameArray;
        Iterator itl = fql.iterator();
        Iterator itr = fqr.iterator();
        int writtenSeqs = 0;

        while (itl.hasNext())
        {
            leftSeqRecord = (FastqRecord) itl.next();
            rightSeqRecord = (FastqRecord) itr.next();
            readName = leftSeqRecord.getReadHeader();
            readNameArray = readName.split(" ");
            readName = new String(readNameArray[0]);

            if (unmappedSeqs.contains(readName))
            {
                leftSeqs.write(leftSeqRecord);
                rightSeqs.write(rightSeqRecord);
                unmappedSeqs.remove(readName);
                writtenSeqs++;
            }

        }
        leftSeqs.close();
        rightSeqs.close();

        System.out.println("Written " + writtenSeqs + " to files");
    }

    public HashSet getListOfUnmappedPairedSeqs(File bamFile)
    {
        HashSet<String> readList = new HashSet<String>();
        final SAMFileReader reader = new SAMFileReader(bamFile);
        reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        SAMRecord samRecord;
        boolean unMapped;
        boolean mateUnmapped;
        // Open an iterator for the particular sequence
        //maybe iterate over each sequence
        SAMRecordIterator iterator = reader.iterator();

        while (iterator.hasNext())
        {
            samRecord = iterator.next();
            unMapped = samRecord.getReadUnmappedFlag();
            mateUnmapped = samRecord.getMateUnmappedFlag();
            if (unMapped && mateUnmapped)
            {
                readList.add(samRecord.getReadName());
            }
        }
        reader.close();
        iterator.close();
        //should be a list of unmapped sequences remaining
        System.out.println("There are " + readList.size() + " unmapped sequences");
        return readList;
    }

//    public static void main(String[] args)
//    {
//        File bam = new File("/Users/ethering/temp/bam/Galaxy5-[Map_with_BWA_for_Illumina_on_data_4_and_data_3__mapped_reads].sam");
//        String str = " ";
//        GetUnmappedSeqs gus = new GetUnmappedSeqs();
//        HashSet<String> hs = new HashSet<String>();
//        hs = gus.getListOfUnmappedPairedSeqs(bam);
//
//    }
}
