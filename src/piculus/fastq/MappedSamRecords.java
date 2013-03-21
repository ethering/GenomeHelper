package piculus.fastq;



/*
 * Takes a SAM File, prints the header and then prints any entries where the
 * read or its mate (or both) are mapped
 */
import java.io.*;
import java.util.HashSet;
import java.util.Iterator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import piculus.fasta.ContigPosition;

/**
 *
 * @author ethering
 */
public class MappedSamRecords
{

    public HashSet getMappedSamRecords(String bamFileString)
    {

        File bamFile = new File(bamFileString);
        HashSet<String> mappedReads = new HashSet<String>();
        int mappedReadCount = 0;

        final SAMFileReader reader = new SAMFileReader(bamFile);
        reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        // Open an iterator for the particular sequence
        SAMRecordIterator iterator = reader.iterator();

        while (iterator.hasNext())
        {
            SAMRecord samRecord = iterator.next();
            //is the read mapped?
            boolean unMapped = samRecord.getReadUnmappedFlag();

            //..if so...
            if (unMapped == false)
            {
                //get the read name
                String currentRead = samRecord.getReadName();
                //convert it from HWI-EAS396:4:1:4:1832#0/1_1 to HWI-EAS396:4:1:4:1832#
                int hashIndex = currentRead.indexOf("#");
                currentRead = currentRead.substring(0, hashIndex);
                mappedReads.add(currentRead);
                mappedReadCount++;

            }
        }
        System.out.println("Found " + mappedReadCount + " mapped reads");
        return mappedReads;
    }

    public HashSet getMappedSamRecordsFromCuffdiff(String bamFileString, String cuffdiffFileString) throws FileNotFoundException, IOException
    {

        File bamFile = new File(bamFileString);
        File cuffDiffFile = new File(cuffdiffFileString);
        HashSet<ContigPosition> highEx = new HashSet<ContigPosition>();
        BufferedReader input = new BufferedReader(new FileReader(cuffDiffFile));
        try
        {
            String line = null; //not declared within while loop
        /*
             * readLine is a bit quirky : it returns the content of a line MINUS
             * the newline. it returns null only for the END of the stream. it
             * returns an empty String if two newlines appear in a row.
             */
            while ((line = input.readLine()) != null)
            {
                
                String[] array = line.split("\\t");
                
                String significance = array[array.length - 1];
                if (significance.equalsIgnoreCase("yes"))
                {
                    String locus = array[3];
                    //System.out.println("Locus = " + locus);
                    String chr = locus.substring(0, locus.indexOf(":"));
                    String coords = locus.substring(locus.indexOf(":") + 1, locus.length());
                    String[] cordsArray = coords.split("-");
                    int start = Integer.parseInt(cordsArray[0]);
                    int end = Integer.parseInt(cordsArray[1]);
                    ContigPosition cp = new ContigPosition(chr, start, end);
                    highEx.add(cp);
                    //System.out.println("Chr = " + chr + " start = " + start + " end = " + end);
                }
            }
        }
        finally
        {
            input.close();
        }

        HashSet<String> mappedReads = new HashSet<String>();
        int mappedReadCount = 0;

        final SAMFileReader reader = new SAMFileReader(bamFile);
        reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        // Open an iterator for the particular sequence
        SAMRecordIterator iterator = reader.iterator();

        while (iterator.hasNext())
        {
            SAMRecord samRecord = iterator.next();
            //is the read mapped?
            boolean unMapped = samRecord.getReadUnmappedFlag();

            //..if so...
            if (unMapped == false)
            {
                //get the read name
                String currentRead = samRecord.getReadName();
                String ref = samRecord.getReferenceName();
                int alStart = samRecord.getAlignmentStart();
                int alEnd = samRecord.getAlignmentEnd();
                Iterator it = highEx.iterator();
                while (it.hasNext())
                {
                    ContigPosition cp = (ContigPosition) it.next();
                    String cpChr = cp.getContigid();
                    int cpStart = cp.getStart();
                    int cpEnd = cp.getEnd();
                    if (ref.equalsIgnoreCase(cpChr) && ((alStart >= cpStart && alStart <= cpEnd) || (alEnd >= cpStart && alEnd <= cpEnd)))
                    {
                        //convert it from HWI-EAS396:4:1:4:1832#0/1_1 to HWI-EAS396:4:1:4:1832#
                        int hashIndex = currentRead.indexOf("#");
                        currentRead = currentRead.substring(0, hashIndex);
                        mappedReads.add(currentRead);
                        mappedReadCount++;
                    }
                }
            }
        }
        System.out.println("Found " + mappedReadCount + " mapped reads");
        return mappedReads;
    }
    
    
}