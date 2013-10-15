package uk.ac.tsl.etherington.piculus.fasta;

import java.util.ArrayList;
import java.util.Iterator;
import net.sf.samtools.SAMRecord;

/**
 *This class represents the co-ordinates of any feature or region of a contig, chromosome, scaffold, etc. 
 * It contains a reference name (contigId) and start and end position and a strand ('+' by default). 
 * The start and end positions are the 1-based co-ordinates of the first and last nucleotide, so the length of a given ContigPosition can be calculated as end-start+1.
 * 
 * @author Graham Etherington
 */
public class ContigPosition 
{

    String contigid;
    int start;
    int end;
    char strand;

    /**
     * Creates a new ContigPosition object with strand set to '+'
     *
     * @param contigid  the name of the reference
     * @param start  the nucleotide start position on the reference
     * @param end  the nucleotide end position on the reference
     *
     */
    public ContigPosition(String contigid, int start, int end)
    {
        this.contigid = contigid;
        this.start = start;
        this.end = end;
        this.strand = '+';
    }

    /**
     * Creates a ContigPosition with a defined strand
     *
     * @param contigid  the name of the reference
     * @param start  the nucleotide start position on the reference
     * @param end  the nucleotide end position on the reference
     * @param strand  the strand on which the contig sits (can only be '+' or
     * '-')
     */
    public ContigPosition(String contigid, int start, int end, char strand)
    {
        boolean goodChar = validateStrand(strand);
        if (goodChar)
        {
            this.contigid = contigid;
            this.start = start;
            this.end = end;
            this.strand = strand;
        }
    }

    /**
     * Creates a ContigPosition from a net.sf.samtools.SAMRecord object
     *
     * @param samRecord  a net.sf.samtools.SAMRecord object
     */
    public ContigPosition(SAMRecord samRecord)
    {
        this.contigid = samRecord.getReferenceName();
        this.start = samRecord.getAlignmentStart();
        this.end = samRecord.getAlignmentEnd();

        boolean onMinusStrand = samRecord.getReadNegativeStrandFlag();

        if (onMinusStrand)
        {
            this.strand = '-';
        }
        else
        {
            this.strand = '+';
        }
    }
    
     /**
     * Validates that the strand char is '+' or '-'
     * @return  returns true if the strand is '+' or '-'. Otherwise returns
     * false
     */
    public final boolean validateStrand(char strand)
    {
        if (strand != '+' || strand != '-')
        {
            System.err.println("Error: strand must be '+' or '-' ");
            return false;
        }
        else
        {
            return true;
        }
    }

    /**
     * Sets the strand of the ContigPosition
     *
     * @param strand  '+' for plus strand or '-' for minus strand
     */
    public void setStrand(char strand)
    {
        boolean goodChar = validateStrand(strand);
        if (goodChar)
        {
            this.strand = strand;
        }      
    }

    /**
     * 
     * @return  returns the strand char ('+' or '-')
     */
    public char getStrand()
    {
        return strand;
    }

    /**
     *
     * @return The reference name
     */
    public String getContigid()
    {
        return contigid;
    }

    /**
     *
     * @return The start position
     */
    public int getStart()
    {
        return start;
    }

    /**
     *
     * @return The end position
     */
    public int getEnd()
    {
        return end;
    }

    /**
     * Sets the name of the reference sequence for a given ContigPosition
     * @param contigid The reference name
     */
    public void setContigid(String contigid)
    {
        this.contigid = contigid;
    }

    /**
     * Sets the start position for a given ContigPosition
     * @param start The start position on the reference
     */
    public void setStart(int start)
    {
        this.start = start;
    }

    /**
     *Sets the end positions for a given ContigPosition
     * @param end The end position on the reference
     */
    public void setEnd(int end)
    {
        this.end = end;
    }

    /**
     *
     * @param contigPosition the ContigPosition object to compare
     * @return true if the two ContigPosition objects are equal, false if not
     */
    @Override
    public boolean equals(Object contigPosition)
    {
        if (contigPosition == null)
        {
            return false;
        }
        if (getClass() != contigPosition.getClass())
        {
            return false;
        }
        final ContigPosition other = (ContigPosition) contigPosition;
        if ((this.contigid == null) ? (other.contigid != null) : !this.contigid.equals(other.contigid))
        {
            return false;
        }
        if (this.start != other.start)
        {
            return false;
        }
        if (this.end != other.end)
        {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode()
    {
        int hash = 7;
        hash = 79 * hash + (this.contigid != null ? this.contigid.hashCode() : 0);
        hash = 79 * hash + this.start;
        hash = 79 * hash + this.end;
        return hash;
    }


    /**
     * Tests to see if two ContigPosition objects overlap on the same strand
     * @param contigPositionToCompare the ContigPosition object to compare
     * @return true if they overlap, false if they don't
     */
    
    public boolean overlaps(ContigPosition contigPositionToCompare)
    {
        if (this.getContigid().equalsIgnoreCase(contigPositionToCompare.getContigid()) && this.getStrand() == contigPositionToCompare.getStrand() && ((this.getStart() >= contigPositionToCompare.getStart() && this.getStart() <= contigPositionToCompare.getEnd()) || (this.getEnd() <= contigPositionToCompare.getEnd() && this.getEnd() >= contigPositionToCompare.getStart()) || (this.getStart() <= contigPositionToCompare.getStart() && this.getEnd() >= contigPositionToCompare.getEnd()) || (this.getStart() >= contigPositionToCompare.getStart() && this.getEnd() <= contigPositionToCompare.getEnd())))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    /**
     * Tests if a ContigPosition object already exists in an ArrayList of ContigPosition objects
     * @param cpArray  an ArrayList of ContigPosition objects
     * @return true if the ContigPosition exists in the ArrayList, false if it doesn't
     */
    
    public boolean exists(ArrayList<ContigPosition> cpArray)
    {
        Iterator it = cpArray.iterator();
        while (it.hasNext())
        {
            if (this.equals(it.next()))
            {
                return true;
            }
        }
        return false;

    }
    /**
     * Creates a human readable version of a ContigPosition
     * @return a human readable version of a ContigPosition
     */

    @Override
    public String toString()
    {
        String str = this.contigid.concat(":").concat(Integer.toString(this.start)).concat("-").concat(Integer.toString(this.end)).concat(String.valueOf(strand));

        return str;
    }
}
