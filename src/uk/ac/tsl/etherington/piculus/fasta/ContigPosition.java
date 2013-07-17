package uk.ac.tsl.etherington.piculus.fasta;

import java.util.ArrayList;
import java.util.Iterator;
import net.sf.samtools.SAMRecord;

/**
 *
 * @author Graham Etherington
 */
public class ContigPosition implements Comparable<ContigPosition>
{

    String contigid;
    int start;
    int end;
    char strand;

    /**
     * Creates a new ContigPosition object where the String is the name of the
     * reference/contig and the int is a position on the reference/contig
     *
     * @param contigid The name of the reference
     * @param start The nucleotide position on the reference
     */
    public ContigPosition(String contigid, int start, int end)
    {
        this.contigid = contigid;
        this.start = start;
        this.end = end;
        this.strand = '+';
    }

    public ContigPosition(String contigid, int start, int end, char strand)
    {
        this.contigid = contigid;
        this.start = start;
        this.end = end;
        this.strand = strand;
    }
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
    

    public void setStrand(char strand)
    {
        this.strand = strand;
    }

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
     * @return The e position
     */
    public int getEnd()
    {
        return end;
    }

    /**
     *
     * @param contigid The reference name
     */
    public void setContigid(String contigid)
    {
        this.contigid = contigid;
    }

    /**
     *
     * @param start The start position on the reference
     */
    public void setStart(int start)
    {
        this.start = start;
    }

    /**
     *
     * @param end The end position on the reference
     */
    public void setEnd(int end)
    {
        this.end = end;
    }

    /**
     *
     * @param obj the object to compare
     * @return true if the objects are equal, false if not
     */
    @Override
    public boolean equals(Object obj)
    {
        if (obj == null)
        {
            return false;
        }
        if (getClass() != obj.getClass())
        {
            return false;
        }
        final ContigPosition other = (ContigPosition) obj;
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

    public int compareTo(ContigPosition t)
    {
        return this.getStart() - t.getStart();
    }

    public boolean overlaps(ContigPosition other)
    {
        if (this.getContigid().equalsIgnoreCase(other.getContigid()) && this.getStrand() == other.getStrand() && ((this.getStart() >= other.getStart() && this.getStart() <= other.getEnd() )|| (this.getEnd() <= other.getEnd() && this.getEnd() >= other.getStart()) || (this.getStart() <= other.getStart() && this.getEnd() >= other.getEnd()) || (this.getStart() >= other.getStart() && this.getEnd() <= other.getEnd())))
        {
            return true;
        }
        else
        {
            return false;
        }


    }

    public boolean exists(ArrayList<ContigPosition> cp)
    {
        Iterator it = cp.iterator();
        while (it.hasNext())
        {
            if (this.equals(it.next()))
            {
                return true;
            }
        }
        return false;

    }
    
    @Override
    public String toString()
    {
        String str = this.contigid.concat(":").concat(Integer.toString(this.start)).concat("-").concat(Integer.toString(this.end)).concat(String.valueOf(strand));
        
        return str;
    }

    

    
}
