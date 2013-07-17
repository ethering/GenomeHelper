/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package uk.ac.tsl.etherington.piculus.fastq;

import net.sf.picard.fastq.FastqRecord;

/**
 *
 * @author ethering
 */

public class PEFastqRead 
{
    String readName;
    FastqRecord leftRead;
    FastqRecord rightRead;

    public PEFastqRead(String readName, FastqRecord leftRead, FastqRecord rightRead)
    {
        this.readName = readName;
        this.leftRead = leftRead;
        this.rightRead = rightRead;
    }

    public FastqRecord getLeftRead()
    {
        return leftRead;
    }

    public void setLeftRead(FastqRecord leftRead)
    {
        this.leftRead = leftRead;
    }

    public String getReadName()
    {
        return readName;
    }

    public void setReadName(String readName)
    {
        this.readName = readName;
    }

    public FastqRecord getRightRead()
    {
        return rightRead;
    }

    public void setRightRead(FastqRecord rightRead)
    {
        this.rightRead = rightRead;
    }
    
    

}
