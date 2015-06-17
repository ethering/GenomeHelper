/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.vcf;

import java.util.HashSet;

/**
 *
 * @author ethering
 */
public class SnpLoci
{

    private String chr;
    private int pos;
    private String ref;
    private HashSet alts;

    

    public SnpLoci(String chr, int pos, String ref, HashSet alts)
    {
        this.chr = chr;
        this.pos = pos;
        this.ref = ref;
        this.alts = alts;

    }

    public String getChr()
    {
        return chr;
    }

    public void setChr(String chr)
    {
        this.chr = chr;
    }

    public int getPos()
    {
        return pos;
    }

    public void setPos(int pos)
    {
        this.pos = pos;
    }

    public String getRef()
    {
        return ref;
    }

    public void setRef(String ref)
    {
        this.ref = ref;
    }

    public HashSet getAlts()
    {
        return alts;
    }

    public void setAlts(HashSet alts)
    {
        this.alts = alts;
    }

    @Override
    public String toString()
    {
        String s = this.chr + "\t"+ this.pos + "\t"+ this.ref + "\t"+ this.alts.toString();
        return s;
    }
}
