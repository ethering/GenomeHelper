package uk.ac.tsl.etherington.piculus.gff;


import uk.ac.tsl.etherington.piculus.fasta.ContigPosition;
import org.biojava3.genome.parsers.gff.Feature;
import org.biojava3.genome.parsers.gff.Location;



/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ethering
 */

public class GFFObject extends org.biojava3.genome.parsers.gff.Feature
{
    public char strand;

    public GFFObject(String seqname, String source, String type, Location location, Double score, char strand, int frame, String attributes)
    {
        super(seqname, source, type, location, score, frame, attributes);
        this.strand = strand;
    }

    public GFFObject(Feature feature, char strand)
    {
        super(feature);
        this.strand = strand;
    }
    
    public GFFObject(Feature feature)
    {
        super(feature);
        strand = this.location().bioStrand();
    }
    

    @Override
    public String toString()
    {
        return (this.seqname()+"\t"+this.source()+"\t"+this.type()+"\t"+this.location().bioStart()+"\t"+this.location().bioEnd()+"\t"+this.score()+"\t"+this.location().bioStrand()+"\t"+this.frame()+"\t"+this.attributes());
    
    }
    
    public ContigPosition toContigPosition()
    {
        
        ContigPosition cp = new ContigPosition(this.seqname(), this.location().bioStart(), this.location().bioEnd(), this.strand);
        return cp;
        
    }

    

}
