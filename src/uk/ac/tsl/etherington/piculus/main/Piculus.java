/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.piculus.main;

/**
 *
 * @author ethering
 */
public class Piculus
{

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args)
    {
         if (args[0].equalsIgnoreCase("-help") || args[0].equalsIgnoreCase("-h"))
        {
            System.out.println("Welcome to Piculus. Available programs:");

            System.out.println("Use: FastaUtils get_longest_subseqs  (String [] fastafiles),  File fastaOut");
            System.out.println("Use: FastaUtils translate_dna_counts  File infile,  File outfile");
            System.out.println("Use: FastaUtils gc_content  String fastafile");
            System.out.println("Use: FastaUtils seq_from_command_line  String fastafile, String outfile, String seqId");
            System.out.println("Use: FastaUtils seq_from_command_line  String fastafile, String outfile, String seqId, int subseqStart, int subseqEnd");
            System.out.println("Use: FastaUtils select_random_seqs  String fastafile, int numberOfSeqs, String outfile");

        }
        else if (args[0].equalsIgnoreCase("get_longest_subseqs"))
        {
            
        }
    }
}
