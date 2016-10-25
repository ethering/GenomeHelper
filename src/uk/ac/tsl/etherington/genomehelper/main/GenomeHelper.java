/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.tsl.etherington.genomehelper.main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.jtr.transliterate.CharacterParseException;
import uk.ac.tsl.etherington.genomehelper.bam.MappedSamRecords;
import uk.ac.tsl.etherington.genomehelper.fasta.FastaFeatures;
import uk.ac.tsl.etherington.genomehelper.fasta.FastaMotifFinder;
import uk.ac.tsl.etherington.genomehelper.fasta.FastaParser;
import uk.ac.tsl.etherington.genomehelper.fasta.FastaSubstrings;
import uk.ac.tsl.etherington.genomehelper.fasta.FastaTranslator;
import uk.ac.tsl.etherington.genomehelper.fasta.RandomFasta;
import uk.ac.tsl.etherington.genomehelper.fastq.FastqCompression;
import uk.ac.tsl.etherington.genomehelper.fastq.FastqInterlacer;
import uk.ac.tsl.etherington.genomehelper.fastq.FastqJoiner;
import uk.ac.tsl.etherington.genomehelper.fastq.FastqMotifFinder;
import uk.ac.tsl.etherington.genomehelper.fastq.FastqParser;
import uk.ac.tsl.etherington.genomehelper.fastq.FastqQC;
import uk.ac.tsl.etherington.genomehelper.gff.GFFFeatureStats;
import uk.ac.tsl.etherington.genomehelper.utils.Interval;
import uk.ac.tsl.etherington.genomehelper.vcf.VCFParser;

/**
 *
 * @author ethering
 */
public class GenomeHelper
{

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     * @throws org.jtr.transliterate.CharacterParseException
     *
     */
    public static String footer = "\nPlease report issues at https://github.com/ethering/GenomeHelper/issues";

    @SuppressWarnings("static-access")
    public static
            void main(String[] args) throws IOException, CharacterParseException, Exception
    {
        if (args.length == 0 || args[0].equalsIgnoreCase("-help") || args[0].equalsIgnoreCase("-h"))
        {
            System.out.println("Welcome to GenomeFinder\nThe following programs are available. Use <program_name> -h for help with each program:\n");
            System.out.println("The following options are used in multiple categories of programs");
            System.out.println("-in             an input file pertinent to the method being used");
            System.out.println("-out            an output file pertinent to the method being used");

            System.out.println("\nFasta-related programs:");
            System.out.println("The following options are specific to this section");
            System.out.println("-fasta          a fasta input file");

            System.out.println("Usage: FastaMotifFinder -fasta -searchMotif -motifCounts -proteinCounts -minCount");
            System.out.println("-searchMotif    a string motif to search for");
            System.out.println("-motifCounts   a tab-delimited file of the occurrence of each motif");
            System.out.println("-proteinCounts a tab-delimited file of the occurrence of each pattern when translated into ammino acids");
            System.out.println("-minCount      the minimum number of times a motif must be found to be included in the results");

            System.out.println("Usage: MultiFastaToSingleFasta -fastain -fastaout");
            System.out.println("-fastain    input fasta file containing multiple fasta sequences");
            System.out.println("-fastaout   output fasta containing single concatenated sequence");

            System.out.println("Usage: FastaGetLengths -fastain ");
            System.out.println("-fastain    input fasta file containing fasta sequence(s)");

            System.out.println("Usage: FastaGetLongestSubstring  <path to files>  -out.");
            System.out.println("Usage: FastaGetGenomeLength  infile");
            System.out.println("Usage: FastaTranslate infile outfile");
            System.out.println("Usage: FastaGetGCContent  fastafile");
            System.out.println("Usage: FastaGetSingleFromMultiFasta  infile outfile  seqId subsequence_start (optional) subsequence_end (optional)");
            System.out.println("Usage: FastaSelectRandomSequences fastaIn numberOfRandomSeqs randomSeqsoutfile");
            System.out.println("Usage: FastaToFastq fastaIn fastqOut");
            System.out.println("Usage: FastaFilterByLength fastaIn fastaOut minLength");
            System.out.println("Usage: GetNStats fastaIn");

            System.out.println("\nFastq-related programs:");
            System.out.println("Usage: FastqCompress fastqIn fastqOut.gz");
            System.out.println("Usage: FastqInterlace leftReads rightReads interlacedFastqFile singlesFile");
            System.out.println("Usage: FastqDeinterlace interlacedFastqFile leftReads rightReads  leftSinglesFile rightSinglesFile");
            System.out.println("Usage: FastqInterlaceKnownPairs leftReads rightReads interlacedFastqFile");
            System.out.println("Usage: FastqJoin leftReads rightReads fastqJoinedFile fastqSinglesFile");
            System.out.println("Usage: FastqSplit joinedFastqFile leftPairdReads rightPairedReads");
            System.out.println("Usage: FastqMotifFinder fastqFile searchMotif motifCountsFile proteinCountsFile minCount");
            System.out.println("Usage: FastqGetPairedEndSequencesWithMotifMatch interlacedFastqFile searchMotif matchingFastqFiles");
            System.out.println("Usage: FastqGetPairedEndSequencesFromFile listFile fastqFileInLeft fastqFileInRight fastqFileOutLeft fastqFileOutRight");
            System.out.println("Usage: FastqGetSingleEndSequencesFromFile listFile fastqFileIn fastqFileOut");
            System.out.println("Usage: FastqToFasta fastqIn fastaOut ");
            System.out.println("Usage: FastqTranslate fastqIn fastaOut includeOriginalDNASequence ('true' or 'false')");
            System.out.println("Usage: FastqCountNucleotides fastqIn");
            System.out.println("Usage: FastqFindKmer fastqIn kmer");

            System.out.println("\nQuality-control programs:");
            System.out.println("Usage: QCPairedReads fastqInLeft fastqInRight fastqLeftOut fastqrRightOut singleReads readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')");
            System.out.println("Usage: QCSingleEndReads fastqIn fastqOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')");
            System.out.println("Usage: RemoveNsFromPairedReads fastqInLeft fastqInRight fastqLeftOut fastqrRightOut");
            System.out.println("Usage: QCInterlacedReads fastqIn fastqOut singleReads readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')");
            System.out.println("Usage: QCInterlacedReadsToPairs fastqIn fastqLeftOut fastqrRightOut singleReads readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')");
            System.out.println("Usage: QCJoinedReads fastqIn fastqLeftOut fastqrRightOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')");
            System.out.println("Usage: QCVerifyReads fastqIn");
            System.out.println("Usage: QCVerifyPairedEndReads fastqLeft fastqRight");
            System.out.println("Usage: QCRemoveKmerPairedReads fastqInLeft fastqInRight fastqLeftOut fastqrRightOut kmerFile");
            System.out.println("Usage: QCRemoveKmerSingleReads fastqIn fastqOut kmerFile");

            System.out.println("\nSAM/BAM-related programs:");
            System.out.println("Usage: BAMGetPairedUnmappedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
            System.out.println("Usage: BAMGetPairedMappedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
            System.out.println("Usage: BAMGetBothPairedUnmappedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
            System.out.println("Usage: BAMListBothPairedUnmappedReads bamfile listFile");
            System.out.println("Usage: BAMListSingleUnmappedPairedReads bamfile listFile");
            System.out.println("Usage: GetReadsFromList listFile, leftFastq, rightFastq, readsOut");

            System.out.println("Usage: BAMGetSingleUnmappedPairedReads bamfile, fastqInLeft, fastqInRight, fastqSingles");
            System.out.println("Usage: BAMGetSingleMappedReads bamfile, fastqInLeft, fastqInRight, fastqSingles");

            System.out.println("\nGFF-related programs:");
            System.out.println("Usage: GFFGetMeanFeatureLengthWithSplicing gffFile featureName refSeq");
            System.out.println("Usage: GFFGetMeanFeatureLength gffFile featureName");
            System.out.println("Usage: GFFGetMeanFeatureLengthOfGeneIDs gffFile featureName fileOfIds");
            System.out.println("Usage: GFFCreateNonCodingGenome gffFile refSeq nonCodingGenome");
            System.out.println("Usage: GFFCreateCodingGenome gffFile featureName refSeq codingGenome");
            System.out.println("Usage: GFFGetMeanIntronLength gffFile featureName refSeq ");
            System.out.println("Usage: GFFGetMeanTargetIntronLength gffFile featureName targets");
            System.out.println("Usage: GFFCalculateCodingRegion gffFile refSeq attribute");
            System.out.println("Usage: GFFGetStats gffFile refSeq attribute");

            System.out.println("\nVCF-related programs");
            System.out.println("The following options are specific to this section");

            System.out.println("Usage: CalculateGATKparams -in -max");
            System.out.println("-in a VCF file");
            System.out.println("-max (Optional) maximum number of VCF entries to analyse (default = all)");

            System.out.println("\nOther Utility programs:");
            System.out.println("Usage: gatkToSamInterval bam gatkInterval");

        }

        else if (args[0].equalsIgnoreCase("FastaMotifFinder"))
        {
//            if (args[1].equalsIgnoreCase("-h"))
//            {
//                System.out.println("Usage: FastaMotifFinder -fasta -searchMotif -motifCounts -proteinCounts -minCount");
//                System.out.println("-searchMotif - the motif to search for. The motif can contain regular expressions, such as ATG[CG]G[AT],"
//                        + " which will search for ATGCGA, ATGCGT, ATGGGA and ATGGGT. Refer to java.util.regex.Pattern for more details.");
//                System.out.println("-motifCounts   a tab-delimited file of the occurrence of each motif");
//                System.out.println("-proteinCounts a tab-delimited file of the occurrence of each pattern when translated into ammino acids");
//                System.out.println("-minCount      the minimum number of times a motif must be found to be included in the results");
//
//            } else
//            {

            // create Options object
            Options options = new Options();
            options.addOption(OptionBuilder.withArgName("file")
                    .hasArg()
                    .isRequired()
                    .withDescription("the fasta file to search")
                    .create('f'));
            options.addOption(OptionBuilder.withArgName("string or regex")
                    .hasArg()
                    .isRequired()
                    .withDescription("the motif to search for. The motif can contain regular expressions, such as ATG[CG]G[AT],"
                                     + " which will search for ATGCGA, ATGCGT, ATGGGA and ATGGGT. Refer to java.util.regex.Pattern for more details.")
                    .create('s'));
            options.addOption(OptionBuilder.withArgName("file")
                    .hasArg()
                    .isRequired()
                    .withDescription("a tab-delimited file of the occurrence of each motif")
                    .create('m'));
            options.addOption(OptionBuilder.withArgName("file")
                    .hasArg()
                    .isRequired()
                    .withDescription("a tab-delimited file of the occurrence of each pattern when translated into ammino acids")
                    .create('p'));
            options.addOption(OptionBuilder.withArgName("int")
                    .hasArg()
                    .isRequired()
                    .withDescription("the minimum number of times a motif must be found to be included in the results")
                    .create('c'));
            String header = "Searches fasta file an input file\n\n";
            String footer = "\nPlease report issues to me";

            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("FastaMotifFinder", header, options, footer, true);

            CommandLineParser parser = new BasicParser();
            CommandLine cmd = parser.parse(options, args);

            File fastaFile = new File(cmd.getOptionValue("fasta"));
            String searchMotif = cmd.getOptionValue("searchMotif");
            File motifCounts = new File(cmd.getOptionValue("motifCounts"));
            File aaMotifCounts = new File(cmd.getOptionValue("proteinCounts"));
            int minCount = Integer.parseInt(cmd.getOptionValue("minCount"));
            FastaMotifFinder fmf = new FastaMotifFinder();
            fmf.findMatches(fastaFile, searchMotif, motifCounts, aaMotifCounts, minCount);

            // }
        }
        else if (args[0].equalsIgnoreCase("CalculateGATKparams"))
        {

            // create Options object
            Options options = new Options();
            options.addOption(OptionBuilder.withArgName("file")
                    .hasArg()
                    .isRequired()
                    .withDescription("the VCF file to search")
                    .create("in"));
            options.addOption(OptionBuilder.withArgName("Max records")
                    .hasArg()
                    .withDescription("Optional: The maximum number of VCF lines to read (default = all). Using a lower number, e.g. 500,000 may speed up analysis for large files "
                                     + "whilst giving very similar results to examining all entries")
                    .create("max"));
            options.addOption(OptionBuilder.withLongOpt("help").create('h'));

            String header = "Calculates the mean and standard lower standard deviations for various vcf fields.\n"
                            + "These can be used for input parameters into the GATK VariantFiltration tool\n";

            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("CalculateGATKparams", header, options, footer, false);
            CommandLineParser parser = new BasicParser();
            CommandLine cmd = parser.parse(options, args);
            File in = new File(cmd.getOptionValue("in"));
            VCFParser vcfParser = new VCFParser();

            if (cmd.hasOption("max"))
            {
                System.out.println("Using -max option");
                int max = Integer.parseInt(cmd.getOptionValue("max"));
                vcfParser.calculateGATKParams(in, max);
            }

            else
            {
                vcfParser.calculateGATKParams(in);
            }
        }

        else if (args[0].equalsIgnoreCase("MultiFastaToSingleFasta"))
        {

            // create Options object
            Options options = new Options();
            options.addOption(OptionBuilder.withArgName("infile")
                    .hasArg()
                    .isRequired()
                    .withDescription("input fasta file containing multiple fasta sequences")
                    .create('f'));
            options.addOption(OptionBuilder.withArgName("outfile")
                    .hasArg()
                    .isRequired()
                    .withDescription("a tab-delimited file of the occurrence of each pattern when translated into ammino acids")
                    .create('o'));
            options.addOption(OptionBuilder.withLongOpt("help").create('h'));

            String header = "Concatenates all fasta sequences in multifasta file to a single fasta sequence\n";

            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("MultiFastaToSingleFasta", header, options, footer, false);
            CommandLineParser parser = new BasicParser();
            CommandLine cmd = parser.parse(options, args);

            File in = new File(cmd.getOptionValue("f"));
            File out = new File(cmd.getOptionValue("o"));

            FastaParser fp = new FastaParser();
            fp.multiFastaToSingleFasta(in, out);

            // }
        }

        else if (args[0].equalsIgnoreCase("FastaGetLengths"))
        {

            // create Options object
            Options options = new Options();
            options.addOption(OptionBuilder.withArgName("infile")
                    .hasArg()
                    .isRequired()
                    .withDescription("input fasta file containing multiple fasta sequences")
                    .create('f'));

            options.addOption(OptionBuilder.withLongOpt("help").create('h'));

            String header = "Calculates lengths of sequences in fasta file\n";

            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("FastaGetLengths", header, options, footer, false);
            CommandLineParser parser = new BasicParser();
            CommandLine cmd = parser.parse(options, args);

            File in = new File(cmd.getOptionValue("f"));

            FastaFeatures fp = new FastaFeatures();
            HashMap<String, Integer> seqLengths = new HashMap<>(FastaFeatures.getSequenceLengths(in));
            for (Map.Entry<String, Integer> entry : seqLengths.entrySet())
            {
                String seqName = entry.getKey();
                int genomeLength = entry.getValue();
                System.out.println("Sequence name: " + seqName);
                System.out.println("Length: " + genomeLength);

            }

            // }
        }

        else if (args[0].equalsIgnoreCase("FastaGetLongestSubstring"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaGetLongestSubstring <fastaFiles> outfile");
                System.out.println("Takes a list of fasta files (more than one) and finds the longest subsequences that are common to all files.");
                System.out.println("<fasta_files> - two or more (multi)-fasta files.");
                System.out.println("outfile - the last file path provided. Will contain the longest subsequences in the fasta files common to all input files.");
            }
            else
            {
                // create Options object
                Options options = new Options();
                options.addOption(OptionBuilder.withArgName("infile")
                        .hasArg()
                        .isRequired()
                        .withDescription("input fasta file containing multiple fasta sequences")
                        .create('f'));
                options.addOption(OptionBuilder.withArgName("outfile")
                        .hasArg()
                        .isRequired()
                        .withDescription("Will contain the longest subsequences in the fasta files common to all input files.")
                        .create('o'));
                options.addOption(OptionBuilder.withLongOpt("help").create('h'));

                String header = "Takes a list of fasta files (more than one) and finds the longest subsequences that are common to all files.\n";

                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("FastaGetLongestSubstring", header, options, footer, false);
                CommandLineParser parser = new BasicParser();
                CommandLine cmd = parser.parse(options, args);
                //CAN BE MULTIPLE FILES - HANDLE THIS
                File in = new File(cmd.getOptionValue("f"));
                File out = new File(cmd.getOptionValue("o"));
                
                ArrayList<String> al = new ArrayList<>();
                for (int i = 1; i < args.length - 1; i++)
                {
                    al.add(args[i]);
                }
                File outfile = new File(args[args.length - 1]);
                System.out.println("Outfile " + outfile.toString());
                FastaSubstrings fs = new FastaSubstrings();
                fs.findLongestCommonSequences(al, outfile);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastaGetGenomeLength"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaGetGenomeLength fastaFile");
                System.out.println("Calculates the combined length, in nucleotides, of a file of fasta sequences");
                System.out.println("fastaFile - fasta file to calculated the length");
            }
            else
            {
                File fastaFile = new File(args[1]);
                double genomeLength = FastaFeatures.getGenomeSize(fastaFile);
                System.out.printf("Geneome length =  %.0f\n", genomeLength);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "fastaToFastq"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaToFastq fastaIn fastqOut");
                System.out.println("Changes fasta to fastq files. All quality scores will be '#");
                System.out.println("fastqIn - fasta file to parse");
                System.out.println("fastqOut - the fastq output file");
            }
            else
            {
                File fastaFile = new File(args[1]);
                File fastqFile = new File(args[2]);
                FastaParser fp = new FastaParser();
                fp.fastaToFastq(fastaFile, fastqFile);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastaTranslate"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaTranslate infile outfile");
                System.out.println("Provides a six-frame translation of a fasta file");
                System.out.println("infile - the fasta file to translate");
                System.out.println("outfile - a fasta file containing six traslated protein sequences (1 for each reading frame) for every DNA sequence in the infile.");
            }
            else
            {
                File infile = new File(args[1]);
                File outfile = new File(args[2]);
                FastaTranslator ft = new FastaTranslator();
                ft.translateMultiFasta(infile, outfile);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastaGetGCContent"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaGetGCContent fastaFile ");
                System.out.println("Calculates the GC content of a fasta file");
                System.out.println("fastaFile - the fasta file from which GC content will be calculated");
            }
            else
            {
                File fastaFile = new File(args[1]);
                FastaFeatures ff = new FastaFeatures();
                ff.getGCContent(fastaFile);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastaSelectRandomSequences"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaSelectRandomSequences infile noRandSeqs outfile");
                System.out.println("Selects a given number of random sequences from a fasta file. Sequences will never be included more than once. This method will fail if the number of random sequences requested is larger than that in the original file.");
                System.out.println("infile - the multi-asta file from which to select random sequences");
                System.out.println("noRandSeqs - the number of random sequences to select");
                System.out.println("outfile - the file containing the random sequences");
            }
            else
            {
                File infile = new File(args[1]);
                int noRandSeqs = Integer.parseInt(args[2]);
                File outfile = new File(args[3]);
                RandomFasta rf = new RandomFasta();
                rf.selectRandomSequences(infile, noRandSeqs, outfile);
            }
        }
        else if (args[0].equalsIgnoreCase("FastaFilterByLength"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaFilterByLength infile outfile minLength");
                System.out.println("Filters a fasta file for sequences over a given length.");
                System.out.println("infile - the multi-asta file from which to select  sequences");
                System.out.println("outfile - the file containing the chosen sequences");
                System.out.println("minLength - the minimum length a sequence must be to remain");
            }
            else
            {
                File infile = new File(args[1]);
                File outfile = new File(args[2]);
                int minLength = Integer.parseInt(args[3]);

                FastaParser fp = new FastaParser();
                fp.filterFastaByLength(infile, outfile, minLength);
            }
        }
        else if (args[0].equalsIgnoreCase("FastaGetSingleFromMultiFasta"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaGetSingleFromMultiFasta multiFastaFile outfile seqid");
                System.out.println("Usage: FastaGetSingleFromMultiFasta multiFastaFile outfile seqid start end");
                System.out.println("Extracts a single fasta sequence from a multi-fasta file.");
                System.out.println("Optionally (with addition of parameters 'start' and 'end') extracts a sub-sequence using the start and end co-ordinates");
                System.out.println("multiFastaFile - a multi-fasta file from which a single fasta sequence is required");
                System.out.println("outfile - the file path which to write the single fasta sequence to");
                System.out.println("seqid - the accession number of sequence id of the fasta sequence to extract");
                System.out.println("(optional) start - the 1-based start co-ordinate for a sub-sequence");
                System.out.println("(optional) end - the 1-based end co-ordinate for a sub-sequence");
            }
            else
            {
                File multiFastaFile = new File(args[1]);
                File outfile = new File(args[2]);
                String seqid = args[3];
                int start;
                int end;
                FastaSubstrings ff = new FastaSubstrings();
                if (args.length == 4)
                {
                    ff.getSequence(multiFastaFile, outfile, seqid);
                }
                else if (args.length == 6)
                {
                    start = Integer.parseInt(args[4]);
                    end = Integer.parseInt(args[5]);
                    ff.getSubSequence(multiFastaFile, outfile, seqid, start, end);
                }
                else
                {
                    System.err.println("Wrong number of paramters try GenomeHelper.jar FastaGetSingleFromMultiFasta -h for help");
                }
            }
        }

//        else if (args[0].equalsIgnoreCase("GetNStats"))
//        {
//            if (args[1].equalsIgnoreCase("-h"))
//            {
//                System.out.println("Usage: GetNStats fastaIn");
//                System.out.println("Calculates N10 - N90 for a fasta assembly");
//                System.out.println("fastaInn - the multi-asta file from which to select  sequences");
//            }
//            else
//            {
//                File infile = new File(args[1]);
//                FastaFeatures ff = new FastaFeatures();
//                ArrayList<Integer> seqlengths = ff.getSequenceAsSortedIntArrayList(infile);
//                ff.getNStats(seqlengths);
//            }
//        }
        else if (args[0].equalsIgnoreCase("GetNStats"))
        {

            // create Options object
            Options options = new Options();
            options.addOption(OptionBuilder.withArgName("file")
                    .hasArg()
                    .isRequired()
                    .withDescription("the fasta file to search")
                    .create('f'));

            options.addOption(OptionBuilder.withArgName("min")
                    .hasArg()
                    .withDescription("minimum contig length to include")
                    .create('m'));
            options.addOption(OptionBuilder.withLongOpt("help").create('h'));

            String header = "Calculates N10 - N90 and longest contig for a fasta assembly.\n";
            String footer = "\nPlease report issues to me";
            HelpFormatter formatter = new HelpFormatter();

            formatter.printHelp("GetNStats", header, options, footer);
            CommandLineParser parser = new BasicParser();
            CommandLine cmd = parser.parse(options, args);

            FastaFeatures ff = new FastaFeatures();
            ArrayList<Integer> seqlengths;
            if (cmd.hasOption("m"))
            {
                File in = new File(cmd.getOptionValue("f"));
                int length = Integer.parseInt(cmd.getOptionValue("m"));
                seqlengths = ff.getSequenceAsSortedIntArrayList(in, length);
                System.out.println(in + " " + length);
            }
            else
            {
                File in = new File(cmd.getOptionValue("f"));
                seqlengths = ff.getSequenceAsSortedIntArrayList(in);
                System.out.println(in);
            }

            ff.getNStats(seqlengths);

        }

        else if (args[0].equalsIgnoreCase(
                "GetReadsFromList"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GetReadsFromList listFile, leftFastq, rightFastq, readsOut");

                System.out.println("Takes a tab-delimited file in the format <readname><pair-end>, where pair-end is either '1',meaning a left-hand (forward) read or '2', meaning a right-hand (reverse) read.");
                System.out.println("listfile - a tab-delimited file containing read-names and pair-end annotation");
                System.out.println("leftFastq - the left-hand reads");
                System.out.println("rightFastq the right-hand reads");
                System.out.println("readsOut  the file to write the reads to");
            }
            else
            {
                File listfile = new File(args[1]);
                File leftFastq = new File(args[2]);
                File rightFastq = new File(args[3]);
                File out = new File(args[4]);
                MappedSamRecords msr = new MappedSamRecords();
                msr.getReadsFromList(listfile, leftFastq, rightFastq, out);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqTranslate"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqTranslate in out includeDNA");
                System.out.println("Provides a six-frame translation of a fastq file. Optionally includes the original DNA sequence in the outfile. Output is fasta-format.");
                System.out.println("infile - the fasta file to translate");
                System.out.println("outfile - a fasta file containing six traslated protein sequences (1 for each reading frame) for every DNA sequence in the infile.");
                System.out.println("includeDNA (optional) - a boolean ('true' or 'false' (default)). If set to 'true' the original DNA sequence will be included in outfile. If set to 'false' (default) no DNA sequence will be included.");
            }
            else
            {
                boolean includeDNA = false;
                File in = new File(args[1]);
                File out = new File(args[2]);
                if (args.length == 4)
                {
                    includeDNA = Boolean.parseBoolean(args[3]);
                }

                FastqParser fp = new FastqParser();
                fp.fastqToFastaSixFrameTranslation(in, out, includeDNA);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqMotifFinder"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqMotifFinder fastaFile searchMotif motifCountsFile proteinCountsFile minCount");
                System.out.println("fastaFile - the path to a (multi)fasta file in which to search for the motif");
                System.out.println("searchMotif - the motif to search for. The motif can contain regular expressions, such as ATG[CG]G[AT],"
                                   + " which will search for ATGCGA, ATGCGT, ATGGGA and ATGGGT. Refer to java.util.regex.Pattern for more details.");
                System.out.println("motifCounts - a tab-delimited file of the occurrence of each pattern");
                System.out.println("proteinCounts - a tab-delimited file of the occurrence of each pattern when translated into ammino acids");
                System.out.println("minCount - the minimum number of times a motif must be found to be included in the results\n");
            }
            else
            {
                File fastqFile = new File(args[1]);
                String searchMotif = args[2];
                File motifCounts = new File(args[3]);
                File aaMotifCounts = new File(args[4]);
                int minCount = Integer.parseInt(args[5]);
                FastqMotifFinder fmf = new FastqMotifFinder();
                fmf.findMatches(fastqFile, searchMotif, motifCounts, aaMotifCounts, minCount);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqToFasta"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqToFasta fastq fasta");
                System.out.println("Transforms a fastq dataset into a fasta-formatted dataset");
                System.out.println("fastq - the fastq which to transform into fasta");
                System.out.println("fasta - the resulting fasta file");
            }
            else
            {
                File fastq = new File(args[1]);
                File fasta = new File(args[2]);
                FastqParser fp = new FastqParser();
                fp.fastqToFastaFile(fastq, fasta);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqGetSingleEndSequencesFromFile"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqGetSingleEndSequencesFromFile readNames in out");
                System.out.println("Extracts a subset of single-end sequences from a list of sequence names conatined in a given file (each sequence name on a seperate line).");
                System.out.println("readNames - a file containing a list of sequence names to extract from a single-end dataset (each sequence name on a seperate line)");
                System.out.println("in - the reads from which to extract the subsetted sequences");
                System.out.println("out - the subsetted reads");
            }
            else
            {
                File listFile = new File(args[1]);
                File in = new File(args[2]);
                File out = new File(args[3]);
                FastqParser fp = new FastqParser();
                fp.getOneSideFastqSeqsFromList(listFile, in, out);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqGetPairedEndSequencesFromFile"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqGetPairedEndSequencesFromFile readNames leftIn rightIn leftOut rightOut");
                System.out.println("Extracts a subset of paired-end sequences from a list of sequence names conatined in a given file (each sequence name on a seperate line).");
                System.out.println("readNames - a file containing a list of sequence names to extract from a paired-end dataset (each sequence name on a seperate line)");
                System.out.println("leftIn - the left-handed reads from which to extract the subseted sequences");
                System.out.println("rightIn - the right-handed reads from which to extract the subseted sequences");
                System.out.println("leftOut - the sub-setted left-handed reads");
                System.out.println("rightOut - the sub-setted right-handed reads");
            }
            else
            {
                File listFile = new File(args[1]);
                File leftIn = new File(args[2]);
                File rightIn = new File(args[3]);
                File leftOut = new File(args[4]);
                File rightOut = new File(args[5]);
                FastqParser fp = new FastqParser();
                HashSet readNames = fp.readNamesToHashSet(listFile);
                fp.getPairedFastqSeqsFromHashSet(readNames, leftIn, rightIn, leftOut, rightOut);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqGetPairedEndSequencesWithMotifMatch"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqGetPairedEndSequencesWithMotifMatch fastqIn searchPattern matchingFastq");
                System.out.println("Finds paired-end sequences where either pair contains a given sequence motif. Reverse complements of each read is also examined.");
                System.out.println("fastqIn - an interlaced file of paired-end reads");
                System.out.println("searchPattern - the sequenced motif to search for");
                System.out.println("matchingFastq - the paired-end fastq files in which at least one of the pair contains the searchPattern");
            }
            else
            {
                File fastqIn = new File(args[1]);
                String searchPattern = args[2];
                File matchingFastq = new File(args[3]);
                FastqMotifFinder fmf = new FastqMotifFinder();
                fmf.getPEFastqReadsFromMotif(fastqIn, searchPattern, matchingFastq);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqCompress"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqCompress fastq compFastqOut");
                System.out.println("Compresses an input fastq file");
                System.out.println("fastq - the fastq file to compress");
                System.out.println("compFastqOut - the compressed fastq file");
            }
            else
            {
                File fastq = new File(args[1]);
                File compFastqOut = new File(args[2]);
                FastqCompression comp = new FastqCompression();
                comp.compressFastq(fastq, compFastqOut);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqInterlace"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqInterlace leftReads rightReads interlacedFastqFile singlesFile");
                System.out.println("Interlaces two fastq files of left-handed and right-handed reads. Read pairs need not be in their corresponding order. Singletons are permitted.");
                System.out.println("interlacedFastqFile - the input interlaced fastq file");
                System.out.println("leftReads - the left-handed pair of a paired-end fastq file");
                System.out.println("rightReads - the right-handed pair of a paired-end fastq file");
                System.out.println("interlacedFastqFile - the output interlaced pairs (in fastq format)");
                System.out.println("singlesFile - singleton fastq files without a pair");
            }
            else
            {
                File leftReads = new File(args[1]);
                File rightReads = new File(args[2]);
                File interlacedFastqFile = new File(args[3]);
                File singlesFile = new File(args[4]);
                FastqInterlacer ft = new FastqInterlacer();
                ft.interlace(leftReads, rightReads, interlacedFastqFile, singlesFile);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqInterlaceKnownPairs"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqInterlaceKnownPairs leftReads rightReads interlacedFastqFile");
                System.out.println("Interlaces two fastq files of left-handed and right-handed reads. Each pair of reads must be the corresponding position of its mate in both files. No singletons are permitted.");
                System.out.println("leftReads - the left-handed pair of a paired-end fastq file");
                System.out.println("rightReads - the right-handed pair of a paired-end fastq file");
                System.out.println("interlacedFastqFile - the output interlaced fastq file");
            }
            else
            {
                File leftReads = new File(args[1]);
                File rightReads = new File(args[2]);
                File interlacedFastqFile = new File(args[3]);
                FastqInterlacer ft = new FastqInterlacer();
                ft.interlaceKnownPairs(leftReads, rightReads, interlacedFastqFile);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqDeinterlace"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqDeinterlace interlacedFastqFile leftReads rightReads leftSinglesFile rightSinglesFile");
                System.out.println("Deinterlaces a file of interlaced fastq files. The file need not be in order and may contain unpaired reads");
                System.out.println("interlacedFastqFile - the input interlaced fastq file");
                System.out.println("leftReads - the left-handed pair of a paired-end fastq file");
                System.out.println("rightReads - the right-handed pair of a paired-end fastq file");
                System.out.println("fastqJoinedFile - the output joined pairs (in fastq format)");
                System.out.println("fastqSinglesFile - singleton fastq files without a pair");
            }
            else
            {
                File interlacedFastqFile = new File(args[1]);
                File leftReads = new File(args[2]);
                File rightReads = new File(args[3]);
                File leftSinglesFile = new File(args[4]);
                File rightSinglesFile = new File(args[5]);

                FastqInterlacer ft = new FastqInterlacer();
                ft.deinterlace(interlacedFastqFile, leftReads, rightReads, leftSinglesFile, rightSinglesFile);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqJoin"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqJoin leftReads rightReads fastqJoinedFile fastqSinglesFile");
                System.out.println("Joins paired-end reads into one read by joining the two pairs together.");
                System.out.println("leftReads - the left-handed pair of a paired-end fastq file");
                System.out.println("rightReads - the right-handed pair of a paired-end fastq file");
                System.out.println("fastqJoinedFile - the output joined pairs (in fastq format)");
                System.out.println("fastqSinglesFile - singleton fastq files without a pair");
            }
            else
            {
                File leftReads = new File(args[1]);
                File rightReads = new File(args[2]);
                File fastqJoinedFile = new File(args[3]);
                File fastqSinglesFile = new File(args[4]);
                FastqJoiner fj = new FastqJoiner();
                fj.join(leftReads, rightReads, fastqJoinedFile, fastqSinglesFile);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqSplit"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqSplit joinedFastqFile leftReads rightReads");
                System.out.println("Splits a single fastq file which contains paired-end reads where each read pair is joined");
                System.out.println("leftReads - the split left-handed reads");
                System.out.println("rightReads - the split right-handed reads");
            }
            else
            {
                File joinedFastqFile = new File(args[1]);
                File leftReads = new File(args[2]);
                File rightReads = new File(args[3]);
                FastqJoiner fj = new FastqJoiner();
                fj.split(joinedFastqFile, leftReads, rightReads);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqFindKmer"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqFindKmer fastqIn kmer");
                System.out.println("Finds reads with a given kmer or sub-sequence");
                System.out.println("fastqin - the reads to search");
                System.out.println("kmer - the kmer or subsequence to search for");
                System.out.println("Prints read names and sequence to STDOUT");
            }
            else
            {
                File fastqin = new File(args[1]);
                String kmer = args[2];
                FastqParser fp = new FastqParser();
                fp.findKmerInReads(fastqin, kmer);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "FastqCountNucleotides"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqCountNucleotides fastqIn");
                System.out.println("Counts the number of reads and combined read lengths for a given fastq file");
                System.out.println("fastqIn - the fastq file to count");
            }
            else
            {
                File fastqIn = new File(args[1]);
                FastqQC fq = new FastqQC();
                fq.getNucleotideCount(fastqIn);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "QCVerifyPairedEndReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCVerifyPairedEndReads fastqInLeft fastqInRight");
                System.out.println("Verifies that two files of left/right paired-end reads are in order and contain no fastq format errors .");
                System.out.println("fastqInLeft - the left-handed reads to verify");
                System.out.println("fastqInRight - the right-handed reads to verify");
            }
            else
            {
                File fastqInLeft = new File(args[1]);
                File fastqInRight = new File(args[2]);
                FastqQC check = new FastqQC();
                check.veryfiyPairedReads(fastqInLeft, fastqInRight);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "QCVerifyReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCVerifyReads fastqIn ");
                System.out.println("Checks a fastq file to verify that all the reads can be parsed into a FastqRecord. Errors should be thrown if any of the reads are not formatted correctly.");
                System.out.println("fastqIn - the fastq file to verify");
            }
            else
            {
                File fastqIn = new File(args[1]);

                FastqQC check = new FastqQC();
                check.veryfiyReads(fastqIn);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "QCJoinedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCJoinedReads fastqIn fastqOutLeft fastqOutRight readLength format writeBadReads");
                System.out.println("Takes an fastq file which contains joined reads and returns two fastq files containing QCd paired-end sequences. Removes reads which contain short/long reads or a read that contain an 'N',"
                                   + " from single-end fastq files. Fastq files with a qulaity-score format of 'illumina' (pre-Illumina 1.8)"
                                   + " will have their quality score format changed to sanger format. Optionally writes bad reads to a file.");
                System.out.println("fastqIn - the reads to QC");
                System.out.println("fastqOut - the QCd reads");
                System.out.println("readLength - the single-end read length for the input files");
                System.out.println("format - 'illumina' or 'sanger'. 'illumina' should be used for sequences with quality-scores that are pre-Illumina 1.8."
                                   + " They will automatically be reformated into sanger-quality scores.");
                System.out.println("writeBadReads (optional) - a boolean ('true' or 'false' (default)). If set to 'true' the paired-reads failing QC"
                                   + " will be written to a 'bad reads' file (bad reads file name will start with 'bad_'");
            }
            else
            {
                File fastqIn = new File(args[1]);
                File fastqOutLeft = new File(args[2]);
                File fastqOutRight = new File(args[3]);
                int readLength = Integer.parseInt(args[4]);
                String format = args[5];
                boolean writeBadReads = Boolean.parseBoolean(args[6]);

                System.out.println("write bad reads = " + writeBadReads);
                FastqQC check = new FastqQC();
                check.qcJoinedReads(fastqIn, fastqOutLeft, fastqOutRight, readLength, format, writeBadReads);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "QCInterlacedReadsToPairs"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCInterlacedReads fastqIn  fastqLeftOut fastqRightOut singleReads readLength format writeBadReads");
                System.out.println("Takes an interlaced fastq file and returns two fastq files containing QCd paired-end sequences. Removes reads which contain short/long reads or a read that contain an 'N',"
                                   + " from single-end fastq files. Fastq files with a qulaity-score format of 'illumina' (pre-Illumina 1.8)"
                                   + " will have their quality score format changed to sanger format. Optionally writes bad reads to a file.");
                System.out.println("fastqIn - the reads to QC");
                System.out.println("fastqOut - the QCd reads");
                System.out.println("singleReads - the QCd reads where the opposite pair failed QC");
                System.out.println("readLength - the single-end read length for the input files");
                System.out.println("format - 'illumina' or 'sanger'. 'illumina' should be used for sequences with quality-scores that are pre-Illumina 1.8."
                                   + " They will automatically be reformated into sanger-quality scores.");
                System.out.println("writeBadReads (optional) - a boolean ('true' or 'false' (default)). If set to 'true' the paired-reads failing QC"
                                   + " will be written to a 'bad reads' file (bad reads file name will start with 'bad_'");
            }
            else
            {
                boolean writeBadReads = false;
                File fastqIn = new File(args[1]);
                File fastqLeftOut = new File(args[2]);
                File fastRightOut = new File(args[3]);
                File singles = new File(args[4]);
                int readLength = Integer.parseInt(args[5]);
                String format = args[6];
                if (args.length == 8)
                {
                    writeBadReads = Boolean.parseBoolean(args[7]);
                }
                System.out.println("write bad reads = " + writeBadReads);
                FastqQC check = new FastqQC();
                check.qcInterlacedReadsToPairs(fastqIn, fastqLeftOut, fastRightOut, singles, readLength, format, writeBadReads);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "QCInterlacedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCInterlacedReads fastqIn fastqOut singleReads readLength format writeBadReads");
                System.out.println("Takes and returns an interlaced fastq file and removes reads which contain short/long reads or a read that contain an 'N',"
                                   + " from single-end fastq files. Fastq files with a qulaity-score format of 'illumina' (pre-Illumina 1.8)"
                                   + " will have their quality score format changed to sanger format. Optionally writes bad reads to a file.");
                System.out.println("fastqIn - the reads to QC");
                System.out.println("fastqOut - the QCd reads");
                System.out.println("singleReads - the QCd reads where the opposite pair failed QC");
                System.out.println("readLength - the single-end read length for the input files");
                System.out.println("format - 'illumina' or 'sanger'. 'illumina' should be used for sequences with quality-scores that are pre-Illumina 1.8."
                                   + " They will automatically be reformated into sanger-quality scores.");
                System.out.println("writeBadReads (optional) - a boolean ('true' or 'false' (default)). If set to 'true' the paired-reads failing QC"
                                   + " will be written to a 'bad reads' file (bad reads file name will start with 'bad_'");
            }
            else
            {
                boolean writeBadReads = false;
                File fastqIn = new File(args[1]);
                File fastqOut = new File(args[2]);
                File singles = new File(args[3]);
                int readLength = Integer.parseInt(args[4]);
                String format = args[5];
                if (args.length == 7)
                {
                    writeBadReads = Boolean.parseBoolean(args[6]);
                }
                System.out.println("write bad reads = " + writeBadReads);
                FastqQC check = new FastqQC();
                check.qcInterlacedReads(fastqIn, fastqOut, singles, readLength, format, writeBadReads);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "RemoveNsFromPairedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: RemoveNsFromPairedReads fastqIn fastqOutLeft fastqOutRight ");
                System.out.println("Takes fastq paired reads and returns  paired-end reads that contain no 'N',"
                                   + " from single-end fastq files. Expects quality scores in fastqsanger");
                System.out.println("fastqInLeft - the left-handed reads to QC");
                System.out.println("fastqInRight - the right-handed reads to QC");
                System.out.println("fastqOutLeft - the QCd left-handed reads");
                System.out.println("fastqOutRight - the QCd right-handed reads");

            }
            else
            {

                File fastqInLeft = new File(args[1]);
                File fastqInRight = new File(args[2]);
                File fastqOutLeft = new File(args[3]);
                File fastqOutRight = new File(args[4]);

                FastqQC check = new FastqQC();
                check.removeNsFromPairedReads(fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "QCSingleEndReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCSingleEndReads fastqIn  fastqOut  readLength format writeBadReads");
                System.out.println("Removes reads which contain short/long reads or a read that contain an 'N',"
                                   + " from single-end fastq files. Fastq files with a qulaity-score format of 'illumina' (pre-Illumina 1.8)"
                                   + " will have their quality score format changed to sanger format. Optionally writes bad reads to a file.");
                System.out.println("fastqIn - the reads to QC");
                System.out.println("fastqOut - the QCd reads");
                System.out.println("readLength - the single-end read length for the input files");
                System.out.println("format - 'illumina' or 'sanger'. 'illumina' should be used for sequences with quality-scores that are pre-Illumina 1.8."
                                   + " They will automatically be reformated into sanger-quality scores.");
                System.out.println("writeBadReads (optional) - a boolean ('true' or 'false' (default)). If set to 'true' the paired-reads failing QC"
                                   + " will be written to a 'bad reads' file (bad reads file name will start with 'bad_'");
            }
            else
            {
                boolean writeBadReads = false;
                File fastqIn = new File(args[1]);
                File fastqOut = new File(args[2]);
                int readLength = Integer.parseInt(args[4]);
                String format = args[5];
                if (args.length == 7)
                {
                    writeBadReads = Boolean.parseBoolean(args[6]);
                }
                System.out.println("write bad reads = " + writeBadReads);
                FastqQC check = new FastqQC();
                check.qcSingleEndReads(fastqIn, fastqOut, readLength, format, writeBadReads);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "QCPairedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCPairedReads fastqInLeft fastqInRight fastqOutLeft fastqOutRight singleReads readLength format writeBadReads");
                System.out.println("Removes pairs of reads where at least one of the pair contains short/long reads or a read that contain an 'N',"
                                   + " from paired-end fastq files. Fastq files with a qulaity-score format of 'illumina' (pre-Illumina 1.8)"
                                   + " will have their quality score format changed to sanger format. Optionally writes bad reads to a file.");
                System.out.println("fastqInLeft - the left-handed reads to QC");
                System.out.println("fastqInRight - the right-handed reads to QC");
                System.out.println("fastqOutLeft - the QCd left-handed reads");
                System.out.println("fastqOutRight - the QCd right-handed reads");
                System.out.println("singleReads - the QCd reads where the opposite pair failed QC");
                System.out.println("readLength - the single-end read length for the input files");
                System.out.println("format - 'illumina' or 'sanger'. 'illumina' should be used for sequences with quality-scores that are pre-Illumina 1.8."
                                   + " They will automatically be reformated into sanger-quality scores.");
                System.out.println("writeBadReads (optional) - a boolean ('true' or 'false' (default)). If set to 'true' the paired-reads failing QC"
                                   + " will be written to a 'bad reads' file (bad reads file name will start with 'bad_'");
            }
            else
            {
                boolean writeBadReads = false;
                File fastqInLeft = new File(args[1]);
                File fastqInRight = new File(args[2]);
                File fastqOutLeft = new File(args[3]);
                File fastqOutRight = new File(args[4]);
                File singles = new File(args[5]);
                int readLength = Integer.parseInt(args[6]);
                String format = args[7];
                if (args.length == 9)
                {
                    writeBadReads = Boolean.parseBoolean(args[8].toLowerCase());
                }

                System.out.println("write bad reads = " + writeBadReads);
                FastqQC check = new FastqQC();
                check.qcPairedReads(fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight, singles, readLength, format, writeBadReads);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "QCRemoveKmerPairedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCRemoveKmerPairedReads fastqInLeft fastqInRight fastqOutLeft fastqOutRight kmerFile");
                System.out.println("Removes pairs of reads where any number of provided kmers is found in either read");
                System.out.println("fastqInLeft - the left-handed reads to filter");
                System.out.println("fastqInRight - the right-handed reads to filter");
                System.out.println("fastqOutLeft - the filtered left-handed reads");
                System.out.println("fastqOutRight - the filtered right-handed reads");
                System.out.println("kmerFile - a file containg the kmers to filter. One kmer per line");

            }
            else
            {
                File fastqInLeft = new File(args[1]);
                File fastqInRight = new File(args[2]);
                File fastqOutLeft = new File(args[3]);
                File fastqOutRight = new File(args[4]);
                File kmerFile = new File(args[5]);
                ArrayList<String> kmers = new ArrayList<>();
                FastqQC check = new FastqQC();
                try (BufferedReader br = new BufferedReader(new FileReader(kmerFile)))
                {
                    String line;
                    while ((line = br.readLine()) != null)
                    {
                        line = line.trim(); // remove leading and trailing whitespace

                        if (!line.equals("")) // don't write out blank lines
                        {
                            kmers.add(line);
                        }
                    }
                }
                check.removePairedReadsWithKmers(fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight, kmers);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "QCRemoveKmerSingleReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCRemoveKmerSingleReads fastqIn  fastqOut  kmerFile");
                System.out.println("Removes single reads where any number of provided kmers is found in either read");
                System.out.println("fastqIn - the  reads to filter");
                System.out.println("fastqOut - the filtered  reads");
                System.out.println("kmerFile - a file containg the kmers to filter. One kmer per line");

            }
            else
            {
                File fastqIn = new File(args[1]);
                File fastqOUt = new File(args[2]);
                File kmerFile = new File(args[3]);
                ArrayList<String> kmers = new ArrayList<>();
                FastqQC check = new FastqQC();
                BufferedReader br = new BufferedReader(new FileReader(kmerFile));
                String line;
                while ((line = br.readLine()) != null)
                {
                    line = line.trim(); // remove leading and trailing whitespace

                    if (!line.equals("")) // don't write out blank lines
                    {
                        kmers.add(line);
                    }
                }

                br.close();
                check.removeSingleReadsWithKmers(fastqIn, fastqOUt, kmers);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "BAMGetPairedUnmappedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetPairedUnmappedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the pairs of sequences where one or both of the pairs are unmapped.");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqInLeft - the left-handed reads that were used in the mapping");
                System.out.println("fastqInRight - the right-handed reads that were used in the mapping");
                System.out.println("fastqOutLeft - The mapped left-handed paired reads");
                System.out.println("fastqOutRight - The mapped right-handed paired reads");
            }
            else
            {
                File bamfile = new File(args[1]);
                File fastqInLeft = new File(args[2]);
                File fastqInRight = new File(args[3]);
                File fastqOutLeft = new File(args[4]);
                File fastqOutRight = new File(args[5]);

                MappedSamRecords msr = new MappedSamRecords();
                HashSet hs = msr.listEitherPairedReadUnmappedFromBam(bamfile);
                msr.writePairedReadsFromHashSet(hs, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "BAMGetPairedMappedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetPairedMappedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the pairs of sequences where one or both of the pairs are mapped.");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqInLeft - the left-handed reads that were used in the mapping");
                System.out.println("fastqInRight - the right-handed reads that were used in the mapping");
                System.out.println("fastqOutLeft - The mapped left-handed paired reads");
                System.out.println("fastqOutRight - The mapped right-handed paired reads");
            }
            else
            {
                File bamfile = new File(args[1]);
                File fastqInLeft = new File(args[2]);
                File fastqInRight = new File(args[3]);
                File fastqOutLeft = new File(args[4]);
                File fastqOutRight = new File(args[5]);

                MappedSamRecords msr = new MappedSamRecords();
                HashSet hs = msr.listEitherPairedReadMappedFromBam(bamfile);
                msr.writePairedReadsFromHashSet(hs, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "BAMGetBothPairedUnmappedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetBothPairedUnmappedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the pairs of sequences where both of the pairs are unmapped.");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqInLeft - the left-handed reads that were used in the mapping");
                System.out.println("fastqInRight - the right-handed reads that were used in the mapping");
                System.out.println("fastqOutLeft - The unmapped left-handed paired reads");
                System.out.println("fastqOutRight - The unmapped right-handed paired reads");
            }
            else
            {
                File bamfile = new File(args[1]);
                File fastqInLeft = new File(args[2]);
                File fastqInRight = new File(args[3]);
                File fastqOutLeft = new File(args[4]);
                File fastqOutRight = new File(args[5]);

                MappedSamRecords msr = new MappedSamRecords();
                HashSet hs = msr.listBothPairedReadsUnmappedFromBam(bamfile);
                msr.writePairedReadsFromHashSet(hs, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "BAMListBothPairedUnmappedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetBothPairedUnmappedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns a text file of the pairs of reads where both of the pairs are unmapped.");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqInLeft - the left-handed reads that were used in the mapping");
                System.out.println("fastqInRight - the right-handed reads that were used in the mapping");
                System.out.println("fastqOutLeft - The unmapped left-handed paired reads");
                System.out.println("fastqOutRight - The unmapped right-handed paired reads");
            }
            else
            {
                File bamfile = new File(args[1]);
                File listFile = new File(args[2]);

                MappedSamRecords msr = new MappedSamRecords();
                HashSet hs = msr.listBothPairedReadsUnmappedFromBam(bamfile);
                msr.hashSetToTextFile(hs, listFile);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "BAMGetSingleUnmappedPairedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetSingleUnmappedPairedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the reads that are paired but only one of the pair is unmapped.");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqInLeft - the left-handed reads that were used in the mapping");
                System.out.println("fastqInRight - the right-handed reads that were used in the mapping");
                System.out.println("fastqSingles - The mapped left-handed paired reads");

            }
            else
            {
                File bamfile = new File(args[1]);
                File fastqInLeft = new File(args[2]);
                File fastqInRight = new File(args[3]);
                File fastqSingles = new File(args[4]);

                MappedSamRecords msr = new MappedSamRecords();
                HashMap hm = msr.listSinglePairedReadUnmappedFromBam(bamfile);
                msr.writeSingleReadsFromHashMap(hm, fastqInLeft, fastqInRight, fastqSingles);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "BAMListSingleUnmappedPairedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMListSingleUnmappedPairedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the reads that are paired but only one of the pair is unmapped.");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("outfile - the unmapped read headers ('1' denontes left-hand reads, '2' denotes right-hand reads");

            }
            else
            {
                File bamfile = new File(args[1]);
                File outfile = new File(args[2]);

                MappedSamRecords msr = new MappedSamRecords();
                HashMap hm = msr.listSinglePairedReadUnmappedFromBam(bamfile);
                msr.hashMapToTextFile(hm, outfile);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "BAMPrintReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetBothMappedPairedRead bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the pairs of sequences where both of the pairs have mapped to a genome.");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqInLeft - the left-handed reads that were used in the mapping");
                System.out.println("fastqInRight - the right-handed reads that were used in the mapping");

            }
            else
            {
                File bamfile = new File(args[1]);
                File fastqInLeft = new File(args[2]);
                File fastqInRight = new File(args[3]);

                MappedSamRecords msr = new MappedSamRecords();
                msr.printReadsFromBamAndFastq(bamfile, fastqInLeft, fastqInRight);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "GFFGetMeanIntronLength"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFGetMeanIntronLength gffFile attribute targets");
                System.out.println("Calculates the mean intron length within a genome of a subset of target genes");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("attribute - the name of the attribute that will make the genes unique (e.g. 'name', 'gene_id', etc))");
                System.out.println("refSeq - the reference sequence for the gff file");
            }
            else
            {
                String gffFile = args[1];
                String attribute = args[2];
                File refSeq = new File(args[3]);

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);
                double genomeSize = FastaFeatures.getGenomeSize(refSeq);
                gffs.getMeanIntronLength(fl, attribute, genomeSize);
            }

        }
        else if (args[0].equalsIgnoreCase(
                "GFFCalculateCodingRegion"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFCalculateCodingRegion gffFile refSeq attribute");
                System.out.println("Calculates proportion of the genome which contains coding features");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("refSeq - the reference sequence for the gff file");
                System.out.println("attribute - the name of the attribute that will make the genes unique (e.g. 'name', 'gene_id', etc))");
            }
            else
            {
                String gffFile = args[1];
                File refSeq = new File(args[2]);
                String attribute = args[3];

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);
                HashMap<String, int[]> genomeMap = new HashMap<>(FastaFeatures.getSequenceAsHashMapIntArray(refSeq));
                double genomeSize = gffs.getGenomeSizeFromIntArrayHashMap(genomeMap);
                gffs.calculateCodingRegion(fl, genomeMap, attribute, genomeSize);
            }

        }
        else if (args[0].equalsIgnoreCase(
                "GFFGetMeanTargetIntronLength"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFGetMeanTargetIntronLength gffFile attribute targets");
                System.out.println("Calculates the mean intron length within a genome of a subset of target genes");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("attribute - the name of the attribute that will make the genes unique (e.g. 'name', 'gene_id', etc))");
                System.out.println("targets - a file of gene names (one per line) for which the mean intron length will be calculated");
            }
            else
            {
                String gffFile = args[1];
                String attribute = args[2];
                File targets = new File(args[3]);

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);

                gffs.getMeanSecretedIntronLength(fl, attribute, targets);
            }

        }
        else if (args[0].equalsIgnoreCase(
                "GFFCreateCodingGenome"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFCreateCodingGenome gffFile refSeq featureName cgenome");
                System.out.println("Writes the DNA sequence of the coding portion of the genome (i.e. exons).");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("refSeq - the reference sequence for the gff file");
                System.out.println("featureName - the gff feature from which to calculate the coding part of the genome (e.g. exon, cds, mRNA)");
                System.out.println("cgenome - the coding genome in fasta format");
            }
            else
            {
                String gffFile = args[1];
                File refSeq = new File(args[2]);
                String feature = args[3];
                File cgenome = new File(args[4]);
                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);

                gffs.createCodingGenome(fl, feature, refSeq, cgenome);
            }

        }
        else if (args[0].equalsIgnoreCase(
                "GFFCreateNonCodingGenome"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFCreateNonCodingGenome gffFile refSeq ncgenome");
                System.out.println("Writes the DNA sequence of the non-coding portion of the genome (i.e. intergenic and intron).");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("refSeq - the reference sequence for the gff file");
                System.out.println("ncgenome - the non-coding genome in fasta format");
            }
            else
            {
                String gffFile = args[1];
                File refSeq = new File(args[2]);
                File ncgenome = new File(args[3]);
                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);

                gffs.createNonCodingGenome(fl, refSeq, ncgenome);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "GFFGetMeanFeatureLength"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFGetMeanFeatureLength gffFile featureName");
                System.out.println("Calculates the mean length of any feature (third column) in a gff file.");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("featureName - the feature to be analysed (e.g. mRNA, exon, CDS, etc.)");
            }
            else
            {
                String gffFile = args[1];
                String featureName = args[2];

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);
                gffs.getMeanFeatureLength(fl, featureName);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "GFFGetMeanFeatureLengthOfGeneIDs"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFGetMeanFeatureLengthOfGeneIDs gffFile refSeq geneIds featureName attribute");
                System.out.println("Calculates the mean length of any feature (third column) in a gff file.");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("refSeq - the reference sequence for the gff file");
                System.out.println("geneIds - a file of gene IDs, one per line. Must matuch up with the value given in the attribute parameter."
                                   + " E.g. gene_id=geneX, 'gene_id' would be the attribute and 'geneX' would be the name of a gene in the geneIds file");
                System.out.println("featureName - the feature to be analysed (e.g. mRNA, exon, CDS, etc.)");
                System.out.println("attribute - the name of the attribute that will identify the genes in the geneIds file (e.g. 'name', 'gene_id', 'ID', etc.))");
            }
            else
            {
                String gffFile = args[1];
                File refSeq = new File(args[2]);
                File geneIds = new File(args[3]);
                String featureName = args[4];
                String attribute = args[5];

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);
                HashMap<String, int[]> genomeMap = new HashMap<>(FastaFeatures.getSequenceAsHashMapIntArray(refSeq));
                double result = gffs.getMeanFeatureLength(fl, genomeMap, featureName, geneIds, attribute);

            }
        }
        else if (args[0].equalsIgnoreCase(
                "GFFGetMeanFeatureLengthWithSplicing"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFGetMeanFeatureLengthWithSplicing gffFile refSeq featureName");
                System.out.println("Calculates the mean length of any feature (third column) in a gff file. Additionally, prints a number of other stats.");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("refSeq - the reference sequence for the gff file");
                System.out.println("featureName - the feature to be analysed (e.g. mRNA, exon, CDS, etc.)");
            }
            else
            {
                String gffFile = args[1];
                String featureName = args[2];
                File refSeq = new File(args[3]);

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);
                HashMap<String, int[]> genomeMap = new HashMap<>(FastaFeatures.getSequenceAsHashMapIntArray(refSeq));
                gffs.getMeanFeatureLength(fl, genomeMap, featureName);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "GFFGetStats"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFGetStats gffFile refSeq featureName");
                System.out.println("Calculates the mean length of any CDS, exons and introns in a gff file.");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("refSeq - the reference sequence for the gff file");
                System.out.println("attribute - the name of the attribute that will make the genes unique (e.g. 'name', 'gene_id', etc))");
            }
            else
            {
                String gffFile = args[1];
                File refSeq = new File(args[2]);
                String attribute = args[3];

                GFFFeatureStats gffs = new GFFFeatureStats();
                gffs.getStats(gffFile, refSeq, attribute);
            }
        }
        else if (args[0].equalsIgnoreCase(
                "gatkToSamInterval"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: gatkToSamInterval bamFile");
                System.out.println("Changes gatk interval format into SAM interval format");
                System.out.println("bam - the bam file for the interval data");
                System.out.println("gatkInterval - the gatk interval file");
                System.out.println("Output file name will be the name of the bam file, with '.interval_list' replacing '.bam'");
            }
            else
            {
                File bamFile = new File(args[1]);
                File interval = new File(args[2]);

                Interval i = new Interval();
                i.gatkToSamInterval(bamFile, interval);
            }
        }
        else if (args[0].equalsIgnoreCase("scrambleGenome"))
        {
            // create Options object
            Options options = new Options();
            options.addOption(OptionBuilder.withArgName("infile")
                    .hasArg()
                    .isRequired()
                    .withDescription("input fasta file containing multiple fasta sequences")
                    .create('f'));
            options.addOption(OptionBuilder.withArgName("no_genomes")
                    .hasArg()
                    .isRequired()
                    .withDescription("the number of shuffled genomes you require")
                    .create('n'));
            options.addOption(OptionBuilder.withArgName("prefix")
                    .hasArg()
                    .isRequired()
                    .withDescription("prefix file names with <string>")
                    .create('p'));

            options.addOption(OptionBuilder.withLongOpt("help").create('h'));

            String header = "Provides shuffled versions of a given fasta file\n";

            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("scrambleGenomes", header, options, footer, false);
            CommandLineParser parser = new BasicParser();
            CommandLine cmd = parser.parse(options, args);

            File fasta = new File(cmd.getOptionValue("f"));

            int noRandSeqs = Integer.parseInt(cmd.getOptionValue("n"));
            String prefix = cmd.getOptionValue("p");
            RandomFasta rand = new RandomFasta();
            rand.scrambleGenome(fasta, noRandSeqs, prefix);

        }

        else
        {
            System.err.println("Unknow program, use GenomeHelper.jar -h for help");
            System.exit(0);
        }
    }
}
