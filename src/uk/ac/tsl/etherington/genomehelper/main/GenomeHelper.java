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

/**
 *
 * @author ethering
 */
public class GenomeHelper
{

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, CharacterParseException, Exception
    {
        if (args.length == 0 || args[0].equalsIgnoreCase("-help") || args[0].equalsIgnoreCase("-h"))
        {
            System.out.println("Welcome to GenomeFinder\nThe following programs are available. Use <program_name> -h for help with each program:\n");
            System.out.println("\nFasta-related programs:");
            System.out.println("Usage: FastaMotifFinder fastaFile searchMotif motifCountsFile proteinCountsFile minCount");
            System.out.println("Usage: FastaGetLongestSubstring  <path to files>  outfile.");
            System.out.println("Usage: FastaGetGenomeLength  infile");
            System.out.println("Usage: FastaTranslate infile outfile");
            System.out.println("Usage: FastaGetGCContent  fastafile");
            System.out.println("Usage: FastaGetSingleFromMultiFasta  infileoutfile  seqId subsequence_start (optional) subsequence_end (optional)");
            System.out.println("Usage: FastaSelectRandomSequences fastaIn numberOfRandomSeqs randomSeqsoutfile");
            System.out.println("Usage: FastaToFastq fastaIn fastqOut");

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
            System.out.println("Usage: QCPairedReads fastqInLeft fastqInRight fastqLeftOut fastqrRightOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')");
            System.out.println("Usage: QCSingleEndReads fastqIn fastqOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')");
            System.out.println("Usage: QCInterlacedReads fastqIn fastqOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')");
            System.out.println("Usage: QCInterlacedReadsToPairs fastqIn fastqLeftOut fastqrRightOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')");
            System.out.println("Usage: QCJoinedReads fastqIn fastqLeftOut fastqrRightOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')");
            System.out.println("Usage: QCVerifyReads fastqIn");
            System.out.println("Usage: QCVerifyPairedEndReads fastqLeft fastqRight");
            System.out.println("Usage: QCRemoveKmerPairedReads fastqInLeft fastqInRight fastqLeftOut fastqrRightOut kmerFile");
            System.out.println("Usage: QCRemoveKmerSingleReads fastqIn fastqOut kmerFile");

            System.out.println("\nSAM/BAM-related programs:");
            System.out.println("Usage: BAMGetMappedPairedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
            System.out.println("Usage: BAMGetUnmappedPairedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
            System.out.println("Usage: BAMGetBothUnmappedPairedReads bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
            System.out.println("Usage: BAMGetBothMappedPairedRead bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
            System.out.println("Usage: BAMGetSingleUnmappedPairedReads bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
            System.out.println("Usage: BAMGetSingleMappedPairedReads bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");

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

        } else if (args[0].equalsIgnoreCase("FastaMotifFinder"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaMotifFinder fastaFile searchMotif motifCountsFile proteinCountsFile minCount");
                System.out.println("fastaFile - the path to a (multi)fasta file in which to search for the motif");
                System.out.println("searchMotif - the motif to search for. The motif can contain regular expressions, such as ATG[CG]G[AT],"
                        + " which will search for ATGCGA, ATGCGT, ATGGGA and ATGGGT. Refer to java.util.regex.Pattern for more details.");
                System.out.println("motifCounts - a tab-delimited file of the occurrence of each pattern");
                System.out.println("proteinCounts - a tab-delimited file of the occurrence of each pattern when translated into ammino acids");
                System.out.println("minCount - the minimum number of times a motif must be found to be included in the results\n");
            } else
            {
                File fastaFile = new File(args[1]);
                String searchMotif = args[2];
                File motifCounts = new File(args[3]);
                File aaMotifCounts = new File(args[4]);
                int minCount = Integer.parseInt(args[5]);
                FastaMotifFinder fmf = new FastaMotifFinder();
                fmf.findMatches(fastaFile, searchMotif, motifCounts, aaMotifCounts, minCount);
            }
        } else if (args[0].equalsIgnoreCase("FastaGetLongestSubstring"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaGetLongestSubstring <fastaFiles> outfile");
                System.out.println("Takes a list of fasta files (more than one) and finds the longest subsequences that are common to all files.");
                System.out.println("<fasta_files> - two or more (multi)-fasta files.");
                System.out.println("outfile - the last file path provided. Will contain the longest subsequences in the fasta files common to all input files.");
            } else
            {
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
        } else if (args[0].equalsIgnoreCase("FastaGetGenomeLength"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaGetGenomeLength fastaFile");
                System.out.println("Calculates the combined length, in nucleotides, of a file of fasta sequences");
                System.out.println("fastaFile - fasta file to calculated the length");
            } else
            {
                File fastaFile = new File(args[1]);
                double genomeLength = FastaFeatures.getGenomeSize(fastaFile);
                System.out.printf("Geneome length =  %.0f\n", genomeLength);
            }
        } else if (args[0].equalsIgnoreCase("fastaToFastq"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaToFastq fastaIn fastqOut");
                System.out.println("Changes fasta to fastq files. All quality scores will be '#");
                System.out.println("fastqIn - fasta file to parse");
                System.out.println("fastqOut - the fastq output file");
            } else
            {
                File fastaFile = new File(args[1]);
                File fastqFile = new File(args[2]);
                FastaParser fp = new FastaParser();
                fp.fastaToFastq(fastaFile, fastqFile);
            }
        } else if (args[0].equalsIgnoreCase("FastaTranslate"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaTranslate infile outfile");
                System.out.println("Provides a six-frame translation of a fasta file");
                System.out.println("infile - the fasta file to translate");
                System.out.println("outfile - a fasta file containing six traslated protein sequences (1 for each reading frame) for every DNA sequence in the infile.");
            } else
            {
                File infile = new File(args[1]);
                File outfile = new File(args[2]);
                FastaTranslator ft = new FastaTranslator();
                ft.translateMultiFasta(infile, outfile);
            }
        } else if (args[0].equalsIgnoreCase("FastaGetGCContent"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaGetGCContent fastaFile ");
                System.out.println("Calculates the GC content of a fasta file");
                System.out.println("fastaFile - the fasta file from which GC content will be calculated");
            } else
            {
                File fastaFile = new File(args[1]);
                FastaFeatures ff = new FastaFeatures();
                ff.getGCContent(fastaFile);
            }
        } else if (args[0].equalsIgnoreCase("FastaSelectRandomSequences"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastaSelectRandomSequences infile noRandSeqs outfile");
                System.out.println("Selects a given number of random sequences from a fasta file. Sequences will never be included more than once. This method will fail if the number of random sequences requested is larger than that in the original file.");
                System.out.println("infile - the multi-asta file from which to select random sequences");
                System.out.println("noRandSeqs - the number of random sequences to select");
                System.out.println("outfile - the file containing the random sequences");
            } else
            {
                File infile = new File(args[1]);
                int noRandSeqs = Integer.parseInt(args[2]);
                File outfile = new File(args[3]);
                RandomFasta rf = new RandomFasta();
                rf.selectRandomSequences(infile, noRandSeqs, outfile);
            }
        } else if (args[0].equalsIgnoreCase("FastaGetSingleFromMultiFasta"))
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
            } else
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
                } else if (args.length == 6)
                {
                    start = Integer.parseInt(args[4]);
                    end = Integer.parseInt(args[5]);
                    ff.getSubSequence(multiFastaFile, outfile, seqid, start, end);
                } else
                {
                    System.err.println("Wrong number of paramters try GenomeHelper.jar FastaGetSingleFromMultiFasta -h for help");
                }
            }
        } else if (args[0].equalsIgnoreCase("FastqTranslate"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqTranslate in out includeDNA");
                System.out.println("Provides a six-frame translation of a fastq file. Optionally includes the original DNA sequence in the outfile. Output is fasta-format.");
                System.out.println("infile - the fasta file to translate");
                System.out.println("outfile - a fasta file containing six traslated protein sequences (1 for each reading frame) for every DNA sequence in the infile.");
                System.out.println("includeDNA (optional) - a boolean ('true' or 'false' (default)). If set to 'true' the original DNA sequence will be included in outfile. If set to 'false' (default) no DNA sequence will be included.");
            } else
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
        } else if (args[0].equalsIgnoreCase("FastqMotifFinder"))
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
            } else
            {
                File fastqFile = new File(args[1]);
                String searchMotif = args[2];
                File motifCounts = new File(args[3]);
                File aaMotifCounts = new File(args[4]);
                int minCount = Integer.parseInt(args[5]);
                FastqMotifFinder fmf = new FastqMotifFinder();
                fmf.findMatches(fastqFile, searchMotif, motifCounts, aaMotifCounts, minCount);
            }
        } else if (args[0].equalsIgnoreCase("FastqToFasta"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqToFasta fastq fasta");
                System.out.println("Transforms a fastq dataset into a fasta-formatted dataset");
                System.out.println("fastq - the fastq which to transform into fasta");
                System.out.println("fasta - the resulting fasta file");
            } else
            {
                File fastq = new File(args[1]);
                File fasta = new File(args[2]);
                FastqParser fp = new FastqParser();
                fp.fastqToFastaFile(fastq, fasta);
            }
        } else if (args[0].equalsIgnoreCase("FastqGetSingleEndSequencesFromFile"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqGetSingleEndSequencesFromFile readNames in out");
                System.out.println("Extracts a subset of single-end sequences from a list of sequence names conatined in a given file (each sequence name on a seperate line).");
                System.out.println("readNames - a file containing a list of sequence names to extract from a single-end dataset (each sequence name on a seperate line)");
                System.out.println("in - the reads from which to extract the subsetted sequences");
                System.out.println("out - the subsetted reads");
            } else
            {
                File listFile = new File(args[1]);
                File in = new File(args[2]);
                File out = new File(args[3]);
                FastqParser fp = new FastqParser();
                fp.getOneSideFastqSeqsFromList(listFile, in, out);
            }
        } else if (args[0].equalsIgnoreCase("FastqGetPairedEndSequencesFromFile"))
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
            } else
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
        } else if (args[0].equalsIgnoreCase("FastqGetPairedEndSequencesWithMotifMatch"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqGetPairedEndSequencesWithMotifMatch fastqIn searchPattern matchingFastq");
                System.out.println("Finds paired-end sequences where either pair contains a given sequence motif. Reverse complements of each read is also examined.");
                System.out.println("fastqIn - an interlaced file of paired-end reads");
                System.out.println("searchPattern - the sequenced motif to search for");
                System.out.println("matchingFastq - the paired-end fastq files in which at least one of the pair contains the searchPattern");
            } else
            {
                File fastqIn = new File(args[1]);
                String searchPattern = args[2];
                File matchingFastq = new File(args[3]);
                FastqMotifFinder fmf = new FastqMotifFinder();
                fmf.getPEFastqReadsFromMotif(fastqIn, searchPattern, matchingFastq);
            }
        } else if (args[0].equalsIgnoreCase("FastqCompress"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqCompress fastq compFastqOut");
                System.out.println("Compresses an input fastq file");
                System.out.println("fastq - the fastq file to compress");
                System.out.println("compFastqOut - the compressed fastq file");
            } else
            {
                File fastq = new File(args[1]);
                File compFastqOut = new File(args[2]);
                FastqCompression comp = new FastqCompression();
                comp.compressFastq(fastq, compFastqOut);
            }
        } else if (args[0].equalsIgnoreCase("FastqInterlace"))
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
            } else
            {
                File leftReads = new File(args[1]);
                File rightReads = new File(args[2]);
                File interlacedFastqFile = new File(args[3]);
                File singlesFile = new File(args[4]);
                FastqInterlacer ft = new FastqInterlacer();
                ft.interlace(leftReads, rightReads, interlacedFastqFile, singlesFile);
            }
        } else if (args[0].equalsIgnoreCase("FastqInterlaceKnownPairs"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqInterlaceKnownPairs leftReads rightReads interlacedFastqFile");
                System.out.println("Interlaces two fastq files of left-handed and right-handed reads. Each pair of reads must be the corresponding position of its mate in both files. No singletons are permitted.");
                System.out.println("leftReads - the left-handed pair of a paired-end fastq file");
                System.out.println("rightReads - the right-handed pair of a paired-end fastq file");
                System.out.println("interlacedFastqFile - the output interlaced fastq file");
            } else
            {
                File leftReads = new File(args[1]);
                File rightReads = new File(args[2]);
                File interlacedFastqFile = new File(args[3]);
                FastqInterlacer ft = new FastqInterlacer();
                ft.interlaceKnownPairs(leftReads, rightReads, interlacedFastqFile);
            }
        } else if (args[0].equalsIgnoreCase("FastqDeinterlace"))
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
            } else
            {
                File interlacedFastqFile = new File(args[1]);
                File leftReads = new File(args[2]);
                File rightReads = new File(args[3]);
                File leftSinglesFile = new File(args[4]);
                File rightSinglesFile = new File(args[5]);

                FastqInterlacer ft = new FastqInterlacer();
                ft.deinterlace(interlacedFastqFile, leftReads, rightReads, leftSinglesFile, rightSinglesFile);
            }
        } else if (args[0].equalsIgnoreCase("FastqJoin"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqJoin leftReads rightReads fastqJoinedFile fastqSinglesFile");
                System.out.println("Joins paired-end reads into one read by joining the two pairs together.");
                System.out.println("leftReads - the left-handed pair of a paired-end fastq file");
                System.out.println("rightReads - the right-handed pair of a paired-end fastq file");
                System.out.println("fastqJoinedFile - the output joined pairs (in fastq format)");
                System.out.println("fastqSinglesFile - singleton fastq files without a pair");
            } else
            {
                File leftReads = new File(args[1]);
                File rightReads = new File(args[2]);
                File fastqJoinedFile = new File(args[3]);
                File fastqSinglesFile = new File(args[4]);
                FastqJoiner fj = new FastqJoiner();
                fj.join(leftReads, rightReads, fastqJoinedFile, fastqSinglesFile);
            }
        } else if (args[0].equalsIgnoreCase("FastqSplit"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqSplit joinedFastqFile leftReads rightReads");
                System.out.println("Splits a single fastq file which contains paired-end reads where each read pair is joined");
                System.out.println("leftReads - the split left-handed reads");
                System.out.println("rightReads - the split right-handed reads");
            } else
            {
                File joinedFastqFile = new File(args[1]);
                File leftReads = new File(args[2]);
                File rightReads = new File(args[3]);
                FastqJoiner fj = new FastqJoiner();
                fj.split(joinedFastqFile, leftReads, rightReads);
            }
        } else if (args[0].equalsIgnoreCase("FastqFindKmer"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqFindKmer fastqIn kmer");
                System.out.println("Finds reads with a given kmer or sub-sequence");
                System.out.println("fastqin - the reads to search");
                System.out.println("kmer - the kmer or subsequence to search for");
                System.out.println("Prints read names and sequence to STDOUT");
            } else
            {
                File fastqin = new File(args[1]);
                String kmer = args[2];
                FastqParser fp = new FastqParser();
                fp.findKmerInReads(fastqin, kmer);
            }
        } else if (args[0].equalsIgnoreCase("FastqCountNucleotides"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: FastqCountNucleotides fastqIn");
                System.out.println("Counts the number of reads and combined read lengths for a given fastq file");
                System.out.println("fastqIn - the fastq file to count");
            } else
            {
                File fastqIn = new File(args[1]);
                FastqQC fq = new FastqQC();
                fq.getNucleotideCount(fastqIn);
            }
        } else if (args[0].equalsIgnoreCase("QCVerifyPairedEndReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCVerifyPairedEndReads fastqInLeft fastqInRight");
                System.out.println("Verifies that two files of left/right paired-end reads are in order and contain no fastq format errors .");
                System.out.println("fastqInLeft - the left-handed reads to verify");
                System.out.println("fastqInRight - the right-handed reads to verify");
            } else
            {
                File fastqInLeft = new File(args[1]);
                File fastqInRight = new File(args[2]);
                FastqQC check = new FastqQC();
                check.veryfiyPairedReads(fastqInLeft, fastqInRight);
            }
        } else if (args[0].equalsIgnoreCase("QCVerifyReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCVerifyReads fastqIn ");
                System.out.println("Checks a fastq file to verify that all the reads can be parsed into a FastqRecord. Errors should be thrown if any of the reads are not formatted correctly.");
                System.out.println("fastqIn - the fastq file to verify");
            } else
            {
                File fastqIn = new File(args[1]);

                FastqQC check = new FastqQC();
                check.veryfiyReads(fastqIn);
            }
        } else if (args[0].equalsIgnoreCase("QCJoinedReads"))
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
            } else
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
        } else if (args[0].equalsIgnoreCase("QCInterlacedReadsToPairs"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCInterlacedReads fastqIn  fastqLeftOut fastqRightOut  readLength format writeBadReads");
                System.out.println("Takes an interlaced fastq file and returns two fastq files containing QCd paired-end sequences. Removes reads which contain short/long reads or a read that contain an 'N',"
                        + " from single-end fastq files. Fastq files with a qulaity-score format of 'illumina' (pre-Illumina 1.8)"
                        + " will have their quality score format changed to sanger format. Optionally writes bad reads to a file.");
                System.out.println("fastqIn - the reads to QC");
                System.out.println("fastqOut - the QCd reads");
                System.out.println("readLength - the single-end read length for the input files");
                System.out.println("format - 'illumina' or 'sanger'. 'illumina' should be used for sequences with quality-scores that are pre-Illumina 1.8."
                        + " They will automatically be reformated into sanger-quality scores.");
                System.out.println("writeBadReads (optional) - a boolean ('true' or 'false' (default)). If set to 'true' the paired-reads failing QC"
                        + " will be written to a 'bad reads' file (bad reads file name will start with 'bad_'");
            } else
            {
                boolean writeBadReads = false;
                File fastqIn = new File(args[1]);
                File fastqLeftOut = new File(args[2]);
                File fastRightOut = new File(args[3]);
                int readLength = Integer.parseInt(args[4]);
                String format = args[5];
                if (args.length == 7)
                {
                    writeBadReads = Boolean.parseBoolean(args[6]);
                }
                System.out.println("write bad reads = " + writeBadReads);
                FastqQC check = new FastqQC();
                check.qcInterlacedReadsToPairs(fastqIn, fastqLeftOut, fastRightOut, readLength, format, writeBadReads);
            }
        } else if (args[0].equalsIgnoreCase("QCInterlacedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCInterlacedReads fastqIn fastqOut readLength format writeBadReads");
                System.out.println("Takes and returns an interlaced fastq file and removes reads which contain short/long reads or a read that contain an 'N',"
                        + " from single-end fastq files. Fastq files with a qulaity-score format of 'illumina' (pre-Illumina 1.8)"
                        + " will have their quality score format changed to sanger format. Optionally writes bad reads to a file.");
                System.out.println("fastqIn - the reads to QC");
                System.out.println("fastqOut - the QCd reads");
                System.out.println("readLength - the single-end read length for the input files");
                System.out.println("format - 'illumina' or 'sanger'. 'illumina' should be used for sequences with quality-scores that are pre-Illumina 1.8."
                        + " They will automatically be reformated into sanger-quality scores.");
                System.out.println("writeBadReads (optional) - a boolean ('true' or 'false' (default)). If set to 'true' the paired-reads failing QC"
                        + " will be written to a 'bad reads' file (bad reads file name will start with 'bad_'");
            } else
            {
                boolean writeBadReads = false;
                File fastqIn = new File(args[1]);
                File fastqOut = new File(args[2]);
                int readLength = Integer.parseInt(args[3]);
                String format = args[4];
                if (args.length == 6)
                {
                    writeBadReads = Boolean.parseBoolean(args[5]);
                }
                System.out.println("write bad reads = " + writeBadReads);
                FastqQC check = new FastqQC();
                check.qcInterlacedReads(fastqIn, fastqOut, readLength, format, writeBadReads);
            }
        } else if (args[0].equalsIgnoreCase("QCSingleEndReads"))
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
            } else
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
        } else if (args[0].equalsIgnoreCase("QCPairedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCPairedReads fastqInLeft fastqInRight fastqOutLeft fastqOutRight readLength format writeBadReads");
                System.out.println("Removes pairs of reads where at least one of the pair contains short/long reads or a read that contain an 'N',"
                        + " from paired-end fastq files. Fastq files with a qulaity-score format of 'illumina' (pre-Illumina 1.8)"
                        + " will have their quality score format changed to sanger format. Optionally writes bad reads to a file.");
                System.out.println("fastqInLeft - the left-handed reads to QC");
                System.out.println("fastqInRight - the right-handed reads to QC");
                System.out.println("fastqOutLeft - the QCd left-handed reads");
                System.out.println("fastqOutRight - the QCd right-handed reads");
                System.out.println("readLength - the single-end read length for the input files");
                System.out.println("format - 'illumina' or 'sanger'. 'illumina' should be used for sequences with quality-scores that are pre-Illumina 1.8."
                        + " They will automatically be reformated into sanger-quality scores.");
                System.out.println("writeBadReads (optional) - a boolean ('true' or 'false' (default)). If set to 'true' the paired-reads failing QC"
                        + " will be written to a 'bad reads' file (bad reads file name will start with 'bad_'");
            } else
            {
                boolean writeBadReads = false;
                File fastqInLeft = new File(args[1]);
                File fastqInRight = new File(args[2]);
                File fastqOutLeft = new File(args[3]);
                File fastqOutRight = new File(args[4]);
                int readLength = Integer.parseInt(args[5]);
                String format = args[6];
                if (args.length == 8)
                {
                    writeBadReads = Boolean.parseBoolean(args[7].toLowerCase());
                }

                System.out.println("write bad reads = " + writeBadReads);
                FastqQC check = new FastqQC();
                check.qcPairedReads(fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight, readLength, format, writeBadReads);
            }
        } else if (args[0].equalsIgnoreCase("QCRemoveKmerPairedReads"))
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

            } else
            {
                File fastqInLeft = new File(args[1]);
                File fastqInRight = new File(args[2]);
                File fastqOutLeft = new File(args[3]);
                File fastqOutRight = new File(args[4]);
                File kmerFile = new File(args[5]);
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
                check.removePairedReadsWithKmers(fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight, kmers);
            }
        } else if (args[0].equalsIgnoreCase("QCRemoveKmerSingleReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: QCRemoveKmerSingleReads fastqIn  fastqOut  kmerFile");
                System.out.println("Removes single reads where any number of provided kmers is found in either read");
                System.out.println("fastqIn - the  reads to filter");
                System.out.println("fastqOut - the filtered  reads");
                System.out.println("kmerFile - a file containg the kmers to filter. One kmer per line");

            } else
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
        } else if (args[0].equalsIgnoreCase("BAMGetBothMappedPairedRead"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetBothMappedPairedRead bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the pairs of sequences where both of the pairs have mapped to a genome.");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqInLeft - the left-handed reads that were used in the mapping");
                System.out.println("fastqInRight - the right-handed reads that were used in the mapping");
                System.out.println("fastqOutLeft - The mapped left-handed paired reads");
                System.out.println("fastqOutRight - The mapped right-handed paired reads");
            } else
            {
                File bamfile = new File(args[1]);
                File fastqInLeft = new File(args[2]);
                File fastqInRight = new File(args[3]);
                File fastqOutLeft = new File(args[4]);
                File fastqOutRight = new File(args[5]);

                MappedSamRecords msr = new MappedSamRecords();
                msr.getOnePairedMappedSamRecords(bamfile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
            }
        } else if (args[0].equalsIgnoreCase("BAMGetBothUnmappedPairedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetBothUnmappedPairedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the pairs of sequences where both of the pairs has been not mapped to a genome.");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqInLeft - the left-handed reads that were used in the mapping");
                System.out.println("fastqInRight - the right-handed reads that were used in the mapping");
                System.out.println("fastqOutLeft - The unmapped left-handed paired reads");
                System.out.println("fastqOutRight - The unmapped right-handed paired reads");
            } else
            {
                File bamfile = new File(args[1]);
                File fastqInLeft = new File(args[2]);
                File fastqInRight = new File(args[3]);
                File fastqOutLeft = new File(args[4]);
                File fastqOutRight = new File(args[5]);

                MappedSamRecords msr = new MappedSamRecords();
                msr.getOnePairedUnmappedSamRecords(bamfile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
            }
        } else if (args[0].equalsIgnoreCase("BAMGetMappedPairedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetMappedPairedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the pairs of sequences where at least one of the pair has been mapped to a genome.");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqInLeft - the left-handed reads that were used in the mapping");
                System.out.println("fastqInRight - the right-handed reads that were used in the mapping");
                System.out.println("fastqOutLeft - The mapped left-handed paired reads");
                System.out.println("fastqOutRight - The mapped right-handed paired reads");
            } else
            {
                File bamfile = new File(args[1]);
                File fastqInLeft = new File(args[2]);
                File fastqInRight = new File(args[3]);
                File fastqOutLeft = new File(args[4]);
                File fastqOutRight = new File(args[5]);

                MappedSamRecords msr = new MappedSamRecords();
                msr.getOnePairedMappedSamRecords(bamfile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
            }
        } else if (args[0].equalsIgnoreCase("BAMGetUnmappedPairedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetUnmappedPairedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the pairs of sequences where at least one of the pair has not been mapped to a genome.");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqInLeft - the left-handed reads that were used in the mapping");
                System.out.println("fastqInRight - the right-handed reads that were used in the mapping");
                System.out.println("fastqOutLeft - The unmapped left-handed paired reads");
                System.out.println("fastqOutRight - The unmapped right-handed paired reads");
            } else
            {
                File bamfile = new File(args[1]);
                File fastqInLeft = new File(args[2]);
                File fastqInRight = new File(args[3]);
                File fastqOutLeft = new File(args[4]);
                File fastqOutRight = new File(args[5]);

                MappedSamRecords msr = new MappedSamRecords();
                msr.getOnePairedUnmappedSamRecords(bamfile, fastqInLeft, fastqInRight, fastqOutLeft, fastqOutRight);
            }
        } else if (args[0].equalsIgnoreCase("BAMGetSingleUnmappedPairedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetSingleUnmappedPairedReads bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the fastq sequences that have been mapped to a genome. All"
                        + " fastq reads must either be left-handed or right-handed");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqIn - the fastq files that were used in the mapping");
                System.out.println("fastqOut - The mapped reads");
                System.out.println("isRightHandedReads - a boolean ('true' or 'false'). If the the reads are left-handed reads, "
                        + "set to 'true'. If they are right-handed, set to 'false'");
            } else
            {
                File bamfile = new File(args[1]);
                File fastqIn = new File(args[2]);
                File fastqOut = new File(args[3]);
                boolean isRightHandedReads = Boolean.parseBoolean(args[4]);

                MappedSamRecords msr = new MappedSamRecords();
                msr.getSingleMappedSamRecords(bamfile, fastqIn, fastqOut, isRightHandedReads);
            }
        } else if (args[0].equalsIgnoreCase("BAMGetSingleMappedPairedReads"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: BAMGetSingleMappedPairedReads bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight");
                System.out.println("Returns the fastq sequences that have not been mapped to a genome. All"
                        + " fastq reads must either be left-handed or right-handed");
                System.out.println("bamfile - the sam or bam file to examin");
                System.out.println("fastqIn - the fastq files that were used in the mapping");
                System.out.println("fastqOut - The unmapped reads");
                System.out.println("isRightHandedReads - a boolean ('true' or 'false'). If the the reads are left-handed reads, "
                        + "set to 'true'. If they are right-handed, set to 'false'");
            } else
            {
                File bamfile = new File(args[1]);
                File fastqIn = new File(args[2]);
                File fastqOut = new File(args[3]);
                boolean isRightHandedReads = Boolean.parseBoolean(args[4]);

                MappedSamRecords msr = new MappedSamRecords();
                msr.getSingleUnmappedSamRecords(bamfile, fastqIn, fastqOut, isRightHandedReads);
            }
        } else if (args[0].equalsIgnoreCase("GFFGetMeanIntronLength"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFGetMeanIntronLength gffFile attribute targets");
                System.out.println("Calculates the mean intron length within a genome of a subset of target genes");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("attribute - the name of the attribute that will make the genes unique (e.g. 'name', 'gene_id', etc))");
                System.out.println("refSeq - the reference sequence for the gff file");
            } else
            {
                String gffFile = args[1];
                String attribute = args[2];
                File refSeq = new File(args[3]);

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);
                double genomeSize = FastaFeatures.getGenomeSize(refSeq);
                gffs.getMeanIntronLength(fl, attribute, genomeSize);
            }

        } else if (args[0].equalsIgnoreCase("GFFCalculateCodingRegion"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFCalculateCodingRegion gffFile refSeq attribute");
                System.out.println("Calculates proportion of the genome which contains coding features");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("refSeq - the reference sequence for the gff file");
                System.out.println("attribute - the name of the attribute that will make the genes unique (e.g. 'name', 'gene_id', etc))");
            } else
            {
                String gffFile = args[1];
                File refSeq = new File(args[2]);
                String attribute = args[3];

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);
                HashMap<String, int[]> genomeMap = new HashMap<>(FastaFeatures.getSequenceAsIntArray(refSeq));
                double genomeSize = gffs.getGenomeSizeFromIntArrayHashMap(genomeMap);
                gffs.calculateCodingRegion(fl, genomeMap, attribute, genomeSize);
            }

        } else if (args[0].equalsIgnoreCase("GFFGetMeanTargetIntronLength"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFGetMeanTargetIntronLength gffFile attribute targets");
                System.out.println("Calculates the mean intron length within a genome of a subset of target genes");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("attribute - the name of the attribute that will make the genes unique (e.g. 'name', 'gene_id', etc))");
                System.out.println("targets - a file of gene names (one per line) for which the mean intron length will be calculated");
            } else
            {
                String gffFile = args[1];
                String attribute = args[2];
                File targets = new File(args[3]);

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);

                gffs.getMeanSecretedIntronLength(fl, attribute, targets);
            }

        } else if (args[0].equalsIgnoreCase("GFFCreateCodingGenome"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFCreateCodingGenome gffFile refSeq featureName cgenome");
                System.out.println("Writes the DNA sequence of the non-coding portion of the genome (i.e. intergenic and intron).");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("refSeq - the reference sequence for the gff file");
                System.out.println("featureName - the gff feature from which to calculate the coding part of the genome (e.g. exon, cds, mRNA)");
                System.out.println("cgenome - the coding genome in fasta format");
            } else
            {
                String gffFile = args[1];
                File refSeq = new File(args[2]);
                String feature = args[3];
                File cgenome = new File(args[4]);
                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);

                gffs.createCodingGenome(fl, feature, refSeq, cgenome);
            }

        } else if (args[0].equalsIgnoreCase("GFFCreateNonCodingGenome"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFCreateNonCodingGenome gffFile refSeq ncgenome");
                System.out.println("Writes the DNA sequence of the non-coding portion of the genome (i.e. intergenic and intron).");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("refSeq - the reference sequence for the gff file");
                System.out.println("ncgenome - the non-coding genome in fasta format");
            } else
            {
                String gffFile = args[1];
                File refSeq = new File(args[2]);
                File ncgenome = new File(args[3]);
                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);

                gffs.createNonCodingGenome(fl, refSeq, ncgenome);
            }
        } else if (args[0].equalsIgnoreCase("GFFGetMeanFeatureLength"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFGetMeanFeatureLength gffFile featureName");
                System.out.println("Calculates the mean length of any feature (third column) in a gff file.");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("featureName - the feature to be analysed (e.g. mRNA, exon, CDS, etc.)");
            } else
            {
                String gffFile = args[1];
                String featureName = args[2];

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);
                gffs.getMeanFeatureLength(fl, featureName);
            }
        } else if (args[0].equalsIgnoreCase("GFFGetMeanFeatureLengthOfGeneIDs"))
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
            } else
            {
                String gffFile = args[1];
                File refSeq = new File(args[2]);
                File geneIds = new File(args[3]);
                String featureName = args[4];
                String attribute = args[5];

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);
                HashMap<String, int[]> genomeMap = new HashMap<>(FastaFeatures.getSequenceAsIntArray(refSeq));
                double result = gffs.getMeanFeatureLength(fl, genomeMap, featureName, geneIds, attribute);

            }
        } else if (args[0].equalsIgnoreCase("GFFGetMeanFeatureLengthWithSplicing"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFGetMeanFeatureLengthWithSplicing gffFile refSeq featureName");
                System.out.println("Calculates the mean length of any feature (third column) in a gff file. Additionally, prints a number of other stats.");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("refSeq - the reference sequence for the gff file");
                System.out.println("featureName - the feature to be analysed (e.g. mRNA, exon, CDS, etc.)");
            } else
            {
                String gffFile = args[1];
                String featureName = args[2];
                File refSeq = new File(args[3]);

                GFFFeatureStats gffs = new GFFFeatureStats();
                FeatureList fl = gffs.getFeatureList(gffFile);
                HashMap<String, int[]> genomeMap = new HashMap<>(FastaFeatures.getSequenceAsIntArray(refSeq));
                gffs.getMeanFeatureLength(fl, genomeMap, featureName);
            }
        } else if (args[0].equalsIgnoreCase("GFFGetStats"))
        {
            if (args[1].equalsIgnoreCase("-h"))
            {
                System.out.println("Usage: GFFGetStats gffFile refSeq featureName");
                System.out.println("Calculates the mean length of any CDS, exons and introns in a gff file.");
                System.out.println("gffFile - the gff or gtf file in which the features are stored");
                System.out.println("refSeq - the reference sequence for the gff file");
                System.out.println("attribute - the name of the attribute that will make the genes unique (e.g. 'name', 'gene_id', etc))");
            } else
            {
                String gffFile = args[1];
                File refSeq = new File(args[2]);
                String attribute = args[3];

                GFFFeatureStats gffs = new GFFFeatureStats();
                gffs.getStats(gffFile, refSeq, attribute);
            }
        } else
        {
            System.err.println("Unknow program, use GenomeHelper.jar -h for help");
            System.exit(0);
        }
    }
}
