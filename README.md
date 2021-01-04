# GenomeHelper - A command Line program and Java API of helper tools for genomics

### Use `java -jar GenomeHelper.jar <program-name> -h` for help with each program:
### The following programs are available. 

   
### Fasta-related programs:   
Usage: FastaMotifFinder fastaFile searchMotif motifCountsFile proteinCountsFile minCount   
Usage: FastaGetLongestSubstring  <path to files>  outfile.   
Usage: FastaGetGenomeLength  infile   
Usage: FastaTranslate infile outfile   
Usage: FastaGetGCContent  fastafile   
Usage: FastaGetSingleFromMultiFasta  infileoutfile  seqId subsequence_start (optional) subsequence_end (optional)   
Usage: FastaSelectRandomSequences fastaIn numberOfRandomSeqs randomSeqsoutfile   
Usage: FastaToFastq fastaIn fastqOut   
   
### Fastq-related programs:   
Usage: FastqCompress fastqIn fastqOut.gz   
Usage: FastqInterlace leftReads rightReads interlacedFastqFile singlesFile   
Usage: FastqDeinterlace interlacedFastqFile leftReads rightReads  leftSinglesFile rightSinglesFile   
Usage: FastqInterlaceKnownPairs leftReads rightReads interlacedFastqFile   
Usage: FastqJoin leftReads rightReads fastqJoinedFile fastqSinglesFile   
Usage: FastqSplit joinedFastqFile leftPairdReads rightPairedReads   
Usage: FastqMotifFinder fastqFile searchMotif motifCountsFile proteinCountsFile minCount   
Usage: FastqGetPairedEndSequencesWithMotifMatch interlacedFastqFile searchMotif matchingFastqFiles   
Usage: FastqGetPairedEndSequencesFromFile listFile fastqFileInLeft fastqFileInRight fastqFileOutLeft fastqFileOutRight   
Usage: FastqGetSingleEndSequencesFromFile listFile fastqFileIn fastqFileOut   
Usage: FastqToFasta fastqIn fastaOut    
Usage: FastqTranslate fastqIn fastaOut includeOriginalDNASequence ('true' or 'false')   
Usage: FastqCountNucleotides fastqIn   
Usage: FastqFindKmer fastqIn kmer   
   
### Quality-control programs:   
Usage: QCPairedReads fastqInLeft fastqInRight fastqLeftOut fastqrRightOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')   
Usage: QCSingleEndReads fastqIn fastqOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')   
Usage: QCInterlacedReads fastqIn fastqOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')   
Usage: QCInterlacedReadsToPairs fastqIn fastqLeftOut fastqrRightOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')   
Usage: QCJoinedReads fastqIn fastqLeftOut fastqrRightOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')   
Usage: QCVerifyReads fastqIn   
Usage: QCVerifyPairedEndReads fastqLeft fastqRight   
Usage: QCRemoveKmerPairedReads fastqInLeft fastqInRight fastqLeftOut fastqrRightOut kmerFile   
Usage: QCRemoveKmerSingleReads fastqIn fastqOut kmerFile   
   
### SAM/BAM-related programs:   
Usage: BAMGetMappedPairedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight   
Usage: BAMGetUnmappedPairedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight   
Usage: BAMGetBothUnmappedPairedReads bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight   
Usage: BAMGetBothMappedPairedRead bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight   
Usage: BAMGetSingleUnmappedPairedReads bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight   
Usage: BAMGetSingleMappedPairedReads bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight   
   
### GFF-related programs:   
Usage: GFFGetMeanFeatureLengthWithSplicing gffFile featureName refSeq   
Usage: GFFGetMeanFeatureLength gffFile featureName   
Usage: GFFGetMeanFeatureLengthOfGeneIDs gffFile featureName fileOfIds   
Usage: GFFCreateNonCodingGenome gffFile refSeq nonCodingGenome   
Usage: GFFCreateCodingGenome gffFile featureName refSeq codingGenome   
Usage: GFFGetMeanIntronLength gffFile featureName refSeq    
Usage: GFFGetMeanTargetIntronLength gffFile featureName targets   
Usage: GFFCalculateCodingRegion gffFile refSeq attribute   
Usage: GFFGetStats gffFile refSeq attribute   
