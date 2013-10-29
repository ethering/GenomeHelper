#GenomeHelper - A command Line program and Java API of helper tools for genomics


###The following programs are available. Use -h for help with each program:

Fasta-related programs:  
Usage: FastaMotifFinder fastaFile searchMotif motifCountsFile proteinCountsFile minCount  
Usage: FastaGetLongestSubstring outfile.  
Usage: FastaGetGenomeLength infile  
Usage: FastaTranslate infile outfile  
Usage: FastaGetGCContent fastafile  
Usage: FastaGetSingleFromMultiFasta infileoutfile seqId subsequence_start (optional) subsequence_end (optional)  
Usage: FastaSelectRandomSequences fastaIn numberOfRandomSeqs randomSeqsoutfile

Fastq-related programs:  
Usage: FastqCompress fastqIn fastqOut.gz  
Usage: FastqInterlace leftReads rightReads interlacedFastqFile singlesFile  
Usage: FastqDeinterlace interlacedFastqFile leftReads rightReads leftSinglesFile rightSinglesFile  
Usage: FastqInterlaceKnownPairs leftReads rightReads interlacedFastqFile  
Usage: FastqJoin leftReads rightReads fastqJoinedFile fastqSinglesFile  
Usage: FastqSplit joinedFastqFile leftPairdReads rightPairedReads  
Usage: FastqMotifFinder fastqFile searchMotif motifCountsFile proteinCountsFile minCount  
Usage: FastqGetPairedEndSequencesWithMotifMatch interlacedFastqFile searchMotif matchingFastqFiles  
Usage: FastqGetPairedEndSequencesFromFile listFile fastqFileInLeft fastqFileInRight fastqFileOutLeft fastqFileOutRight  
Usage: FastqGetSingleEndSequencesFromFile listFile fastqFileIn fastqFileOut  
Usage: FastqToFasta fastqIn fastaOut  
Usage: FastqTranslate fastqIn fastaOut includeOriginalDNASequence ('true' or 'false')

Quality-control programs:  
Usage: QCPairedReads fastqInLeft fastqInRight fastqLeftOut fastqrRightOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')  
Usage: QCSingleEndReads fastqIn fastqOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')  
Usage: QCInterlacedReads fastqIn fastqOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')  
Usage: QCInterlacedReadsToPairs fastqIn fastqLeftOut fastqrRightOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false') 
Usage: QCJoinedReads fastqIn fastqLeftOut fastqrRightOut readLength format('sanger' or 'illumina') writeBadReads ('true' or 'false')  
Usage: QCVerifyReads fastqIn  
Usage: QCVerifyPairedEndReads fastqLeft fastqRight  

SAM/BAM-related programs:  
Usage: BAMGetMappedPairedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight  
Usage: BAMGetUnmappedPairedReads bamfile fastqInLeft fastqInRight fastqOutLeft fastqOutRight  
Usage: BAMGetBothUnmappedPairedReads bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight  
Usage: BAMGetBothMappedPairedRead bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight  
Usage: BAMGetSingleUnmappedPairedReads bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight  
Usage: BAMGetSingleMappedPairedReads bamFile fastqInLeft fastqInRight fastqOutLeft fastqOutRight  

GFF-related programs:  
Usage: GFFGetMeanFeatureLengthWithSplicing gffFile featureName refSeq  
Usage: GFFGetMeanFeatureLength gffFile featureName  
Usage: GFFCreateNonCodingGenome gffFile refSeq nonCodingGeneom  
Usage: GFFCreateCodingGenome gffFile featureName refSeq codingGeneom  
Usage: GFFGetMeanIntronLength gffFile featureName refSeq  
Usage: GFFGetMeanTargetIntronLength gffFile featureName targets  
Usage: GFFCalculateCodingRegion gffFile refSeq attribute