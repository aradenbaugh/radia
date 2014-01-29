=======
RADIA
=======

RADIA:  RNA and DNA Integrated Analysis for Somatic Mutation Detection

RADIA identifies RNA and DNA variants in BAM files.  The program is designed
to take in 4 BAM files:  DNA Normal, RNA Normal, DNA Tumor, and RNA Tumor.  For the 
normal DNA, the program outputs any differences compared to the reference which could 
be potential Germline mutations.  For the normal RNA, the program outputs any differences 
compared to the reference and the normal DNA which could be potential normal RNA-Editing
events.  For the tumor DNA, the program outputs any difference compared to the reference, 
normal DNA and normal RNA which could be potential Somatic mutations.  For the tumor
RNA, the program outputs any difference compared to the reference, normal DNA, normal 
RNA and tumor DNA which could be potential RNA-editing events.

RADIA is typically run on 3 BAM triplets consisting of the Normal DNA, Tumor DNA and 
Tumor RNA.  If no RNA is available from the tumor, then it is run on the normal/tumor
pairs.


=====================
PREREQUISITES
=====================

1) python (version 2.7)
The RADIA code is written and compiled in python.  In order to run the commands, 
you'll need python (version 2.7).

2) samtools (version 0.1.18)
RADIA uses samtools (version 0.1.18 or higher) to examine pileups of reads across
each sample in parallel.  You must install samtools prior to running RADIA.

3) BLAT
RADIA uses BLAT to check the mapping of reads for all Triple BAM calls.

4) snpEff
RADIA uses snpEff to annotate passing variants and to filter out calls from the 
Triple BAM method that land in genes with high sequence similarity.


=====================
DATA PREPARATION
=====================

1) BAM files
The BAM files need to be indexed with the samtools index command and located in
the same directory as the BAM file itself.

2) FASTA files
The fasta files need to be indexed with the samtools faidx command and located in
the same directory as the fasta file itself.

When running RADIA, you need to specify the appropriate fasta file - typically the
one that was used during the alignment which is usually specified in the BAM header.
Some fasta files use the "chr" prefix and some do not.  Some BAM files use the "chr"
prefix and some do not.  If a BAM file uses the "chr" prefix, then the fasta file
that is specified must also use the "chr" prefix and vice versa.

There are multiple ways to specify the fasta files when running RADIA.  You can use
the -f parameter for a fasta file that can be used for multiple BAM files.  Typically, 
the fasta file for the normal and tumor DNA BAMs are the same.  You can also overwrite
the default fasta file specified with the -f parameter with the following BAM specific
fasta files:

--dnaNormalFasta
--rnaNormalFasta
--dnaTumorFasta
--rnaTumorFasta

If the "chr" prefix is neeeded, then add the corresponding flag:
--dnaNormalUseChr
--rnaNormalUseChr
--dnaTumorUseChr
--rnaTumorUseChr


=======================
TEST SAMTOOLS COMMAND
=======================

In order to see if RADIA will be able to execute the samtools command without any errors,
test the command prior to running RADIA.  Here is an example command:

If the "chr" prefix should be used:
samtools mpileup -f fastaFilename.fa -E -r chr7:55248979-55249079 normalBamFilename.bam

If no "chr" prefix is needed:
samtools mpileup -f fastaFilename.fa -E -r 7:55248979-55249079 normalBamFilename.bam


===========================
RUN RADIA INITIAL COMMAND
===========================

RADIA is run on each chromosome separately.  You need to specify the BAM files and the
corresponding FASTA files.  You can specify the output filename where the VCF files will
be output, otherwise it will be sent to STDOUT.

1) Run RADIA on 3 BAM triplets:
python radia.pyc patientId 1 -n normalDnaBamFilename.bam -t tumorDnaBamFilename.bam -r tumorRnaBamFilename.bam -f hg19.fa --rnaTumorUseChr --rnaTumorFasta=hg19_w_chr_prefix.fa -o patientId_chr1.vcf -e hg19 -u url_to_fasta.fa

2) Run RADIA on 2 BAM files:
python radia.pyc patientId 1 -n normalDnaBamFilename.bam -t tumorDnaBamFilename.bam -f hg19.fa -o patientId_chr1.vcf -e hg19 -u url_to_fasta.fa

For the full list of optional parameters, type:
python radia.pyc -h


