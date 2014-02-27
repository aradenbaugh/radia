=======
RADIA
=======

RADIA:  RNA and DNA Integrated Analysis for Somatic Mutation Detection

RADIA identifies RNA and DNA variants in BAM files.  RADIA is typically run on 3 BAM 
files consisting of the Normal DNA, Tumor DNA and Tumor RNA.  If no RNA is available 
from the tumor, then it is run on the normal/tumor pairs.  For the normal DNA, RADIA 
outputs any differences compared to the reference which could be potential Germline 
mutations.  For the tumor DNA, RADIA outputs any differences compared to the reference 
and the normal DNA which could be potential Somatic mutations.  RADIA combines the 
tumor DNA and tumor RNA to augment the somatic mutation calls.  It also uses the 
tumor RNA to identify potential RNA editing events.

The DNA Only Method (DOM) uses just the tumor/normal pairs of DNA (ignoring the RNA), 
while the Triple BAM Method (TBM) uses all three datasets from the same patient to detect 
somatic mutations.  The mutations from the TBM are further categorized into 2 sub-groups:  RNA 
Confirmation and RNA Rescue calls.  RNA Confirmation calls are those that are made by 
both the DOM and the TBM due to the strong read support in both the DNA and RNA.  RNA Rescue 
calls are those that had very little DNA support, hence not called by the DOM, but strong 
RNA support, and thus called by the TBM.  RNA Rescue calls are typically missed by 
traditional methods that only interrogate the DNA.

=====================
PREREQUISITES
=====================

1) python (version 2.7)<br>
The RADIA code is written and compiled in python.  In order to run the commands, 
you'll need python (version 2.7).

2) samtools (version 0.1.18)<br>
RADIA uses samtools (version 0.1.18 or higher) to examine pileups of reads across
each sample in parallel.  You must install samtools prior to running RADIA.

3) pysam API<br>
RADIA uses the pysam API during the filtering process.

4) BLAT<br>
RADIA uses BLAT to check the mapping of reads for all Triple BAM calls.

5) SnpEff<br>
RADIA uses SnpEff to annotate passing variants and to filter out calls from the 
Triple BAM method that land in genes with high sequence similarity.


=====================
DATA PREPARATION
=====================

1) BAM files<br>
The BAM files need to be indexed with the samtools index command and located in
the same directory as the BAM file itself.

2) FASTA files<br>
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

--dnaNormalFasta<br>
--rnaNormalFasta<br>
--dnaTumorFasta<br>
--rnaTumorFasta<br>

If the "chr" prefix is neeeded, then add the corresponding flag:<br>
--dnaNormalUseChr<br>
--rnaNormalUseChr<br>
--dnaTumorUseChr<br>
--rnaTumorUseChr<br>


=======================
TEST SAMTOOLS COMMAND
=======================

In order to see if RADIA will be able to execute the samtools command without any errors,
test the command prior to running RADIA.  Here is an example command:

If the "chr" prefix should be used:<br>
samtools mpileup -f fastaFilename.fa -E -r chr7:55248979-55249079 normalBamFilename.bam

If no "chr" prefix is needed:<br>
samtools mpileup -f fastaFilename.fa -E -r 7:55248979-55249079 normalBamFilename.bam


===========================
RUN RADIA INITIAL COMMAND
===========================

RADIA is run on each chromosome separately.  You need to specify the BAM files and the
corresponding FASTA files.  You can specify the output filename where the VCF files will
be output, otherwise it will be sent to STDOUT.  If the filename ends with ".gz", the VCF
will be gzipped.

1) Run RADIA on 3 BAM files:<br>
python radia.pyc patientId chromId -n normalDnaBamFilename.bam -t tumorDnaBamFilename.bam -r tumorRnaBamFilename.bam -f hg19.fa --rnaTumorUseChr --rnaTumorFasta=hg19_w_chr_prefix.fa -o patientId_chr1.vcf -i hg19 -u http://url_to_fasta.fa

2) Run RADIA on 2 BAM files:<br>
python radia.pyc patientId chromId -n normalDnaBamFilename.bam -t tumorDnaBamFilename.bam -f hg19.fa -o patientId_chr1.vcf -i hg19 -u http://url_to_fasta.fa

For the full list of optional parameters, type:<br>
python radia.pyc -h


===========================
RUN RADIA FILTER COMMAND 
===========================

There is one filter script that can generally be used for all filtering needs.  By default
all of the following filters are applied:<br>

- 1000 Genomes Blacklists<br>
- Flag dbSnp135<br>
- Flag pseudogenes<br>
- Flag cosmic sites<br>
- mpileup support<br>
- read support<br>

For calls that originate in the RNA, there are further filters:
- Filter dbSnp135, dbSnp132, dbSnp130<br>
- Filter pseudogenes<br>
- Filter by BLAT<br>
- Filter by positional bias<br>
- Annotate with SnpEff<br>
- Filter RNA gene and gene family blacklists<br>

If you only have DNA pairs, use the --dnaOnly flag<br>
If you only want the calls from the Triple BAM method, use the --rnaOnly flag<br>

If you want to exclude a particular filter, there are flags such as --noBlat 
to exclude the BLAT filter.

Many of the filters rely on data that is provided in the radia/data/ directory.  Other
dependencies are on the pysam API and external programs such as BLAT and SnpEff.

Here is an example filtering command:<br>
python filterRadia.pyc patientId 22 /inputDir/patientId_chr22.vcf /outputDir/ /radiaDir/scripts/ -b /radiaDir/data/hg19/blacklists/1000Genomes/phase1/ -x /radiaDir/data/hg19/snp135 -y /radiaDir/data/hg19/snp132/ -z /radiaDir/data/hg19/snp130/ -r /radiaDir/data/hg19/retroGenes/ -p /radiaDir/data/hg19/pseudoGenes/ -c /radiaDir/data/hg19/cosmic/ -t /radiaDir/data/hg19/gaf/2_1/ -d /snpEffDir/ --rnaGeneBlckFile ../data/rnaGeneBlacklist.tab --rnaGeneFamilyBlckFile ../data/rnaGeneFamilyBlacklist.tab

Some default parameters to watch out for:<br>
- The default SnpEff genome is set to "GRCh37.69".<br>
- BLAT FASTA filename:  by default the fasta file specified in the BAM header will be used.  You
can overwrite it with the -f parameter.  We recommend that you use a fasta file that includes
the chrUn_gl000 contigs, chr_gl000_random contigs and the hap contigs (e.g. chr6_apd_hap1, 
chr6_cox_hap2, etc).  Often times, the fasta files that are used during the alignment process of
the bams exclude these contigs.<br>
- By default, the calls are filtered by the GAF 2.1 target regions.  If you don't want to filter
by target regions, then use the --noTargets flag.<br>

For the full list of optional parameters, type:<br>
python filterRadia.pyc -h
