RADIA
=======

RADIA:  RNA and DNA Integrated Analysis for Somatic Mutation Detection

RADIA identifies RNA and DNA variants in BAM files.  RADIA is typically run on 3 BAM 
files consisting of the Normal DNA, Tumor DNA and Tumor RNA.  If no RNA is available 
from the tumor, then it is run on the normal/tumor DNA pairs.  For the normal DNA, RADIA 
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


PREREQUISITES
=====================

1) python (version 2.7)<br>
The RADIA code is written and compiled in python.  In order to run the commands, 
you'll need python (version 2.7).

2) samtools (tested on version 0.1.18 and 0.1.19)<br>
RADIA uses samtools to examine pileups of reads across each sample in parallel.  
You must install samtools prior to running RADIA.

3) pysam API (version 0.8.1 and higher)<br>
RADIA uses the pysam API during the filtering process.

4) BLAT (optional)<br>
RADIA uses BLAT to check the mapping of reads for all Triple BAM calls.

5) SnpEff (optional - tested on version 3.3 and 4.3) <br>
RADIA uses SnpEff to annotate passing variants and to filter out calls from the 
Triple BAM method that land in genes with high sequence similarity.


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
the fasta file for the normal and tumor DNA BAMs are the same.  If the fasta file is 
not the same for all BAM files, you can overwrite the default fasta file specified 
with the -f parameter with the following BAM specific fasta files:

--dnaNormalFasta<br>
--rnaNormalFasta<br>
--dnaTumorFasta<br>
--rnaTumorFasta<br>

If the "chr" prefix is neeeded, then add the corresponding flag:<br>
--dnaNormalUseChr<br>
--rnaNormalUseChr<br>
--dnaTumorUseChr<br>
--rnaTumorUseChr<br>


TEST SAMTOOLS COMMAND
=======================

In order to see if RADIA will be able to execute the samtools command without any errors,
test the command prior to running RADIA.  Here is an example command:

If the "chr" prefix should be used:<br>
samtools mpileup -f fastaFilename.fa -E -r chr7:55248979-55249079 normalBamFilename.bam

If no "chr" prefix is needed:<br>
samtools mpileup -f fastaFilename.fa -E -r 7:55248979-55249079 normalBamFilename.bam


RUN RADIA INITIAL COMMAND
===========================

RADIA is run on each chromosome separately.  You need to specify the BAM files and the
corresponding FASTA files.  You can specify the output filename where the VCF files will
be output, otherwise it will be sent to STDOUT.  If the filename ends with ".gz", the VCF
will be gzipped.

1) Run RADIA on 3 BAM files:<br>
python radia.py patientId chromId -n normalDnaBamFilename.bam -t tumorDnaBamFilename.bam -r tumorRnaBamFilename.bam -f hg19.fa --rnaTumorUseChr --rnaTumorFasta=hg19_w_chr_prefix.fa -o /radia/raw/patientId_chr1.vcf.gz -i hg19 -u http://url_to_fasta.fa

2) Run RADIA on 2 BAM files:<br>
python radia.py patientId chromId -n normalDnaBamFilename.bam -t tumorDnaBamFilename.bam -f hg19.fa -o /radia/raw/patientId_chr1.vcf.gz -i hg19 -u http://url_to_fasta.fa

For the full list of optional parameters, type:<br>
python radia.py -h


RUN RADIA FILTER COMMAND 
===========================

There is one filter script that can generally be used for all filtering needs.  Note:  there 
is a difference between "flagging" and "filtering".  A site that is "flagged" will have the 
information added to the INFO column of the VCF, a site that is "filtered" will have the
filter added to the FILTER column of the VCF.  By default all of the following filters are applied:<br>

- 1000 Genomes Blacklists<br>
- Flag dbSNP<br>
- Flag Pseudogenes<br>
- Flag Retrogenes<br>
- Flag COSMIC sites<br>
- Flag RADAR sites<br>
- Flag DARNED sites<br>
- MPileup support<br>
- Read support<br>

For calls that originate in the RNA, there are further filters:
- Filter dbSNP<br>
- Filter Pseudogenes<br>
- Filter Retrogenes<br>
- Filter by BLAT (optional)<br>
- Annotate with SnpEff (optional)<br>
- Filter RNA gene and gene families<br>

If you only have DNA pairs, use the --dnaOnly flag<br>
If you only want the calls from the Triple BAM method, use the --rnaOnly flag<br>

If you want to exclude a particular filter, there are flags such as --noBlat 
to exclude the BLAT filter.

Many of the filters rely on data that is provided in the radia/data/ directory.  Other
dependencies are on the pysam API and external programs such as BLAT and SnpEff.

Here is an example filtering command:<br>
python filterRadia.py patientId 22 /radia/raw/patientId_chr22.vcf /radia/filtered/ /radiaDir/scripts/ -b /radiaDir/data/hg19/blacklists/1000Genomes/phase3/ -d /radiaDir/data/hg19/snp150/ -r /radiaDir/data/hg19/retroGenes/ -p /radiaDir/data/hg19/pseudoGenes/ -c /radiaDir/data/hg19/cosmic/ -t /radiaDir/data/hg19/gencode/basic/ -a /radiaDir/data/hg19/radar/ -n /radiaDir/data/hg19/darned/ -s /snpEffDir/ --rnaGeneBlckFile ../data/rnaGeneBlacklist.tab --rnaGeneFamilyBlckFile ../data/rnaGeneFamilyBlacklist.tab

Some default parameters to watch out for:<br>
- The default SnpEff genome is set to "GRCh37.75".  If you are using a different version, 
be sure to upate the --snpEffGenome parameter.<br>
- BLAT FASTA filename:  by default the fasta file specified in the BAM header will be used.  You
can overwrite it with the -f parameter.  We recommend that you use a fasta file that includes
the chrUn_gl000 contigs, chr_gl000_random contigs and the hap contigs (e.g. chr6_apd_hap1, 
chr6_cox_hap2, etc).  Often times, the fasta files that are used during the alignment process of
the bams exclude these contigs.<br>
- By default, the calls are filtered by the GENCODE basic gene regions.  By specifying the 
--targetsInfo flag, the calls will be flagged (tagged in the INFO column) instead of filtered. 
If you don't want to filter by target regions at all, then use the --noTargets flag.<br>

For the full list of optional parameters, type:<br>
python filterRadia.py -h


RUN RADIA MERGE COMMAND 
===========================

To merge all of the filtered chromosome files into one VCF file for the patient, execute the following command:<br>
python mergeChroms.py patientId /radia/filtered/ /radia/filtered/ --gzip

This will merge all of the files with the names: patientId_chr\*.vcf or patientId_chr\*.vcf.gz into one file called patientId.vcf or patientId.vcf.gz (if you specify the --gzip parameter).


CITATION
===========
If you use RADIA, please cite the method:<br>
Radenbaugh AJ, Ma S, Ewing A, Stuart JM, Collisson EA, Zhu J, Haussler D. (2014) RADIA: RNA and DNA Integrated Analysis for Somatic Mutation Detection. PLoS ONE 9(11): e111516. doi:10.1371/journal.pone.0111516


LICENSE
=========

RNA and DNA Integrated Analysis (RADIA) identifies RNA and DNA variants in NGS data.
Copyright (C) 2010  Amie J. Radenbaugh, Ph.D.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

A copy of the GNU Affero General Public License has been provided 
along with this program.  Otherwise, see <http://www.gnu.org/licenses/>.
