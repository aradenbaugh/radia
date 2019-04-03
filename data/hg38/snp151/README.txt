# This is a switch from using the entire snp151 database to the snp151Common database.
# The snp151Common database is a subset of the snp151 database and only contains snps
# that are mapped to a single location in the reference genome and have a minor allele 
# frequency of at least 1%.

# In the future, these large DB dump files may move to something like Git Large File Storage (LFS).

# To download the data yourself, follow these steps:

# 1) download the UCSC hg38 dbSNP 151 Common table file:
#        On the UCSC genome browser page, click on the 'Common SNPs(151) Track Settings'.
#        At the bottom of the page, there will be information about 'Data Access':
#        http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/

wget --timestamping 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz' -O snp151Common.txt.gz

# 2) Parse out the chrom, chromStart, chromEnd and name fields from the file with commands like these:

zcat snp151Common.txt.gz | awk '$2 == "chr1" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr1.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr2" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr2.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr3" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr3.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr4" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr4.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr5" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr5.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr6" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr6.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr7" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr7.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr8" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr8.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr9" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr9.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr10" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr10.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr11" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr11.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr12" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr12.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr13" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr13.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr14" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr14.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr15" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr15.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr16" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr16.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr17" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr17.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr18" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr18.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr19" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr19.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr20" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr20.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr21" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr21.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chr21" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chr22.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chrX" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chrX.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chrY" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chrY.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chrM" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chrM.bed.gz
zcat snp151Common.txt.gz | awk '$2 == "chrM" {print}' | cut -f 2,3,4,5 | gzip > /path/to/radia/data/hg38/snp151/chrMT.bed.gz
