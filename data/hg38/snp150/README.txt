# This is a switch from using the entire snp150 database to the snp150Common database.
# The snp150Common database is a subset of the snp150 database and only contains snps
# that are mapped to a single location in the reference genome and have a minor allele 
# frequency of at least 1%.

# In the future, these large DB dump files may move to something like Git Large File Storage (LFS).

# To download the data yourself, follow these steps:

# 1) connect to the UCSC hg38 database and download the dbSNP 150 Common table with commands like these:

"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr1'" > /path/to/radia/data/hg38/snp150/chr1_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr2'" > /path/to/radia/data/hg38/snp150/chr2_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr3'" > /path/to/radia/data/hg38/snp150/chr3_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr4'" > /path/to/radia/data/hg38/snp150/chr4_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr5'" > /path/to/radia/data/hg38/snp150/chr5_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr6'" > /path/to/radia/data/hg38/snp150/chr6_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr7'" > /path/to/radia/data/hg38/snp150/chr7_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr8'" > /path/to/radia/data/hg38/snp150/chr8_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr9'" > /path/to/radia/data/hg38/snp150/chr9_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr10'" > /path/to/radia/data/hg38/snp150/chr10_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr11'" > /path/to/radia/data/hg38/snp150/chr11_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr12'" > /path/to/radia/data/hg38/snp150/chr12_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr13'" > /path/to/radia/data/hg38/snp150/chr13_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr14'" > /path/to/radia/data/hg38/snp150/chr14_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr15'" > /path/to/radia/data/hg38/snp150/chr15_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr16'" > /path/to/radia/data/hg38/snp150/chr16_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr17'" > /path/to/radia/data/hg38/snp150/chr17_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr18'" > /path/to/radia/data/hg38/snp150/chr18_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr19'" > /path/to/radia/data/hg38/snp150/chr19_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr20'" > /path/to/radia/data/hg38/snp150/chr20_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr21'" > /path/to/radia/data/hg38/snp150/chr21_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chr22'" > /path/to/radia/data/hg38/snp150/chr22_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chrX'" > /path/to/radia/data/hg38/snp150/chrX_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chrY'" > /path/to/radia/data/hg38/snp150/chrY_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chrM'" > /path/to/radia/data/hg38/snp150/chrM_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp150Common where chrom = 'chrM'" > /path/to/radia/data/hg38/snp150/chrMT_unsorted.bed

# 2) delete the header line consisting of chrom chromStart chromEnd name from each file

# 3) sort and gzip the files with commands like these:

sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr1_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr1.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr2_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr2.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr3_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr3.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr4_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr4.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr5_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr5.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr6_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr6.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr7_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr7.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr8_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr8.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr9_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr9.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr10_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr10.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr11_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr11.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr12_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr12.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr13_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr13.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr14_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr14.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr15_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr15.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr16_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr16.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr17_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr17.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr18_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr18.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr19_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr19.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr20_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr20.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr21_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr21.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chr22_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chr22.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chrX_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chrX.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chrY_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chrY.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chrM_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chrM.bed.gz
sort -t $'\t' +1n -3 /radia/data/hg38/snp150/chrMT_unsorted.bed | gzip > /path/to/radia/data/hg38/snp150/chrMT.bed.gz

# 4) remove the *_unsorted.bed files
