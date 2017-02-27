# Some of the hg38 dbSNP 147 files (chroms 1-4) were too large (> 100 MB) to upload to git, therefore they have been broken into 2 parts.  
# After you download the data, please cat the 2 parts into 1 file before filtering.  For example, you can use a command like this:

# cat chr1_part1.bed.gz chr1_part2.bed.gz > chr1.bed.gz

# In the future, these large DB dump files may move to something like Git Large File Storage (LFS).


# To download the data yourself:
# connect to the UCSC hg38 database
# download the dbSNP 147 tables with commands like these:
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr1'" > /radia/data/hg38/snp147/chr1_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr2'" > /radia/data/hg38/snp147/chr2_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr3'" > /radia/data/hg38/snp147/chr3_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr4'" > /radia/data/hg38/snp147/chr4_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr5'" > /radia/data/hg38/snp147/chr5_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr6'" > /radia/data/hg38/snp147/chr6_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr7'" > /radia/data/hg38/snp147/chr7_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr8'" > /radia/data/hg38/snp147/chr8_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr9'" > /radia/data/hg38/snp147/chr9_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr10'" > /radia/data/hg38/snp147/chr10_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr11'" > /radia/data/hg38/snp147/chr11_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr12'" > /radia/data/hg38/snp147/chr12_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr13'" > /radia/data/hg38/snp147/chr13_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr14'" > /radia/data/hg38/snp147/chr14_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr15'" > /radia/data/hg38/snp147/chr15_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr16'" > /radia/data/hg38/snp147/chr16_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr17'" > /radia/data/hg38/snp147/chr17_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr18'" > /radia/data/hg38/snp147/chr18_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr19'" > /radia/data/hg38/snp147/chr19_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr20'" > /radia/data/hg38/snp147/chr20_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr21'" > /radia/data/hg38/snp147/chr21_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chr22'" > /radia/data/hg38/snp147/chr22_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chrX'" > /radia/data/hg38/snp147/chrX_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chrY'" > /radia/data/hg38/snp147/chrY_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chrM'" > /radia/data/hg38/snp147/chrM_unsorted.bed
"select chrom,chromStart,chromEnd,name from snp147 where chrom = 'chrMT'" > /radia/data/hg38/snp147/chrMT_unsorted.bed

# be sure to sort them:
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr1_unsorted.bed > /radia/data/hg38/snp147/chr1.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr2_unsorted.bed > /radia/data/hg38/snp147/chr2.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr3_unsorted.bed > /radia/data/hg38/snp147/chr3.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr4_unsorted.bed > /radia/data/hg38/snp147/chr4.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr5_unsorted.bed > /radia/data/hg38/snp147/chr5.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr6_unsorted.bed > /radia/data/hg38/snp147/chr6.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr7_unsorted.bed > /radia/data/hg38/snp147/chr7.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr8_unsorted.bed > /radia/data/hg38/snp147/chr8.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr9_unsorted.bed > /radia/data/hg38/snp147/chr9.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr10_unsorted.bed > /radia/data/hg38/snp147/chr10.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr11_unsorted.bed > /radia/data/hg38/snp147/chr11.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr12_unsorted.bed > /radia/data/hg38/snp147/chr12.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr13_unsorted.bed > /radia/data/hg38/snp147/chr13.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr14_unsorted.bed > /radia/data/hg38/snp147/chr14.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr15_unsorted.bed > /radia/data/hg38/snp147/chr15.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr16_unsorted.bed > /radia/data/hg38/snp147/chr16.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr17_unsorted.bed > /radia/data/hg38/snp147/chr17.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr18_unsorted.bed > /radia/data/hg38/snp147/chr18.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr19_unsorted.bed > /radia/data/hg38/snp147/chr19.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr20_unsorted.bed > /radia/data/hg38/snp147/chr20.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr21_unsorted.bed > /radia/data/hg38/snp147/chr21.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chr22_unsorted.bed > /radia/data/hg38/snp147/chr22.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chrX_unsorted.bed > /radia/data/hg38/snp147/chrX.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chrY_unsorted.bed > /radia/data/hg38/snp147/chrY.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chrM_unsorted.bed > /radia/data/hg38/snp147/chrM.bed
sort -t $'\t' +1n -3 /radia/data/hg38/snp147/chrMT_unsorted.bed > /radia/data/hg38/snp147/chrMT.bed
