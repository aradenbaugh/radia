# The GENCODE basic gene set is replacing the old GAF (Generic Annotation File) data sets.
# If you want to filter out mutations that are not in the GENCODE basic gene set,
# these are the files that should be specified with the --targetDir option when filtering.
# If you don't want to filter out mutations based on the GENCODE basic gene set, 
# use the --noTargets flag when filtering.

# To download the data yourself, follow these steps:

# 1) connect to the UCSC hg19 database and download the GENCODE basic gene set table with commands like these:

echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr1'" > /path/to/radia/data/hg19/gencode/basic/chr1_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr2'" > /path/to/radia/data/hg19/gencode/basic/chr2_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr3'" > /path/to/radia/data/hg19/gencode/basic/chr3_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr4'" > /path/to/radia/data/hg19/gencode/basic/chr4_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr5'" > /path/to/radia/data/hg19/gencode/basic/chr5_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr6'" > /path/to/radia/data/hg19/gencode/basic/chr6_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr7'" > /path/to/radia/data/hg19/gencode/basic/chr7_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr8'" > /path/to/radia/data/hg19/gencode/basic/chr8_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr9'" > /path/to/radia/data/hg19/gencode/basic/chr9_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr10'" > /path/to/radia/data/hg19/gencode/basic/chr10_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr11'" > /path/to/radia/data/hg19/gencode/basic/chr11_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr12'" > /path/to/radia/data/hg19/gencode/basic/chr12_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr13'" > /path/to/radia/data/hg19/gencode/basic/chr13_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr14'" > /path/to/radia/data/hg19/gencode/basic/chr14_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr15'" > /path/to/radia/data/hg19/gencode/basic/chr15_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr16'" > /path/to/radia/data/hg19/gencode/basic/chr16_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr17'" > /path/to/radia/data/hg19/gencode/basic/chr17_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr18'" > /path/to/radia/data/hg19/gencode/basic/chr18_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr19'" > /path/to/radia/data/hg19/gencode/basic/chr19_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr20'" > /path/to/radia/data/hg19/gencode/basic/chr20_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr21'" > /path/to/radia/data/hg19/gencode/basic/chr21_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chr22'" > /path/to/radia/data/hg19/gencode/basic/chr22_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chrX'" > /path/to/radia/data/hg19/gencode/basic/chrX_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chrY'" > /path/to/radia/data/hg19/gencode/basic/chrY_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chrM'" > /path/to/radia/data/hg19/gencode/basic/chrM_unsorted.bed
echo "select chrom,txStart,txEnd from wgEncodeGencodeBasicV27lift37 where chrom = 'chrM'" > /path/to/radia/data/hg19/gencode/basic/chrMT_unsorted.bed

# 2) delete the header line consisting of 'chrom txStart txEnd' from each file

# 3) sort and gzip the files with commands like these:

sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr1_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr1.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr2_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr2.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr3_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr3.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr4_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr4.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr5_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr5.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr6_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr6.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr7_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr7.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr8_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr8.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr9_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr9.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr10_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr10.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr11_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr11.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr12_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr12.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr13_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr13.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr14_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr14.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr15_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr15.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr16_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr16.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr17_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr17.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr18_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr18.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr19_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr19.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr20_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr20.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr21_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr21.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chr22_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chr22.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chrX_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chrX.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chrY_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chrY.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chrM_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chrM.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/gencode/basic/chrMT_unsorted.bed | gzip > /path/to/radia/data/hg19/gencode/basic/chrMT.bed.gz

# 4) remove the *_unsorted.bed files
