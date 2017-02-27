# The GENCODE basic gene set is replacing the old GAF (Generic Annotation File) data sets.
# If you want to filter out mutations that are not in the GENCODE basic gene set,
# these are the files that should be specified with the --targetDir option when filtering.
# If you don't want to filter out mutations based on the GENCODE basic gene set, 
# use the --noTargets flag when filtering.

# connect to the UCSC hg38 database
# download the GENCODE basic gene set with commands like these:
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr1'" > /radia/data/hg38/gencode/basic/chr1_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr2'" > /radia/data/hg38/gencode/basic/chr2_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr3'" > /radia/data/hg38/gencode/basic/chr3_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr4'" > /radia/data/hg38/gencode/basic/chr4_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr5'" > /radia/data/hg38/gencode/basic/chr5_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr6'" > /radia/data/hg38/gencode/basic/chr6_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr7'" > /radia/data/hg38/gencode/basic/chr7_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr8'" > /radia/data/hg38/gencode/basic/chr8_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr9'" > /radia/data/hg38/gencode/basic/chr9_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr10'" > /radia/data/hg38/gencode/basic/chr10_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr11'" > /radia/data/hg38/gencode/basic/chr11_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr12'" > /radia/data/hg38/gencode/basic/chr12_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr13'" > /radia/data/hg38/gencode/basic/chr13_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr14'" > /radia/data/hg38/gencode/basic/chr14_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr15'" > /radia/data/hg38/gencode/basic/chr15_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr16'" > /radia/data/hg38/gencode/basic/chr16_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr17'" > /radia/data/hg38/gencode/basic/chr17_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr18'" > /radia/data/hg38/gencode/basic/chr18_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr19'" > /radia/data/hg38/gencode/basic/chr19_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr20'" > /radia/data/hg38/gencode/basic/chr20_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr21'" > /radia/data/hg38/gencode/basic/chr21_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chr22'" > /radia/data/hg38/gencode/basic/chr22_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chrX'" > /radia/data/hg38/gencode/basic/chrX_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chrY'" > /radia/data/hg38/gencode/basic/chrY_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chrM'" > /radia/data/hg38/gencode/basic/chrM_unsorted.bed
echo "select chrom,txStart,txEnd from knownGene where chrom = 'chrM'" > /radia/data/hg38/gencode/basic/chrMT_unsorted.bed

# sort them with commands like these:
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr1_unsorted.bed > /radia/data/hg38/gencode/basic/chr1.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr2_unsorted.bed > /radia/data/hg38/gencode/basic/chr2.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr3_unsorted.bed > /radia/data/hg38/gencode/basic/chr3.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr4_unsorted.bed > /radia/data/hg38/gencode/basic/chr4.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr5_unsorted.bed > /radia/data/hg38/gencode/basic/chr5.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr6_unsorted.bed > /radia/data/hg38/gencode/basic/chr6.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr7_unsorted.bed > /radia/data/hg38/gencode/basic/chr7.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr8_unsorted.bed > /radia/data/hg38/gencode/basic/chr8.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr9_unsorted.bed > /radia/data/hg38/gencode/basic/chr9.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr10_unsorted.bed > /radia/data/hg38/gencode/basic/chr10.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr11_unsorted.bed > /radia/data/hg38/gencode/basic/chr11.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr12_unsorted.bed > /radia/data/hg38/gencode/basic/chr12.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr13_unsorted.bed > /radia/data/hg38/gencode/basic/chr13.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr14_unsorted.bed > /radia/data/hg38/gencode/basic/chr14.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr15_unsorted.bed > /radia/data/hg38/gencode/basic/chr15.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr16_unsorted.bed > /radia/data/hg38/gencode/basic/chr16.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr17_unsorted.bed > /radia/data/hg38/gencode/basic/chr17.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr18_unsorted.bed > /radia/data/hg38/gencode/basic/chr18.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr19_unsorted.bed > /radia/data/hg38/gencode/basic/chr19.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr20_unsorted.bed > /radia/data/hg38/gencode/basic/chr20.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr21_unsorted.bed > /radia/data/hg38/gencode/basic/chr21.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chr22_unsorted.bed > /radia/data/hg38/gencode/basic/chr22.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chrX_unsorted.bed > /radia/data/hg38/gencode/basic/chrX.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chrY_unsorted.bed > /radia/data/hg38/gencode/basic/chrY.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chrM_unsorted.bed > /radia/data/hg38/gencode/basic/chrM.bed
sort -t $'\t' +1n -3 /radia/data/hg38/gencode/basic/chrMT_unsorted.bed > /radia/data/hg38/gencode/basic/chrMT.bed
