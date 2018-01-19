# Use the Encode Gencode pseudo genes annotation  

# To download the data yourself, follow these steps:

# 1) connect to the UCSC hg19 database and download the Encode Gencode pseudo genes with commands like these:

"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr1'" > /path/to/radia/data/hg19/pseudoGenes/chr1_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr2'" > /path/to/radia/data/hg19/pseudoGenes/chr2_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr3'" > /path/to/radia/data/hg19/pseudoGenes/chr3_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr4'" > /path/to/radia/data/hg19/pseudoGenes/chr4_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr5'" > /path/to/radia/data/hg19/pseudoGenes/chr5_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr6'" > /path/to/radia/data/hg19/pseudoGenes/chr6_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr7'" > /path/to/radia/data/hg19/pseudoGenes/chr7_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr8'" > /path/to/radia/data/hg19/pseudoGenes/chr8_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr9'" > /path/to/radia/data/hg19/pseudoGenes/chr9_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr10'" > /path/to/radia/data/hg19/pseudoGenes/chr10_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr11'" > /path/to/radia/data/hg19/pseudoGenes/chr11_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr12'" > /path/to/radia/data/hg19/pseudoGenes/chr12_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr13'" > /path/to/radia/data/hg19/pseudoGenes/chr13_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr14'" > /path/to/radia/data/hg19/pseudoGenes/chr14_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr15'" > /path/to/radia/data/hg19/pseudoGenes/chr15_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr16'" > /path/to/radia/data/hg19/pseudoGenes/chr16_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr17'" > /path/to/radia/data/hg19/pseudoGenes/chr17_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr18'" > /path/to/radia/data/hg19/pseudoGenes/chr18_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr19'" > /path/to/radia/data/hg19/pseudoGenes/chr19_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr20'" > /path/to/radia/data/hg19/pseudoGenes/chr20_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr21'" > /path/to/radia/data/hg19/pseudoGenes/chr21_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chr22'" > /path/to/radia/data/hg19/pseudoGenes/chr22_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chrX'" > /path/to/radia/data/hg19/pseudoGenes/chrX_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chrY'" > /path/to/radia/data/hg19/pseudoGenes/chrY_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chrM'" > /path/to/radia/data/hg19/pseudoGenes/chrM_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV27lift37 where chrom = 'chrM'" > /path/to/radia/data/hg19/pseudoGenes/chrMT_unsorted.bed

# 2) delete the header line consisting of 'chrom txStart txEnd name' from each file

# 3) sort and gzip the files with commands like these:

sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr1_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr1.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr2_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr2.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr3_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr3.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr4_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr4.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr5_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr5.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr6_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr6.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr7_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr7.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr8_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr8.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr9_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr9.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr10_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr10.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr11_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr11.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr12_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr12.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr13_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr13.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr14_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr14.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr15_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr15.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr16_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr16.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr17_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr17.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr18_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr18.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr19_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr19.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr20_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr20.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr21_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr21.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chr22_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chr22.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chrX_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chrX.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chrY_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chrY.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chrM_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chrM.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/pseudoGenes/chrMT_unsorted.bed | gzip > /path/to/radia/data/hg19/pseudoGenes/chrMT.bed.gz

# 4) remove the *_unsorted.bed files
