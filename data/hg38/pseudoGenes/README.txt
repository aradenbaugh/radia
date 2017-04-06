# connect to the UCSC hg38 database:
# download the Endcode Gencode pseudo genes with commands like these:
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr1'" > /radia/data/hg38/pseudoGenes/chr1_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr2'" > /radia/data/hg38/pseudoGenes/chr2_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr3'" > /radia/data/hg38/pseudoGenes/chr3_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr4'" > /radia/data/hg38/pseudoGenes/chr4_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr5'" > /radia/data/hg38/pseudoGenes/chr5_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr6'" > /radia/data/hg38/pseudoGenes/chr6_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr7'" > /radia/data/hg38/pseudoGenes/chr7_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr8'" > /radia/data/hg38/pseudoGenes/chr8_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr9'" > /radia/data/hg38/pseudoGenes/chr9_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr10'" > /radia/data/hg38/pseudoGenes/chr10_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr11'" > /radia/data/hg38/pseudoGenes/chr11_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr12'" > /radia/data/hg38/pseudoGenes/chr12_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr13'" > /radia/data/hg38/pseudoGenes/chr13_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr14'" > /radia/data/hg38/pseudoGenes/chr14_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr15'" > /radia/data/hg38/pseudoGenes/chr15_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr16'" > /radia/data/hg38/pseudoGenes/chr16_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr17'" > /radia/data/hg38/pseudoGenes/chr17_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr18'" > /radia/data/hg38/pseudoGenes/chr18_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr19'" > /radia/data/hg38/pseudoGenes/chr19_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr20'" > /radia/data/hg38/pseudoGenes/chr20_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr21'" > /radia/data/hg38/pseudoGenes/chr21_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chr22'" > /radia/data/hg38/pseudoGenes/chr22_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chrX'" > /radia/data/hg38/pseudoGenes/chrX_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chrY'" > /radia/data/hg38/pseudoGenes/chrY_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chrM'" > /radia/data/hg38/pseudoGenes/chrM_unsorted.bed
"select chrom,txStart,txEnd from wgEncodeGencodePseudoGeneV24 where chrom = 'chrM'" > /radia/data/hg38/pseudoGenes/chrMT_unsorted.bed

# be sure to sort them:
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr1_unsorted.bed > /radia/data/hg38/pseudoGenes/chr1.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr2_unsorted.bed > /radia/data/hg38/pseudoGenes/chr2.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr3_unsorted.bed > /radia/data/hg38/pseudoGenes/chr3.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr4_unsorted.bed > /radia/data/hg38/pseudoGenes/chr4.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr5_unsorted.bed > /radia/data/hg38/pseudoGenes/chr5.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr6_unsorted.bed > /radia/data/hg38/pseudoGenes/chr6.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr7_unsorted.bed > /radia/data/hg38/pseudoGenes/chr7.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr8_unsorted.bed > /radia/data/hg38/pseudoGenes/chr8.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr9_unsorted.bed > /radia/data/hg38/pseudoGenes/chr9.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr10_unsorted.bed > /radia/data/hg38/pseudoGenes/chr10.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr11_unsorted.bed > /radia/data/hg38/pseudoGenes/chr11.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr12_unsorted.bed > /radia/data/hg38/pseudoGenes/chr12.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr13_unsorted.bed > /radia/data/hg38/pseudoGenes/chr13.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr14_unsorted.bed > /radia/data/hg38/pseudoGenes/chr14.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr15_unsorted.bed > /radia/data/hg38/pseudoGenes/chr15.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr16_unsorted.bed > /radia/data/hg38/pseudoGenes/chr16.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr17_unsorted.bed > /radia/data/hg38/pseudoGenes/chr17.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr18_unsorted.bed > /radia/data/hg38/pseudoGenes/chr18.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr19_unsorted.bed > /radia/data/hg38/pseudoGenes/chr19.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr20_unsorted.bed > /radia/data/hg38/pseudoGenes/chr20.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr21_unsorted.bed > /radia/data/hg38/pseudoGenes/chr21.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chr22_unsorted.bed > /radia/data/hg38/pseudoGenes/chr22.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chrX_unsorted.bed > /radia/data/hg38/pseudoGenes/chrX.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chrY_unsorted.bed > /radia/data/hg38/pseudoGenes/chrY.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chrM_unsorted.bed > /radia/data/hg38/pseudoGenes/chrM.bed
sort -t $'\t' +1n -3 /radia/data/hg38/pseudoGenes/chrMT_unsorted.bed > /radia/data/hg38/pseudoGenes/chrMT.bed
