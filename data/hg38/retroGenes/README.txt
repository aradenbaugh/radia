# Download the Retroposed Genes 

# To download the data yourself, follow these steps:

# 1) connect to the UCSC hg38 database and download the Retroposed Genes table with commands like these:

"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr1'" > /path/to/radia/data/hg38/retroGenes/chr1_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr2'" > /path/to/radia/data/hg38/retroGenes/chr2_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr3'" > /path/to/radia/data/hg38/retroGenes/chr3_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr4'" > /path/to/radia/data/hg38/retroGenes/chr4_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr5'" > /path/to/radia/data/hg38/retroGenes/chr5_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr6'" > /path/to/radia/data/hg38/retroGenes/chr6_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr7'" > /path/to/radia/data/hg38/retroGenes/chr7_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr8'" > /path/to/radia/data/hg38/retroGenes/chr8_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr9'" > /path/to/radia/data/hg38/retroGenes/chr9_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr10'" > /path/to/radia/data/hg38/retroGenes/chr10_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr11'" > /path/to/radia/data/hg38/retroGenes/chr11_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr12'" > /path/to/radia/data/hg38/retroGenes/chr12_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr13'" > /path/to/radia/data/hg38/retroGenes/chr13_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr14'" > /path/to/radia/data/hg38/retroGenes/chr14_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr15'" > /path/to/radia/data/hg38/retroGenes/chr15_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr16'" > /path/to/radia/data/hg38/retroGenes/chr16_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr17'" > /path/to/radia/data/hg38/retroGenes/chr17_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr18'" > /path/to/radia/data/hg38/retroGenes/chr18_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr19'" > /path/to/radia/data/hg38/retroGenes/chr19_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr20'" > /path/to/radia/data/hg38/retroGenes/chr20_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr21'" > /path/to/radia/data/hg38/retroGenes/chr21_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr22'" > /path/to/radia/data/hg38/retroGenes/chr22_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chrX'" > /path/to/radia/data/hg38/retroGenes/chrX_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chrY'" > /path/to/radia/data/hg38/retroGenes/chrY_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chrM'" > /path/to/radia/data/hg38/retroGenes/chrM_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chrMT'" > /path/to/radia/data/hg38/retroGenes/chrMT_unsorted.bed

# 2) delete the header line consisting of 'chrom tStart tEnd name' from each file

# 3) sort and gzip the files with commands like these:

sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr1_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr1.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr2_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr2.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr3_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr3.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr4_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr4.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr5_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr5.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr6_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr6.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr7_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr7.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr8_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr8.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr9_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr9.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr10_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr10.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr11_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr11.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr12_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr12.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr13_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr13.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr14_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr14.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr15_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr15.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr16_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr16.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr17_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr17.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr18_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr18.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr19_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr19.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr20_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr20.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr21_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr21.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chr22_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chr22.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chrX_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chrX.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chrY_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chrY.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chrM_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chrM.bed
sort -t $'\t' +1n -3 /path/to/radia/data/hg38/retroGenes/chrMT_unsorted.bed > /path/to/radia/data/hg38/retroGenes/chrMT.bed

# 4) remove the *_unsorted.bed files
