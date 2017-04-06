# connect to the hg38 database
# download the retro genes with commands like these:
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr1'" > /radia/data/hg38/retroGenes/chr1_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr2'" > /radia/data/hg38/retroGenes/chr2_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr3'" > /radia/data/hg38/retroGenes/chr3_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr4'" > /radia/data/hg38/retroGenes/chr4_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr5'" > /radia/data/hg38/retroGenes/chr5_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr6'" > /radia/data/hg38/retroGenes/chr6_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr7'" > /radia/data/hg38/retroGenes/chr7_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr8'" > /radia/data/hg38/retroGenes/chr8_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr9'" > /radia/data/hg38/retroGenes/chr9_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr10'" > /radia/data/hg38/retroGenes/chr10_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr11'" > /radia/data/hg38/retroGenes/chr11_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr12'" > /radia/data/hg38/retroGenes/chr12_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr13'" > /radia/data/hg38/retroGenes/chr13_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr14'" > /radia/data/hg38/retroGenes/chr14_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr15'" > /radia/data/hg38/retroGenes/chr15_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr16'" > /radia/data/hg38/retroGenes/chr16_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr17'" > /radia/data/hg38/retroGenes/chr17_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr18'" > /radia/data/hg38/retroGenes/chr18_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr19'" > /radia/data/hg38/retroGenes/chr19_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr20'" > /radia/data/hg38/retroGenes/chr20_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr21'" > /radia/data/hg38/retroGenes/chr21_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chr22'" > /radia/data/hg38/retroGenes/chr22_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chrX'" > /radia/data/hg38/retroGenes/chrX_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chrY'" > /radia/data/hg38/retroGenes/chrY_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chrM'" > /radia/data/hg38/retroGenes/chrM_unsorted.bed
"select tName,tStart,tEnd from ucscRetroAli9 where tName = 'chrMT'" > /radia/data/hg38/retroGenes/chrMT_unsorted.bed

# be sure to sort them
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr1_unsorted.bed > /radia/data/hg38/retroGenes/chr1.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr2_unsorted.bed > /radia/data/hg38/retroGenes/chr2.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr3_unsorted.bed > /radia/data/hg38/retroGenes/chr3.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr4_unsorted.bed > /radia/data/hg38/retroGenes/chr4.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr5_unsorted.bed > /radia/data/hg38/retroGenes/chr5.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr6_unsorted.bed > /radia/data/hg38/retroGenes/chr6.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr7_unsorted.bed > /radia/data/hg38/retroGenes/chr7.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr8_unsorted.bed > /radia/data/hg38/retroGenes/chr8.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr9_unsorted.bed > /radia/data/hg38/retroGenes/chr9.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr10_unsorted.bed > /radia/data/hg38/retroGenes/chr10.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr11_unsorted.bed > /radia/data/hg38/retroGenes/chr11.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr12_unsorted.bed > /radia/data/hg38/retroGenes/chr12.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr13_unsorted.bed > /radia/data/hg38/retroGenes/chr13.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr14_unsorted.bed > /radia/data/hg38/retroGenes/chr14.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr15_unsorted.bed > /radia/data/hg38/retroGenes/chr15.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr16_unsorted.bed > /radia/data/hg38/retroGenes/chr16.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr17_unsorted.bed > /radia/data/hg38/retroGenes/chr17.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr18_unsorted.bed > /radia/data/hg38/retroGenes/chr18.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr19_unsorted.bed > /radia/data/hg38/retroGenes/chr19.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr20_unsorted.bed > /radia/data/hg38/retroGenes/chr20.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr21_unsorted.bed > /radia/data/hg38/retroGenes/chr21.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chr22_unsorted.bed > /radia/data/hg38/retroGenes/chr22.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chrX_unsorted.bed > /radia/data/hg38/retroGenes/chrX.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chrY_unsorted.bed > /radia/data/hg38/retroGenes/chrY.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chrM_unsorted.bed > /radia/data/hg38/retroGenes/chrM.bed
sort -t $'\t' +1n -3 /radia/data/hg38/retroGenes/chrMT_unsorted.bed > /radia/data/hg38/retroGenes/chrMT.bed
