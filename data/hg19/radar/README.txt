# RNA-Editing events are annotated by version 2 of RADAR (Rigorously Annotated Database of A-to-I RNA Editing).  
# If you want to annotate variants with the RADAR database, 
# specify these files with the --radarDir option when filtering.
# If you don't want to annotate variants with the RADAR database, 
# use the --noRadar flag when filtering.

# To download the data yourself, follow these steps:

# 1) Download the RADAR database (v2) file for 'All Sites' for Human (hg19) from http://rnaedit.com/download/

# 2) Parse the data into chromosome specific .bed files with commands like these:  

awk '{ if ($1=="chr1") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr1_unsorted.bed
awk '{ if ($1=="chr2") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr2_unsorted.bed
awk '{ if ($1=="chr3") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr3_unsorted.bed
awk '{ if ($1=="chr4") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr4_unsorted.bed
awk '{ if ($1=="chr5") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr5_unsorted.bed
awk '{ if ($1=="chr6") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr6_unsorted.bed
awk '{ if ($1=="chr7") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr7_unsorted.bed
awk '{ if ($1=="chr8") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr8_unsorted.bed
awk '{ if ($1=="chr9") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr9_unsorted.bed
awk '{ if ($1=="chr10") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr10_unsorted.bed
awk '{ if ($1=="chr11") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr11_unsorted.bed
awk '{ if ($1=="chr12") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr12_unsorted.bed
awk '{ if ($1=="chr13") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr13_unsorted.bed
awk '{ if ($1=="chr14") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr14_unsorted.bed
awk '{ if ($1=="chr15") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr15_unsorted.bed
awk '{ if ($1=="chr16") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr16_unsorted.bed
awk '{ if ($1=="chr17") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr17_unsorted.bed
awk '{ if ($1=="chr18") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr18_unsorted.bed
awk '{ if ($1=="chr19") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr19_unsorted.bed
awk '{ if ($1=="chr20") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr20_unsorted.bed
awk '{ if ($1=="chr21") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr21_unsorted.bed
awk '{ if ($1=="chr22") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chr22_unsorted.bed
awk '{ if ($1=="chrX") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chrX_unsorted.bed
awk '{ if ($1=="chrY") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chrY_unsorted.bed
awk '{ if ($1=="chrM") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chrM_unsorted.bed
awk '{ if ($1=="chrM") {if ($7=="no") {a="AN"} else {a="AY"}; if ($8=="no") {a=a",NN"} else {a=a",NY"}; if ($9=="N") {a=a",CN"} else {a=a",CY"}; if ($10=="N") {a=a",RN"} else {a=a",RY"}; if ($11=="N") {a=a",MN"} else {a=a",MY"}; {print $1"\t"$2-1"\t"$2"\t"a}}}' /path/to/radia/data/hg19/radar/Human_AG_all_hg19_v2.txt > /path/to/radia/data/hg19/radar/chrMT_unsorted.bed

# 3) sort and gzip the files with commands like these:

sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr1_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr1.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr2_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr2.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr3_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr3.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr4_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr4.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr5_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr5.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr6_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr6.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr7_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr7.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr8_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr8.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr9_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr9.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr10_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr10.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr11_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr11.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr12_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr12.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr13_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr13.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr14_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr14.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr15_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr15.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr16_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr16.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr17_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr17.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr18_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr18.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr19_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr19.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr20_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr20.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr21_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr21.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chr22_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chr22.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chrX_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chrX.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chrY_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chrY.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chrM_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chrM.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/radar/chrMT_unsorted.bed | gzip > /path/to/radia/data/hg19/radar/chrMT.bed.gz

# 4) remove the *_unsorted.bed files
