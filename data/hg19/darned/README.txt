# RNA-Editing events are annotated by the Database of RNA Editing (DARNED). 
# If you want to annotate variants with the DARNED database, 
# specify these files with the --darnedDir option when filtering.
# If you don't want to annotate variants with the DARNED database, 
# use the --noDarned flag when filtering.

# To download the data yourself, follow these steps:

# 1) Download the Database of RNA Editing (DARNED) file for NCBI37/hg19 from https://darned.ucc.ie/static/downloads/

# 2) Map the long DARNED tissue source text to an abbreviate version using the provided mapping file: darned_tissue_mapping.txt with a command like this:
#    This command first loads the tissue mappings into 'a', then splits the
#    DARNED tissue column by ',' and stores them in 'b'. The leading and 
#    trailing spaces are removed from the elements in 'b' and the abbreviation
#    for the long text is looked up in 'a'.  

awk -F "\t" 'NR==FNR{a[$2]=$1;next}{split($9,b,",")}{gsub(/^[ ]+|[ ]+$/,"",b[1]); text=a[b[1]]; for (i=2; i<=length(b); i++) {gsub(/^[ ]+|[ ]+$/,"",b[i]); text=text","a[b[i]];}}{print $1"\t"$2"\t"text}' darned_tissue_mapping.txt darned_hg19.txt > darned_hg19_mapped.txt

# 3) Parse the data into chromosome specific .bed files with commands like these:  

awk -F '\t' '($1=="1") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr1_unsorted.bed
awk -F '\t' '($1=="2") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr2_unsorted.bed
awk -F '\t' '($1=="3") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr3_unsorted.bed
awk -F '\t' '($1=="4") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr4_unsorted.bed
awk -F '\t' '($1=="5") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr5_unsorted.bed
awk -F '\t' '($1=="6") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr6_unsorted.bed
awk -F '\t' '($1=="7") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr7_unsorted.bed
awk -F '\t' '($1=="8") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr8_unsorted.bed
awk -F '\t' '($1=="9") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr9_unsorted.bed
awk -F '\t' '($1=="10") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr10_unsorted.bed
awk -F '\t' '($1=="11") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr11_unsorted.bed
awk -F '\t' '($1=="12") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr12_unsorted.bed
awk -F '\t' '($1=="13") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr13_unsorted.bed
awk -F '\t' '($1=="14") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr14_unsorted.bed
awk -F '\t' '($1=="15") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr15_unsorted.bed
awk -F '\t' '($1=="16") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr16_unsorted.bed
awk -F '\t' '($1=="17") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr17_unsorted.bed
awk -F '\t' '($1=="18") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr18_unsorted.bed
awk -F '\t' '($1=="19") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr19_unsorted.bed
awk -F '\t' '($1=="20") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr20_unsorted.bed
awk -F '\t' '($1=="21") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr21_unsorted.bed
awk -F '\t' '($1=="22") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chr22_unsorted.bed
awk -F '\t' '($1=="X") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chrX_unsorted.bed
awk -F '\t' '($1=="Y") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chrY_unsorted.bed
awk -F '\t' '($1=="M") {print $1"\t"$2-1"\t"$2"\t"$3}' /path/to/radia/data/hg19/darned/darned_hg19_mapped.txt  > /path/to/radia/data/hg19/darned/chrM_unsorted.bed

# 3) sort and gzip the files with commands like these:

sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr1_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr1.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr2_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr2.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr3_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr3.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr4_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr4.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr5_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr5.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr6_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr6.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr7_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr7.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr8_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr8.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr9_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr9.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr10_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr10.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr11_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr11.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr12_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr12.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr13_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr13.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr14_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr14.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr15_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr15.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr16_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr16.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr17_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr17.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr18_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr18.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr19_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr19.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr20_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr20.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr21_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr21.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chr22_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chr22.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chrX_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chrX.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chrY_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chrY.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chrM_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chrM.bed.gz
sort -t $'\t' +1n -3 /path/to/radia/data/hg19/darned/chrMT_unsorted.bed | gzip > /path/to/radia/data/hg19/darned/chrMT.bed.gz

# 4) remove the *_unsorted.bed files
