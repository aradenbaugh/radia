# For more information about the 1000 Genomes project accessibility masks see here:
# http://www.internationalgenome.org/announcements/genome-accessibility-masks/

# Download the phase 3 all chromosome pilot mask for GRCh38 here:
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/PilotMask/

# Run these commands:
egrep chr1[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr1.bed
egrep chr2[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr2.bed
egrep chr3[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr3.bed
egrep chr4[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr4.bed
egrep chr5[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr5.bed
egrep chr6[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr6.bed
egrep chr7[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr7.bed
egrep chr8[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr8.bed
egrep chr9[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr9.bed
egrep chr10[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr10.bed
egrep chr11[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr11.bed
egrep chr12[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr12.bed
egrep chr13[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr13.bed
egrep chr14[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr14.bed
egrep chr15[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr15.bed
egrep chr16[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr16.bed
egrep chr17[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr17.bed
egrep chr18[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr18.bed
egrep chr19[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr19.bed
egrep chr20[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr20.bed
egrep chr21[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr21.bed
egrep chr22[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chr22.bed
egrep chrX[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chrX.bed
egrep chrY[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chrY.bed
egrep chrM[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chrM.bed
egrep chrMT[[:space:]] 20160622.allChr.pilot_mask.bed | awk '{print $1"\t"$2"\t"$3}' | sort -t $'\t' +1n -3 > /radia/data/hg38/blacklists/1000Genomes/phase3/chrMT.bed
