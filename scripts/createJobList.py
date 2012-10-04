#!/usr/bin/env python2.7

import sys                          # system module
import re                           # regex module
import os                           # os module for checking pathnames
import rnaEditingUtil               # utility functions for rna editing
from optparse import OptionParser   # used for parsing command line arguments
import glob


'''
'    Amie Radenbaugh - 01/08/2011
'    UCSC - RNA Editing  
'    Program name: "createJobList.py"
'   
'    This program can be used to generate job lists that can be used with a batch system such as pyrasol to format different files that 
'    are necessary in the RNA-Editing pipeline.  This program currently supports the following formatting options:
'    
'    Preliminary Steps:
'    1) getRefGenesTx                    Download the transcription start and stop coordinates
'    2) sortRefGenesTx                   Sort the transcription start and stop coordinates
'    3) getRefGenesCDS                   Download the CDS start and stop coordinates
'    4) sortRefGenesCDS                  Sort the CDS start and stop coordinates
'    5) uniqRefGenesCDS                  Get the unique CDS start and stop coordinates
'    6) getRefGenesExons                 Download the exon start and stop coordinates
'    7) prepRefGenesExons                Prep the exon start and stop coordinates
'    8) sortRefGenesExons                Sort the exon start and stop coordinates
'    9) extractBlacklists                Extract the chromosomal blacklist regions from the main blacklist file
'    10) flankBlacklists                 Flank the blacklist regions by subtracting 200 from every start and adding 200 to every stop coordinate
'    11) downloadDbSnp                   Download the dbSNP tables
'    12) sortDbSnp                       Sort the dbSNP files
' 
'    Primary Steps:
'    1) selectMUTs                       Select all of the "MUT" lines out of a .bb file	
'    2) bbToBed                          Converting a bambam file to bed format
'    3) liftOver                         Doing a liftOver of .bed from hg18 to hg19
'    4) mergeHG18AndHG19                 Merge the original hg18 file with the new hg19 coordinates
'    5) validateLiftOver                 Validate that the liftover was done correctly	
'    6) rnaEditing                       Run rnaEditing.py on an RNA-Seq and wg or exome file
'    7) filterBlacklists                 Filter out editing events in blacklist regions  
'    8) filterTranscripts                Filter out editing events in non-transcript regions
'    9) filterCDS                        Filter out editing events in non-CDS regions
'    10) filterExons                     Filter out editing events in non-exon regions
'
'    Cheung Steps:
'    1) downloadCheungDNA                Download all of the DNA files for the 27 Cheung individuals
'
'    Testing/Debugging/Deprecated:
'    1) samtools                         Executing a samtools command on a .bam file	
'    2) ptpn6                            Select for chr 12 PTPN6 region out of .bb files
'    3) sortBed                          Sort a liftover .bed file according to the start and stop coordinates
'    4) sortBB                           Sort a .bb file according to the start and stop coordinates
'
'    The user must specify one of the formatting options described above. The user specifies the directory where the input data can be found 
'    and the directory where the converted files should be output.  For an RNA-Seq input directory, the user specifies the root RNA-Seq 
'    directory where the sub-dirs are named according to the RNA-Seq Ids followed by the file named "all_reads.bam".  The user also specifies 
'    the name of the job.list file.
'
'    The user may also need to specify other parameters depending on the formatting option that is selected.  The user must specify an Id filename 
'    containing all of the Ids for which we have both RNA-Seq and bambam data with the -i option for everything except "sortDbSNP".  An Id filename
'    can be obtained by using the createIdList.py.  The user must specify an hg18 to hg19 over chain file with the -c option when using the "liftOver" 
'    command.  The user must specify an hg19 fasta file with the -g option when using the "rnaEditing" command.  For the "rnaEditing" command, the 
'    user may also specify the output format (BED or GFF, default=BED) using the -f option, the minimum base quality (default=20) with the -b option, 
'    and the minimum mapping quality (default=20) with the -m option.
'
'    Here are some sample commands for the program:
'
'    Get the transcription, CDS and exon coordinates from the Db	
'    python2.7 createJobList.py getRefGeneTx /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/getHg19RefGenesTxJob.list -a hg19
'    python2.7 createJobList.py sortRefGeneTx /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/sortHg19RefGenesTxJob.list -a hg19 
'    python2.7 createJobList.py uniqRefGeneTx /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/uniqHg19RefGenesTxJob.list -a hg19
'
'    python2.7 createJobList.py getRefGeneCDS /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/getHg19RefGenesCDSJob.list -a hg19
'    python2.7 createJobList.py sortRefGeneCDS /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/sortHg19RefGenesCDSJob.list -a hg19  
'    python2.7 createJobList.py uniqRefGeneCDS /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/uniqHg19RefGenesCDSJob.list -a hg19
'
'    python2.7 createJobList.py getRefGeneExons /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/getHg19RefGenesExonsJob.list -a hg19 
'    python2.7 createJobList.py prepRefGeneExons /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/prepHg19RefGenesExonsJob.list -a hg19 
'    python2.7 createJobList.py sortRefGeneExons /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/sortHg19RefGenesExonsJob.list -a hg19 
'    python2.7 createJobList.py uniqRefGeneExons /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/uniqHg19RefGenesExonsJob.list -a hg19
'
'    python2.7 createJobList.py getRefGeneTx /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/getHg18RefGenesTxJob.list -a hg18
'    python2.7 createJobList.py sortRefGeneTx /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/sortHg18RefGenesTxJob.list -a hg18
'    python2.7 createJobList.py uniqRefGeneTx /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/uniqHg18RefGenesTxJob.list -a hg18
'
'    python2.7 createJobList.py getRefGeneCDS /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/getHg18RefGenesCDSJob.list -a hg18
'    python2.7 createJobList.py sortRefGeneCDS /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/sortHg18RefGenesCDSJob.list -a hg18  
'    python2.7 createJobList.py uniqRefGeneCDS /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/uniqHg18RefGenesCDSJob.list -a hg18
'
'    python2.7 createJobList.py getRefGeneExons /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/getHg18RefGenesExonsJob.list -a hg18 
'    python2.7 createJobList.py prepRefGeneExons /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/prepHg18RefGenesExonsJob.list -a hg18 
'    python2.7 createJobList.py sortRefGeneExons /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/sortHg18RefGenesExonsJob.list -a hg18
'    python2.7 createJobList.py uniqRefGeneExons /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/data/hg18/refGenes/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/uniqHg18RefGenesExonsJob.list -a hg18
'
'    Download and sort the dnSNP tables
'    python2.7 createJobList.py downloadDbSnp /inside/depot/ann/hg19/snp135/ /inside/depot/ann/hg19/snp135/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/downloadDbSNP135Job.list -a hg19
'    python2.7 createJobList.py sortDbSnp /inside/depot/ann/hg19/snp135/ /inside/depot/ann/hg19/snp135/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/sortDbSNP135Job.list -a hg19
'    python2.7 createJobList.py downloadDbSnp /inside/depot/ann/hg19/snp132/ /inside/depot/ann/hg19/snp132/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/downloadDbSNP132Job.list -a hg19
'    python2.7 createJobList.py sortDbSnp /inside/depot/ann/hg19/snp132/ /inside/depot/ann/hg19/snp132/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/sortDbSNP132Job.list -a hg19
'    python2.7 createJobList.py sortDbSnp /inside/depot/ann/hg18/snp130/ /inside/depot/ann/hg18/snp130/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/sortDbSNP130Job.list -a hg18
'    python2.7 createJobList.py downloadDbSnp /inside/depot/ann/hg18/snp129/ /inside/depot/ann/hg18/snp129/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/downloadDbSNP129Job.list -a hg18
'    python2.7 createJobList.py sortDbSnp /inside/depot/ann/hg18/snp129/ /inside/depot/ann/hg18/snp129/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/sortDbSNP129Job.list -a hg18
'
'    Extract and sort chromosomal regions of master blacklist file
'    python2.7 createJobList.py extractBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/ /inside/home/aradenba/rnaEditing/data/hg19/blacklists/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/hg19ExtractBlacklistsJob.list
'    python2.7 createJobList.py extractBlacklists /inside/home/aradenba/rnaEditing/data/hg18/blacklists/ /inside/home/aradenba/rnaEditing/data/hg18/blacklists/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/hg18ExtractBlacklistsJob.list
'
'    python2.7 createJobList.py extractBlacklists /inside/depot/blacklists/hg18/ /inside/home/aradenba/rnaEditing/data/hg18/blacklists/1kGenomes/depth/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/hg18ExtractBlacklistsJob.list
'    python2.7 createJobList.py extractBlacklists /inside/depot/blacklists/hg18/ /inside/home/aradenba/rnaEditing/data/hg18/blacklists/1kGenomes/mapQual/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/hg18ExtractBlacklistsJob.list
'    python2.7 createJobList.py extractBlacklists /inside/depot/blacklists/hg19/ /inside/home/aradenba/rnaEditing/data/hg19/blacklists/1kGenomes/depth/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/hg19ExtractBlacklistsJob.list
'    python2.7 createJobList.py extractBlacklists /inside/depot/blacklists/hg19/ /inside/home/aradenba/rnaEditing/data/hg19/blacklists/1kGenomes/mapQual/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/hg19ExtractBlacklistsJob.list
'
'
'    Add flanking regions of +/- 200 to the start and stop coordinates in the blacklist
'    python2.7 createJobList.py flankBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/ /inside/home/aradenba/rnaEditing/data/hg19/blacklists/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/hg19FlankBlacklistsJob.list -s /inside/home/aradenba/rnaEditing/data/hg19/hg19_chromSizes.tab
'    python2.7 createJobList.py flankBlacklists /inside/home/aradenba/rnaEditing/data/hg18/blacklists/ /inside/home/aradenba/rnaEditing/data/hg18/blacklists/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/hg18FlankBlacklistsJob.list -s /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes.tab
'
'    python2.7 createJobList.py getMutClassData /inside/depot/users/aradenba/rnaEditing/extracts/ /inside/depot/users/aradenba/rnaEditing/extracts/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgGetMutClassDataJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt  
'    python2.7 createJobList.py getMutClassData /inside/depot/users/aradenba/rnaEditing/extracts/ /inside/depot/users/aradenba/rnaEditing/extracts/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeGetMutClassDataJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt  
'    python2.7 createJobList.py runMutClass /inside/depot/users/aradenba/rnaEditing/extracts/mutClass/ /inside/depot/users/aradenba/rnaEditing/extracts/mutClass/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgRunMutClassJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt  
'    python2.7 createJobList.py runMutClass /inside/depot/users/aradenba/rnaEditing/extracts/mutClass/ /inside/depot/users/aradenba/rnaEditing/extracts/mutClass/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeRunMutClassJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt    
'
'    Select all of the "MUT" lines out of a .bb file  
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/aml/wg/ /inside/grotto/users/aradenba/data/hg18/aml/wg/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgGetMUTsHG18BBJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt  
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/aml/exome/ /inside/grotto/users/aradenba/data/hg18/aml/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeGetMUTsHG18BBJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/aml/wg/ /inside/grotto/users/aradenba/data/hg18/aml/wg/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgGetNCNWCsHG18BBJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt  
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/aml/exome/ /inside/grotto/users/aradenba/data/hg18/aml/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeGetNCNWCsHG18BBJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/aml/wg/ /inside/grotto/users/aradenba/data/hg18/aml/wg/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgGetNCsHG18BBJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt  
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/aml/exome/ /inside/grotto/users/aradenba/data/hg18/aml/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeGetNCsHG18BBJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/cheung/hg19/exome/ /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheung454GetMUTsBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheung454DnaIds.txt
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/cheung/hg19/exome/ /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungIlluminaGetMUTsBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIlluminaDnaIds.txt
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/cheung/hg19/exome/ /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheung454GetUSMsBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheung454DnaIds.txt
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/cheung/hg19/exome/ /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungIlluminaGetUSMsBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIlluminaDnaIds.txt
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/cheung/hg19/exome/ /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheung454GetNCsBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheung454DnaIds.txt
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/cheung/hg19/exome/ /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungIlluminaGetNCsBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIlluminaDnaIds.txt
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/cheung/hg19/exome/ /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheung454GetNWCsBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheung454DnaIds.txt
'    python2.7 createJobList.py selectMUTs /inside/grotto/bambam/cheung/hg19/exome/ /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungIlluminaGetNWCsBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIlluminaDnaIds.txt
'
'    Convert bambam to bed format 
'    python2.7 createJobList.py bbToBed /inside/grotto/users/aradenba/data/hg18/aml/wg/bambam/ /inside/grotto/users/aradenba/data/hg18/aml/wg/bed/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgBBToBedJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt  
'    python2.7 createJobList.py bbToBed /inside/grotto/users/aradenba/data/hg18/aml/exome/bambam/ /inside/grotto/users/aradenba/data/hg18/aml/exome/bed/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeBBToBedJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt
'    python2.7 createJobList.py bbToBed /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/ /inside/grotto/users/aradenba/data/hg19/cheung/exome/bed/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheung454BBToBedJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheung454DnaIds.txt  
'    python2.7 createJobList.py bbToBed /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/ /inside/grotto/users/aradenba/data/hg19/cheung/exome/bed/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungIlluminaBBToBedJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIlluminaDnaIds.txt
'    python2.7 createJobList.py radToBed /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptPatientSupport/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/radToBed/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/radToBedJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py radToBed /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/radToBed/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/radToBedJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'
'    Do liftover of .bed from hg18 to hg19
'    python2.7 createJobList.py liftOver /inside/grotto/users/aradenba/data/hg18/aml/wg/bed/ /inside/grotto/users/aradenba/data/hg19/aml/wg/bed/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgBedLiftOverJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt -c /inside/home/aradenba/rnaEditing/data/hg18ToHg19.over.chain  
'    python2.7 createJobList.py liftOver /inside/grotto/users/aradenba/data/hg18/aml/exome/bed/ /inside/grotto/users/aradenba/data/hg19/aml/exome/bed/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeBedLiftOverJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt -c /inside/home/aradenba/rnaEditing/data/hg18ToHg19.over.chain
'    python2.7 createJobList.py liftOver /inside/grotto/users/aradenba/data/hg19/cheung/exome/bed/ /inside/grotto/users/aradenba/data/hg18/cheung/exome/bed/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheung454BedLiftOverJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheung454DnaIds.txt -c /inside/home/aradenba/rnaEditing/data/hg19ToHg18.over.chain  
'    python2.7 createJobList.py liftOver /inside/grotto/users/aradenba/data/hg19/cheung/exome/bed/ /inside/grotto/users/aradenba/data/hg18/cheung/exome/bed/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungIlluminaBedLiftOverJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIlluminaDnaIds.txt -c /inside/home/aradenba/rnaEditing/data/hg19ToHg18.over.chain
'    python2.7 createJobList.py liftOver /inside/grotto/users/aradenba/data/hg18/aml/wg/bed/ /inside/grotto/users/aradenba/data/hg19/aml/wg/bed/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgBedLiftOverJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt -c /inside/home/aradenba/rnaEditing/data/hg18ToHg19.over.chain
'    python2.7 createJobList.py liftOver /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/radToBed/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/liftover/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/radLiftOverJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt -c /inside/home/aradenba/rnaEditing/data/hg19ToHg18.over.chain
'    python2.7 createJobList.py liftOver /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/radToBed/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/liftover/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/radLiftOverJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt -c /inside/home/aradenba/rnaEditing/data/hg18ToHg19.over.chain
'
'    Merge the original hg18 file with the new hg19 coordinates     
'    python2.7 createJobList.py mergeHG18AndHG19 /inside/grotto/users/aradenba/data/hg18/aml/wg/bambam/,/inside/grotto/users/aradenba/data/hg19/aml/wg/bed/ /inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgMergeHG18AndHG19Job.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt
'    python2.7 createJobList.py mergeHG18AndHG19 /inside/grotto/users/aradenba/data/hg18/aml/exome/bambam/,/inside/grotto/users/aradenba/data/hg19/aml/exome/bed/ /inside/grotto/users/aradenba/data/hg19/aml/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeMergeHG18AndHG19Job.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt
'    python2.7 createJobList.py mergeHG18AndHG19 /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/,/inside/grotto/users/aradenba/data/hg18/cheung/exome/bed/ /inside/grotto/users/aradenba/data/hg18/cheung/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheung454MergeHG18AndHG19Job.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheung454DnaIds.txt
'    python2.7 createJobList.py mergeHG18AndHG19 /inside/grotto/users/aradenba/data/hg19/cheung/exome/bambam/,/inside/grotto/users/aradenba/data/hg18/cheung/exome/bed/ /inside/grotto/users/aradenba/data/hg18/cheung/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungIlluminaMergeHG18AndHG19Job.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIlluminaDnaIds.txt
'    python2.7 createJobList.py mergeHG18AndHG19 /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptPatientSupport/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/liftover/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_transcriptPatientSupport/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/mergeRadLiftOverJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py mergeHG18AndHG19 /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_final/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/liftover/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/mergeRadLiftOverJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin10/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin20/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin50/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin75/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin150/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin250/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin350/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin500/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/rnaMinSupport5/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/rnaMinSupport10/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/rnaMin20/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterDnaReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg19_final/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dna100rna40/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDnaReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin10/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/dnaMin10/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin20/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/dnaMin20/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin50/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/dnaMin50/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin75/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/dnaMin75/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin150/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/dnaMin150/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin250/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/dnaMin250/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin350/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/dnaMin350/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dnaMin500/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/dnaMin500/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/rnaMinSupport5/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/rnaMinSupport5/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/rnaMinSupport10/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/rnaMinSupport10/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/rnaMin20/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/rnaMin20/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/dna100rna40/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/dna100rna40/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateDnaMinEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'
'    Validate that the liftover was done correctly
'    python2.7 createJobList.py validateLiftOver /inside/grotto/users/aradenba/data/hg18/aml/wg/bambam/,/inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/,/inside/depot/ann/hg18/snp130/,/inside/depot/ann/hg19/snp132/ /inside/grotto/users/aradenba/data/hg19/aml/validateLiftOver/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgValidateLiftOverJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt
'    python2.7 createJobList.py validateLiftOver /inside/grotto/users/aradenba/data/hg18/aml/exome/bambam/,/inside/grotto/users/aradenba/data/hg19/aml/exome/bambam/,/inside/depot/ann/hg18/snp130/,/inside/depot/ann/hg19/snp132/ /inside/grotto/users/aradenba/data/hg19/aml/validateLiftOver/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeValidateLiftOverJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt
'
'    Do some analysis on the BB MUT calls
'    python2.7 createJobList.py bbAnalysis /inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/ /inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgBBAnalysisJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt  
'    python2.7 createJobList.py bbAnalysis /inside/grotto/users/aradenba/data/hg19/aml/exome/bambam/ /inside/grotto/users/aradenba/data/hg19/aml/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeBBAnalysisJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt  
' 
'    Remove all the temporary files
'    python2.7 createJobList.py removeTmpFiles 
'    python2.7 createJobList.py removeTmpFiles
'
'    Dirs:
'    /inside/grotto/users/aradenba/data/hg18/aml/wg/bambam/
'    /inside/grotto/users/aradenba/data/hg18/aml/wg/bed/
'    /inside/grotto/users/aradenba/data/hg19/aml/wg/bed/
'
'    /inside/grotto/users/aradenba/data/hg18/aml/exome/bambam/
'    /inside/grotto/users/aradenba/data/hg18/aml/exome/bed/
'    /inside/grotto/users/aradenba/data/hg19/aml/exome/bed/
'
'
'    Run rnaEditing.py on an RNA-Seq file and either a WG or exome file
'    python2.7 createJobList.py rnaEditing /inside/depot/wlee/bccagsc_070610/,/inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgRNAEditingJob.list -f RAD -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt -g /inside/depot/fa/hg19.fasta -s /inside/home/aradenba/rnaEditing/data/hg19/hg19_chromSizes.tab -Q 25 -q 20
'    python2.7 createJobList.py rnaEditing /inside/depot/wlee/bccagsc_070610/,/inside/grotto/users/aradenba/data/hg19/aml/exome/bambam/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeRNAEditingJob.list -f RAD -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt -g /inside/depot/fa/hg19.fasta -s /inside/home/aradenba/rnaEditing/data/hg19/hg19_chromSizes.tab -Q 25 -q 20
'
'    First filter out all the blacklisted regions using the filterByCoordinate script 
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/blacklist/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterBlacklistsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt  
'
'    Second filter out all coordinates found in dbSNP using the filterByCoordinate script 
'    Note:  This could be more advanced -> just filter out the ones that make sense
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg19/snp132/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/blacklist/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/snp132/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterDbSnp132Job.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt  
'
'    Filter out all non-transcripts using the filterByCoordinate script 
'    python2.7 createJobList.py filterTranscripts /inside/home/aradenba/rnaEditing/data/hg19/refGenes/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/snp132/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcript/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterTranscriptsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt -a hg19  
'
'    Filter out all non-CDS regions using the filterByCoordinate script 
'    python2.7 createJobList.py filterCDS /inside/home/aradenba/rnaEditing/data/hg19/refGenes/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/snp132/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/cds/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterCDSJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt -a hg19
'
'    Filter out all non-exons using the filterByCoordinate script 
'    python2.7 createJobList.py filterExons /inside/home/aradenba/rnaEditing/data/hg19/refGenes/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcript/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exon/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterExonsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt -a hg19
'    python2.7 createJobList.py filterExons /inside/home/aradenba/rnaEditing/data/hg19/refGenes/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/snp132/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exon/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterExonsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt -a hg19  
'
'    Filter out events that don't meet all the read support requirements
'    python2.7 createJobList.py filterTranscriptReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcript/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterTranscriptReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterExonReadSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exon/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterExonReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'
'    python2.7 createJobList.py filterPatientSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptReadSupport/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptPatientSupport/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterTranscriptPatientSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterPatientSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptReadSupport/bed/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptPatientSupport/bed/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterTranscriptPatientSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterPatientSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptMutClass/input/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptMutClass/filtered/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterTranscriptPatientSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterPatientSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exonReadSupport/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exonPatientSupport/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterExonPatientSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterPatientSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exonReadSupport/bed/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exonPatientSupport/bed/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterExonPatientSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py filterPatientSupport /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exonMutClass/input/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exonMutClass/filtered/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterExonPatientSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'
'    python2.7 createJobList.py getDNAReads /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_transcriptPatientSupport/,/inside/depot/aml/wg/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_addDNA/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgGetDNAJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt
'    python2.7 createJobList.py getDNAReads /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_removeNCs/,/inside/depot/aml/exome/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_addDNA/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeGetDNAJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt
'    python2.7 createJobList.py convertStrand /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_addDNA/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_final/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/convertStrandJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'
'    python2.7 createJobList.py filterBBMuts /inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptPatientSupport/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/removeBBMuts/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgRemoveMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt
'    python2.7 createJobList.py filterBBMuts /inside/grotto/users/aradenba/data/hg19/aml/exome/bambam/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptPatientSupport/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/removeBBMuts/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeRemoveMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt
'    python2.7 createJobList.py filterBBNCs /inside/grotto/users/aradenba/data/hg18/aml/wg/bambam/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_addDNA_tmp/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_addDNA/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgRemoveMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt
'    python2.7 createJobList.py filterBBNCs /inside/grotto/users/aradenba/data/hg18/aml/exome/bambam/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_transcriptPatientSupport/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/hg18_removeNCs/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeRemoveMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt
'    python2.7 createJobList.py filterBBNCs /inside/grotto/users/aradenba/data/hg18/cheung/wg/bambam/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNCs/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNWCs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/wgRemoveMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheung454DnaIds.txt
'    python2.7 createJobList.py filterBBNCs /inside/grotto/users/aradenba/data/hg18/cheung/exome/bambam/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNCs/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNWCs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/exomeRemoveMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIlluminaDnaIds.txt
'
'    Generate stats on the filtering process
'    python2.7 createJobList.py generateFilterStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/,/inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateFilterStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'
'    Remove all the temporary Dirs: 
'    /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/blacklist/
'    /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/snp132/
'    /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcript/
'
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptPatientSupport/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/transcripts/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateTranscriptEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/removeBBMuts/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/removeBBMuts/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateRemoveMutsStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exonPatientSupport/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/v2/exons/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/generateExonEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'
'    python2.7 createJobList.py runMutClass /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptMutClass/filtered/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/transcriptMutClass/output/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/runMutClassJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt -a hg19
'    python2.7 createJobList.py runMutClass /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exonMutClass/filtered/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/exonMutClass/output/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/runMutClassJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt -a hg19
'
'    python2.7 createJobList.py ptpn6Rad /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/ptpn6/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/getPTPN6RadJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt
'
'    Download the Cheung DNA data
'    python2.7 createJobList.py downloadCheungDNA . . /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/downloadDNAJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py sortBamsCheung /inside/depot/cheung/rna/ /inside/depot/cheung/rna/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/sortBamsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py indexBamsCheung /inside/depot/cheung/rna/ /inside/depot/cheung/rna/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/indexBamsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'
'    python2.7 createJobList.py rnaEditingCheung /inside/depot/cheung/rna/,/inside/depot/cheung/dna/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungRNAEditingJob.list -f RAD -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -g /inside/depot/fa/hsap36.1-hg18.fasta -s /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes.tab -Q 25 -q 20
'    python2.7 createJobList.py rnaEditingCheung /inside/depot/cheung/rna/,/inside/depot/cheung/dna/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/test/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungRNAEditingJob.list -f RAD -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -g /inside/depot/fa/hsap36.1-hg18.fasta -s /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes.tab -Q 25 -q 20
'    python2.7 createJobList.py rnaEditingCheung /inside/depot/cheung/rna/,/inside/depot/cheung/dna/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/minOne/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungRNAEditingJob.list -f RAD -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -g /inside/depot/fa/hsap36.1-hg18.fasta -s /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes.tab -Q 25 -q 20
'    python2.7 createJobList.py rnaEditingCheung /inside/depot/cheung/rna/,/inside/depot/cheung/dna/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungRNAEditingJob.list -f RAD -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -g /inside/depot/fa/hsap36.1-hg18.fasta -s /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes.tab -Q 25 -q 20
'    python2.7 createJobList.py samtoolsSelection /inside/depot/cheung/rna/,/inside/depot/cheung/dna/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungSamtoolsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg18/blacklists/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/blacklist/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterBlacklistsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg18/snp130/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/blacklist/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/snp130/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterDbSnp130Job.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterTranscripts /inside/home/aradenba/rnaEditing/data/hg18/refGenes/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/snp130/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/transcript/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterTranscriptsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -a hg18
'    python2.7 createJobList.py filterCDS /inside/home/aradenba/rnaEditing/data/hg18/refGenes/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/snp130/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/cds/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterCDSJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -a hg18
'    python2.7 createJobList.py filterExons /inside/home/aradenba/rnaEditing/data/hg18/refGenes/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/snp130/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exon/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterExonsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -a hg18
'    python2.7 createJobList.py filterPatientSupport /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exon/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exonPatientSupport/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterExonPatientSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt   
'    python2.7 createJobList.py filterExonReadSupport /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exonPatientSupport/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterExonReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterPatientSupport /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/transcript/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/transcriptPatientSupport/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterTranscriptPatientSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterTranscriptReadSupport /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/transcriptPatientSupport/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterTranscriptReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'
'    python2.7 createJobList.py filterTranscriptReadSupport /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/transcript/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterTranscriptReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterCDSReadSupport /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/cds/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterCDSReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterExonReadSupport /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exon/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterExonReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterBBMuts /inside/grotto/users/aradenba/data/hg18/cheung/exome/bed/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exonReadSupport/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeBBMuts/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filter454BBMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheung454DnaIds.txt
'    python2.7 createJobList.py filterBBMuts /inside/grotto/users/aradenba/data/hg18/cheung/exome/bed/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exonReadSupport/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeBBMuts/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterIlluminaBBMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIlluminaDnaIds.txt
'    python2.7 createJobList.py filterBBMuts /inside/grotto/users/aradenba/data/hg18/cheung/exome/bed/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exonReadSupport/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeUSMs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filter454USMsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheung454DnaIds.txt
'    python2.7 createJobList.py filterBBMuts /inside/grotto/users/aradenba/data/hg18/cheung/exome/bed/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exonReadSupport/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeUSMs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterIlluminaUSMsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIlluminaDnaIds.txt
'    python2.7 createJobList.py filterBBMuts /inside/grotto/users/aradenba/data/hg18/cheung/exome/bed/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeUSMs/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNCs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filter454NCsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheung454DnaIds.txt
'    python2.7 createJobList.py filterBBMuts /inside/grotto/users/aradenba/data/hg18/cheung/exome/bed/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeUSMs/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNCs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterIlluminaNCsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIlluminaDnaIds.txt
'    python2.7 createJobList.py catBBMuts /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeBBMuts/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeBBMuts/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/catBBMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIntersectionDnaIds.txt
'    python2.7 createJobList.py removeCatFiles /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeBBMuts/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeBBMuts/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/removeCatFilesJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIntersectionDnaIds.txt
'    python2.7 createJobList.py catBBMuts /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeUSMs/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeUSMs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/catBBMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIntersectionDnaIds.txt
'    python2.7 createJobList.py removeCatFiles /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeUSMs/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeUSMs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/removeCatFilesJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIntersectionDnaIds.txt
'    python2.7 createJobList.py catBBMuts /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNCs/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNCs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/catBBMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIntersectionDnaIds.txt
'    python2.7 createJobList.py removeCatFiles /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNCs/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNCs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/removeCatFilesJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIntersectionDnaIds.txt
'    python2.7 createJobList.py catBBMuts /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNWCs/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNWCs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/catBBMutsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIntersectionDnaIds.txt
'
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg18/blacklists/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/blacklist/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterBlacklistsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg18/snp129/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/snp129/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterDbSnp129Job.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg18/snp130/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/snp130/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterDbSnp130Job.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg18/blacklists/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/blacklist/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterBlacklistsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg18/snp130/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/blacklist/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/snp130/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterDbSnp130Job.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterTranscripts /inside/home/aradenba/rnaEditing/data/hg18/refGenes/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/snp130/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/transcript/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterTranscriptsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -a hg18
'    python2.7 createJobList.py filterCDS /inside/home/aradenba/rnaEditing/data/hg18/refGenes/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/snp130/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/cds/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterCDSJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -a hg18
'    python2.7 createJobList.py filterExons /inside/home/aradenba/rnaEditing/data/hg18/refGenes/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/snp130/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/exon/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterExonsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -a hg18
'    python2.7 createJobList.py filterExons /inside/home/aradenba/rnaEditing/data/hg18/refGenes/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/transcript/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/exon/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterExonsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -a hg18
'
'    python2.7 createJobList.py filterPatientSupport /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/transcript/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/transcriptPatientSupport/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterTranscriptPatientSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterPatientSupport /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/exon/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/exonPatientSupport/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterExonPatientSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'
'    python2.7 createJobList.py filterTranscriptReadSupport /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/transcriptPatientSupport/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterTranscriptReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py filterExonReadSupport /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/exonPatientSupport/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/filterExonReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'
'    python2.7 createJobList.py generateFilterStats /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/,/inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/ /inside/home/aradenba/rnaEditing/cancer/cheung/stats/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/generateFilterStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py runMutClass /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/transcriptMutClass/input/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/transcriptMutClass/output/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/runMutClassJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -a hg18
'    python2.7 createJobList.py runMutClass /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exonMutClass/input/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exonMutClass/output/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/runMutClassJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt -a hg18 
'
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/transcriptReadSupport/ /inside/home/aradenba/rnaEditing/cancer/cheung/stats/transcripts/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/generateTranscriptEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/cdsReadSupport/ /inside/home/aradenba/rnaEditing/cancer/cheung/stats/cds/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/generateCDSEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exonReadSupport/ /inside/home/aradenba/rnaEditing/cancer/cheung/stats/exons/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/generateExonEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/exonReadSupport/ /inside/home/aradenba/rnaEditing/cancer/cheung/stats/exons2/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/generateExonEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeBBMuts/ /inside/home/aradenba/rnaEditing/cancer/cheung/stats/removeBBMuts/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/generateExonEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungUnionDnaIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/filters/removeNCs/ /inside/home/aradenba/rnaEditing/cancer/cheung/stats/removeNCs/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/generateExonEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungUnionDnaIds.txt
'
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/transcriptReadSupport/ /inside/home/aradenba/rnaEditing/cancer/cheung/stats/transcripts/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/generateTranscriptEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'    python2.7 createJobList.py generateEventStats /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/exonReadSupport/ /inside/home/aradenba/rnaEditing/cancer/cheung/stats/exons/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/generateExonEventStatsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'
'    python2.7 createJobList.py runMutClass /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/transcriptMutClass/input/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/transcriptMutClass/output/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/runMutClassJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt -a hg18
'    python2.7 createJobList.py runMutClass /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/exonMutClass/input/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/exonMutClass/output/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/runMutClassJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/allOverlappingIds.txt -a hg18
'    python2.7 createJobList.py mutNullModel /inside/home/aradenba/rnaEditing/data/hg19/refGenes/ /inside/grotto/users/aradenba/data/hg19/mutNullModel/input/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/mutNullModelJob.list
'    python2.7 createJobList.py runNullMutClass /inside/grotto/users/aradenba/data/hg19/mutNullModel/input/ /inside/grotto/users/aradenba/data/hg19/mutNullModel/output/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/runMutClassJob.list -a hg19
'
'    python2.7 createJobList.py samtoolsSelection /inside/depot/cheung/rna/ /inside/home/aradenba/rnaEditing/cancer/cheung/stats/ /inside/home/aradenba/rnaEditing/cancer/cheung/jobLists/cheungSamtoolsJob.list -i /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt
'
'    Old - No longer needed?
'    Filter events and generate stat files for plotting 
'    python2.7 createJobList.py filterEvents /inside/depot/users/aradenba/dbSNP/,/inside/home/aradenba/rnaEditing/data/refGenes/,/inside/depot/users/aradenba/rnaEditing/ /inside/depot/users/aradenba/rnaEditing/extracts/,/inside/depot/users/aradenba/stats/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/filterEventsJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt,/inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt  
'
'    Sort a .bb file according to the start and stop coordinates 
'    python2.7 createJobList.py sortBB /inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/ /inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgSortHG19BBJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt  
'    python2.7 createJobList.py sortBB /inside/grotto/users/aradenba/data/hg19/aml/exome/bambam/ /inside/grotto/users/aradenba/data/hg19/aml/exome/bambam/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeSortHG19BBJob.list -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt  
'
'    Using the samtools command, select the chr12 PTPN6 region from the RNA-Seq files
'    python2.7 createJobList.py samtools /inside/depot/wlee/bccagsc_070610/ /inside/grotto/users/aradenba/data/hg19/aml/rnaSeq/sam/chr12_ptpn6/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/rnaSeqJob.list -d _chr12_ptpn6 -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/rnaEditIds.txt  
'
'    Select for chr12_ptpn6 out of .bb files
'    python2.7 createJobList.py ptpn6 /inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/ /inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/chr12_ptpn6/ /inside/home/aradenba/rnaEditing/cancer/aml/jobLists/wgBBSelectionJob.list -d _chr12_ptpn6 -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt  
'    python2.7 createJobList.py ptpn6 /inside/grotto/users/aradenba/data/hg19/aml/exome/bambam/ /inside/grotto/users/aradenba/data/hg19/aml/exome/bambam/chr12_ptpn6/ -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt/inside/home/aradenba/rnaEditing/cancer/aml/jobLists/exomeBBSelectionJob.list -d _chr12_ptpn6 -i /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt
'
'    python2.7 createJobList.py rnaDnaDiffs /inside/grotto/astrocytoma/hg18Files/ /inside/grotto/astrocytoma/tripleBam/ /inside/home/aradenba/rnaEditing/cancer/astrocytoma/jobLists/rnaDnaDiffsJob.list -g /inside/depot/fa/all_sequences.fasta -s /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes_sorted.tab -Q 20 -q 20
'
'
'
'    Run radia.py on benchmark files
'    Extract and sort chromosomal regions of master TCGA target file
'    python2.7 createJobList.py extractTcgaTargets /inside/home/cwilks/cancer/tcga/vcf/data/ /inside/home/aradenba/rnaEditing/data/hg19/tcgaTargets/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/hg19ExtractTcgaTargetsJob.list
'    python2.7 createJobList.py extractTcgaTargets /inside/home/cwilks/cancer/tcga/vcf/data/ /inside/home/aradenba/rnaEditing/data/hg18/tcgaTargets/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/hg18ExtractTcgaTargetsJob.list
'    python2.7 createJobList.py extractTcgaTargets /inside/home/singer/data/bed/ /inside/home/aradenba/rnaEditing/data/hg19/washuTargets/ /inside/home/aradenba/rnaEditing/cancer/common/jobLists/hg19ExtractWashuTargetsJob.list
'
'    Create grade3 astrocytoma symbolic links
'    cut -f 1,3 /inside/depot2/grade3_astrocytoma/bams/id_to_bam_map.txt.sorted | awk '{print "ln -s /inside/depot2/grade3_astrocytoma/bams/"$1"*.bam /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/exome/"$2".bam"}' > /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/symbolicLinksDNA.list
'    cut -f 1,3 /inside/depot2/grade3_astrocytoma/bams/id_to_bam_map.txt.sorted | awk '{print "cp /inside/depot2/grade3_astrocytoma/bams/"$1"*.bam.bai /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/exome/"$2".bam.bai"}' > /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/cpBaiDNA.list
'    cut -f 1,3 /inside/depot2/grade3_astrocytoma/rna/id_to_rna_bam_map.txt | awk '{print "ln -s /inside/depot2/grade3_astrocytoma/rna/"$1"*.bam /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/rna/"$2"-RNA-seq.bam"}' > /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/symbolicLinksRNA.list
'    cut -f 1,3 /inside/depot2/grade3_astrocytoma/rna/id_to_rna_bam_map.txt | awk '{print "cp /inside/depot2/grade3_astrocytoma/rna/"$1"*.bam.bai /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/rna/"$2"-RNA-seq.bam.bai"}' > /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/cpBaiRNA.list
'
'    Create grade3 astrocytoma ids and then manually edit
'    cut -f 3 /inside/depot2/grade3_astrocytoma/rna/id_to_rna_bam_map.txt > /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/dna.ids
'    cut -f 3 /inside/depot2/grade3_astrocytoma/rna/id_to_rna_bam_map.txt > /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/rna.ids
'    cut -f 3 /inside/depot2/grade3_astrocytoma/rna/id_to_rna_bam_map.txt > /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/overlapping.ids
'
'    Get non-hypermutated UCEC ids from Broad .maf
'    Download An_UCEC_194.aggregated.tcga.maf from TCGA website, then:
'    grep "TCGA-" An_UCEC_194.aggregated.tcga.maf | cut -f 16 | sort -u | cut -c 1-12 > ucec_non_hypermutated_dna.ids
'
'    python2.7 createJobList.py samtoolsIndex /inside/grotto/users/aradenba/data/hg19/eric/rna/ /inside/grotto/users/aradenba/data/hg19/eric/rna/ /inside/home/aradenba/rnaEditing/cancer/eric/jobLists/samtoolsIndex.list -i /inside/home/aradenba/rnaEditing/cancer/eric/idLists/rna.ids
'
'    python2.7 createJobList.py symbolicLinks /inside/home/aradenba/rnaEditing/cancer/brca/dirList/ /inside/grotto/users/aradenba/data/hg19/brca/ /inside/home/aradenba/rnaEditing/cancer/cancer/brca/jobLists/symbolicLink.list
'
'    python2.7 createJobList.py samtoolsSelection /inside/grotto/users/aradenba/data/hg18/gbm/exome/,/inside/grotto/users/aradenba/data/hg18/gbm/rna/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/samtoolsTestHg18Job.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py samtoolsSelection /inside/grotto/users/aradenba/data/hg19/gbm/exome/,/inside/grotto/users/aradenba/data/hg19/gbm/rna/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/samtoolsTestHg19Job.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py samtoolsSelection /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/exome/,/inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/rna/ /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/ /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/samtoolsTestGrade3Job.list -i /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/overlapping.ids
'    python2.7 createJobList.py samtoolsSelection /inside/grotto/users/aradenba/data/hg18/grade2_astrocytoma/exome/,/inside/grotto/users/aradenba/data/hg18/grade2_astrocytoma/rna/ /inside/grotto/users/aradenba/data/hg18/grade2_astrocytoma/radia/ /inside/home/aradenba/rnaEditing/cancer/grade2_astrocytoma/jobLists/samtoolsTestJob.list -i /inside/home/aradenba/rnaEditing/cancer/grade2_astrocytoma/idLists/overlapping.ids
'
'    python2.7 createJobList.py radia /inside/home/cwilks/bb_pipeline/runs/mutcomp_ex3_usrerun/bams/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt -g /inside/depot/fa/all_sequences.fasta -s /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes.tab -Q 10 -q 10
'    python2.7 createJobList.py radia /inside/grotto/users/aradenba/data/hg18/gbm/exome/,/inside/grotto/users/aradenba/data/hg18/gbm/rna/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaHg18Job.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids -g /inside/depot/fa/NCBI36-HG18_Broad_variant.fasta -s /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes.tab -Q 10 -q 10
'    python2.7 createJobList.py radia /inside/grotto/users/aradenba/data/hg19/gbm/exome/,/inside/grotto/users/aradenba/data/hg19/gbm/rna/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaHg19Job.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids -g /inside/depot/fa/Homo_sapiens_assembly19.fasta -s /inside/home/aradenba/rnaEditing/data/hg19/hg19_chromSizes.tab -Q 10 -q 10
'    python2.7 createJobList.py radia /inside/grotto/users/aradenba/data/hg18/gbm/wg/,/inside/grotto/users/aradenba/data/hg18/gbm/rna/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/wg/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaHg18WGJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18_wg.ids -g /inside/depot/fa/NCBI36-HG18_Broad_variant.fasta -s /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes.tab -Q 10 -q 10
'    python2.7 createJobList.py radia /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/exome/,/inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/rna/ /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/ /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/radiaJob.list -i /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/overlapping.ids -g /inside/depot/fa/all_sequences.fasta -s /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes.tab -Q 10 -q 10
'    python2.7 createJobList.py radia /inside/grotto/users/aradenba/data/hg18/grade2_astrocytoma/exome/,/inside/grotto/users/aradenba/data/hg18/grade2_astrocytoma/rna/ /inside/grotto/users/aradenba/data/hg18/grade2_astrocytoma/radia/ /inside/home/aradenba/rnaEditing/cancer/grade2_astrocytoma/jobLists/radiaJob.list -i /inside/home/aradenba/rnaEditing/cancer/grade2_astrocytoma/idLists/overlapping.ids -g /inside/depot/fa/all_sequences.fasta -s /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes.tab -Q 10 -q 10
'    python2.7 createJobList.py radia /inside/home/singer/data/ucecfreeze/bams1/,/inside/home/singer/data/ucecfreeze/bams2/,/inside/home/singer/data/ucecfreeze/bams3/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/radiaJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids -g /inside/depot/fa/GRCh37-lite.fa -s /inside/home/aradenba/rnaEditing/data/hg19/hg19_chromSizes.tab -Q 10 -q 10
'
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg18/blacklists/,/inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/blacklist/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterBlacklistsJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg18/blacklists/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/blacklist/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterHg18BlacklistsJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/blacklist/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterHg19BlacklistsJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg18/blacklists/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/wg/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/wg/filters/blacklist/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterHg18WGBlacklistsJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18_wg.ids
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg18/blacklists/,/inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/ /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/blacklist/ /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/filterBlacklistsJob.list -i /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/overlapping.ids
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/qsubfilterBlacklistsJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/1000Genomes/depth/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_depth/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/qsubFilterBlacklistsDepthJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/1000Genomes/mapQual/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_depth/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_mapQual/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/qsubFilterBlacklistsMapQJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/1000Genomes/depth/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_depth_fixed/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/filterBlacklistsDepthJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/1000Genomes/mapQual/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_depth_fixed/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_mapQual_fixed/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/filterBlacklistsMapQJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/1000Genomes/depth/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_depth_pybed/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/filterBlacklistsDepthPybedJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/1000Genomes/mapQual/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_depth_pybed/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_mapQual_pybed/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/filterBlacklistsMapQPybedJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterBlacklists /inside/home/aradenba/rnaEditing/data/hg19/blacklists/,/inside/depot4/radia/brca/hg19/raw/ /inside/depot4/radia/brca/hg19/exome/filters/blacklist/ /inside/home/aradenba/rnaEditing/cancer/brca/jobLists/qsubFilterBlacklistsJob.list -i /inside/home/aradenba/rnaEditing/cancer/brca/idLists/overlappingIds.txt
'
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg18/snp130/,/inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/blacklist/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/snp130/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterDbSnp130Job.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg18/snp130/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/blacklist/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/snp130/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterDbSnp130Job.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg19/snp135/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/blacklist/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/snp135/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterDbSnp135Job.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg18/snp130/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/wg/filters/blacklist/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/wg/filters/snp130/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterDbSnp130WGJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18_wg.ids
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg18/snp130/,/inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/blacklist/ /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/snp130/ /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/filterDbSnp130Job.list -i /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/overlapping.ids
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg19/snp135/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/snp135/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/filterDbSnp135Job.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg19/snp130/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/snp130/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/filterDbSnp130Job.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterDbSnp /inside/depot/ann/hg19/snp135/,/inside/depot4/radia/brca/hg19/exome/filters/blacklist/ /inside/depot4/radia/brca/hg19/exome/filters/snp135/ /inside/home/aradenba/rnaEditing/cancer/brca/jobLists/filterDbSnp135Job.list -i /inside/home/aradenba/rnaEditing/cancer/brca/idLists/overlappingIds.txt
'
'    python2.7 createJobList.py filterTargets /inside/home/aradenba/rnaEditing/data/hg18/tcgaTargets/,/inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/snp130/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterTcgaTargetsJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterTargets /inside/home/aradenba/rnaEditing/data/hg18/tcgaTargets/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/snp130/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/tcgaTargets/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterTcgaTargetsJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py filterTargets /inside/home/aradenba/rnaEditing/data/hg19/tcgaTargets/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/snp135/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/tcgaTargets/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterTcgaTargetsJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterTargets /inside/home/aradenba/rnaEditing/data/hg18/tcgaTargets/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/wg/filters/snp130/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/wg/filters/tcgaTargets/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterTcgaTargetsJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18_wg.ids
'    python2.7 createJobList.py filterTargets /inside/home/aradenba/rnaEditing/data/hg19/tcgaTargets/,/inside/grotto/users/aradenba/data/hg19/gbm/newBamBam/ /inside/grotto/users/aradenba/data/hg19/gbm/newBamBam/filters/tcgaTargets/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterTcgaTargetsJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterTargets /inside/home/aradenba/rnaEditing/data/hg19/tcgaTargets/,/inside/grotto/users/aradenba/data/hg19/gbm/oldBamBam/ /inside/grotto/users/aradenba/data/hg19/gbm/oldBamBam/filters/tcgaTargets/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterTcgaTargetsJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterTargets /inside/home/aradenba/rnaEditing/data/hg18/tcgaTargets/,/inside/grotto/users/aradenba/data/hg18/benchmark/newBamBam/ /inside/grotto/users/aradenba/data/hg18/benchmark/newBamBam/filters/tcgaTargets/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterTcgaTargetsJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterTargets /inside/home/aradenba/rnaEditing/data/hg18/tcgaTargets/,/inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/snp130/ /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/tcgaTargets/ /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/filterTcgaTargetsJob.list -i /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/overlapping.ids
'    python2.7 createJobList.py filterTargets /inside/home/aradenba/rnaEditing/data/hg19/tcgaTargets/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/snp135/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/tcgaTargets/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/filterTcgaTargetsJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterTargets /inside/home/aradenba/rnaEditing/data/hg19/washuTargets/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/snp135/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/washuTargets/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/filterWashuTargetsJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterTargets /inside/home/aradenba/rnaEditing/data/hg19/washuTargets/,/inside/depot4/radia/brca/hg19/exome/filters/snp135/ /inside/depot4/radia/brca/hg19/exome/filters/washuTargets/ /inside/home/aradenba/rnaEditing/cancer/brca/jobLists/filterWashuTargetsJob.list -i /inside/home/aradenba/rnaEditing/cancer/brca/idLists/overlappingIds.txt
'
'    python2.7 createJobList.py createMutClassInput /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportMAB4Err01/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/mutClassInput/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/createMutInput.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py runMutClass /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/mutClassInput/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/mutClassOutput/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/runMutClass.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt -a hg18
'    python2.7 createJobList.py mergeVCFAndMutClassOutput /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportMAB4Err01/,/inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/mutClassOutput/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/mutClassMerged/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/mergeMutClass.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'
'    python2.7 createJobList.py createMutClassInput /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassInput/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/createMutInput.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py runMutClass /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassInput/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassOutput/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/runMutClass.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids -a hg19
'    python2.7 createJobList.py mergeVCFAndMutClassOutput /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/tcgaTargets/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassOutput/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/mergeMutClass.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'
'    python2.7 createJobList.py createMutClassInput /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/mutClassInput/ /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/createMutInput.list -i /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/overlapping.ids
'    python2.7 createJobList.py runMutClass /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/mutClassInput/ /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/mutClassOutput/ /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/runMutClass.list -i /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/overlapping.ids -a hg18
'    python2.7 createJobList.py mergeVCFAndMutClassOutput /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/tcgaTargets/,/inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/mutClassOutput/ /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/mutClassMerged/ /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/mergeMutClass.list -i /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/overlapping.ids
' 
'    python2.7 createJobList.py createMutClassInput /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/washuTargets/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/mutClassInput/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/createMutInput.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py runMutClass /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/mutClassInput/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/mutClassOutput/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/runMutClass.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids -a hg19
'    python2.7 createJobList.py runMutClass /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/mutClassInput/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/mutClassOutput/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/runMutClassPK.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids -a hg19
'    python2.7 createJobList.py mergeVCFAndMutClassOutput /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/tcgaTargets/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/mutClassOutput/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/mutClassMerged/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/mergeMutClass.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'
'    python2.7 createJobList.py createMutClassInput /inside/depot4/radia/brca/hg19/exome/filters/washuTargets/ /inside/depot4/radia/brca/hg19/exome/filters/mutClassInput/ /inside/home/aradenba/rnaEditing/cancer/brca/jobLists/createMutInput.list -i /inside/home/aradenba/rnaEditing/cancer/brca/idLists/overlappingIds.txt
'    python2.7 createJobList.py runMutClass /inside/depot4/radia/brca/hg19/exome/filters/mutClassInput/ /inside/depot4/radia/brca/hg19/exome/filters/mutClassOutput/ /inside/home/aradenba/rnaEditing/cancer/brca/jobLists/runMutClass.list -i /inside/home/aradenba/rnaEditing/cancer/brca/idLists/overlappingIds.txt -a hg19
'    python2.7 createJobList.py mergeVCFAndMutClassOutput /inside/depot4/radia/brca/hg19/exome/filters/tcgaTargets/,/inside/depot4/radia/brca/hg19/exome/filters/mutClassOutput/ /inside/depot4/radia/brca/hg19/exome/filters/mutClassMerged/ /inside/home/aradenba/rnaEditing/cancer/brca/jobLists/mergeMutClass.list -i /inside/home/aradenba/rnaEditing/cancer/brca/idLists/overlappingIds.txt
'
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportMAB2/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt    
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportMAB4/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportMAB2Err01/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt    
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportMAB4Err01/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBB/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT4/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT6/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT8/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT10/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT10MP3/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT10MP5/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT10MP8/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT10MP10/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/rs10_4_10_01/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/rs10_4_10_01_01/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNAStrictRNANone/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrict/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrict55/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrict1002/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrict1001/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrictAB4/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrictAll/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrictReallyAll/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNALoose/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNANone/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNALooseRNAStrict/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10MA1/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT20/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT20MA1/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAStrict/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT30/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT30MA1/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01_plus/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsRNAStrict_10_8/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10_strbias/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10MA1_strbias/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAStrict_strbias/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids 
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNALoose_strbias/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNALooseRNAStrict_strbias/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/rs10_4_10_01_01/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/washuTargets/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/rs10_4_10_01_01/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py filterReadSupport /inside/depot4/radia/brca/hg19/exome/filters/washuTargets/ /inside/depot4/radia/brca/hg19/exome/filters/rs10_4_10_01_01/ /inside/home/aradenba/rnaEditing/cancer/brca/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/brca/idLists/overlappingIds.txt
'    
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/mutClassMerged/ /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/radia/filters/rs10_4_10_01_01/ /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/grade3_astrocytoma/idLists/overlapping.ids
'
'    python2.7 createJobList.py extractInfo /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/depths/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterExtractInfoJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py extractInfo /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/blacklist/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/depths/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterExtractInfoJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py extractInfo /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/depths/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterExtractInfoJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py extractInfo /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/strbias/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/extractSBJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py extractInfo /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/strbias/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/extractSBJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py extractInfo /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/strbias/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/extractSBJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py extractInfo /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/strbias/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/extractSBJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py extractInfo /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNALooseRNAStrict/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/strbias/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/extractSBJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py extractInfo /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsRNAStrict_10_8/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/strbias/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/extractSBJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNAStrictRNALoose2/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py filterReadSupport /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrict2/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/filterReadSupportJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'
'    python2.7 createJobList.py extractVCF /inside/home/singer/data/crc/UCSCbambamtargeted/ /inside/grotto/users/aradenba/data/hg18/benchmark/bambam/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/extractBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py extractVCF /inside/home/singer/data/crc/BaylorIllumina/ /inside/grotto/users/aradenba/data/hg18/benchmark/baylor/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/extractBaylorJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py extractVCF /inside/home/singer/data/crc/Broad/ /inside/grotto/users/aradenba/data/hg18/benchmark/solid/broad/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/extractBroadJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py extractVCF /inside/grotto/bambam/coad_read/exome/testsinger/illumina40bambam3qual40/ /inside/grotto/users/aradenba/data/hg18/benchmark/newBamBam/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/extractBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'
'    python2.7 createJobList.py extractVCF /inside/grotto/bambam/gbm/exome/ /inside/grotto/users/aradenba/data/hg19/gbm/newBamBam/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/extractBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py extractVCF /inside/grotto/bambam/gbm/exomeold/ /inside/grotto/users/aradenba/data/hg19/gbm/oldBamBam/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/extractBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'
'    python2.7 createJobList.py extractVCF /inside/depot3/bambam/endometrial/washutarget/ /inside/grotto/users/aradenba/data/hg19/ucec/bambam/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/extractBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py extractVCF /inside/depot3/bambam/endometrial/alltarget/ /inside/grotto/users/aradenba/data/hg19/ucec/bambam/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/extractBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py extractMAF /inside/home/singer/data/maf/endometrial/ /inside/grotto/users/aradenba/data/hg19/ucec/maf/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/extractMafJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py extractMAF /inside/home/singer/data/maf/endometrial/ /inside/grotto/users/aradenba/data/hg19/ucec/mafWithAllBroad/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/extractMafJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py catChromVCFFiles /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/rs10_4_10_01_01/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/patientVCFs/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/catFiles.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py extractMAF /inside/home/singer/data/maf/brca/ /inside/depot4/radia/brca/hg19/maf/ /inside/home/aradenba/rnaEditing/cancer/brca/jobLists/extractMafJob.list -i /inside/home/aradenba/rnaEditing/cancer/brca/idLists/overlappingIds.txt
'    python2.7 createJobList.py extractVCF /inside/depot3/bambam/brca/targetted/ /inside/depot4/radia/brca/hg19/oldBamBam/ /inside/home/aradenba/rnaEditing/cancer/brca/jobLists/extractBBJob.list -i /inside/home/aradenba/rnaEditing/cancer/brca/idLists/overlappingIds.txt
'
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/rs10_4_10_01_01/,/inside/grotto/users/aradenba/data/hg19/ucec/bambam/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/stats/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/rs10_4_10_01_01/overlapsBB/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/rs10_4_10_01_01/nonOvlpBB/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/rs10_4_10_01_01/,/inside/grotto/users/aradenba/data/hg19/ucec/bambam/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/stats/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/rs10_4_10_01_01/overlapsBBwDB/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/rs10_4_10_01_01/nonOvlpBBwDB/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_mapQual/,/inside/grotto/users/aradenba/data/hg19/ucec/bambam/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/stats/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_mapQual/overlapsBB/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_mapQual/nonOvlpBB/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py radiaCompareMaf /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/patientVCFs/,/inside/grotto/users/aradenba/data/hg19/ucec/mafWithAllBroad/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/stats/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/patientVCFs/overlapsMafWithBroad/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/patientVCFs/nonOvlpMafWithBroad/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py radiaCompareMaf /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/patientVCFs/,/inside/grotto/users/aradenba/data/hg19/ucec/maf/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/stats/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/patientVCFs/overlapsMaf/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/patientVCFs/nonOvlpMaf/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_depth_fixed/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_depth_pybed/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/stats/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_depth_fixed/overlapsPybed/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_depth_fixed/nonOvlpPybed/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_mapQual_fixed/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_mapQual_pybed/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/stats/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_mapQual_fixed/overlapsPybed/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/blacklist_mapQual_fixed/nonOvlpPybed/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/snp130/,/inside/grotto/users/aradenba/data/hg19/ucec/bambam/ /inside/grotto/users/aradenba/data/hg19/ucec/radia/stats/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/snp130/overlapsBB/,/inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/snp130/nonOvlpBB/ /inside/home/aradenba/rnaEditing/cancer/ucec/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/ucec/idLists/ucec_non_hypermutated_dna.ids
'
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportMAB4Err01/,/inside/grotto/users/aradenba/data/hg18/benchmark/bambam/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBB/,/inside/grotto/users/aradenba/data/hg18/benchmark/bambam/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT4/,/inside/grotto/users/aradenba/data/hg18/benchmark/bambam/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT10MP3/,/inside/grotto/users/aradenba/data/hg18/benchmark/bambam/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBB/,/inside/grotto/users/aradenba/data/hg18/benchmark/baylor/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT10/,/inside/grotto/users/aradenba/data/hg18/benchmark/baylor/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT10MP3/,/inside/grotto/users/aradenba/data/hg18/benchmark/baylor/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBBMT10MP10/,/inside/grotto/users/aradenba/data/hg18/benchmark/baylor/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/rs10_4_10_01/,/inside/grotto/users/aradenba/data/hg18/benchmark/oldBamBam/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/,/inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/rs10_4_10_01/overlapsBB/,/inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/rs10_4_10_01/nonOvlpBB/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/rs10_4_10_01/,/inside/grotto/users/aradenba/data/hg18/benchmark/newBamBam/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/,/inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/rs10_4_10_01/overlapsNewBB/,/inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/rs10_4_10_01/nonOvlpNewBB/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/rs10_4_10_01_01/,/inside/grotto/users/aradenba/data/hg18/benchmark/newBamBam/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/,/inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/rs10_4_10_01_01/overlapsNewBB/,/inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/rs10_4_10_01_01/nonOvlpNewBB/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/benchmark/oldBamBam/,/inside/grotto/users/aradenba/data/hg18/benchmark/newBamBam/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/,/inside/grotto/users/aradenba/data/hg18/benchmark/oldBamBam/overlapsNewBB/,/inside/grotto/users/aradenba/data/hg18/benchmark/oldBamBam/nonOvlpNewBB/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNAStrictRNANone/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrict/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/stats/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNAStrictRNANone/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrict55/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/stats/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNAStrictRNANone/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrict1002/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/stats/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNAStrictRNANone/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrict1001/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/stats/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNAStrictRNANone/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrictAB4/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/stats/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNAStrictRNANone/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrictAll/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/stats/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNAStrictRNANone2/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrict10_4_10_01/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/stats/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNAStrictRNANone2/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrictReallyAll/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/stats/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNAStrictRNALoose/,/inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/filters/rsDNALooseRNAStrict1001/ /inside/grotto/users/aradenba/data/hg18/gbm/radia/exome/stats/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg18.ids
'
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNALoose/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNALooseRNAStrict/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNALoose/overlapsStrictRNA/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNALoose/nonOvlpStrictRNA/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNANone/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNALooseRNAStrict/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNANone/overlapsStrictRNA/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNANone/nonOvlpStrictRNA/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10MA1/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10/overlapsMT10MA1/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10/nonOvlpMT10MA1/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/,/inside/grotto/users/aradenba/data/hg19/gbm/newBamBam/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/overlapsBB/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/nonOvlpBB/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/,/inside/grotto/users/aradenba/data/hg19/gbm/oldBamBam/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/overlapsBB/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/nonOvlpBB/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/,/inside/grotto/users/aradenba/data/hg19/gbm/oldBamBam/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/overlapsOldBB/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01/nonOvlpOldBB/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/oldBamBam/filters/tcgaTargets/,/inside/grotto/users/aradenba/data/hg19/gbm/newBamBam/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/oldBamBam/filters/tcgaTargets/overlapsOldBB/,/inside/grotto/users/aradenba/data/hg19/gbm/oldBamBam/filters/tcgaTargets/nonOvlpOldBB/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01/,/inside/grotto/users/aradenba/data/hg19/gbm/newBamBam/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01/overlapsNewBB/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01/nonOvlpNewBB/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01_plus/,/inside/grotto/users/aradenba/data/hg19/gbm/newBamBam/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01_plus/overlapsNewBB/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01_plus/nonOvlpNewBB/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01_plus/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAStrict_strbias/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01_plus/overlapsStrictRNA/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rs10_4_10_01_01_plus/nonOvlpStrictRNA/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10_strbias/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10MA1_strbias/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10_strbias/overlapsMT10MA1/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAMT10_strbias/nonOvlpMT10MA1/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'    python2.7 createJobList.py radiaCompare /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNALoose_strbias/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNALooseRNAStrict_strbias/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/stats/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNALoose_strbias/overlapsStrictRNA/,/inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNALoose_strbias/nonOvlpStrictRNA/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/radiaCompareJob.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids
'
'    python2.7 createJobList.py vcf2Excel /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAStrict2/ /inside/grotto/users/aradenba/data/hg19/gbm/radia/exome/filters/rsDNAStrictRNAStrict2/excel/ /inside/home/aradenba/rnaEditing/cancer/gbm/jobLists/vcf2Excel.list -i /inside/home/aradenba/rnaEditing/cancer/gbm/idLists/hg19.ids -a hg19
'
'    # get the pk jobs that failed:  
'    find . -maxdepth 1 -type f -size +1000c | grep "chr14"
'
'    # get the pk jobs that are supposedly running
'    find . -maxdepth 1 -type f -size -55c | grep "chr14"
'
'    find . -maxdepth 1 -type f | grep "chr11" | egrep -e "TCGA-(\\w){2}-(\\w){4}" -o | sort -u > successfulChr11.list
'    find . -maxdepth 1 -type f -size 53c | grep "chr11" | egrep -e "TCGA.*" -o | awk '{print $1}'
'    qstat -u "aradenba" | grep "Eqw" | awk '{print $1}' | xargs -i echo qmod -cj {}
'    find . -maxdepth 1 -type f -size 53c | grep "chr11" | egrep -e "TCGA-(\\w){2}-(\\w){4}" -o | awk '{print "rnaEditing/cancer/brca/jobLists/qsub/"$1"_chr11.sh"}' | xargs -i echo qsub {}
'    find . -maxdepth 1 -type f -size 53c | grep "chr14" | egrep -e "TCGA-(\\w){2}-(\\w){4}_chr14.sh.e(\\w){7}" -o | awk '{print $1}' | xargs -i rm {}
'
'    grep "GERM" *06-0879* | grep "PASS" | grep -v "LOH" | wc
'    grep "GERM" *06-0879* | grep "PASS" | grep -v "LOH" | grep "DB" | wc
'    grep "SOM" *02-0033*.vcf | grep "PASS" | grep -v "LOH" | grep -v "DB" | grep "MIS" | wc
'    awk '$1 == "TCGA-02-0033" {sum += $18} END {print sum}' cmpStrictLoose.tab
'
'    aradenba@cruncher /inside/grotto/users/aradenba/data/hg19/ucec/maf  $ grep "SNP" *.maf | wc
'    38372 1902304 19236705
'    aradenba@cruncher /inside/grotto/users/aradenba/data/hg19/ucec/maf  $ awk '$10 == "SNP"' *.maf | wc
'    38372 1902304 18584381
'    aradenba@cruncher /inside/grotto/users/aradenba/data/hg19/ucec/maf  $ awk '$10 == "SNP" && $3 == "genome.wustl.edu;ucsc.edu;broad.mit.edu" {print}' *.maf | wc
'    23780 1192627 11847606
'    aradenba@cruncher /inside/grotto/users/aradenba/data/hg19/ucec/bambam  $ grep "PASS" *.vcf | grep "SS=Somatic" | grep "SNP" | wc
'    32948  362428 6950937
'    aradenba@cruncher /inside/grotto/users/aradenba/data/hg19/ucec/maf  $ awk '$10 == "SNP" {print}' *.maf | grep "ucsc" | wc
'    30653 1533191 14987837
'    aradenba@pk ~  $ awk '{sum += $12} END {print sum}' ucecRadiaMafCmp.log 
'    29019
'    aradenba@pk ~  $ awk '{sum += $13} END {print sum}' ucecRadiaMafCmp.log 
'    38372
'    aradenba@pk ~  $ awk '{sum += $14} END {print sum}' ucecRadiaMafCmp.log 
'    23258
'    aradenba@pk ~  $ egrep "^TCGA" *.sh.o* > ucecRadiaMafCmp2.log
'    aradenba@pk ~  $ grep "no radia call" *.sh.o* | wc
'    3301  173374 1685325
'    aradenba@pk ~  $ grep "new radia call" *.sh.o* | wc
'    5757   80598 1750731
'    aradenba@pk ~  $ grep "found call" *.sh.o* | wc
'    23258 1210336 12273654
'    aradenba@pk ~  $ grep "found but no radia pass" *.sh.o* | wc
'    11164  178624 3819416
'    aradenba@pk ~  $ grep "overlap but not same type" *.sh.o* | wc
'    937   14992  323495
'    aradenba@cruncher ~  $ grep "overlap but not same type" *.sh.o* | grep "PASS" | grep "GERM" | wc
'    31     496   10028
'    aradenba@cruncher ~  $ grep "overlap but not same type" *.sh.o* | grep "PASS" | grep "LOH" | wc
'    72    1152   23264
'    aradenba@cruncher ~  $ grep "overlap but not same type" *.sh.o* | grep "LOH" | wc
'    180    2880   62046
'    aradenba@cruncher ~  $ grep "overlap but not same type" *.sh.o* | grep "GERM" | wc
'    757   12112  261449
'    aradenba@cruncher ~  $ grep "found but no radia pass" *.sh.o* | wc
'    11164  178624 3819416
'    aradenba@cruncher ~  $ grep "found but no radia pass" *.sh.o* | grep "blck" | wc
'    2224   35584  730420
'    aradenba@cruncher ~  $ grep "found but no radia pass" *.sh.o* | grep "tb" | wc
'    2150   34400  709083
'    aradenba@pk ~  $ grep "filter:" *.sh.o* | grep "blck" | awk '{sum += $3} END {print sum}'
'    2224
'    aradenba@pk ~  $ grep "filter:" *.sh.o* | grep "multi" | awk '{sum += $3} END {print sum}'
'    5804
'    aradenba@pk ~  $ grep "filter:" *.sh.o* | grep "bias" | awk '{sum += $3} END {print sum}'
'    1584
'    aradenba@pk ~  $ grep "filter:" *.sh.o* | grep "mntb" | awk '{sum += $3} END {print sum}'
'    984
'    aradenba@pk ~  $ grep "filter:" *.sh.o* | grep "mnab" | awk '{sum += $3} END {print sum}'
'    1569
'    aradenba@pk ~  $ grep "filter:" *.sh.o* | grep "mnap" | awk '{sum += $3} END {print sum}'
'    3779
'    aradenba@pk ~  $ grep "filter:" *.sh.o* | grep "bias" | awk '{sum += $3} END {print sum}'
'    1584
'    aradenba@pk ~  $ grep "filter:" *.sh.o* | grep "mxerr" | awk '{sum += $3} END {print sum}'
'    7531
'    aradenba@cruncher ~  $ grep "found call" *.sh.o* | awk '$5 == "genome.wustl.edu;ucsc.edu;broad.mit.edu" {print}' | wc
'    18443  962037 9862534
'    aradenba@cruncher /inside/grotto/users/aradenba/data/hg19/ucec/radia/filters/patientVCFs  $ grep "PASS" *.vcf | grep "SOM" | wc
'    29019  319209 8220989
'    
'    grep "found but no radia pass" *.sh.o* | awk '$19 == "genome.wustl.edu;ucsc.edu;broad.mit.edu" {print}' | grep "mxerr" | less
'
'    grep "no radia call" *.sh.o* | awk '$6 == "ucsc.edu" {print}' | less
'    grep "no radia call" *.sh.o* | awk '$6 == "ucsc.edu" {print$8"\t"$9"\t"$17"\t"$18}' | less
'    aradenba@pk ~  $ grep "no radia call" *.sh.o* | awk '$6 == "genome.wustl.edu;ucsc.edu;broad.mit.edu" {print "samtools mpileup -f /inside/depot/fa/GRCh37-lite.fa -r "$8":"$9"-"$9" /inside/home/singer/data/ucecfreeze/bams1/"$18"*.bam\nsamtools mpileup -f /inside/depot/fa/GRCh37-lite.fa -r "$8":"$9"-"$9" /inside/home/singer/data/ucecfreeze/bams1/"$19"*.bam\nsamtools mpileup -f /inside/depot/fa/GRCh37-lite.fa -r "$8":"$9"-"$9" /inside/home/singer/data/ucecfreeze/bams2/"$18"*.bam\nsamtools mpileup -f /inside/depot/fa/GRCh37-lite.fa -r "$8":"$9"-"$9" /inside/home/singer/data/ucecfreeze/bams2/"$19"*.bam\n"}' | less
'
'    Create an inspectPileups list
'    grep "no radia call" *.sh.o* | awk -F "\t" '$3 == "genome.wustl.edu;ucsc.edu;broad.mit.edu" {print $16"\t"$5"\t"$6"\t"$7"\t/inside/home/singer/data/ucecfreeze/bams1/"$16"*.bam\t/inside/depot/fa/GRCh37-lite.fa\t0\t0\tFalse\n"$17"\t"$5"\t"$6"\t"$7"\t/inside/home/singer/data/ucecfreeze/bams1/"$17"*.bam\t/inside/depot/fa/GRCh37-lite.fa\t0\t0\tFalse"}'
'
'    python2.7 createJobList.py filterAmbiguous /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/ /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/ambiguous/ /inside/home/aradenba/rnaEditing/cancer/benchmark/jobLists/makeAmbCallsJob.list -i /inside/home/aradenba/rnaEditing/cancer/benchmark/idLists/ids.txt
'
'''


def get_ids(anIdFilename):
    '''
    '  Return all of the Ids found in the given file.
    '  The file is formatted with one line of comma separated Ids.
    '  Returns a set of Ids with the format "TCGA-XX-YYYY" where X 
    '  is a TCGA cancer type such as "AB" and Y is alphanumeric.  
    '  For example, "TCGA-AB-1234".
    '''
    
    i_fileHandler = open(anIdFilename, "r")
    ids = list()
    
    for line in i_fileHandler:
        line = line.rstrip("\r\n")
        ids += line.split(",")
    
    i_fileHandler.close()
    return ids


def get_chrom_size(aChrom, anInputStream, anIsDebug):
    '''
    ' This function reads from a .tab file that contains the sizes for each chromosome.
    ' The file has two columns:  the first column has the chromosome identifier, the
    ' second column has the size of the chromosome.  The columns are separated by tabs.
    ' Here is an example:
    ' chr1    249250621
    ' chr2    243199373
    ' ...       ...
    '
    ' aChrom: The chrom size to return
    ' anInputStream: The input stream for the .sam file
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
     
    for line in anInputStream:
          
        # if it is an empty line, then just continue
        if (line.isspace() or line.startswith("#")):
            continue;

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        #if (anIsDebug):
        #    print >> sys.stderr, line    
        
        # split the line on the tab
        splitLine = re.split("\t", line)

        # the coordinate is the second element
        chr = splitLine[0]
        size = int(splitLine[1])
        
        # the chromosome from the file has the following format "chr1"
        # the aChrom variable only has the integer 1
        if (chr == "chr" + str(aChrom)):
            return size
    return


def get_dir_list(aDirFilename):
    '''
    '  Return a list of directories that are found in the directory file that is given. 
    '''
    dirList = list()
    fileHandler = open(aDirFilename, "r")
    for line in fileHandler:
        line = line.rstrip("\r\n")
        dirList.append(line)   
    fileHandler.close()
    return dirList


def main():
    '''
    ' This code is divided into 3 main sections: 
    '    1) Initial error-handling and variable initiation
    '    2) Job lists that are created by only looping over all chromosomes.
    '    3) Job lists that are created by looping over all chromosomes and all ids.
    '
    ' This is the start of section 1 which takes care of initial error-handling and variable initiation:
    '''
    i_debug = True
    
    if (i_debug):
        print >> sys.stderr, len(sys.argv), sys.argv
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(5,22,2)
    i_argLength = str(len(sys.argv))
    
    if (i_argLength in i_possibleArgLengths):
        print "usage: python2.7 createJobList.py formattingOption inputDir(s) outputDir(s) job.list -i id.lists -f outputFormat -c overChain -a hgAssembly -g hg19FastaFile -s chromSizesFile -q mappingQual -Q baseQual"
        sys.exit(1)
        
    # add the optional parameters    
    i_cmdLineParser = OptionParser()
    i_cmdLineParser.add_option("-i", "--idFilename", help="a file with comma separated Ids", dest="idFilename")
    i_cmdLineParser.add_option("-c", "--hg18ToHg19OverChain", help="the over chain file for hg18 to hg19", dest="overChainFilename")
    i_cmdLineParser.add_option("-g", "--hg19FastaFilename", help="the fasta file for hg19", dest="hg19FastaFilename")
    i_cmdLineParser.add_option("-a", "--hgAssembly", help="hg18 or hg19", dest="hgAssembly")
    i_cmdLineParser.add_option("-s", "--chromSizeFilename", help="the file with the chrom sizes", dest="chromSizesFilename")
    i_cmdLineParser.add_option("-f", "--rnaEditingOutputFormat", help="the format that the RNA-Editing events should be output in {RAD, BED, GFF}", dest="outputFormat", default="RAD")
    i_cmdLineParser.add_option("-q", "--minMappingQual", help="minimum mapping quality", dest="minMappingQual", default="20")
    i_cmdLineParser.add_option("-Q", "--minBaseQual", help="minimum base quality", dest="minBaseQual", default="20")
    
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_formattingOption = str(i_cmdLineArgs[0])
    i_inputDir = str(i_cmdLineArgs[1])
    i_outputDir = str(i_cmdLineArgs[2])
    i_joblistFilename = str(i_cmdLineArgs[3])
    
    # get the optional params with default values
    i_outputFormat = str(i_cmdLineOptions.outputFormat)
    i_minMappingQual = str(i_cmdLineOptions.minMappingQual)
    i_minBaseQual = str(i_cmdLineOptions.minBaseQual)
    
    # initialize some variables to None
    # this is needed for error-handling later
    i_idFilename = None
    i_overChainFilename = None
    i_hg19FastaFilename = None
    i_chromSizesFilename = None
    i_hgAssembly = None
    
    # try to get any optional parameters with no defaults
    if (i_cmdLineOptions.idFilename != None):
        i_idFilename = str(i_cmdLineOptions.idFilename)
        
    if (i_cmdLineOptions.overChainFilename != None):
        i_overChainFilename = str(i_cmdLineOptions.overChainFilename)
        
    if (i_cmdLineOptions.hg19FastaFilename != None):
        i_hg19FastaFilename = str(i_cmdLineOptions.hg19FastaFilename)
    
    if (i_cmdLineOptions.chromSizesFilename != None):
        i_chromSizesFilename = str(i_cmdLineOptions.chromSizesFilename)        
    
    if (i_cmdLineOptions.hgAssembly != None):
        i_hgAssembly = str(i_cmdLineOptions.hgAssembly)   
        
    # output some debug info
    if (i_debug):
        print >> sys.stderr, "i_formattingOption =", i_formattingOption
        print >> sys.stderr, "i_idFilename =", i_idFilename
        print >> sys.stderr, "i_overChainFilename =", i_overChainFilename
        print >> sys.stderr, "i_hg19FastaFilename =", i_hg19FastaFilename
        print >> sys.stderr, "i_hgAssembly =", i_hgAssembly
        print >> sys.stderr, "i_chromSizesFilename =", i_chromSizesFilename
        print >> sys.stderr, "i_inputDir =", i_inputDir
        print >> sys.stderr, "i_outputDir =", i_outputDir
        print >> sys.stderr, "i_joblistFilename =", i_joblistFilename
        print >> sys.stderr, "i_minMappingQual =", i_minMappingQual
        print >> sys.stderr, "i_minBaseQual =", i_minBaseQual
    
    # some of the jobs need multiple input and/or output dirs
    i_inputDirs = i_inputDir.split(",")
    i_outputDirs = i_outputDir.split(",")
    joblistDir = os.path.split(i_joblistFilename)[0]
    i_qsubDir = os.path.join(joblistDir, "qsub/")
    i_dirList = i_inputDirs + i_outputDirs + [i_qsubDir]
    
    # gather all the files to be read from
    i_readFilenameList = list()
    if (i_idFilename != None):
        # some of the jobs need multiple Id lists
        i_idFilenames = i_idFilename.split(",")
        i_readFilenameList = i_idFilenames
    if (i_overChainFilename != None):
        i_readFilenameList += [i_overChainFilename]
    if (i_hg19FastaFilename != None):
        i_readFilenameList += [i_hg19FastaFilename]    
    if (i_chromSizesFilename != None):
        i_readFilenameList += [i_chromSizesFilename]
    
    # we're always creating a job list
    i_writeFilenameList = [i_joblistFilename]

    if (i_debug):
        print >> sys.stderr, "i_readFilenameList =", i_readFilenameList
        print >> sys.stderr, "i_writeFilenameList =", i_writeFilenameList
        print >> sys.stderr, "i_dirList =", i_dirList

    # check to make sure all dirs and files exist
    if (not rnaEditingUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)
    
    # more specific error-handling
    if (i_idFilename == None and i_formattingOption != "symbolicLinks"  and i_formattingOption != "extractTcgaTargets" and i_formattingOption != "rnaDnaDiffs" and i_formattingOption != "runNullMutClass" and i_formattingOption != "mutNullModel" and i_formattingOption != "uniqRefGeneTx" and i_formattingOption != "uniqRefGeneExons" and i_formattingOption != "flankBlacklists" and i_formattingOption != "extractBlacklists" and i_formattingOption != "downloadCheungDNA" and i_formattingOption != "generateEventStats" and i_formattingOption != "generateFilterStats" and i_formattingOption != "sortDbSnp" and i_formattingOption != "getRefGeneTx" and i_formattingOption != "getRefGeneCDS" and i_formattingOption != "getRefGeneExons" and i_formattingOption != "sortRefGeneTx" and i_formattingOption != "sortRefGeneCDS" and i_formattingOption != "sortRefGeneExons" and i_formattingOption != "prepRefGeneExons" and i_formattingOption != "uniqRefGeneCDS" and i_formattingOption != "flankBlacklists" and i_formattingOption != "downloadDbSnp"):
        print >> sys.stderr, "You must specify an Id filename for the", i_formattingOption, "command."
        sys.exit(1)
    elif (i_formattingOption == "liftOver" and i_overChainFilename == None):
        print >> sys.stderr, "You must specify an hg18 to hg19 over chain file when doing a liftOver."
        sys.exit(1)
    elif (i_formattingOption == "rnaEditing" or i_formattingOption == "rnaDnaDiffs" and (i_hg19FastaFilename == None or i_chromSizesFilename == None)):
        print >> sys.stderr, "You must specify a fasta file and a file with the chrom sizes when calculating RNA-editing events."
        sys.exit(1) 
    elif (i_formattingOption == "flankBlacklists" and i_chromSizesFilename == None):
        print >> sys.stderr, "You must specify a file with the chrom sizes when flanking the blacklist regions."
        sys.exit(1)
    elif (i_formattingOption == "vcf2Excel" and i_formattingOption == "filterTranscripts" and i_formattingOption == "filterCDS" and i_formattingOption == "filterExons" and i_formattingOption == "getRefGeneTx" and i_formattingOption == "sortRefGeneTx" and i_formattingOption == "getRefGeneCDS" and i_formattingOption == "sortRefGeneCDS" and i_formattingOption == "uniqRefGeneCDS" and i_formattingOption == "getRefGeneExons" and i_formattingOption == "prepRefGeneExons" and i_formattingOption == "sortRefGeneExons" and i_formattingOption == "runMutClass" and (i_hgAssembly == None)):
        print >> sys.stderr, "You must specify either the hg18 or hg19 assembly using the -a option."
        sys.exit(1)
            
    # create the job.list
    i_fileHandler = open(i_joblistFilename, "w")
    i_fileHandler.write("> formattingOption: " + i_formattingOption + "\n")
   
    # create the chromId list
    
    #i_chromIdsOne = range(2, 6)
    #i_chromIdsTwo = range(11, 23)
    #i_chromIds = i_chromIdsOne + i_chromIdsTwo
    i_chromIds = range(9, 13)
    i_chromIds = map(str, i_chromIds)
    #i_chromIds.append("X")
    #i_chromIds.append("Y")
    #i_chromIds.append("M")
    
    #i_chromIds = ["M"]
    
    i_numOfJobsCounter = 0
    i_numOfExpectedJobs = 1
    if (i_formattingOption == "generateFilterStats"):
        # python2.7 generateFilterStats.py /inside/home/aradenba/rnaEditing/cancer/idLists/wgOverlappingIds.txt /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/filters/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/
        i_numOfJobsCounter += 1
        script = "/inside/home/aradenba/rnaEditing/scripts/generateFilterStats.py"
        i_fileHandler.write("python2.7 " + script + " " + i_idFilename + " " + i_inputDirs[0] + " " + i_inputDirs[1] + " " + i_outputDir + "\n")
    elif (i_formattingOption == "generateEventStats"):
        # python2.7 generateEventStats.py /inside/home/aradenba/rnaEditing/cancer/idLists/wgOverlappingIds.txt /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/v2/transcriptReadSupport/ /inside/home/aradenba/rnaEditing/cancer/aml/stats/
        i_numOfJobsCounter += 1
        script = "/inside/home/aradenba/rnaEditing/scripts/generateEventStats.py"
        i_fileHandler.write("python2.7 " + script + " " + i_idFilename + " " + i_inputDir + " " + i_outputDir + "\n")
    elif (i_formattingOption == "sortBamsCheung"): 
        # get Ids
        i_ids = get_ids(i_idFilename)
        i_numOfExpectedJobs = len(i_ids)

        for id in i_ids:
            unsortedBam = i_inputDir + "GSE25840_*_GM" + id + "*.bam"
            sortedBam = i_outputDir + "GSE25840_GM" + id + "_gencode_spliced_sorted"
            
            # keep track of the number of jobs that are generated for validation
            i_numOfJobsCounter += 1
            i_fileHandler.write("samtools sort " + unsortedBam + " " + sortedBam + "\n") 
        
    elif (i_formattingOption == "symbolicLinks"): 
        
        # get all the input directories to search through
        # the user can specify a comma-separated list which is stored in inputDirs
        # or they specify the filename, and we extract the list of dirs here
        if (i_inputDirs == None):
            i_inputFilename = i_inputDir + "allDirList.txt"
            i_inputDirs = get_dir_list(i_inputFilename)
            
        cancerTypesInput = "A1,A2,A7,A8,AN,AO,AQ,AR,B6,BH,C8,CG,D8,E2"
        cancerTypes = cancerTypesInput.split(",")
        
        skip = ['ncbi_enc', 'md5']
               
        i_numOfExpectedJobs = 0
        for dir in i_inputDirs:
            for file in os.listdir(dir):
                #print "file", file
                
                # if there is no bam in the filename, then skip this one
                if file.find('bam') == -1:
                    continue
                
                # we are only looking for .bam and .bai files,
                # so skip any other files that have bam in them
                toskip = False
                for s in skip:
                    if file.find(s) != -1:
                        toskip = True
                        break
                if toskip:
                    continue
        
                # look for the barcode
                if 'TCGA' in file:
                    i = file.index('TCGA')
                    #pname = file[i:i+12]             # TCGA-AB-2973-03A-01D-0739-09_whole.bam
                    ctype = file[i+5:i+7]            # cancer type = AB  
                    #tid   = file[i+13:i+15]          # sample tissue type = 03
                    sid   = file[i+19]               # D=DNA, R=RNA, T=Total_RNA, W/G/X=Whole_Genome
                else:
                    continue
                #print pname, ctype, tid, sid
                
                # if this isn't the study we are interested in,
                # then skip it
                if ctype not in cancerTypes:
                    continue
                
                # find the correct sub-directory
                subDir = ""
                if (sid == "W" or file.find("_whole") != -1):
                    subDir = "wg"
                elif (sid == "D" or file.find("_exome") != -1 or file.find("_capture") != -1):
                    subDir = "exome"
                elif (sid == "R"):
                    subDir = "rna"
                
                # copy over the .bai files
                if (file.find("bam") != -1 and file.find("bai") != -1):
                    #cp /inside/depot2/grade3_astrocytoma/bams/"$1"*.bam.bai /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/exome/"$2".bam.bai
                    i_fileHandler.write("cp " + os.path.join(dir, file) + os.path.join(i_outputDir, subDir, file) + "\n")
                    # keep track of the number of jobs that are generated for validation
                    i_numOfJobsCounter += 1
                    i_numOfExpectedJobs += 1
                # symbolic link the .bam files
                elif (file.find("bam") != -1):
                    #ln -s /inside/depot2/grade3_astrocytoma/bams/"$1"*.bam /inside/grotto/users/aradenba/data/hg18/grade3_astrocytoma/exome/"$2".bam
                    i_fileHandler.write("ln -s " + os.path.join(dir, file) + os.path.join(i_outputDir, subDir, file) + "\n")
                    # keep track of the number of jobs that are generated for validation
                    i_numOfJobsCounter += 1
                    i_numOfExpectedJobs += 1

    elif (i_formattingOption == "samtoolsIndex"): 
        # get Ids
        i_ids = get_ids(i_idFilename)
        i_numOfExpectedJobs = len(i_ids)
        
        for id in i_ids:
            #sortedBam = i_inputDir + "GSE25840_GM" + id + "*_sorted.bam"
            sortedBam = i_inputDir + id + "*.sorted.bam"
            
            # keep track of the number of jobs that are generated for validation
            i_numOfJobsCounter += 1
            i_fileHandler.write("samtools index " + sortedBam + "\n") 
    elif (i_formattingOption == "samtoolsSelection"):
        # get Ids
        i_ids = get_ids(i_idFilename)
        i_numOfExpectedJobs = len(i_ids) * 3
        
        for id in i_ids:
            #dnaNormalFilename = i_inputDirs[0] + id + "-10*.bam"
            #dnaTumorFilename = i_inputDirs[0] + id + "-01*.bam"
            #rnaTumorFilename = i_inputDirs[1] + id + "*.bam"
            
            dnaNormalFilename = i_inputDirs[0] + id + "N*.bam"
            dnaTumorFilename = i_inputDirs[0] + id + "T*.bam"
            rnaTumorFilename = i_inputDirs[1] + id + "T*.bam"
            
            # keep track of the number of jobs that are generated for validation
            i_numOfJobsCounter += 3
            # samtools mpileup -f /inside/depot/fa/hsap36.1-hg18.fasta -Q 25 -q 20 -r chr3:0-1000000 /inside/depot/cheung/rna/GSE25840_GM06985_gencode_spliced_sorted.bam
            # /inside/depot/fa/hsap36.1-hg18.fasta
            # /inside/depot/fa/all_sequences.fasta
            # /inside/depot/fa/Homo_sapiens_assembly18.fasta
            # /inside/depot/fa/NCBI36-HG18_Broad_variant.fasta
            i_fileHandler.write("samtools mpileup -f /inside/depot/fa/all_sequences.fasta -Q 20 -q 20 -r 22:49366846-49366849 " + dnaNormalFilename + " > " + id + "_dnaNormal_samtoolsTest.log\n")
            i_fileHandler.write("samtools mpileup -f /inside/depot/fa/all_sequences.fasta -Q 20 -q 20 -r 22:49366846-49366849 " + dnaTumorFilename + " > " + id + "_dnaTumor_samtoolsTest.log\n")
            i_fileHandler.write("samtools mpileup -f /inside/depot/fa/all_sequences.fasta -Q 20 -q 20 -r 22:49366846-49366849 " + rnaTumorFilename + " > " + id + "_rnaTumor_samtoolsTest.log\n")
    
    
    elif (i_formattingOption == "catChromVCFFiles"):
        # this is dumb
        i_ids = get_ids(i_idFilename)
        i_numOfExpectedJobs = len(i_ids)
        
        # cat chroms 2-22, X, and Y files together and remove all the headers, and save it in a temp file
        # then cat the chrom 1 file and the temp file and save it to a file labeled with just the id
        # remove the temp file
        for id in i_ids:
            catList = "cat "
            for chromId in i_chromIds:
                
                # skip chrom 1
                if (chromId == "1"):
                    continue;
                
                catList += i_inputDir + id + "_chr" + chromId + ".vcf "
                
            # keep track of the number of jobs that are generated for validation
            i_numOfJobsCounter += 1        
            patientVCF = i_outputDir + id + ".vcf"
            patientTmpVCF = i_outputDir + id + "_tmp.vcf"
            
            i_fileHandler.write("qsub " + i_qsubDir + id + ".sh\n")
                        
            tempFileHandler = open(i_qsubDir + id + ".sh", "w")
            tempFileHandler.write("#!/bin/bash\n")
            #tempFileHandler.write("#$ -M aradenba@soe.ucsc.edu\n")
            #tempFileHandler.write("#$ -m be\n")
            tempFileHandler.write("source ~/.bashrc\n")
            # cat chroms 2-22, X, and Y files together and remove all the headers, and save it in a temp file
            tempFileHandler.write(catList + "| sed '/^#/d' > " + patientTmpVCF + "\n")
            # then cat the chrom 1 file and the temp file and save it to a file labeled with just the id
            tempFileHandler.write("cat " + i_inputDir + id + "_chr1.vcf " + patientTmpVCF + " > " + patientVCF + "\n")
            # remove the temp file
            tempFileHandler.write("rm " + patientTmpVCF + "\n")
            tempFileHandler.close()    
            
    elif (i_formattingOption == "extractMAF"):
        i_ids = get_ids(i_idFilename)
        i_numOfExpectedJobs = len(i_ids)
        
        for id in i_ids:
                
            # keep track of the number of jobs that are generated for validation
            i_numOfJobsCounter += 1        
            #mafFile = i_inputDir + "ucec_3_center_222_cases_vars_within_capture_targets.sorted.maf"
            #mafFile = i_inputDir + "all222UCEC.sorted.maf"
            mafFile = i_inputDir + "genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf"
            outputFile = i_outputDir + id + ".maf"
            
            #i_fileHandler.write("grep \"" + id + "\" " + mafFile + " | awk '$5 != \"MT\" && $5 != \"M\" {print}' > " + outputFile + "\n")
            #i_fileHandler.write("grep \"" + id + "\" " + mafFile + " | awk '$10 == \"SNP\" && $5 != \"MT\" && $5 != \"M\" {print}' > " + outputFile + "\n")
            i_fileHandler.write("grep \"" + id + "\" " + mafFile + " | grep -v \"11A-\" | awk '$10 == \"SNP\" && $5 != \"MT\" && $5 != \"M\" {print}' > " + outputFile + "\n")
    elif (i_formattingOption == "radiaCompareMaf"):
        
        script = "/inside/home/aradenba/rnaEditing/scripts/radiaCompare.py"
        i_ids = get_ids(i_idFilename)
        i_numOfExpectedJobs = len(i_ids)
        
        for id in i_ids:
            radFilename = i_inputDirs[0] + id + ".vcf"
            cmpFilename = i_inputDirs[1] + id + ".maf"
            outputFilename = i_outputDirs[0] + "cmpStrictLoose.tab"
            overlapFilename = i_outputDirs[1] + id + "_ovlp.vcf"
            nonOverlapFilename = i_outputDirs[2] + id + "_nonovlp.vcf"
            
            # keep track of the number of jobs that are generated for validation
            i_numOfJobsCounter += 1
            i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + radFilename + " " + cmpFilename + " -c \"GERM=Germline;SOM=Somatic\" -s " + outputFilename + " -o " + overlapFilename + " -n " + nonOverlapFilename +"\n")
            
            '''
            i_fileHandler.write("qsub " + i_qsubDir + id + ".sh\n")
            
            tempFileHandler = open(i_qsubDir + id + ".sh", "w")
            tempFileHandler.write("#!/bin/bash\n")
            #tempFileHandler.write("#$ -M aradenba@soe.ucsc.edu\n")
            #tempFileHandler.write("#$ -m be\n")
            tempFileHandler.write("source ~/.bashrc\n")
            tempFileHandler.write("python2.7 " + script + " " + id + " Z " + radFilename + " " + cmpFilename + " -c \"SOM=Somatic\" -s " + outputFilename + " -o " + overlapFilename + " -n " + nonOverlapFilename +"\n")
            tempFileHandler.close()
            '''    
    # Select ptpn6 region out of hg19 bambam files 
    elif (i_formattingOption == "ptpn6Rad"):
        # get Ids
        i_ids = get_ids(i_idFilename)
        i_numOfExpectedJobs = len(i_ids)
        
        for id in i_ids:
            radInputFilename = i_inputDir + id + "_chr12.rad"
            radOutputFile = i_outputDir + id + "_chr12_ptpn6.rad"
            
            # keep track of the number of jobs that are generated for validation
            i_numOfJobsCounter += 1
            i_fileHandler.write("awk '$2 >= 7055740 && $3 <= 7070479 {print$1\"\\t\"$2\"\\t\"$3\"\\t\"$4}' " + radInputFilename + " > " + radOutputFile + "\n")
            # awk '$1 == "TCGA-02-0033" {sum += $18} END {print sum}' cmpStrictLoose.tab
    else:
        # loop over each chromosome 
        i_numOfJobsCounter = 0
        for chromId in i_chromIds:
            
            '''
            ' This code is divided into 4 main sections: 
            '    1) Initial error-handling and variable initiation
            '    2) Job lists that are require no looping or customized looping
            '    2) Job lists that are created by only looping over all chromosomes.
            '    3) Job lists that are created by looping over all chromosomes and all ids.
            '
            ' This is the beginning of section 3 which loops over all chromosomes but not all ids:
            '''
    
            # initialize expected jobs
            i_numOfExpectedJobs = len(i_chromIds)
    
            # these formatting options loop over each chrom but not over each sample Id
            if (i_formattingOption == "getRefGeneTx"):
                i_numOfJobsCounter += 1
                outputFile = i_outputDir + i_hgAssembly + "_chr" + chromId + "_tx.bed"
                i_fileHandler.write("echo \"select chrom,txStart,txEnd,name2,strand from refGene where chrom = \'chr" + chromId + "\'\" | hgsql " + i_hgAssembly + " | uniq > " + outputFile + "\n")
            elif (i_formattingOption == "getRefGeneCDS"):
                i_numOfJobsCounter += 1
                outputFile = i_outputDir + i_hgAssembly + "_chr" + chromId + "_cds.bed"
                i_fileHandler.write("echo \"select chrom,cdsStart,cdsEnd,name2,strand from refGene where chrom = \'chr" + chromId + "\'\" | hgsql " + i_hgAssembly + " | uniq > " + outputFile + "\n")
            elif (i_formattingOption == "getRefGeneExons"):
                i_numOfJobsCounter += 1
                outputFile = i_outputDir + i_hgAssembly + "_chr" + chromId + "_exons.bed"
                i_fileHandler.write("echo \"select chrom,exonStarts,exonEnds,name2,strand from refGene where chrom = \'chr" + chromId + "\'\" | hgsql " + i_hgAssembly + " | uniq > " + outputFile + "\n")
            elif (i_formattingOption == "prepRefGeneExons"):
                i_numOfJobsCounter += 1
                inputFile = i_inputDir + i_hgAssembly + "_chr" + chromId + "_exons.bed"
                outputFile = i_outputDir + i_hgAssembly + "_chr" + chromId + "_exons_extracted.bed"
                script = "/inside/home/aradenba/rnaEditing/scripts/prepExons.py"
                i_fileHandler.write("python2.7 " + script + " " + inputFile + " " + outputFile + "\n")    
            elif (i_formattingOption == "sortRefGeneTx"):
                i_numOfJobsCounter += 1
                inputFile = i_inputDir + i_hgAssembly + "_chr" + chromId + "_tx.bed"
                outputFile = i_outputDir + i_hgAssembly + "_chr" + chromId + "_tx_sorted.bed"
                # sort numerically on columns 1 and 2 with tab as separator
                # .bed file has columns: 0=chr, 1=startCoordinate, 2=stopCoordinate, etc.           
                i_fileHandler.write("sort -t $'\\t' +1n -3 " + inputFile + " > " + outputFile + "\n")
            elif (i_formattingOption == "sortRefGeneCDS"):
                i_numOfJobsCounter += 1
                inputFile = i_inputDir + i_hgAssembly + "_chr" + chromId + "_cds.bed"
                outputFile = i_outputDir + i_hgAssembly + "_chr" + chromId + "_cds_sorted.bed"
                # sort numerically on columns 1 and 2 with tab as separator
                # .bed file has columns: 0=chr, 1=startCoordinate, 2=stopCoordinate, etc.           
                i_fileHandler.write("sort -t $'\\t' +1n -3 " + inputFile + " > " + outputFile + "\n")
            elif (i_formattingOption == "sortRefGeneExons"):
                i_numOfJobsCounter += 1
                inputFile = i_inputDir + i_hgAssembly + "_chr" + chromId + "_exons_extracted.bed"
                outputFile = i_outputDir + i_hgAssembly + "_chr" + chromId + "_exons_sorted.bed"
                # sort numerically on columns 1 and 2 with tab as separator
                # .bed file has columns: 0=chr, 1=startCoordinate, 2=stopCoordinate, etc.           
                i_fileHandler.write("sort -t $'\\t' +1n -3 " + inputFile + " > " + outputFile + "\n")
            elif (i_formattingOption == "uniqRefGeneTx"):
                i_numOfJobsCounter += 1
                inputFile = i_inputDir + i_hgAssembly + "_chr" + chromId + "_tx_sorted.bed"
                outputFile = i_outputDir + i_hgAssembly + "_chr" + chromId + "_tx_sorted_uniq.bed"           
                i_fileHandler.write("uniq " + inputFile + " > " + outputFile + "\n")
            elif (i_formattingOption == "uniqRefGeneCDS"):
                i_numOfJobsCounter += 1
                inputFile = i_inputDir + i_hgAssembly + "_chr" + chromId + "_cds_sorted.bed"
                outputFile = i_outputDir + i_hgAssembly + "_chr" + chromId + "_cds_sorted_uniq.bed"           
                i_fileHandler.write("uniq " + inputFile + " > " + outputFile + "\n")
            elif (i_formattingOption == "uniqRefGeneExons"):
                i_numOfJobsCounter += 1
                inputFile = i_inputDir + i_hgAssembly + "_chr" + chromId + "_exons_sorted.bed"
                outputFile = i_outputDir + i_hgAssembly + "_chr" + chromId + "_exons_sorted_uniq.bed"           
                i_fileHandler.write("uniq " + inputFile + " > " + outputFile + "\n")
            elif (i_formattingOption == "mutNullModel"):
                i_numOfJobsCounter += 1
                # hg19_chr17_tx_sorted_uniq.bed
                # python2.7 mutationNullModel.py chr txFile outputFilename
                script = "/inside/home/aradenba/rnaEditing/scripts/mutationNullModel.py"
                inputFile = i_inputDir + "hg19_chr" + chromId + "_tx_sorted_uniq.bed"
                outputFile = i_outputDir + "chr" + chromId + ".bed"           
                i_fileHandler.write("python2.7 " + script + " " + chromId + " " + inputFile + " " + outputFile + "\n")
            elif (i_formattingOption == "runNullMutClass"):
                script = "mutClass -lazyLoading -skipRegulatory -oneBased " + i_hgAssembly
                inputFile = i_inputDir + "chr" + chromId + ".bed"
                outputFile = i_outputDir + "chr" + chromId + ".bed"
                i_fileHandler.write(script + " " + inputFile + " > " + outputFile + "\n")
            elif (i_formattingOption == "downloadDbSnp"):
                i_numOfJobsCounter += 1
                outputFile = i_outputDir + "chr" + chromId + "_unsorted.bed"
                i_fileHandler.write("echo \"select chrom,chromStart,chromEnd,name from snp135 where chrom = \'chr" + chromId + "\'\" | hgsql " + i_hgAssembly + " > " + outputFile + "\n")
            elif (i_formattingOption == "sortDbSnp"):
                i_numOfJobsCounter += 1
                inputFile = i_inputDir + "chr" + chromId + "_unsorted.bed"
                outputFile = i_outputDir + "chr" + chromId + ".bed"
                # sort numerically on columns 1 and 2 with tab as separator
                # .bed file has columns: 0=chr, 1=startCoordinate, 2=stopCoordinate, etc.           
                i_fileHandler.write("sort -t $'\\t' +1n -3 " + inputFile + " > " + outputFile + "\n")
            elif (i_formattingOption == "extractBlacklists"):
                i_numOfJobsCounter += 1
                #inputFile = i_inputDir + "original_blacklist.bed"
                #inputFile = i_inputDir + "covMask1kGPilotLowCovUnionDepth.bed"
                inputFile = i_inputDir + "covMask1kGPilotLowCovUnionMapQ.bed"
                outputFile = i_outputDir + "chr" + chromId + "_sorted.bed"
                i_fileHandler.write("egrep chr" + chromId + "[[:space:]] " + inputFile + " | sort -t $'\\t' +1n -3 > " + outputFile + "\n")
            elif (i_formattingOption == "extractTcgaTargets"):
                i_numOfJobsCounter += 1
                #inputFile = i_inputDir + "gaf_1.0_exon_coordinates_hg19.uniq.bed.sorted.merged"
                inputFile = i_inputDir + "washu_tier1_annotation_space_on_grch37_lite.bed"
                outputFile = i_outputDir + "chr" + chromId + "_sorted.bed"
                i_fileHandler.write("awk '$1 == \"" + chromId + "\" {print$1\"\\t\"$2\"\\t\"$3}' " + inputFile + " | sort -t $'\\t' +1n -3 > " + outputFile + "\n")    
            elif (i_formattingOption == "flankBlacklists"):
                i_numOfJobsCounter += 1
                # get the start and end coordinates
                i_chromSizeFileHandler = open(i_chromSizesFilename, "r")
                i_stopCoordinate = get_chrom_size(chromId, i_chromSizeFileHandler, i_debug)
                i_chromSizeFileHandler.close()
                inputFile = i_inputDir + "chr" + chromId + "_sorted.bed"
                outputFile = i_outputDir + "chr" + chromId + "_sorted_flanked.bed"
                i_fileHandler.write("awk '{if ($2 <= 200) $2 = 1; else $2=$2-200} {if ($3+200 > " + str(i_stopCoordinate) + ") $3 = " + str(i_stopCoordinate) + "; else $3=$3+200} {print $1\"\\t\"$2\"\\t\"$3}' " + inputFile + " > " + outputFile + "\n")    
            elif (i_formattingOption == "filterPatientSupport"):
                # python2.7 filterByPatientSupport.py 1 /inside/home/aradenba/rnaEditing/cancer/cheung/idLists/cheungIds.txt /inside/home/aradenba/rnaEditing/data/validation/rddSites/rdd_sites_expanded_hg18.bed /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/transcriptReadSupport/ /inside/grotto/users/aradenba/data/hg18/cheung/rnaEditing/noBaq/filters/transcriptPatientSupport/ -s 2
                script = "/inside/home/aradenba/rnaEditing/scripts/filterByPatientSupport.py"
                #rddFilename = "/inside/home/aradenba/rnaEditing/data/validation/rddSites/rdd_sites_expanded_hg18.bed"
                rddFilename = "/inside/home/aradenba/rnaEditing/data/validation/rddSites/rdd_sites_merged_hg19.bed"
                i_fileHandler.write("python2.7 " + script + " " + chromId + " " + i_idFilename + " " + rddFilename + " " + i_inputDir + " " + i_outputDir + " -s 2 \n")
            # Filter events and generate stat files for plotting 
            elif (i_formattingOption == "filterEvents"): 
                i_numOfJobsCounter += 1
                #python2.7 filterEvents.py /inside/home/aradenba/rnaEditing/cancer/aml/idLists/wgOverlappingIds.txt /inside/home/aradenba/rnaEditing/cancer/aml/idLists/exomeOverlappingIds.txt /inside/grotto/users/aradenba/data/hg19/dbSNP/ /inside/home/aradenba/rnaEditing/data/refGenes/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/ /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/extracts/ /inside/grotto/users/aradenba/data/hg19/aml/stats/  
        
                script = "/inside/home/aradenba/rnaEditing/scripts/filterEvents.py"
                wgIdFilename = i_idFilenames[0]
                exomeIdFilename = i_idFilenames[1]
                dbSnpDir = i_inputDirs[0]
                refGeneDir = i_inputDirs[1]
                inputDir = i_inputDirs[2] 
                outputDir = i_outputDirs[0] 
                outputStatsDir = i_outputDirs[1] 
                print >> sys.stderr, script, wgIdFilename, exomeIdFilename, dbSnpDir, refGeneDir, inputDir, outputDir, outputStatsDir
                i_fileHandler.write("python2.7 " + script + " " + wgIdFilename + " " + exomeIdFilename + " " + dbSnpDir + " " + refGeneDir + " " + inputDir + " " + outputDir + " " + outputStatsDir + "\n")
            elif (i_formattingOption == "rnaDnaDiffs"):
                i_numOfJobsCounter += 1
                # get the start and end coordinates
                # python2.7 rnaDNAPipeline.py SF4454 X -n /inside/grotto/astrocytoma/hg18Files/A02404_81MJ7ABXX_2.bam -t /inside/grotto/astrocytoma/hg18Files/A02397_81MH7ABXX_5.bam -r /inside/grotto/astrocytoma/hg18Files/A02406_819T0ABXX_2.bam -f /inside/depot/fa/all_sequences.fasta -c /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes_sorted.tab -Q 20 -q 20 -o /inside/grotto/astrocytoma/tripleBam/SF4454_chrX.rad
                script = "/inside/home/aradenba/rnaEditing/scripts/rnaDNAPipeline.py"
                normalDnaFile = i_inputDir + "A02404_81MJ7ABXX_2.bam"
                tumorDnaFile = i_inputDir + "A02397_81MH7ABXX_5.bam"
                rnaSeqFile = i_inputDir + "A02406_819T0ABXX_2.bam"
                outputFile = i_outputDir + "SF4454_chr" + chromId + ".rad"
                i_fileHandler.write("python2.7 " + script + " SF4454 " + chromId + " -n " + normalDnaFile + " -t " + tumorDnaFile + " -r " + rnaSeqFile + " -f " + i_hg19FastaFilename + " -c " + i_chromSizesFilename + " -Q 20 -q 20 -o " + outputFile + "\n")    
            else:
                
                # some formatting options need to do some things per chromosome and other things per chrom and id
                # before starting to loop over the ids, this section just loops over the chrom 
                # initialize the validation file
                if (i_formattingOption == "validateLiftOver"):
                    # the output file is only generated per chrom instead of per chrom and id, so it needs
                    # to be initialized only once here and then each validation step will append to it
                    outputFile = i_outputDir + "chr" + chromId + "_liftOver_validation.log"
                    outputFileHandler = open(outputFile, "w")
                    outputFileHandler.write("Chr\tId\tHG18_%\tHG19_%\n")
                    outputFileHandler.close()
                
                '''
                ' This code is divided into 3 main sections: 
                '    1) Initial error-handling and variable initiation
                '    2) Job lists that are created by only looping over all chromosomes.
                '    3) Job lists that are created by looping over all chromosomes and all ids.
                '
                ' This is the start of section 3 which loops over all chroms and all ids:
                '''
                
                # get Ids
                i_ids = get_ids(i_idFilename)
                
                # initialize expected jobs
                i_numOfExpectedJobs = len(i_chromIds) * len(i_ids)
                
                for id in i_ids:
                
                    # select MUT lines from the original BB files
                    if (i_formattingOption == "selectMUTs"):
                        # the whole genome files have names like this:  TCGA-AB-2968_D_whole_12.bb
                        # the exome files have names like this:  TCGA-AB-2820_chr12.bb or TCGA-AB-2820_W_capture_chr12.bb
                        # rarely, some TCGA AML Ids have 2 exome files with names: TCGA-AB-2820_chr12.bb and TCGA-AB-2820_W_capture_chr12.bb 
                        # we will use the ones with TCGA-AB-2820_chr20.bb when both files exist
                       
                        # for exome, look for these first: TCGA-AB-2820_chr20.bb 
                        if (os.path.isfile(i_inputDir + id + "_chr" + chromId + ".bb")): 
                            inputFile = i_inputDir + id + "_chr" + chromId + ".bb"
                        # for exome, look for these second: TCGA-AB-2820_W_capture_chr12.bb 
                        elif (os.path.isfile(i_inputDir + id + "_W_capture_chr" + chromId + ".bb")):
                            inputFile = i_inputDir + id + "_W_capture_chr" + chromId + ".bb"
                        # for whole genome, look for these: TCGA-AB-2968_D_whole_12.bb 
                        elif (os.path.isfile(i_inputDir + id + "_D_whole_" + chromId + ".bb")):	
                            inputFile = i_inputDir + id + "_D_whole_" + chromId + ".bb"         
                                               
                        # for cheung data from 454, look for these: NA12043.chrom5.LS454.ssaha2.CEU.exon_targetted.20100311_5.bb  
                        elif (os.path.isfile(i_inputDir + "NA" + id + ".chrom" + chromId + ".LS454.ssaha2.CEU.exon_targetted.20100311_" + chromId + ".bb")):    
                            inputFile = i_inputDir + "NA" + id + ".chrom" + chromId + ".LS454.ssaha2.CEU.exon_targetted.20100311_" + chromId + ".bb"
                        # for cheung data from Illumina, look for these: NA12043.chrom5.ILLUMINA.bwa.CEU.exon_targetted.20100311_5.bb  
                        elif (os.path.isfile(i_inputDir + "NA" + id + ".chrom" + chromId + ".ILLUMINA.bwa.CEU.exon_targetted.20100311_" + chromId + ".bb")):    
                            inputFile = i_inputDir + "NA" + id + ".chrom" + chromId + ".ILLUMINA.bwa.CEU.exon_targetted.20100311_" + chromId + ".bb"
                        
                        else:
                            inputFile = None
                            print >> sys.stderr, "Couldn't find a BamBam file for chrom", chromId, "and Id", id, "in the", i_inputDir, "directory." 
                            
                        if (inputFile != None):  
                            # keep track of the number of jobs that are generated for validation
                            i_numOfJobsCounter += 1
                              
                            # grep all the lines that have "MUT" in them
                            if (inputFile.find("LS454") != -1):
                                outputFile = i_outputDir + id + "_LS454_nwc_chr" + chromId + ".bb"
                            elif (inputFile.find("ILLUMINA") != -1):
                                outputFile = i_outputDir + id + "_ILLUMINA_nwc_chr" + chromId + ".bb"
                            else:
                                outputFile = i_outputDir + id + "_ncnwc_chr" + chromId + ".bb"
                                
                            #i_fileHandler.write("grep \"MUT\" " + inputFile + " > " + outputFile + "\n")
                            #i_fileHandler.write("grep \"USM\" " + inputFile + " > " + outputFile + "\n")
                            #i_fileHandler.write("grep \"NC\\|NWC\" " + inputFile + " > " + outputFile + "\n")
                            #i_fileHandler.write("grep \"NC\" " + inputFile + " > " + outputFile + "\n")
                            i_fileHandler.write("grep \"NWC\" " + inputFile + " > " + outputFile + "\n")
                                        
                    # get the DNA Reads
                    elif (i_formattingOption == "getDNAReads"):
                        # the whole genome files have names like this:  
                        #    TCGA-AB-2907-03A-01D-0739-09_whole.bam and TCGA-AB-2907-11A-01D-0739-09_whole.bam
                        # the exome files have names like this:
                        #    TCGA-AB-2804-03B-01W-0728-08.bam and TCGA-AB-2804-11B-01W-0728-08.bam
                        #    TCGA-AB-2805-03B-01W-0728-08_IlluminaGA-DNASeq_capture.bam and TCGA-AB-2805-11B-01W-0728-08_IlluminaGA-DNASeq_capture.bam
                        #    and more...
                        
                        #print "tumor", os.system("ls " + i_inputDirs[1] + id + "-03*.bam")
                        #print "normal", os.system("ls " + i_inputDirs[1] + id + "-11*.bam")
                        
                        #print "tumor", glob.glob(i_inputDirs[1] + id + "-03*.bam")
                        #print "normal", glob.glob(i_inputDirs[1] + id + "-11*.bam")
                        
                        # for exome: 
                        if (len(glob.glob(i_inputDirs[1] + id + "-03*.bam")) > 0): 
                            tumorDNAFile = glob.glob(i_inputDirs[1] + id + "-03*.bam")[0]
                        # for whole genome, look for tumor:  TCGA-AB-2907-03A-01D-0739-09_whole.bam
                        # for whole genome, look for normal:  TCGA-AB-2907-11A-01D-0739-09_whole.bam 
                        elif (os.path.isfile(i_inputDirs[1] + id + "-03A-01D-0739-09_whole.bam")):    
                            tumorDNAFile = i_inputDirs[1] + id + "-03A-01D-0739-09_whole.bam"   
                        else:
                            tumorDNAFile = None
                            print >> sys.stderr, "Couldn't find a tumor DNA Bam file for chrom", chromId, "and Id", id, "in the", i_inputDirs[1], "directory." 
                        
                        # for exome: 
                        if (len(glob.glob(i_inputDirs[1] + id + "-11*.bam")) > 0): 
                            normalDNAFile = glob.glob(i_inputDirs[1] + id + "-11*.bam")[0]
                        # for whole genome, look for normal:  TCGA-AB-2907-11A-01D-0739-09_whole.bam
                        elif (os.path.isfile(i_inputDirs[1] + id + "-11A-01D-0739-09_whole.bam")):    
                            normalDNAFile = i_inputDirs[1] + id + "-11A-01D-0739-09_whole.bam"         
                        else:
                            normalDNAFile = None
                            print >> sys.stderr, "Couldn't find a matched-normal DNA Bam file for chrom", chromId, "and Id", id, "in the", i_inputDirs[1], "directory." 
                            
                        if (tumorDNAFile != None and normalDNAFile != None):  
                            # keep track of the number of jobs that are generated for validation
                            i_numOfJobsCounter += 1
                            script = "/inside/home/aradenba/rnaEditing/scripts/getDNAReads.py"
                              
                            inputFile = i_inputDirs[0] + id + "_chr" + chromId + ".rad"
                            outputFile = i_outputDir + id + "_chr" + chromId + ".rad"
                            
                            i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + inputFile + " " + normalDNAFile + " " + tumorDNAFile + " " + outputFile + "\n")
                        
                    # convert .bb to .bed, minimalistically 
                    elif (i_formattingOption == "bbToBed"): 
                        #inputFile = i_inputDir + id + "_muts_chr" + chromId + ".bb"
                        #outputFile = i_outputDir + id + "_muts_chr" + chromId + ".bed"
                        
                        inputFile = i_inputDir + id + "_LS454_nwc_chr" + chromId + ".bb"
                        outputFile = i_outputDir + id + "_LS454_nwc_chr" + chromId + ".bed"
                        
                        #inputFile = i_inputDir + id + "_ILLUMINA_nwc_chr" + chromId + ".bb"
                        #outputFile = i_outputDir + id + "_ILLUMINA_nwc_chr" + chromId + ".bed"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        
                        # extract the 3rd and 4th columns from the .bb files which correspond to the hg18 start and stop coordinates
                        # create a new .bed file with the chrId, startCoordinate, and stopCoordinate, and for later validation:
                        # a string consisting of the start and stop coordinates concatenated with an underscore, and the original .bb line number
                        i_fileHandler.write("awk '{print \"chr" + chromId + "\\t\"$3\"\\t\"$4\"\\t\"$3\"_\"$4\"\\t\"FNR}' " + inputFile + " > " + outputFile + "\n")
                    
                    # convert .rad to .bed, minimalistically 
                    elif (i_formattingOption == "radToBed"): 
                        inputFile = i_inputDir + id + "_chr" + chromId + ".rad"
                        outputFile = i_outputDir + id + "_chr" + chromId + ".bed"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        
                        # extract the 2nd and 3rd columns from the .bb files which correspond to the hg19 start and stop coordinates
                        # create a new .bed file with the chrId, startCoordinate, and stopCoordinate, and for later validation:
                        # a string consisting of the start and stop coordinates concatenated with an underscore, and the original .rad line number
                        i_fileHandler.write("awk '{print \"chr" + chromId + "\\t\"$2\"\\t\"$3\"\\t\"$2\"_\"$3\"\\t\"FNR}' " + inputFile + " > " + outputFile + "\n")
                    
                    # do the liftover from hg18 .bed to hg19 .bed
                    elif (i_formattingOption == "liftOver"):
                        # liftOver test.bed hg18ToHg19.over.chain test_liftover.bed unMapped
                        #inputFile = i_inputDir + id + "_muts_chr" + chromId + ".bed"
                        #outputFile = i_outputDir + id + "_muts_chr" + chromId + ".bed"
                        
                        #inputFile = i_inputDir + id + "_LS454_nwc_chr" + chromId + ".bed"
                        #outputFile = i_outputDir + id + "_LS454_nwc_chr" + chromId + ".bed"
                        
                        inputFile = i_inputDir + id + "_ILLUMINA_nwc_chr" + chromId + ".bed"
                        outputFile = i_outputDir + id + "_ILLUMINA_nwc_chr" + chromId + ".bed"
                        
                        #inputFile = i_inputDir + id + "_chr" + chromId + ".bed"
                        #outputFile = i_outputDir + id + "_chr" + chromId + ".bed"
                        
                        logFile = i_outputDir + "unmapped/" + id + "_chr" + chromId + "_liftOver_unmapped.log"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("liftOver " +  inputFile + " " + i_overChainFilename + " " + outputFile + " " + logFile + "\n")
                
                    # merge the original hg18 .bb file with the new hg19 .bed coordinates file
                    elif (i_formattingOption == "mergeHG18AndHG19"):
                        # python2.7 mergeBedAndBB.py TCGA-AB-2805 12 /inside/home/aradenba/rnaEditing/data/hg18_TCGA-AB-2805_muts_chr12.bb /inside/home/aradenba/rnaEditing/data/hg19_TCGA-AB-2805_muts_chr12.bed /inside/home/aradenba/rnaEditing/data/hg19_TCGA-AB-2805_muts_chr12.bb /inside/home/aradenba/rnaEditing/data/unmapped/TCGA-AB-2805_chr12_merged_unmapped.log
                        script = "/inside/home/aradenba/rnaEditing/scripts/mergeBedAndBB.py"
                        '''
                        originalBB = i_inputDirs[0] + id + "_muts_chr" + chromId + ".bb"
                        liftoverBed = i_inputDirs[1] + id + "_muts_chr"+ chromId + ".bed"
                        outputBBFile = i_outputDir + id + "_muts_chr" + chromId + ".bb"
                        outputUnmappedLog = i_outputDir + "unmapped/" + id + "_chr" + chromId + "_merged_unmapped.log"                        
                        '''
                        '''
                        originalBB = i_inputDirs[0] + id + "_LS454_nwc_chr" + chromId + ".bb"
                        liftoverBed = i_inputDirs[1] + id + "_LS454_nwc_chr"+ chromId + ".bed"
                        outputBBFile = i_outputDir + id + "_LS454_nwc_chr" + chromId + ".bb"
                        outputUnmappedLog = i_outputDir + "unmapped/" + id + "_chr" + chromId + "_LS454_merged_unmapped.log"
                        '''
                        
                        originalBB = i_inputDirs[0] + id + "_ILLUMINA_nwc_chr" + chromId + ".bb"
                        liftoverBed = i_inputDirs[1] + id + "_ILLUMINA_nwc_chr"+ chromId + ".bed"
                        outputBBFile = i_outputDir + id + "_ILLUMINA_nwc_chr" + chromId + ".bb"
                        outputUnmappedLog = i_outputDir + "unmapped/" + id + "_chr" + chromId + "_ILLUMINA_merged_unmapped.log"
                        
                        '''
                        originalBB = i_inputDirs[0] + id + "_chr" + chromId + ".rad"
                        liftoverBed = i_inputDirs[1] + id + "_chr"+ chromId + ".bed"
                        outputBBFile = i_outputDir + id + "_chr" + chromId + ".rad"
                        outputUnmappedLog = i_outputDir + "unmapped/" + id + "_chr" + chromId + "_merged_unmapped.log"
                        '''
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + originalBB + " " + liftoverBed + " " + outputBBFile + " " + outputUnmappedLog + "\n")
                        
                    # validate the liftOver by comparing the overlap of SNPs with Db SNP in both the HG18 and HG19 files
                    elif (i_formattingOption == "validateLiftOver"):
                        # python2.7 validateLiftOver.py TCGA-AB-1234 12 /inside/home/aradenba/rnaEditing/data/hg18.bb /inside/home/aradenba/rnaEditing/data/hg18DbSnp.bed /inside/home/aradenba/rnaEditing/data/hg19.bb /inside/home/aradenba/rnaEditing/data/hg19DbSnp.bed /inside/home/aradenba/rnaEditing/data/chr12_mergeValidation.log
                        script = "/inside/home/aradenba/rnaEditing/scripts/validateLiftOver.py"
                        hg18BB = i_inputDirs[0] + id + "_muts_chr" + chromId + ".bb"
                        hg19BB = i_inputDirs[1] + id + "_muts_chr" + chromId + ".bb"
                        hg18DbSNP = i_inputDirs[2] + "chr" + chromId + ".bed"
                        hg19DbSNP = i_inputDirs[3] + "chr" + chromId + ".bed"
                        outputFile = i_outputDir + "chr" + chromId + "_liftOver_validation.log"       
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + hg18BB + " " + hg18DbSNP + " " + hg19BB + " " + hg19DbSNP + " " + outputFile + "\n")
                    
                    # convert strands
                    elif (i_formattingOption == "convertStrand"): 
                        script = "/inside/home/aradenba/rnaEditing/scripts/convertStrand.py"
                        inputFilename = i_inputDir + id + "_chr" + chromId + ".rad"
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".rad"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + inputFilename + " " + outputFilename + "\n")                                        
                     
                    elif (i_formattingOption == "extractVCF"):
                        inputFile = i_inputDir + id + "*.vcf"
                        # for each file that matches the pattern
                        for file in glob.glob(inputFile): 
                            
                            tempFileHandler = open(file, "r")
                            for line in tempFileHandler:
                                if ("ID=NORMAL" in line and "-11A-" not in line):
                                    i_numOfJobsCounter += 1
                                    outputFile = i_outputDir + id + "_chr" + chromId + ".vcf"
                                    i_fileHandler.write("awk '$1 == \"" + chromId + "\" {print}' " + file + " > " + outputFile + "\n")
                                   
                    # Run RADIA on the data that we have
                    elif (i_formattingOption == "radia"): 
                        #python2.7 radia.py TCGA-AG-A008 22 -n /inside/home/cwilks/bb_pipeline/runs/mutcomp_ex3_usrerun/bams/TCGA-AG-A008-10A-01W-A00K-09_IlluminaGA-DNASeq_exome.bam_HOLD_QC_PENDING.bam
                        #                                   -t /inside/home/cwilks/bb_pipeline/runs/mutcomp_ex3_usrerun/bams/TCGA-AG-A008-01A-01W-A00K-09_IlluminaGA-DNASeq_exome.bam
                        #                                   -f /inside/depot/fa/all_sequences.fasta 
                        #                                   -c /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes_sorted.tab 
                        #                                   -o /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/
                        #                                   -e hg18 -u /inside/depot/fa/all_sequences.fasta
                        #                                   -s /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/radia.py"
                        #dnaNormalFilename = i_inputDirs[0] + id + "-10*.bam"
                        #dnaTumorFilename = i_inputDirs[0] + id + "-01*.bam"
                        #rnaTumorFilename = i_inputDirs[1] + id + "-01*.bam"
                        
                        #dnaNormalFilename = i_inputDirs[0] + id + "N*.bam"
                        #dnaTumorFilename = i_inputDirs[0] + id + "T*.bam"
                        #rnaTumorFilename = i_inputDirs[1] + id + "T*.bam"
            
                        # if we have multiple dirs to search through
                        for directory in i_inputDirs:
                            if (len(glob.glob(directory + id + "-01*.bam")) > 0): 
                                dnaTumorFilename = directory + id + "-01*.bam"
                                if (len(glob.glob(directory + id + "-10*.bam")) > 0):
                                    dnaNormalFilename = directory + id + "-10*.bam"
                                elif (len(glob.glob(directory + id + "-11*.bam")) > 0):
                                    dnaNormalFilename = directory + id + "-11*.bam"
                                break;
                                
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".vcf"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " --batchSize 250000 -n " + dnaNormalFilename + " -t " + dnaTumorFilename + " -r " + rnaTumorFilename + " -f " + i_hg19FastaFilename + " -c " + i_chromSizesFilename + " -o " + outputFilename + " -s /inside/grotto/users/aradenba/data/hg18/grade2_astrocytoma/stats/ -e hg18 -u " + i_hg19FastaFilename + "\n")
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " --batchSize 1000000 -n " + dnaNormalFilename + " -t " + dnaTumorFilename + " -f " + i_hg19FastaFilename + " -c " + i_chromSizesFilename + " -o " + outputFilename + " -s /inside/grotto/users/aradenba/data/hg19/ucec/radia/stats/ -e hg19 -u " + i_hg19FastaFilename + "\n")
                        
                        '''
                        i_fileHandler.write("qsub "  + i_qsubDir + id + "_chr" + chromId + ".sh\n")
                        
                        tempFileHandler = open(i_qsubDir  + id + "_chr" + chromId + ".sh", "w")
                        tempFileHandler.write("#!/bin/bash\n")
                        #tempFileHandler.write("#$ -M aradenba@soe.ucsc.edu\n")
                        #tempFileHandler.write("#$ -m be\n")
                        tempFileHandler.write("source ~/.bashrc\n")
                        tempFileHandler.write("python2.7 " + script + " " + id + " " + chromId + " --batchSize 100000 -n " + dnaNormalFilename + " -t " + dnaTumorFilename + " -f " + i_hg19FastaFilename + " -c " + i_chromSizesFilename + " -o " + outputFilename + " -s /inside/grotto/users/aradenba/data/hg19/ucec/radia/stats/ -e hg19 -u " + i_hg19FastaFilename + "\n")
                        tempFileHandler.close()
                        '''                        
                    elif (i_formattingOption == "extractInfo"): 
                        # python2.7 extractInfo.py TCGA-AB-2995 12 /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995_chr12.vcf ../stats/TCGA-AB-2995_12_blck_DP.tab -f "PASS" -m "DP" --log=DEBUG
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/extractInfo.py"
                        
                        
                        inputFilename = i_inputDir + id + "_chr" + chromId + ".vcf"
                        outputFilename = i_outputDir + id + "_" + chromId + "_rs_RnaMoreStrict_ref_SB.tab"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + inputFilename + " " + outputFilename + " -f PASS -m SB"  + "\n")
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + inputFilename + " " + outputFilename + " -f PASS,GERM,!SOM,!LOH,!EDIT -m SB"  + "\n")
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + inputFilename + " " + outputFilename + " -f PASS,!GERM,SOM,!LOH,!EDIT -m SB"  + "\n")
                    
                    # Run RNA-Editing on the RNA-Seq and wg and exome files
                    elif (i_formattingOption == "rnaEditing"): 
                        #python2.7 rnaEditing.py TCGA-AB-2995 12 /inside/depot/wlee/bccagsc_070610/TCGA-AB-2995-03A-01T-0735-13/all_reads.bam /inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/TCGA-AB-2995_muts_chr12.bb -r /inside/depot/fa/hg19.fasta -c /inside/home/aradenba/rnaEditing/data/chromSizes.tab -Q 25 -q 20 -o /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/TCGA-AB-2995_chr12_test.rad 
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/rnaEditing.py"
                        rnaSeqFilename = i_inputDirs[0] + id + "*/all_reads.bam"
                        bbFilename = i_inputDirs[1] + id + "_muts_chr" + chromId + ".bb"
                        outputFilename = i_outputDir + id + "_chr" + chromId + "." + i_outputFormat.lower()
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + rnaSeqFilename + " " + bbFilename + " -r " + i_hg19FastaFilename + " -c " + i_chromSizesFilename + " -Q " + i_minBaseQual + " -q " + i_minMappingQual + " -o " + outputFilename + " -f " + i_outputFormat.upper() + "\n")
                    
                    # Run RNA-Editing on the RNA-Seq and wg and exome files
                    elif (i_formattingOption == "rnaEditingCheung"): 
                        #python2.7 rnaEditing.py TCGA-AB-2995 12 /inside/depot/wlee/bccagsc_070610/TCGA-AB-2995-03A-01T-0735-13/all_reads.bam /inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/TCGA-AB-2995_muts_chr12.bb -r /inside/depot/fa/hg19.fasta -c /inside/home/aradenba/rnaEditing/data/chromSizes.tab -Q 25 -q 20 -o /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/TCGA-AB-2995_chr12_test.rad 
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/rnaEditing.py"
                        rnaSeqFilename = i_inputDirs[0] + "GSE25840_GM" + id + "*_sorted.bam"
                        bbFilename = i_inputDirs[0] + "GSE25840_GM" + id + "*_sorted.bam"
                        outputFilename = i_outputDir + id + "_chr" + chromId + "." + i_outputFormat.lower()
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + rnaSeqFilename + " " + bbFilename + " -r " + i_hg19FastaFilename + " -c " + i_chromSizesFilename + " -Q " + i_minBaseQual + " -q " + i_minMappingQual + " -m 1 -o " + outputFilename + " -f " + i_outputFormat.upper() + "\n")
                    
                    # First filter out all the blacklisted regions using the filterByCoordinate script
                    elif (i_formattingOption == "filterBlacklists"):
                        #python2.7 filterByCoordinate.py TCGA-AB-2995 12 /inside/home/aradenba/rnaEditing/data/test/filterBlacklist.bed /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995.vcf blck --includeFilterName -f "##FILTER=<ID=blck,Description=\"Position overlaps 1000 Genomes Project blacklist\">"
                        #python2.7 filterByCoordinateBruteForce.py TCGA-AB-2995 12 /inside/home/aradenba/rnaEditing/data/test/filterBlacklist.bed /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995.vcf blck --includeOverlaps --includeFilterName --log=DEBUG -f "##FILTER=<ID=blck,Description=\"Position overlaps 1000 Genomes Project blacklist\">"
                        # ##FILTER=<ID=blq,Description="Position overlaps 1000 Genomes Project mapping quality blacklist">
                        # ##FILTER=<ID=bldp,Description="Position overlap 1000 Genomes Project depth blacklist">
                        
                        #script = "/inside/home/aradenba/rnaEditing/scripts/filterByCoordinateBruteForceLists.py"
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByPybed.py"
                        blacklistFilename = i_inputDirs[0] + "chr" + chromId + "_sorted.bed"
                        #radFilename = i_inputDirs[1] + id + "_chr" + chromId + ".rad"
                        #outputFilename = i_outputDir + id + "_chr" + chromId + ".rad"
                        vcfFilename = i_inputDirs[1] + id + "_chr" + chromId + ".vcf"
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".vcf"
                        #includeOverlaps = "False"
                        #includeFilteringName = "False"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + blacklistFilename + " " + vcfFilename + " blck --includeOverlaps --includeFilterName -f \"##FILTER=<ID=blck,Description=\\\"Position overlaps 1000 Genomes Project blacklist\\\">\" -o " + outputFilename + "\n")
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + blacklistFilename + " " + vcfFilename + " bldp --includeOverlaps --includeFilterName -f \"##FILTER=<ID=bldp,Description=\\\"Position overlaps 1000 Genomes Project depth blacklist\\\">\" -o " + outputFilename + "\n")
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + blacklistFilename + " " + vcfFilename + " blq --includeOverlaps --includeFilterName -f \"##FILTER=<ID=blck,Description=\\\"Position overlaps 1000 Genomes Project mapping quality blacklist\\\">\" -o " + outputFilename + "\n")
                
                        
                        i_fileHandler.write("qsub "  + i_qsubDir + id + "_chr" + chromId + ".sh\n")
                        
                        tempFileHandler = open(i_qsubDir  + id + "_chr" + chromId + ".sh", "w")
                        tempFileHandler.write("#!/bin/bash\n")
                        #tempFileHandler.write("#$ -M aradenba@soe.ucsc.edu\n")
                        #tempFileHandler.write("#$ -m be\n")
                        tempFileHandler.write("source ~/.bashrc\n")
                        tempFileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + blacklistFilename + " " + vcfFilename + " blck --includeOverlaps --includeFilterName -f \"##FILTER=<ID=blck,Description=\\\"Position overlaps 1000 Genomes Project blacklist\\\">\" -o " + outputFilename + "\n")
                        #tempFileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + blacklistFilename + " " + vcfFilename + " bldp --includeOverlaps --includeFilterName -f \"##FILTER=<ID=bldp,Description=\\\"Position overlaps 1000 Genomes Project depth blacklist\\\">\" -o " + outputFilename + "\n")
                        #tempFileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + blacklistFilename + " " + vcfFilename + " blq --includeOverlaps --includeFilterName -f \"##FILTER=<ID=blq,Description=\\\"Position overlaps 1000 Genomes Project mapping quality blacklist\\\">\" -o " + outputFilename + "\n")
                        tempFileHandler.close()
                        
                        
                    # First filter out all the blacklisted regions using the filterByCoordinate script
                    elif (i_formattingOption == "filterTargets"):
                        #python2.7 filterByCoordinate.py TCGA-AB-2995 12 /inside/home/aradenba/rnaEditing/data/test/filterBlacklist.bed /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995.vcf blck --includeFilterName -f "##FILTER=<ID=blck,Description=\"Position overlaps 1000 Genomes Project blacklist\">"
                        #python2.7 filterByCoordinateBruteForce.py TCGA-AB-2995 12 /inside/home/aradenba/rnaEditing/data/test/filterBlacklist.bed /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995.vcf blck --includeOverlaps --includeFilterName --log=DEBUG -f "##FILTER=<ID=blck,Description=\"Position overlaps 1000 Genomes Project blacklist\">"
                        
                        #script = "/inside/home/aradenba/rnaEditing/scripts/filterByCoordinateBruteForceLists.py"
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByPybed.py"
                        targetFilename = i_inputDirs[0] + "chr" + chromId + "_sorted.bed"
                        #radFilename = i_inputDirs[1] + id + "_chr" + chromId + ".rad"
                        #outputFilename = i_outputDir + id + "_chr" + chromId + ".rad"
                        vcfFilename = i_inputDirs[1] + id + "_chr" + chromId + ".vcf"
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".vcf"
                        #includeOverlaps = "False"
                        #includeFilteringName = "False"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + targetFilename + " " + vcfFilename + " ntr --includeFilterName -f \"##FILTER=<ID=ntr,Description=\\\"Position does not overlap with a TCGA target region\\\">\" -o " + outputFilename + "\n")
                        
                        i_fileHandler.write("qsub "  + i_qsubDir + id + "_chr" + chromId + ".sh\n")
                        
                        tempFileHandler = open(i_qsubDir  + id + "_chr" + chromId + ".sh", "w")
                        tempFileHandler.write("#!/bin/bash\n")
                        #tempFileHandler.write("#$ -M aradenba@soe.ucsc.edu\n")
                        #tempFileHandler.write("#$ -m be\n")
                        tempFileHandler.write("source ~/.bashrc\n")
                        tempFileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + targetFilename + " " + vcfFilename + " ntr --includeFilterName -f \"##FILTER=<ID=ntr,Description=\\\"Position does not overlap with a TCGA target region\\\">\" -o " + outputFilename + "\n")
                        tempFileHandler.close()
                        
                    # First filter out all the blacklisted regions using the filterByCoordinate script
                    elif (i_formattingOption == "filterAmbiguous"):
                        #python2.7 filterByCoordinate.py TCGA-AB-2995 12 /inside/home/aradenba/rnaEditing/data/test/filterBlacklist.bed /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995.vcf blck --includeFilterName -f "##FILTER=<ID=blck,Description=\"Position overlaps 1000 Genomes Project blacklist\">"
                        #python2.7 filterByCoordinateBruteForce.py TCGA-AB-2995 12 /inside/home/aradenba/rnaEditing/data/test/filterBlacklist.bed /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995.vcf blck --includeOverlaps --includeFilterName --log=DEBUG -f "##FILTER=<ID=blck,Description=\"Position overlaps 1000 Genomes Project blacklist\">"
                        #python2.7 makeAmbiguousCalls.py TCGA-AG-A016 Y /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/TCGA-AG-A016_chrY.vcf -o /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/ambiguous/TCGA-AG-A016_chrY.vcf
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/makeAmbiguousCalls.py"
                        vcfFilename = i_inputDir + id + "_chr" + chromId + ".vcf"
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".vcf"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + vcfFilename + " -o " + outputFilename + "\n")
                
                    # 
                    elif (i_formattingOption == "filterReadSupport"):
                        #python2.7 filterByReadSuppportVCF.py TCGA-AB-2995 12 /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995.vcf
                        #python2.7 filterByReadSuppportVCF.py TCGA-AG-A016 Y /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995.vcf /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/tcgaTargets/TCGA-AG-A016_chrY.vcf -o /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupport/TCGA-AG-A016_chrY.vcf
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByReadSupportVCF.py"
                        vcfFilename = i_inputDir + id + "_chr" + chromId + ".vcf"
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".vcf"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + vcfFilename + " -o " + outputFilename + "\n")
                        
                        i_fileHandler.write("qsub " + i_qsubDir + id + "_chr" + chromId + ".sh\n")
                        
                        tempFileHandler = open(i_qsubDir  + id + "_chr" + chromId + ".sh", "w")
                        tempFileHandler.write("#!/bin/bash\n")
                        #tempFileHandler.write("#$ -M aradenba@soe.ucsc.edu\n")
                        #tempFileHandler.write("#$ -m be\n")
                        tempFileHandler.write("source ~/.bashrc\n")
                        tempFileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + vcfFilename + " -o " + outputFilename + "\n")
                        tempFileHandler.close()
                        
                    elif (i_formattingOption == "vcf2Excel"):
                        # python2.7 vcf2Excel.py TCGA-02-0033 7 hg19 /inside/home/aradenba/rnaEditing/data/test/TCGA-02-0033_chr7.vcf /inside/home/aradenba/rnaEditing/data/test/TCGA-02-0033_chr7.tab
                        script = "/inside/home/aradenba/rnaEditing/scripts/vcf2Excel.py"
                        
                        vcfFilename = i_inputDir + id + "_chr" + chromId + ".vcf"
                        tabFilename = i_outputDir + id + "_chr" + chromId + ".tab"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + i_hgAssembly + " " + vcfFilename + " " + tabFilename + "\n")
                
                    elif (i_formattingOption == "radiaCompare"):
                        #python2.7 radiaCompare.py TCGA-AB-2995 12 /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995_radia.vcf /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995_bb.vcf -c blck=blck,DB=DB,GERM=Germline,SOM=Somatic -s /inside/home/aradenba/rnaEditing/data/test/stats/stats.tab --log=DEBUG
                        #python2.7 radiaCompare.py TCGA-AG-A016 Y /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/filters/readSupportLikeBB/TCGA-AG-A016_chrY.vcf /inside/grotto/users/aradenba/data/hg18/benchmark/bambam/TCGA-AG-A016_chrY.vcf -c blck=blck,DB=DB,GERM=Germline,SOM=Somatic -s /inside/grotto/users/aradenba/data/hg18/benchmark/radia/illumina/stats/cmpRadBB.tab
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/radiaCompare.py"
                        radFilename = i_inputDirs[0] + id + "_chr" + chromId + ".vcf"
                        cmpFilename = i_inputDirs[1] + id + "_chr" + chromId + ".vcf"
                        outputFilename = i_outputDirs[0] + "cmpStrictLoose.tab"
                        overlapFilename = i_outputDirs[1] + id + "_chr" + chromId + "_ovlp.vcf"
                        nonOverlapFilename = i_outputDirs[2] + id + "_chr" + chromId + "_nonovlp.vcf"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + radFilename + " " + cmpFilename + " -c \"GERM=Germline;SOM=Somatic\" -s " + outputFilename + " -o " + overlapFilename + " -n " + nonOverlapFilename +"\n")
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + radFilename + " " + cmpFilename + " -c \"DB=DB\" -s " + outputFilename + " -o " + overlapFilename + " -n " + nonOverlapFilename +"\n")
                        
                        '''
                        i_fileHandler.write("qsub " + i_qsubDir + id + "_chr" + chromId + ".sh\n")
                        
                        tempFileHandler = open(i_qsubDir  + id + "_chr" + chromId + ".sh", "w")
                        tempFileHandler.write("#!/bin/bash\n")
                        #tempFileHandler.write("#$ -M aradenba@soe.ucsc.edu\n")
                        #tempFileHandler.write("#$ -m be\n")
                        tempFileHandler.write("source ~/.bashrc\n")
                        tempFileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + radFilename + " " + cmpFilename + " -c \"blq=blq;bldp=bldp\" -s " + outputFilename + " -o " + overlapFilename + " -n " + nonOverlapFilename +"\n")
                        tempFileHandler.close()
                        '''
                        
                    # Second filter out all the dbSNP regions using the filterByCoordinate script
                    elif (i_formattingOption == "filterDbSnp"):
                        #python2.7 filterByCoordinate.py TCGA-AB-2995 12 /inside/home/aradenba/rnaEditing/data/test/filterSNP.bed /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995.vcf DB -d INFO --includeFilterName --includeIdName -f "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 135\">"
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByCoordinateBruteForceDict.py"
                        dbSNPFilename = i_inputDirs[0] + "chr" + chromId + ".bed"
                        #radFilename = i_inputDirs[1] + id + "_chr" + chromId + ".rad"
                        #outputFilename = i_outputDir + id + "_chr" + chromId + ".rad"
                        vcfFilename = i_inputDirs[1] + id + "_chr" + chromId + ".vcf"
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".vcf"
                        #includeOverlaps = "False"
                        #includeFilteringName = "False"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + dbSNPFilename + " " + vcfFilename + " DB --includeOverlaps --includeFilterName --includeIdName -d INFO -f \"##INFO=<ID=DB,Number=0,Type=Flag,Description=\\\"dbSNP membership, build 130\\\">\" -o " + outputFilename + "\n")
                    
                        i_fileHandler.write("qsub " + i_qsubDir  + id + "_chr" + chromId + ".sh\n")
                        
                        tempFileHandler = open(i_qsubDir  + id + "_chr" + chromId + ".sh", "w")
                        tempFileHandler.write("#!/bin/bash\n")
                        #tempFileHandler.write("#$ -M aradenba@soe.ucsc.edu\n")
                        #tempFileHandler.write("#$ -m be\n")
                        tempFileHandler.write("source ~/.bashrc\n")
                        tempFileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + dbSNPFilename + " " + vcfFilename + " DB --includeOverlaps --includeFilterName --includeIdName -d INFO -f \"##INFO=<ID=DB,Number=0,Type=Flag,Description=\\\"dbSNP membership, build 135\\\">\" -o " + outputFilename + "\n")
                        tempFileHandler.close()
                        
                        
                    elif (i_formattingOption == "createMutClassInput"):
                        # python2.7 createMutClassInput.py TCGA-AB-2995 12 /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995.vcf /inside/home/aradenba/rnaEditing/data/test/TCGA-AB-2995_mutClass_input.bed
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/createMutClassInput.py"
                        vcfFilename = i_inputDir + id + "_chr" + chromId + ".vcf"
                        mutClassInputFilename = i_outputDir + id + "_chr" + chromId + ".bed"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + vcfFilename + " " + mutClassInputFilename + "\n")
                        
                        i_fileHandler.write("qsub " + i_qsubDir  + id + "_chr" + chromId + ".sh\n")
                        
                        tempFileHandler = open(i_qsubDir  + id + "_chr" + chromId + ".sh", "w")
                        tempFileHandler.write("#!/bin/bash\n")
                        #tempFileHandler.write("#$ -M aradenba@soe.ucsc.edu\n")
                        #tempFileHandler.write("#$ -m be\n")
                        tempFileHandler.write("source ~/.bashrc\n")
                        tempFileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + vcfFilename + " " + mutClassInputFilename + "\n")
                        tempFileHandler.close()
                    
                    # Filter out all the BB MUT calls in the Cheung data
                    elif (i_formattingOption == "filterBBMuts"):
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByCoordinate.py"
                        #bbMutFilename = i_inputDirs[0] + id + "_LS454_nc_chr" + chromId + ".bed"
                        #outputFilename = i_outputDir + id + "_chr" + chromId + "_454.rad"
                        #bbMutFilename = i_inputDirs[0] + id + "_ILLUMINA_nc_chr" + chromId + ".bed"
                        #outputFilename = i_outputDir + id + "_chr" + chromId + "_ILL.rad"
                        
                        bbMutFilename = bbFilename = i_inputDirs[0] + id + "_sorted_muts_chr" + chromId + ".bb"
                        radFilename = i_inputDirs[1] + id + "_chr" + chromId + ".rad"
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".rad"
                        includeOverlaps = "False"
                        includeFilteringName = "False"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + bbMutFilename + " " + radFilename + " " + includeOverlaps + " -i " + includeFilteringName + " -o " + outputFilename + "\n")
                    
                    # Filter out all the not covered (NC) regions using the filterByCoordinate script
                    elif (i_formattingOption == "filterBBNCs"):
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByCoordinate.py"
                        bbFilename = i_inputDirs[0] + id + "_nwc_chr" + chromId + ".bb"
                        radFilename = i_inputDirs[1] + id + "_chr" + chromId + ".rad"
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".rad"
                        includeOverlaps = "False"
                        includeFilteringName = "False"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + bbFilename + " " + radFilename + " " + includeOverlaps + " -i " + includeFilteringName + " -o " + outputFilename + "\n")
                
                    elif (i_formattingOption == "catBBMuts"): 
                        file454 = i_inputDir + id + "_chr" + chromId + "_454.rad"
                        fileIllumina = i_inputDir + id + "_chr" + chromId + "_ILL.rad"
                        outputFile = i_outputDir + id + "_chr" + chromId + ".rad"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        
                        i_fileHandler.write("cat " + file454 + " " + fileIllumina + " | sort -u +1n -2 > " + outputFile + "\n")                        
                
                    elif (i_formattingOption == "removeCatFiles"): 
                        file454 = i_inputDir + id + "_chr" + chromId + "_454.rad"
                        fileIllumina = i_inputDir + id + "_chr" + chromId + "_ILL.rad"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        
                        i_fileHandler.write("rm " + file454 + "\n")
                        i_fileHandler.write("rm " + fileIllumina + "\n")
                     
                    # Filter out all non-transcript regions using the filterByCoordinate script
                    elif (i_formattingOption == "filterTranscripts"):
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByCoordinate.py"
                        refGeneFilename = i_inputDirs[0] + i_hgAssembly + "_chr" + chromId + "_tx_sorted_uniq.bed"
                        radFilename = i_inputDirs[1] + id + "_chr" + chromId + ".rad"
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".rad"
                        includeOverlaps = "True"
                        includeFilteringName = "True"
                        includeStrand = "True"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + refGeneFilename + " " + radFilename + " " + includeOverlaps + " -i " + includeFilteringName + " -s " + includeStrand + " -o " + outputFilename + "\n")
                    
                    # Filter out all non-CDS regions using the filterByCoordinate script
                    elif (i_formattingOption == "filterCDS"):
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByCoordinate.py"
                        refGeneFilename = i_inputDirs[0] + i_hgAssembly + "_chr" + chromId + "_cds_sorted_uniq.bed"
                        radFilename = i_inputDirs[1] + id + "_chr" + chromId + ".rad"
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".rad"
                        includeOverlaps = "True"
                        includeFilteringName = "True"
                        includeStrand = "True"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + refGeneFilename + " " + radFilename + " " + includeOverlaps + " -i " + includeFilteringName + " -s " + includeStrand + " -o " + outputFilename + "\n")
                        
                    # Filter out all non-exon regions using the filterByCoordinate script
                    elif (i_formattingOption == "filterExons"):
                        
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByCoordinate.py"
                        refGeneFilename = i_inputDirs[0] + i_hgAssembly + "_chr" + chromId + "_exons_sorted_uniq.bed"
                        radFilename = i_inputDirs[1] + id + "_chr" + chromId + ".rad"
                        outputFilename = i_outputDir + id + "_chr" + chromId + ".rad"
                        includeOverlaps = "True"
                        includeFilteringName = "True"
                        includeStrand = "True"
                        #includeStrand = "False"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + refGeneFilename + " " + radFilename + " " + includeOverlaps + " -i " + includeFilteringName + " -s " + includeStrand + " -o " + outputFilename + "\n")
                    
                    # Filter out all events that don't meet the minimum read support requirements using the filterByReadSupport script
                    elif (i_formattingOption == "filterTranscriptReadSupport"):
                        # python2.7 filterByReadSupport.py radFile minTotalBBReads minTotalRNASeqReads minSupportingRNASeqReads minPercentReadsSupportingEvent minPercentTotalReads -o outputFile -m mutOutputFile
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByReadSupport.py"
                        radFilename = i_inputDir + id + "_chr" + chromId + ".rad"
                        radOutputFile = i_outputDir + "transcriptReadSupport/" + id + "_chr" + chromId + ".rad"
                        bedOutputFile = i_outputDir + "transcriptReadSupport/bed/" + id + "_chr" + chromId + ".bed"
                        mutOutputFile = i_outputDir + "transcriptMutClass/input/" + id + "_chr" + chromId + ".bed"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + radFilename + " 4 10 2 10 90 -r " + radOutputFile + " -b " + bedOutputFile + " -m " + mutOutputFile + "\n")
                    
                    # Filter out all events that don't meet the minimum read support requirements using the filterByReadSupport script
                    elif (i_formattingOption == "filterCDSReadSupport"):
                        # python2.7 filterByReadSupport.py radFile minTotalBBReads minTotalRNASeqReads minSupportingRNASeqReads minPercentReadsSupportingEvent minPercentTotalReads -o outputFile -m mutOutputFile
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByReadSupport.py"
                        radFilename = i_inputDir + id + "_chr" + chromId + ".rad"
                        radOutputFile = i_outputDir + "cdsReadSupport/" + id + "_chr" + chromId + ".rad"
                        bedOutputFile = i_outputDir + "cdsReadSupport/bed/" + id + "_chr" + chromId + ".bed"
                        mutOutputFile = i_outputDir + "cdsMutClass/input/" + id + "_chr" + chromId + ".bed"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + radFilename + " 4 10 2 10 90 -r " + radOutputFile + " -b " + bedOutputFile + " -m " + mutOutputFile + "\n")
                    
                    # Filter out all events that don't meet the minimum read support requirements using the filterByReadSupport script
                    elif (i_formattingOption == "filterExonReadSupport"):
                        # python2.7 filterByReadSupport.py radFile minTotalBBReads minTotalRNASeqReads minSupportingRNASeqReads minPercentReadsSupportingEvent minPercentTotalReads -o outputFile -m mutOutputFile
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByReadSupport.py"
                        radFilename = i_inputDir + id + "_chr" + chromId + ".rad"
                        radOutputFile = i_outputDir + "exonReadSupport/" + id + "_chr" + chromId + ".rad"
                        bedOutputFile = i_outputDir + "exonReadSupport/bed/" + id + "_chr" + chromId + ".bed"
                        mutOutputFile = i_outputDir + "exonMutClass/input/" + id + "_chr" + chromId + ".bed"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + radFilename + " 4 10 2 10 90 -r " + radOutputFile + " -b " + bedOutputFile + " -m " + mutOutputFile + "\n")
                    
                    # Filter out all events that don't meet the minimum read support requirements using the filterByReadSupport script
                    elif (i_formattingOption == "filterDnaReadSupport"):
                        # python2.7 filterByReadSupport.py radFile minTotalBBReads minTotalRNASeqReads minSupportingRNASeqReads minPercentReadsSupportingEvent minPercentTotalReads -o outputFile -m mutOutputFile
                        script = "/inside/home/aradenba/rnaEditing/scripts/filterByReadSupport.py"
                        radFilename = i_inputDir + id + "_chr" + chromId + ".rad"
                        radOutputFile = i_outputDir + id + "_chr" + chromId + ".rad"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + radFilename + " 100 40 10 10 90 -r " + radOutputFile + "\n")
                                        
                    # analyze the BB MUT calls
                    elif (i_formattingOption == "bbAnalysis"):
                        script = "/inside/home/aradenba/rnaEditing/scripts/bbAnalysis.py"
                        inputFile = i_inputDir + id + "_muts_chr" + chromId + ".bb"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        
                        i_fileHandler.write("python2.7 " + script + " " + inputFile + "\n")
                        
                    # Sort a .bb file according to the start and stop coordinates
                    elif (i_formattingOption == "sortBB"):
                        inputFile = i_inputDir + id + "_muts_chr" + chromId + ".bb"
                        outputFile = i_outputDir + id + "_sorted_muts_chr" + chromId + ".bb"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        # sort numerically on columns 2 and 3 with tab as separator
                        # .bb file has columns: 0=bbCall, 1=chr, 2=startCoordinate, 3=stopCoordinate, etc.
                        i_fileHandler.write("sort -t $'\\t' +2n -4 " + inputFile + " > " + outputFile + "\n")    
                    
                    # Sort a .bed file according to the start and stop coordinates
                    elif (i_formattingOption == "sortBed"):
                        inputFile = i_inputDir + id + "_muts_chr" + chromId + ".bed"
                        outputFile = i_outputDir + id + "_sorted.bed"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        # sort numerically on columns 1 and 2 with tab as separator
                        # .bed file has columns: 0=chr, 1=startCoordinate, 2=stopCoordinate, etc.                    
                        i_fileHandler.write("sort -t $'\\t' +1n -3 " + inputFile + " > " + outputFile + "\n")
                    
                    elif (i_formattingOption == "getMutClassData"):
                        script = "/inside/home/aradenba/rnaEditing/scripts/extractReadSupport.py"
                        inputFile = i_inputDir + id + "_chr" + chromId + "_exon.rad"
                        outputFile = i_outputDir + id + "_chr" + chromId + "_readSupport.bed"
                        mutOutputFile = i_outputDir + "mutClass/" + id + "_chr" + chromId + "_mutClass.bed"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("python2.7 " + script + " " + inputFile + " -o " + outputFile + " -m " + mutOutputFile + "\n")
                    
                    elif (i_formattingOption == "runMutClass"):
                        script = "mutClass -lazyLoading -skipRegulatory " + i_hgAssembly
                        inputFile = i_inputDir + id + "_chr" + chromId + ".bed"
                        outputFile = i_outputDir + id + "_chr" + chromId + ".bed"
                        i_numOfJobsCounter += 1
                        #i_fileHandler.write(script + " " + inputFile + " > " + outputFile + "\n")
                        
                        i_fileHandler.write("qsub " + i_qsubDir + id + "_chr" + chromId + ".sh\n")
                        
                        tempFileHandler = open(i_qsubDir + id + "_chr" + chromId + ".sh", "w")
                        tempFileHandler.write("#!/bin/bash\n")
                        #tempFileHandler.write("#$ -M aradenba@soe.ucsc.edu\n")
                        #tempFileHandler.write("#$ -m be\n")
                        tempFileHandler.write("source ~/.bashrc\n")
                        tempFileHandler.write(script + " " + inputFile + " > " + outputFile + "\n")
                        tempFileHandler.close()
                        
                    elif (i_formattingOption == "mergeVCFAndMutClassOutput"):
                        script = "/inside/home/aradenba/rnaEditing/scripts/mergeVCFAndMutClassOutput.py"
                        vcfInputFile = i_inputDirs[0] + id + "_chr" + chromId + ".vcf"
                        mutClassFile = i_inputDirs[1] + id + "_chr" + chromId + ".bed"
                        vcfMergedFile = i_outputDir + id + "_chr" + chromId + ".vcf"
                        i_numOfJobsCounter += 1
                        #i_fileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + vcfInputFile + " " + mutClassFile + " -o " + vcfMergedFile + "\n")
                        
                        i_fileHandler.write("qsub " + i_qsubDir  + id + "_chr" + chromId + ".sh\n")
                        
                        tempFileHandler = open(i_qsubDir  + id + "_chr" + chromId + ".sh", "w")
                        tempFileHandler.write("#!/bin/bash\n")
                        #tempFileHandler.write("#$ -M aradenba@soe.ucsc.edu\n")
                        #tempFileHandler.write("#$ -m be\n")
                        tempFileHandler.write("source ~/.bashrc\n")
                        tempFileHandler.write("python2.7 " + script + " " + id + " " + chromId + " " + vcfInputFile + " " + mutClassFile + " -o " + vcfMergedFile + "\n")
                        tempFileHandler.close()
                        
                    elif (i_formattingOption == "downloadCheungDNA"):
                        # ascp -i asperaweb_id_dsa.putty -QTr -l1000m -k1 anonftp@ftp-trace.ncbi.nlm.nih.gov:/1000genomes/ftp/data/NA06985/alignment/NA06985.chrom22.LS454.ssaha2.CEU.exon_targetted.20100311.bam .
                        
                        # download either from NCBI or EBI
                        #script = "ascp -i asperaweb_id_dsa.putty -QTr -l1000m -k1 anonftp@ftp-trace.ncbi.nlm.nih.gov:/1000genomes/ftp/data/"
                        #script = "ascp -i asperaweb_id_dsa.putty -QTr -l1000m -k1 fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/data/"
                        script = "ascp -i asperaweb_id_dsa.putty -QTr -l1000m -k1 anonftp@ftp-trace.ncbi.nlm.nih.gov:/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/Primary_Assembly/assembled_chromosomes/FASTA/chrY.fa.gz"
                        # download either data aligned with either ssaha2 or bwa 
                        #inputFile = "NA" + id + "/alignment/NA" + id + ".chrom" + chromId + ".LS454.ssaha2.CEU.exon_targetted.20100311.bam"
                        inputFile = "NA" + id + "/alignment/NA" + id + ".chrom" + chromId + ".ILLUMINA.bwa.CEU.exon_targetted.20100311.bam"
                        
                        outputDir = "cheung_data/"
                        i_fileHandler.write(script + inputFile + " " + outputDir + "\n")
                        
                    # select a certain part of ch12 and convert from .bam to .sam
                    elif (i_formattingOption == "samtools"):
                        inputFile = i_inputDir + id + "*/all_reads.bam"
                        outputFile = i_outputDir + id + ".sam"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("samtools view -bh " + inputFile + " chr12:7,055,740-7,070,479 | samtools pileup - > " + outputFile + "\n")
                
                    # Select ptpn6 region out of hg19 bambam files 
                    elif (i_formattingOption == "ptpn6BamBam"):
                        inputFile = i_inputDir + id + "*.bb"
                        outputFile = i_outputDir + id + "_ptpn6.bb"
                        
                        # keep track of the number of jobs that are generated for validation
                        i_numOfJobsCounter += 1
                        i_fileHandler.write("awk '$3 >= 7055740 && $4 <= 7070479' " + inputFile + " > " + outputFile + "\n")
                        
                    else:
                        print >> sys.stderr, "Invalid formatting option specified.  Please check documentation."
       
    # This number can be compared to the number of lines in the job list file
    if (i_debug):         
        print >> sys.stderr, "Number of Jobs:", i_numOfJobsCounter
        
    # This helps catch when expected files are missing.            
    if (i_numOfJobsCounter != 0 and i_numOfJobsCounter != i_numOfExpectedJobs):
        print >> sys.stderr, "Warning:  The number of jobs that were created (" + str(i_numOfJobsCounter) + ") doesn't match the number of expected jobs (" + str(i_numOfExpectedJobs) + ")."
            
    i_fileHandler.close()   
    return

main()
sys.exit(0)
