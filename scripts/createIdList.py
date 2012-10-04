#!/usr/bin/env python2.7

import sys                          # system module
import os                           # os module for checking pathnames
import re                           # re module for matching ids
import rnaEditingUtil               # utility functions for rna editing
from optparse import OptionParser   # used for parsing command line arguments

# this regular expression will match full TCGA sample Ids, e.g. TCGA-AG-A016-01A-01R or TCGA-37-4133-10A-01D
i_tcgaNameRegEx = re.compile("TCGA-(\\w){2}-(\\w){4}-(\\w){3}-(\\w){3}")

'''
'    Amie Radenbaugh - 11/08/2010
'    UCSC - RNA Editing  
'    Program name: "createIdList.py"
'   
'    This program finds all of the TCGA Ids for which we have both RNA-Seq data and either whole-genome or exome-capture data.  The user
'    must specify the type(s) of TCGA Ids to look for as the first argument to the program.  For example, when looking for AML data, specify
'    the TCGA code for AML which is "AB" as the first argument.  
'
'    This program looks for the RNA-Seq data in the directory or directories specified in the first argument.  It looks for the whole-genome 
'    and exome-capture data in the directories specified with the -w and -e options respectively.  Multiple directories can be specified with 
'    a comma-separated list.  A comma-separated list of all of the unique over-lapping Ids will be output to standard out or to the file 
'    provided with the -o option.
'
'    Here are some sample commands for the program:
'    python2.7 createIdList.py AB /inside/depot/wlee/bccagsc_070610/ -w /inside/grotto/bambam/aml/wg/,/inside/grotto/users/aradenba/data/hg19/aml/wg/bambam/ -o ../cancer/aml/idLists/wgOverlappingIds.txt
'    python2.7 createIdList.py AB /inside/depot/wlee/bccagsc_070610/ -e /inside/grotto/bambam/aml/exome/,/inside/grotto/users/aradenba/data/hg19/aml/exome/bambam/ -o ../cancer/aml/idLists/exomeOverlappingIds.txt
'    python2.7 createIdList.py AB /inside/depot/wlee/bccagsc_070610/ -w /inside/grotto/bambam/aml/wg/ -e /inside/grotto/bambam/aml/exome/ -o ../cancer/aml/idLists/allOverlappingIds.txt
'
'    These are the commands that I use:
'
'    For AML: 
'    python2.7 createIdList.py AB /inside/depot/wlee/bccagsc_070610/ -w /inside/grotto/bambam/aml/wg/ -o ../cancer/aml/idLists/wgOverlappingIds.txt
'    python2.7 createIdList.py AB /inside/depot/wlee/bccagsc_070610/ -e /inside/grotto/bambam/aml/exome/ -o ../cancer/aml/idLists/exomeOverlappingIds.txt
'    python2.7 createIdList.py AB /inside/depot/wlee/bccagsc_070610/ -w /inside/grotto/bambam/aml/wg/ -e /inside/grotto/bambam/aml/exome/ -o ../cancer/aml/idLists/allOverlappingIds.txt -d ../cancer/aml/idLists/dnaIds.txt -r ../cancer/aml/idLists/rnaIds.txt
'
'    For Colorectal:
'    python2.7 createIdList.py AA,AG /inside/depot/coad/rnaseq/ -w /inside/grotto/bambam/coad_read_freeze/ -o ../cancer/coad/idLists/wgOverlappingIds.txt
'    python2.7 createIdList.py AA,AG /inside/depot/coad/rnaseq/ -e /inside/grotto/bambam/coad_read_freeze/ -o ../cancer/coad/idLists/exomeOverlappingIds.txt
'
'
'    Significant changes made to code - 08/15/2012
'    Above commands may no longer work
'
'    For BRCA data that was downloaded pre-cgHub, there is a manifest file.  Here are some useful commands to make id lists:
'    grep "Breast" /inside/home/cwilks/cghub/TCGA_phs000178_tranche_Q-3.csv.complete2.added.uuids | grep "RNA" | grep -v "mirna" | egrep -o "TCGA-(\\w){2}-(\\w){4}" | sort -u > rnaIdsFromManifest.txt
'    grep "Breast" /inside/home/cwilks/cghub/TCGA_phs000178_tranche_Q-3.csv.complete2.added.uuids | grep "DNA" | grep "WGS" | egrep -o "TCGA-(\\w){2}-(\\w){4}" | sort -u > wgIdsFromManifest.txt
'    grep "Breast" /inside/home/cwilks/cghub/TCGA_phs000178_tranche_Q-3.csv.complete2.added.uuids | grep "DNA" | grep "exome" | egrep -o "TCGA-(\\w){2}-(\\w){4}" | sort -u > exomeIdsFromManifest.txt
'
'    TCGA Codes for Breast = A1,A2,A7,A8,AN,AO,AQ,AR,B6,BH,C8,CG,D8,E2
'    For each id file, do the following:
'    egrep -o "TCGA-(\\w){2}" exomeIdsFromManifest.txt | egrep -o "A(\\w){1}" | sort -u
'    egrep -o "TCGA-(\\w){2}" exomeIdsFromManifest.txt | egrep -o "B(\\w){1}" | sort -u
'
'    Directories for data:
'    grep "Breast" /inside/home/cwilks/cghub/TCGA_phs000178_tranche_Q-3.csv.complete2.added.uuids | grep "RNA" | grep -v "mirna" | awk '{FS=","} ; {print $18}' | egrep -o "/inside/(.*)/" | sort -u > rnaDirList.txt
'    grep "Breast" /inside/home/cwilks/cghub/TCGA_phs000178_tranche_Q-3.csv.complete2.added.uuids | grep "DNA" | grep "WGS" | awk '{FS=","} ; {print $18}' | egrep -o "/inside/(.*)/" | sort -u > dnaWGDirList.txt
'    grep "Breast" /inside/home/cwilks/cghub/TCGA_phs000178_tranche_Q-3.csv.complete2.added.uuids | grep "DNA" | grep "exome" | awk '{FS=","} ; {print $18}' | egrep -o "/inside/(.*)/" | sort -u > dnaExomeDirList.txt
'
'    For BRCA:
'    python2.7 createIdList.py -t A1,A2,A7,A8,AN,AO,AQ,AR,B6,BH,C8,CG,D8,E2 -w ../cancer/brca/idLists/dnaWGDirList.txt -e ../cancer/brca/idLists/dnaExomeDirList.txt -s ../cancer/brca/idLists/rnaDirList.txt -o ../cancer/brca/idLists/overlappingIds.txt -d ../cancer/brca/idLists/dnaIds.txt -r ../cancer/brca/idLists/rnaIds.txt
'    
'''


def get_ids_from_dirs(aCancerTypesList, aDirList):
    '''
    '  Return all of the unique Ids from the sub-directories found in the list of directories given.
    '  Returns a set of Ids with the format "TCGA-XX-YYYY" where X is a 
    '  TCGA cancer type such as "AB" and Y is alphanumeric.  For example, "TCGA-AB-1234".
    '''
    
    ids = set()
    for directory in aDirList:
        dirList = os.listdir(directory)

        # some files are in directories labeled with the TCGA Id
        for dirName in dirList:
            # get an iterator of match objects for the TCGA id
            iterator = i_tcgaNameRegEx.finditer(dirName)
            
            # for each match object in the iterator
            for match in iterator:
                # get the pattern that matched the reg ex, i.e. TCGA-AB-1234
                tcgaId = match.group()
                # split the id
                idList = tcgaId.split("-")
                # if there are no cancer types, then add the id
                # if there are cancer types, then make sure the type is in the list
                if (len(aCancerTypesList) == 0 or idList[1] in aCancerTypesList):
                    id = idList[0] + "-" + idList[1] + "-" + idList[2]
                    ids.add(id) 
    
    return ids


def get_ids_from_filenames(aCancerTypesList, aDirList, aUniqueKey):
    '''
    '  Return all of the unique Ids from the files found in the list of directories given.
    '  Returns a set of whole-genome or exome-capture Ids with the format "TCGA-XX-YYYY" where X is a 
    '  TCGA cancer type such as "AB" and Y is alphanumeric.  For example, "TCGA-AB-1234".
    '''

    ids = set()
    for directory in aDirList:
        dirList = os.listdir(directory)
    
        # for each file in the directory 
        for fileName in dirList:
            
            if (aUniqueKey in fileName):
                # get an iterator of match objects for the TCGA id
                iterator = i_tcgaNameRegEx.finditer(fileName)
                
                # for each match object in the iterator
                for match in iterator:
                    # get the pattern that matched the reg ex, i.e. TCGA-AB-1234
                    tcgaId = match.group()
                    
                    # split the id 
                    idList = tcgaId.split("-")
                    # if there are no cancer types, then add the id
                    # if there are cancer types, then make sure the type is in the list
                    if (len(aCancerTypesList) == 0 or idList[1] in aCancerTypesList):
                        id = idList[0] + "-" + idList[1] + "-" + idList[2]
                        ids.add(id) 	
                    #else:
                        #print "not correct type", tcgaId
            #else:
                #print "no bam or no uniqueKey", fileName
    return ids


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
    print >> sys.stderr, len(sys.argv), sys.argv
    
    # create the usage statement
    usage = "usage: python2.7 %prog [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="name of file where overlapping Ids should be written, STDOUT by default")
    i_cmdLineParser.add_option("-d", "--dnaFilename", dest="dnaFilename", metavar="DNA_ID_FILE", help="name of file where all DNA Ids should be written")
    i_cmdLineParser.add_option("-r", "--rnaFilename", dest="rnaFilename", metavar="RNA_ID_FILE", help="name of file where all RNA Ids should be written")
    i_cmdLineParser.add_option("-w", "--wholeGenomeDirs", dest="wgDirs", metavar="WG_DIR_LIST_FILE", help="name of file with a list of dirs where whole-genome files are located")
    i_cmdLineParser.add_option("-e", "--exomeCaptureDirs", dest="exomeDirs", metavar="EXOME_DIR_LIST_FILE", help="name of file with a list of dirs where exome-capture files are located")
    i_cmdLineParser.add_option("-s", "--rnaDirs",  dest="rnaDirs", metavar="RNA_DIR_LIST_FILE", help="name of file with a list of dirs where rna-seq files are located")
    i_cmdLineParser.add_option("-t", "--TCGACancerTypes", dest="cancerTypes", metavar="CANCER_TYPES", help="comma-separated list of specific types of cancer that should be selected from dirs")
        
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(3,16,2)
    i_argLength = len(sys.argv)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        # print >> sys.stderr, "usage: python2.7 createIdList.py -t TCGACancerType(s) -s rnaSeqDir(s) -w wgDir(s) -e exomeDir(s) -o out.ids -d dna.ids -r rna.ids"
        i_cmdLineParser.print_help()
        sys.exit(1)
    
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    
    i_writeFilenameList = list()
    i_readFilenameList = list()
    i_outputFilename = None
    i_wgDirFilename = None
    i_exomeDirFilename = None
    i_rnaDirFilename = None
    i_dnaFilename = None
    i_rnaFilename = None
    i_cancerTypes = list()
    
    # try to get any optional parameters with no defaults
    if (i_cmdLineOptions.wgDirs != None):
        i_wgDirFilename = i_cmdLineOptions.wgDirs
        i_readFilenameList.append(i_wgDirFilename)
    if (i_cmdLineOptions.exomeDirs != None):
        i_exomeDirFilename = i_cmdLineOptions.exomeDirs
        i_readFilenameList.append(i_exomeDirFilename)
    if (i_cmdLineOptions.rnaDirs != None):
        i_rnaDirFilename = i_cmdLineOptions.rnaDirs
        i_readFilenameList.append(i_rnaDirFilename)
    if (i_cmdLineOptions.outputFilename != None):
        i_outputFilename = i_cmdLineOptions.outputFilename
        i_writeFilenameList.append(i_outputFilename)
    if (i_cmdLineOptions.dnaFilename != None):
        i_dnaFilename = i_cmdLineOptions.dnaFilename
        i_writeFilenameList.append(i_dnaFilename)
    if (i_cmdLineOptions.rnaFilename != None):
        i_rnaFilename = i_cmdLineOptions.rnaFilename
        i_writeFilenameList.append(i_rnaFilename)
    if (i_cmdLineOptions.cancerTypes != None):
        i_cancerTypes = str(i_cmdLineOptions.cancerTypes).split(",")

    # check for any errors        
    if (not rnaEditingUtil.check_for_argv_errors(None, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)    

    # get the directories
    i_dirList = list()
    i_exomeDirList = list()
    i_wgDirList = list()
    i_rnaDirList = list()
    if (i_rnaDirFilename != None):
        i_rnaDirList = get_dir_list(i_rnaDirFilename)
        i_dirList += i_rnaDirList
    if (i_wgDirFilename != None):
        i_wgDirList = get_dir_list(i_wgDirFilename)
        i_dirList += i_wgDirList
    if (i_exomeDirFilename != None):
        i_exomeDirList = get_dir_list(i_exomeDirFilename)
        i_dirList += i_exomeDirList
        
    print "# of RNA Dirs", len(i_rnaDirList)
    print "# of WG DNA Dirs", len(i_wgDirList)
    print "# of Exome DNA Dirs", len(i_exomeDirList)
    
    # check for any errors        
    if (not rnaEditingUtil.check_for_argv_errors(i_dirList, None, None)):
        sys.exit(1)    
        
    if (len(i_dirList) == 0):
        print >> sys.stderr, "Error:  You must specify at least one directory."
        sys.exit(1)
    
    # initialize some vars
    i_wgIds = list()
    i_exomeIds = list()
    i_rnaIds = list()
    
    # get RNA-Seq Ids
    if (len(i_rnaDirList) > 0):
        i_rnaIds = get_ids_from_filenames(i_cancerTypes, i_rnaDirList, "trimmed.annotated.translated_to_genomic")
        print "RNAIds:", len(i_rnaIds)
    
    # get whole-genome or exome-capture Ids
    if (len(i_wgDirList) != 0):
        i_wgIds = get_ids_from_filenames(i_cancerTypes, i_wgDirList, "DNASeq_whole")
        print "WGIds:", len(i_wgIds)
    
    if (len(i_exomeDirList) != 0):
        i_exomeIds = get_ids_from_filenames(i_cancerTypes, i_exomeDirList, "DNASeq_exome")
        print "ExomeIds:", len(i_exomeIds)
    
    # write the common Ids to the output file
    i_rnaWGIntersection = i_rnaIds.intersection(i_wgIds)
    i_rnaExomeIntersection = i_rnaIds.intersection(i_exomeIds)
    i_rnaDNAUnion = i_rnaWGIntersection.union(i_rnaExomeIntersection)
    
    if (len(i_wgDirList) != 0):
        print >> sys.stderr, "Intersection of RNA-Seq and whole-genome:", len(i_rnaWGIntersection)
    
    if (len(i_exomeDirList) != 0):
        print >> sys.stderr, "Intersection of RNA-Seq and exome-capture:", len(i_rnaExomeIntersection)
    
    print >> sys.stderr, "Union of RNA-Seq and DNA:", len(i_rnaDNAUnion)
    
    if (i_outputFilename != None):
        fileHandler = open(i_outputFilename, "w")
        fileHandler.write(",".join(i_rnaDNAUnion))   
        fileHandler.close()
    elif (len(i_rnaDNAUnion) != 0):
        print >> sys.stdout, ",".join(i_rnaDNAUnion)
        
    # print out all rna-seq ids
    if (i_rnaFilename != None):
        fileHandler = open(i_rnaFilename, "w")
        fileHandler.write(",".join(i_rnaIds))   
        fileHandler.close()

    # print out all wg and exome ids
    if (i_dnaFilename != None):
        # if both sets aren't empty
        if (i_wgIds and i_exomeIds):
            i_dnaIds = i_wgIds.union(i_exomeIds)
        # if wg is not empty
        elif (i_wgIds):
            i_dnaIds = i_wgIds
        else:
            i_dnaIds = i_exomeIds
            
        fileHandler = open(i_dnaFilename, "w")
        fileHandler.write(",".join(i_dnaIds))   
        fileHandler.close()
        
    return

main()
sys.exit(0)
