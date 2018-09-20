#!/usr/bin/env python

import sys
import time
from optparse import OptionParser
import radiaUtil
import os
import logging
import gzip
import collections


'''
'    RNA and DNA Integrated Analysis (RADIA) identifies RNA and DNA variants in NGS data.
'    Copyright (C) 2010-2018  Amie Radenbaugh
'
'    This program is free software: you can redistribute it and/or modify
'    it under the terms of the GNU Affero General Public License as
'    published by the Free Software Foundation, either version 3 of the
'    License, or (at your option) any later version.
'
'    This program is distributed in the hope that it will be useful,
'    but WITHOUT ANY WARRANTY; without even the implied warranty of
'    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
'    GNU Affero General Public License for more details.
'
'    You should have received a copy of the GNU Affero General Public License
'    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''


def get_read_fileHandler(aFilename):
    '''
    ' Open aFilename for reading and return
    ' the file handler.  The file can be 
    ' gzipped or not.
    '''
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename,'rb')
    else:
        return open(aFilename,'r')


def get_write_fileHandler(aFilename):
    '''
    ' Open aFilename for writing and return
    ' the file handler.  The file can be 
    ' gzipped or not.
    '''
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename,'wb')
    else:
        return open(aFilename,'w')
    
    
def merge_filters(anRnaFilterColumn, aDnaFilterColumn):
    # merge the filters for the FILTER column
    # get the rna filters
    rnaFilterSet = set(anRnaFilterColumn.split(";"))
    # get the original filters
    dnaFilterSet = set(aDnaFilterColumn.split(";"))
    # merge the filters
    finalFilterSet = dnaFilterSet.union(rnaFilterSet)
    return ";".join(finalFilterSet)
                            

def merge_mod_filters(anRnaInfoColumn, aDnaInfoColumn):
    
    # merge the mod filters and filter types in the INFO column
    # get the RNA filters and filter types
    rnaInfoList = anRnaInfoColumn.split(";")
    rnaModFilters = []
    rnaModFilterTypes = []
    for info in rnaInfoList:
        keyValueList = info.split("=")
        # some keys are just singular without a value (e.g. DB, SOMATIC, etc.)
        if (len(keyValueList) == 1):
            continue
        elif keyValueList[0] == "MF":
            rnaModFilters = keyValueList[1].split(",")
        elif keyValueList[0] == "MFT":
            rnaModFilterTypes = keyValueList[1].split(",")
        
        # if we found the RNA filters and types, then break
        if rnaModFilters and rnaModFilterTypes:
            break
    
    # get the DNA filters and filter types
    dnaInfoList = aDnaInfoColumn.split(";")
    dnaInfoDict = collections.defaultdict(list)
    for info in dnaInfoList:
        keyValueList = info.split("=")
        # some keys are just singular without a value (e.g. DB, SOMATIC, etc.)
        if (len(keyValueList) == 1):
            dnaInfoDict[keyValueList[0]] = ["True"]
        else:
            # the value can be a comma separated list
            dnaInfoDict[keyValueList[0]] = keyValueList[1].split(",")
            
    # add the RNA filters and filter types to the DNA info
    dnaInfoDict["MF"] += rnaModFilters
    dnaInfoDict["MFT"] += rnaModFilterTypes
    
    # create the INFO field string
    infoField = ""
    for key in sorted(dnaInfoDict.iterkeys()):
        if (len(dnaInfoDict[key]) == 0):
            continue
        elif ("True" in dnaInfoDict[key]):
            infoField += key + ";"
        else:    
            infoField += key + "=" + ",".join(dnaInfoDict[key]) + ";"

    return infoField.rstrip(";")
    
    
def set_sst_field(anInfoField):
    if "GERM" in anInfoField:
        anInfoField = anInfoField.replace("START", "SST=RADIAGerm;START")
    elif "SOM" in anInfoField:
        if "ORIGIN=DNA,RNA" in anInfoField:
            anInfoField = anInfoField.replace("START", "SST=RADIASomRNAConf;START")
        elif "ORIGIN=DNA" in anInfoField:
            anInfoField = anInfoField.replace("START", "SST=RADIASomDNA;START")
        elif "ORIGIN=RNA" in anInfoField:
            anInfoField = anInfoField.replace("START", "SST=RADIASomRNARescue;START")
    elif "TUM_EDIT" in anInfoField:
        anInfoField = anInfoField.replace("START", "SST=RADIATumEdit;START")
    elif "RNA_TUM_VAR" in anInfoField:
        anInfoField = anInfoField.replace("START", "SST=RADIATumRNAVar;START")
    elif "NOR_EDIT" in anInfoField:
        anInfoField = anInfoField.replace("START", "SST=RADIANorEdit;START")
    elif "RNA_NOR_VAR" in anInfoField:
        anInfoField = anInfoField.replace("START", "SST=RADIANorRNAVar;START")
    
    return anInfoField


def merge_vcf_data(aDnaFile, anRnaFile, anOverlapsFile, aNonOverlapsFile, aDnaHeaderOnlyFlag, anIsDebug):
    
    # open the header file
    dnaFileHandler = get_read_fileHandler(aDnaFile)
    rnaFileHandler = get_read_fileHandler(anRnaFile)
    overlapsFileHandler = get_read_fileHandler(anOverlapsFile)
    if (os.path.isfile(aNonOverlapsFile)):
        nonOverlapsFileHandler = get_read_fileHandler(aNonOverlapsFile)
    
    headerList = list()
    coordinateDict = dict()
    
    # the dna file has the results from the dna mpileup filter
    # the rna file has the results from the rna mpileup filter
    # calls that pass in both the DNA and RNA will be in the overlaps file
    # calls that don't pass in the DNA but pass in the RNA are in the non-overlaps file - these are the RNA Rescue and RNA Editing calls
    
    # process all of the calls from the DNA mpileup filter
    for line in dnaFileHandler:
          
        # here are some examples of .vcf lines:
        # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31;DP=53;FA=0.04;INS=0;DEL=0;;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP
        # GT:DP:AD:AF:INS:DEL:START:STOP:BQ:SB      0/0:2:2,0:1.0,0.0:0:0:0:0::36,0:0.0,0.0      0/0:1:1,0:1.0,0.0:0:0:0:0:39,0:1.0,0.0      0/1:50:48,2:0.96,0.04:0:0:2:0:32,18:0.75,0.5
        
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug and not line.startswith("#")):
            logging.debug("DNA mpileup Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # if it is a header line, then add it to the header list
        elif (line.startswith("#")): 
            
            # keep all the header lines
            headerList.append(line + "\n")
                
        # now we are to the data
        # if we only want the header, then break     
        elif (aDnaHeaderOnlyFlag):
            break
        # if we want the DNA data then process it
        else:
            # split the line on the tab
            splitLine = line.split("\t")
    
            # the coordinate is the second element
            stopCoordinate = splitLine[1]
            coordinateDict[stopCoordinate] = line + "\n"
            
    # these are all the calls that pass in both the DNA and RNA   
    for line in overlapsFileHandler:
          
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug and not line.startswith("#")):
            logging.debug("Overlaps file Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # if it is a header line, then continue
        elif (line.startswith("#")): 
            continue;
        
        # now we are to the data    
        else:
            # split the line on the tab
            splitLine = line.split("\t")
    
            # the coordinate is the second element
            stopCoordinate = splitLine[1]
            
            # if the call passed in both the RNA and DNA, then adjust the origin
            if (stopCoordinate in coordinateDict):
                dnaLine = coordinateDict[stopCoordinate]
                if (anIsDebug):
                    logging.debug("passed in both RNA and DNA (from overlaps file) changing the origin to (DNA,RNA) \nDNALine: %s RNALine: %s\n", dnaLine, line)
                dnaLine = dnaLine.replace("ORIGIN=DNA", "ORIGIN=DNA,RNA")
                coordinateDict[stopCoordinate] = dnaLine
            else:
                coordinateDict[stopCoordinate] = line + "\n"
        
      
    # loop through the RNA mpileup filtered calls
    # create 2 dictionaries:  one for passing, one for non-passing
    #
    # if an RNA Rescue or RNA Editing call passes in the anRnaNonOverlapsFile below, 
    # then we want to use the original RNA mpileup passing call to overwrite the DNA call 
    # the non-overlaps file is really the RNA mpileup passing calls that are
    # first filtered by DNA, then grep, and then blat. the filtered by DNA part doesn't
    # select one modType when no call passes, so the final passing call has
    # more than one modType which causes problems in the next filter, therefore
    # use the RNA mpileup passing call.
    #
    # mpileup_rna_origin:
    # 9       17464495        .       G       A       0.0     PASS    
    #    AC=5;AF=0.1;AN=2;BQ=39;DP=49;FA=0.1;INS=0;DEL=0;;MC=G>A;MT=TUM_EDIT;NS=3;ORIGIN=RNA;SB=0.73;SS=5;START=0;STOP=5;VT=SNP
    #    GT:DP:AD:AF:INS:DEL:START:STOP:BQ:SB
    #    0/0:22:22,0:1.0,0.0:0:0:0:1:32,0:0.68,0.0
    #    0/0:10:9,1:0.9,0.1:0:0:0:0:29,13:0.67,1.0
    #    0/1:17:13,4:0.76,0.24:0:0:0:4:60,32:0.92,0.5
    # vs.
    # mpileup_rna_origin->dnaFiltered->blat:
    # 9       17464495        .       G       A       0.0     PASS    
    #    AC=5;AF=0.1;AN=2;BQ=39;DP=49;FA=0.1;INS=0;DEL=0;;MC=G>A,G>A;MF=rnacall,dtmnab_dtmnbq;
    #    MFT=DNA_TUM_EDIT_G>A,DNA_SOM_G>A;MT=TUM_EDIT,SOM;NS=3;ORIGIN=RNA;SB=0.73;SS=5;START=0;STOP=5;VT=SNP     
    #    GT:DP:AD:AF:INS:DEL:START:STOP:BQ:SB
    #    0/0:22:22,0:1.0,0.0:0:0:0:1:32,0:0.68,0.0
    #    0/0:10:9,1:0.9,0.1:0:0:0:0:29,13:0.67,1.0
    #    0/1:17:13,4:0.76,0.24:0:0:0:4:60,32:0.92,0.5
    #
    # when merging a call that passes in the non-overlaps (dnaFiltered or blat) file,
    # replace the DNA call, with the original RNA mpileup passing call
    # 
    # the non-passing dictionary will be used below to help merge filtered calls 
    # when a call gets filtered by both the RNA and DNA
    
    rnaMpileupPassingDict = {}
    rnaMpileupNonpassingDict = {}
    for rnaLine in rnaFileHandler:
          
        # strip the carriage return and newline characters
        rnaLine = rnaLine.rstrip("\r\n")
        
        if (anIsDebug and not rnaLine.startswith("#")):
            logging.debug("RNA mpileup Line: %s", rnaLine)
            
        # if it is an empty line, then just continue
        if (rnaLine.isspace()):
            continue;
        
        # if it is a header line, then just continue
        elif (rnaLine.startswith("#")): 
            continue;
        
        # now we are to the data    
        else:
            # split the line on the tab
            rnaLineSplit = rnaLine.split("\t")
    
            # the coordinate is the second element
            stopCoordinate = rnaLineSplit[1]
            
            # put the call in the right dict
            if "PASS" in rnaLineSplit[6]:
                rnaMpileupPassingDict[stopCoordinate] = rnaLine
            else:
                rnaMpileupNonpassingDict[stopCoordinate] = rnaLine
                
    # these are the RNA Rescue and RNA Editing calls after the initial filtering but before filterByReadSupport.py
    if (os.path.isfile(aNonOverlapsFile)):
        for line in nonOverlapsFileHandler:
          
            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")
        
            if (anIsDebug and not line.startswith("#")):
                logging.debug("Non-overlaps Line: %s", line)
            
            # if it is an empty line, then just continue
            if (line.isspace()):
                continue;
        
            # if it is a header line, then just continue
            elif (line.startswith("#")): 
                continue;
        
            # now we are to the data    
            else:
                # split the line on the tab
                splitLine = line.split("\t")
    
                # the coordinate is the second element
                stopCoordinate = splitLine[1]
                
                # if this call passed in the RNA, then overwrite the DNA call that didn't pass
                if ("PASS" in splitLine[6]):
                    # if this call existed in the DNA
                    if (stopCoordinate in coordinateDict):
                        dnaLine = coordinateDict[stopCoordinate]
                        # get the RNA line from the RNA mpileups passing dict
                        rnaLine = rnaMpileupPassingDict[stopCoordinate]
                        # if it didn't pass in the DNA
                        if ("PASS" not in dnaLine):
                            if (anIsDebug):
                                logging.debug("Overwriting non-passing DNA call with passing RNA Rescue calls \nDNALine: %sRNALineNonOverlaps: %s\nRNALineMpileup: %s\n", dnaLine, line, rnaLine)
                            coordinateDict[stopCoordinate] = rnaLine + "\n"
                        else:
                            # this call passed in both
                            logging.warning("Unusual call in non-overlaps file passed in both the RNA and DNA but they probably don't have the same modType! \nDNALine: %sRNALineNonOverlaps: %s\nRNALineMpileup: %s\n", dnaLine, line, rnaLine)
                            # at this point, there are multiple events that pass all the filters
                            # in this case, pick the passing event in the following order:  GERM, NOR_EDIT, SOM, TUM_EDIT, RNA_TUM_VAR, LOH
                            if ("GERM" in dnaLine or "SOM" in dnaLine):
                                coordinateDict[stopCoordinate] = dnaLine
                            else:
                                coordinateDict[stopCoordinate] = rnaLine
                    # this call didn't exist in the DNA
                    else:
                        logging.warning("Call didn't exist in DNA? RNALine: %s\n", line)
                        coordinateDict[stopCoordinate] = line + "\n"
                # this call didn't pass in the RNA
                else:
                    if (anIsDebug):
                        logging.debug("Call didn't pass in RNA: RNALine: %s\n", line)
                    
                    # if this call existed in the DNA
                    if (stopCoordinate in coordinateDict):
                        dnaLine = coordinateDict[stopCoordinate]
                        # if it didn't pass in the DNA
                        if ("PASS" not in dnaLine):
                            if (anIsDebug):
                                logging.debug("RNANoPass:  Didn't pass in both, so change origin and merge filters \nDNALine: %s RNALine: %s\n", dnaLine, line)
                            # change origin
                            if ("ORIGIN=DNA,RNA" not in dnaLine):
                                dnaLine = dnaLine.replace("ORIGIN=DNA", "ORIGIN=DNA,RNA")
                            dnaLine = dnaLine.rstrip("\r\n")
                            dnaLineSplit = dnaLine.split("\t")
                            
                            # merge the filters for the FILTER column
                            dnaLineSplit[6] = merge_filters(splitLine[6], dnaLineSplit[6])
                            
                            # merge the mod filters and filter types in the INFO column
                            dnaLineSplit[7] = merge_mod_filters(splitLine[7], dnaLineSplit[7])
                            
                            coordinateDict[stopCoordinate] = "\t".join(dnaLineSplit) + "\n"
                            if (anIsDebug):
                                logging.debug("RNANoPass:  After change origin and merge filters \nFinalLine: %s\n", "\t".join(dnaLineSplit))
                        else:
                            # this call passed in both:
                            # DNALine: 17 4857042 .   T   A,G,C   0.0 PASS    
                            #    AB=A,G,C;AC=10,5,8211;AF=0.0,0.0,0.98;AN=4;BQ=31;DP=8379;FA=0.98;INS=0;DEL=0;;MC=T>A;MT=GERM;NS=3;
                            #    ORG_ISO_AD=16_2_1_2615,18_3_1_2791,18_1_2_2805;ORIGIN=DNA;RS_GEN_POS=17:4854383-4860426,17:4854383-4860426,17:4854383-4860426;
                            #    RS_NAME=NM_001193503,NM_001976,NM_053013;RS_ORG_POS=313,484,442;RS_STRAND=+,+,+;SB=0.74;SS=1;START=1;STOP=0;VT=SNP    
                            #    GT:DP:AD:AF:INS:DEL:START:STOP:BQ:SB
                            #    0/1:36:31,4,1,0:0.86,0.11,0.03,0.0:0:0:1:0:29,28,3,0:0.39,0.5,1.0,0.0
                            #    0/0:70:70,0,0,0:1.0,0.0,0.0,0.0:0:0:0:0:31,0,0,0:0.56,0.0,0.0,0.0
                            #    3/3:8273:52,6,4,8211:0.01,0.0,0.0,0.99:0:0:0:0:45,12,12,58:0.94,1.0,1.0,0.97
                            # RNALine: 17    4857042 .   T   A,G,C   0.0 PASS    
                            #    AB=A,G,C;AC=10,5,8211;AF=0.0,0.0,0.98;AN=4;BQ=31;DP=8379;FA=0.98;INS=0;DEL=0;;MC=T>C;MT=TUM_EDIT;NS=3;
                            #    ORG_ISO_AD=16_2_1_2615,18_3_1_2791,18_1_2_2805;ORIGIN=RNA;RS_GEN_POS=17:4854383-4860426,17:4854383-4860426,17:4854383-4860426;
                            #    RS_NAME=NM_001193503,NM_001976,NM_053013;RS_ORG_POS=313,484,442;RS_STRAND=+,+,+;SB=0.74;SS=5;START=1;STOP=0;VT=SNP    
                            #    GT:DP:AD:AF:INS:DEL:START:STOP:BQ:SB
                            #    0/1:36:31,4,1,0:0.86,0.11,0.03,0.0:0:0:1:0:29,28,3,0:0.39,0.5,1.0,0.0
                            #    0/0:70:70,0,0,0:1.0,0.0,0.0,0.0:0:0:0:0:31,0,0,0:0.56,0.0,0.0,0.0
                            #    3/3:8273:52,6,4,8211:0.01,0.0,0.0,0.99:0:0:0:0:45,12,12,58:0.94,1.0,1.0,0.97
                            logging.warning("RNANoPass:  Call passed in both RNA and DNA but they probably don't have the same modType \nDNALine: %s RNALine: %s\n", dnaLine, line)
                            # at this point, there are multiple events that pass all the filters
                            # in this case, pick the passing event in the following order:  GERM, NOR_EDIT, SOM, TUM_EDIT, RNA_TUM_VAR, LOH
                            if ("GERM" in dnaLine or "SOM" in dnaLine):
                                coordinateDict[stopCoordinate] = dnaLine
                            else:
                                coordinateDict[stopCoordinate] = line
                    # this call didn't exist in the DNA
                    else:
                        logging.warning("RNANoPass:  Call didn't exist in DNA? RNALine: %s\n", line)
                        coordinateDict[stopCoordinate] = line + "\n"
       
    # these are needed for merging the RNA mpileup filters
    for (rnaStopCoordinate, rnaLine) in rnaMpileupNonpassingDict.iteritems():
          
        if (anIsDebug and not rnaLine.startswith("#")):
            logging.debug("RNA mpileup non-passing Line: %s", rnaLine)
            
        # split the line on the tab
        rnaLineSplit = rnaLine.split("\t")
        
        # get the original line
        dnaLine = coordinateDict[rnaStopCoordinate]
        dnaLine = dnaLine.rstrip("\r\n")
        dnaLineSplit = dnaLine.split("\t")    
        
        # if the call didn't pass in the RNA or DNA, we want to merge the filters
        if "PASS" not in dnaLineSplit[6]:
            if (anIsDebug):
                logging.debug("Merging filters for \nDNALine: %s \nRNALine: %s", dnaLine, rnaLine)
            
            # merge the filters for the FILTER column
            dnaLineSplit[6] = merge_filters(rnaLineSplit[6], dnaLineSplit[6])
                        
            # merge the mod filters and filter types in the INFO column
            dnaLineSplit[7] = merge_mod_filters(rnaLineSplit[7], dnaLineSplit[7])            
            
            finalLine = "\t".join(dnaLineSplit)
            if ("ORIGIN=DNA,RNA" not in finalLine):
                finalLine = finalLine.replace("ORIGIN=DNA", "ORIGIN=DNA,RNA")
                
            coordinateDict[rnaStopCoordinate] = finalLine + "\n"
            if (anIsDebug):
                logging.debug("Merged filters \nFinalLine: %s", finalLine)
    
    dnaFileHandler.close()
    rnaFileHandler.close()
    overlapsFileHandler.close()
    if (os.path.isfile(aNonOverlapsFile)):
        nonOverlapsFileHandler.close()
    
    return (headerList, coordinateDict)

               
def main():
    
    #python mergeRnaAndDnaFiles.py TCGA-AB-2995 5 ../data/test/TCGA-AB-2995_dnaFile.vcf ../data/test/TCGA-AB-2995_rnaFile.vcf ../data/test/TCGA-AB-2995_rnaFile.vcf ../data/test/
       
    startTime = time.time()
    
    # create the usage statement
    usage = "usage: python %prog id chrom dnaFile rnaFile rnaOverlapsFile rnaNonOverlapsFile outputFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    i_cmdLineParser.add_option("-d", "--dnaHeaderOnly", action="store_true", default=False, dest="dnaHeaderOnly", help="include this argument if only the DNA header should be included")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(6,15,1)
    i_argLength = len(sys.argv)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)
    
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_id = i_cmdLineArgs[0]
    i_chrom = i_cmdLineArgs[1]
    i_dnaFilename = i_cmdLineArgs[2]
    i_rnaFilename = i_cmdLineArgs[3]
    i_overlapsFilename = i_cmdLineArgs[4]
    i_nonOverlapsFilename = i_cmdLineArgs[5]
    i_outputFilename = i_cmdLineArgs[6]
    
    # get the optional params with default values
    i_logLevel = i_cmdLineOptions.logLevel
    i_dnaHeaderOnly = i_cmdLineOptions.dnaHeaderOnly
    
    i_logFilename = None
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
    
    # assuming loglevel is bound to the string value obtained from the
    # command line argument. Convert to upper case to allow the user to
    # specify --log=DEBUG or --log=debug
    i_numericLogLevel = getattr(logging, i_logLevel.upper(), None)
    if not isinstance(i_numericLogLevel, int):
        raise ValueError("Invalid log level: '%s' must be one of the following:  DEBUG, INFO, WARNING, ERROR, CRITICAL", i_logLevel)
    
    # set up the logging
    if (i_logFilename != None):
        logging.basicConfig(level=i_numericLogLevel, filename=i_logFilename, filemode='w', format='%(asctime)s\t%(levelname)s\t%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=i_numericLogLevel, format='%(asctime)s\t%(levelname)s\t%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
        
    # set the debug    
    i_debug = (i_numericLogLevel == logging.DEBUG)
    
    # output some debug info
    if (i_debug):
        logging.debug("id=%s", i_id)
        logging.debug("chrom=%s", i_chrom)
        logging.debug("dnaFilename=%s", i_dnaFilename)
        logging.debug("rnaFilename=%s", i_rnaFilename)
        logging.debug("overlapsFilename=%s", i_overlapsFilename)
        logging.debug("nonOverlapsFilename=%s", i_nonOverlapsFilename)
        logging.debug("outputFilename=%s", i_outputFilename)
                    
    # check for any errors
    i_readFilenameList = [i_dnaFilename, i_rnaFilename, i_overlapsFilename]
    i_writeFilenameList = [i_outputFilename]
    i_dirList = None
    
    if (not radiaUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)
                
    # get the VCF generator
    (headerList, coordinateDict) = merge_vcf_data(i_dnaFilename, i_rnaFilename, i_overlapsFilename, i_nonOverlapsFilename, i_dnaHeaderOnly, i_debug)
    
    outputFileHandler = get_write_fileHandler(i_outputFilename)
    
    for headerLine in headerList:
        outputFileHandler.write(headerLine)

    numericKeys = coordinateDict.keys()
    numericKeys.sort(key=int)
    for coordinate in numericKeys:
        line = coordinateDict[coordinate]
        line = line.rstrip("\r\n")
        # split the line on the tab
        splitLine = line.split("\t")

        # set the SST field in the INFO
        splitLine[7] = set_sst_field(splitLine[7])
        outputFileHandler.write("\t".join(splitLine) + "\n")
    
    stopTime = time.time() 
    logging.info("Total time for Id %s: Total time=%s hrs, %s mins, %s secs", i_id, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime))    
    # close the files 
    outputFileHandler.close()
    
    return
 

main()    
sys.exit(0)
