#!/usr/bin/env python

import sys                          # system module
from optparse import OptionParser   # used for parsing command line arguments
import radiaUtil                    # utility functions for rna editing
import logging
import time
import collections
import gzip


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
    

def get_append_fileHandler(aFilename):
    '''
    ' Open aFilename for appending and return
    ' the file handler.  The file can be 
    ' gzipped or not.
    '''
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename,'ab')
    else:
        return open(aFilename,'a')


def get_vcf_data(anInputFilename, aStatsDict, aCompareDict, aPrefix, anIsDebug):
    '''
    ' The .vcf files must have at least 10 fields:  chromosome, coordinate, id
    ' references, alts, quality score, filters, infos, format, and summary info
    ' for at least one .bam file.
    '
    ' anInputFileHandler: The input stream for the file
    ' aStatsDict: A dictionary holding stats about all the comparisons
    ' aCompareDict: The key,value pairs that should be used in the comparison
    ' aPrefix: "rad" for the RADIA files, "cmp" for the compare files
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    inputFileHandler = get_read_fileHandler(anInputFilename)
    outputDict = {}
     
    for line in inputFileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        #if (anIsDebug):
        #    logging.debug("VCF Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
                
        # these lines are from previous scripts in the pipeline, so output them    
        elif (line.startswith("#")):
            continue;
        
        # now we are to the data
        else:    
            
            # split the line on the tab
            splitLine = line.split("\t")
            
            # get the fields to yield
            #columnHeaders = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            chrom = splitLine[0]
            stopCoordinate = splitLine[1]
             
            #outputDict[chrom + "_" + stopCoordinate] = line
                    
            # keep track of the number of total events per file 
            #aStatsDict[aPrefix + "_events"] += 1
            
            #if ("PASS" in line):
            #    aStatsDict[aPrefix + "_pass_events"] += 1
                    
            # the thing that is being compared to has to have the smaller/limited amount
            # i.e. only the passing som events, otherwise everything will be found
            if (aPrefix == "rad" and "PASS" in line and "SNP" in line and ("SOM" in line or "EDIT" in line or "RNA_TUM_VAR" in line or "RNA_NOR_VAR" in line)):
            #if (aPrefix == "rad" and "SNP" in line):
            #if (aPrefix == "rad"):
                # add the coordinate to the output
                outputDict[chrom + "_" + stopCoordinate] = line
            #elif (aPrefix == "cmp" and "PASS" in line):
            elif (aPrefix == "cmp" and "PASS" in line and "SNP" in line and ("SOM" in line or "EDIT" in line or "RNA_TUM_VAR" in line or "RNA_NOR_VAR" in line)):
            #elif (aPrefix == "cmp" and "SNP" in line):
                outputDict[chrom + "_" + stopCoordinate] = line 
            
            #if ("PASS" in line and "Somatic" in line and "SNP" in line):
            #if ("SOM" in line):
            #    outputDict[chrom + "_" + stopCoordinate] = line
                    
            #    # keep track of the number of total events per file 
            #    aStatsDict[aPrefix + "_events"] += 1
                
            #    if ("PASS" in line):
            #        aStatsDict[aPrefix + "_pass_events"] += 1
            
            #if ("PASS" in line and "SOM" in line):
            #if ("SOM" in line):
            #    outputDict[chrom + "_" + stopCoordinate] = line 
            #    aStatsDict[aPrefix + "_events"] += 1
            #    aStatsDict[aPrefix + "_pass_events"] += 1
            
            # keep track of the total number of comparison events (blck, dnSnp, etc.) per file
            # there can be multiple keys for one filter such as blq and bldp for blacklists
            for (radKeyString, cmpKeyString) in aCompareDict.iteritems():
                if (aPrefix == "rad"):
                    # break up the string to get the individual keys
                    radKeyList = radKeyString.split(",")
                    # search for each one of them
                    for radKey in radKeyList:
                        # if we find one
                        if (radKey in line):
                            # count it using the keyString
                            aStatsDict[aPrefix + "_" + radKeyString] += 1
                            # if this is a passing line, call it using the keyString
                            #if ("PASS" in line and ((radKey == "GERM") or ("DB" not in line))):
                            #if ("PASS" in line and ((radKey == "Germline") or ("DB" not in line))):
                            if ("PASS" in line and "SNP" in line and radKey in line):
                            #if ("PASS" in line):
                            #if (radKey in line):
                            #if ("PASS" in line and "SNP" in line):
                                aStatsDict[aPrefix + "_pass_" + radKeyString] += 1
                            # only count it once
                            break;
                        
                elif (aPrefix == "cmp"):
                    # break up the string to get the individual keys
                    cmpKeyList = cmpKeyString.split(",")
                    # search for each one of them
                    for cmpKey in cmpKeyList:
                        # if we find one
                        if (cmpKey in line):
                            # count it using the keyString
                            #if ("SNP" in line):
                            #    aStatsDict[aPrefix + "_" + cmpKeyString] += 1
                            aStatsDict[aPrefix + "_" + cmpKeyString] += 1
                            # if this is a passing line, call it using the keyString
                            #if ("PASS" in line and ((cmpKey == "GERM") or ("DB" not in line))):
                            #if ("PASS" in line and "SNP" in line and ((cmpKey == "Germline") or ("DB" not in line))):
                            #if ("PASS" in line and "SNP" in line and "SS=2" in line and cmpKey in line):
                            if ("PASS" in line and "SNP" in line and cmpKey in line):
                            #if ("PASS" in line):
                                aStatsDict[aPrefix + "_pass_" + cmpKeyString] += 1
                            # only count it once
                            break;
                        
    inputFileHandler.close()
        
    return (outputDict, aStatsDict)


def get_maf_data(anInputFilename, aStatsDict, aCompareDict, aPrefix, anIsDebug):
    '''
    ' The .vcf files must have at least 10 fields:  chromosome, coordinate, id
    ' references, alts, quality score, filters, infos, format, and summary info
    ' for at least one .bam file.
    '
    ' anInputFileHandler: The input stream for the file
    ' aStatsDict: A dictionary holding stats about all the comparisons
    ' aCompareDict: The key,value pairs that should be used in the comparison
    ' aPrefix: "rad" for the RADIA files, "cmp" for the compare files
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    inputFileHandler = get_read_fileHandler(anInputFilename)
    outputDict = {}
     
    for line in inputFileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        #if (anIsDebug):
        #    logging.debug("MAF Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
                
        # these lines are from previous scripts in the pipeline, so output them    
        elif (line.startswith("#")):
            continue;
        
        # now we are to the data
        else:    
            
            # split the line on the tab
            splitLine = line.split("\t")
            
            # get the fields to yield
            #center = splitLine[2]
            chrom = splitLine[4]
            #startCoordinate = splitLine[5]
            stopCoordinate = splitLine[6]
            #variantType = splitLine[9]
            #dbSnp = splitLine[13]
            
            #if ("Somatic" in line and "SNP" in line):
            if (True):
                #if (chrom + "_" + stopCoordinate) in outputDict:
                #    logging.debug(line + outputDict[chrom + "_" + stopCoordinate])
                
                # add the coordinate to the output
                outputDict[chrom + "_" + stopCoordinate] = line
                
                # keep track of the number of total events per file 
                aStatsDict[aPrefix + "_events"] += 1
                
                # all events are considered passing events
                aStatsDict[aPrefix + "_pass_events"] += 1
                
            # keep track of the total number of comparison events (blck, dnSnp, etc.) per file
            # their can be multiple keys for one filter such as blq and bldp for blacklists
            for (radKeyString, cmpKeyString) in aCompareDict.iteritems():
                if (aPrefix == "rad"):
                    # break up the string to get the individual keys
                    radKeyList = radKeyString.split(",")
                    # search for each one of them
                    for radKey in radKeyList:
                        # if we find one
                        if (radKey in line):
                            # count it using the keyString
                            aStatsDict[aPrefix + "_" + radKeyString] += 1
                            # if this is a passing line, call it using the keyString
                            #if ("PASS" in line and ((radKey == "GERM") or ("DB" not in line))):
                            #if ("PASS" in line and ((radKey == "Germline") or ("DB" not in line))):
                            #if ("PASS" in line and radKey in line):
                            #if ("SNP" in line):
                            if ("SOMATIC" in line):
                                aStatsDict[aPrefix + "_pass_" + radKeyString] += 1
                            # only count it once
                            break;
                        
                elif (aPrefix == "cmp"):
                    # break up the string to get the individual keys
                    cmpKeyList = cmpKeyString.split(",")
                    # search for each one of them
                    for cmpKey in cmpKeyList:
                        # if we find one
                        if (cmpKey in line):
                            # count it using the keyString
                            aStatsDict[aPrefix + "_" + cmpKeyString] += 1
                            # if this is a passing line, call it using the keyString
                            #if ("PASS" in line and ((cmpKey == "GERM") or ("DB" not in line))):
                            #if ("PASS" in line and "SNP" in line and ((cmpKey == "Germline") or ("DB" not in line))):
                            #if ("SNP" in line):
                            #if ("PASS" in line):
                            if ("SOMATIC" in line):
                                aStatsDict[aPrefix + "_pass_" + cmpKeyString] += 1
                            # only count it once
                            break;
                        
    inputFileHandler.close()
        
    return (outputDict, aStatsDict)


def get_validation_data(anInputFilename, aStatsDict, aCompareDict, aPrefix, anIsDebug):
    '''
    ' The validation files must have at least 10 fields:  chromosome, coordinate, id
    ' references, alts, quality score, filters, infos, format, and summary info
    ' for at least one .bam file.
    '
    ' anInputFileHandler: The input stream for the file
    ' aStatsDict: A dictionary holding stats about all the comparisons
    ' aCompareDict: The key,value pairs that should be used in the comparison
    ' aPrefix: "rad" for the RADIA files, "cmp" for the compare files
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    inputFileHandler = get_read_fileHandler(anInputFilename)
    outputDict = {}
     
    for line in inputFileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        #if (anIsDebug):
        #    logging.debug("Validation Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
                
        # these lines are from previous scripts in the pipeline, so skip them    
        elif (line.startswith("#")):
            continue;

        # this is a header line, so skip it   
        elif (line.startswith("chrom")):
            continue;
        
        # now we are to the data
        else:    
            
            # split the line on the tab
            splitLine = line.split("\t")
            
            # get the fields to yield
            # columnHeaders = ["chrom", "chr_start", "chr_stop", "ref", "var", "source", "val_result"]
            # these are 0-based
            chrom = splitLine[0]
            #startCoordinate = splitLine[1]
            stopCoordinate = splitLine[2]
            #ref = splitLine[3]
            #variantAllele = splitLine[4]
            #center = splitLine[5]
            #valResult = splitLine[6]
            
            # add the coordinate to the output
            outputDict[chrom + "_" + stopCoordinate] = line
            
            # keep track of the number of total events per file 
            aStatsDict[aPrefix + "_events"] += 1
            
            # all events are considered passing events
            aStatsDict[aPrefix + "_pass_events"] += 1
            
            # keep track of the total number of comparison events (blck, dnSnp, etc.) per file
            # their can be multiple keys for one filter such as blq and bldp for blacklists
            for (radKeyString, cmpKeyString) in aCompareDict.iteritems():
                if (aPrefix == "cmp"):
                    # break up the string to get the individual keys
                    cmpKeyList = cmpKeyString.split(",")
                    # search for each one of them
                    for cmpKey in cmpKeyList:
                        # if we find one
                        if (cmpKey in line):
                            # count it using the keyString
                            aStatsDict[aPrefix + "_" + cmpKeyString] += 1
                            # if this is a passing line, call it using the keyString
                            if (cmpKey in line):
                                aStatsDict[aPrefix + "_pass_" + cmpKeyString] += 1
                            # only count it once
                            break;
                        
    inputFileHandler.close()
        
    return (outputDict, aStatsDict)


def get_simulation_data(anInputFilename, aStatsDict, aCompareDict, aPrefix, anIsDebug):
    '''
    ' The simulation files have 11 fields:  mutation type, chromosome, start, end, target AF, 
    ' mutation position, base change, coverage in, coverage out, actual AF, highest AF of anything 
    ' linked. The useful ones for comparing against RADIA are chromosome, mutation position, and base change.
    '
    ' anInputFileHandler: The input stream for the file
    ' aStatsDict: A dictionary holding stats about all the comparisons
    ' aCompareDict: The key,value pairs that should be used in the comparison
    ' aPrefix: "rad" for the RADIA files, "cmp" for the compare files
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    inputFileHandler = get_read_fileHandler(anInputFilename)
    outputDict = {}
     
    for line in inputFileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        #if (anIsDebug):
        #    logging.debug("Simulation Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
                
        # these lines are from previous scripts in the pipeline, so skip them    
        elif (line.startswith("#")):
            continue;
        
        # now we are to the data
        else:    
            
            # split the line on the tab
            splitLine = line.split("\t")
            
            #mutType = splitLine[0]
            chrom = splitLine[1]
            #startCoordinate = splitLine[2]
            #stopCoordinate = splitLine[3]
            #targetAF = splitLine[4]
            mutPosition = splitLine[5]
            #baseChange = splitLine[6]
            #coverageIn = splitLine[7]
            #coverageOut = splitLine[8]
            #actualAF = splitLine[9]
            #highestAF = splitLine[10]
            
            if (chrom + "_" + mutPosition) in outputDict:
                logging.debug(line + outputDict[chrom + "_" + mutPosition])
                
            # add the coordinate to the output
            outputDict[chrom + "_" + mutPosition] = line
            
            # keep track of the number of total events per file 
            aStatsDict[aPrefix + "_events"] += 1
            
            # all events are considered passing events
            aStatsDict[aPrefix + "_pass_events"] += 1
            
            # keep track of the total number of comparison events (blck, dnSnp, etc.) per file
            # their can be multiple keys for one filter such as blq and bldp for blacklists
            for (radKeyString, cmpKeyString) in aCompareDict.iteritems():
                if (aPrefix == "cmp"):
                    # break up the string to get the individual keys
                    cmpKeyList = cmpKeyString.split(",")
                    # search for each one of them
                    for cmpKey in cmpKeyList:
                        # if we find one
                        if (cmpKey in line):
                            # count it using the keyString
                            aStatsDict[aPrefix + "_" + cmpKeyString] += 1
                            # if this is a passing line, call it using the keyString
                            aStatsDict[aPrefix + "_pass_" + cmpKeyString] += 1
                            # only count it once
                            break;
                        
    inputFileHandler.close()
        
    return (outputDict, aStatsDict)
        
    
def compare_events(aTCGAId, aChrom, aRadiaFilename, aCompareFilename, aStatsFilename, anOverlapFilename, aNonOverlapFilename, aCompareDict, anIsDebug):
    '''
    ' The function compares variants in one file with variants in another file.  This can be used to compare variants from
    ' different methods, MAF files, or validation files.  At a minimum, the coordinates are compared.  The user can also
    ' specify additional comparisons that should be done such as comparing if the call was classified as somatic in both
    ' methods (e.g. SOM=Somatic).  The keys and values can be comma-separated lists.  For example, a call may be labeled
    ' as blacklisted in one file with "blck" and in another file with "blq" or "bldp", then the comparison string would
    ' be blck=blq,bldp.
    '
    ' aTCGAId: The TCGA Id for this sample
    ' aChrom: The chromosome being filtered
    ' aRadiaFilename: A .vcf file from RADIA
    ' aCompareFilename: A file to compare to
    ' aStatsFilename: A stats file
    ' anOverlapFilename: A file where all the overlaps are output
    ' aNonOverlapFilename: A file where all the non-overlaps are output
    ' aCompareDict: A dictionary of key=value to be compare (coordinate is always compared)
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    # create the generators for the filter and vcf files
    i_statsDict = collections.defaultdict(int)
    i_filterDict = collections.defaultdict(int)
    (i_radDict, i_statsDict) = get_vcf_data(aRadiaFilename, i_statsDict, aCompareDict, "rad", anIsDebug)
    (i_cmpDict, i_statsDict) = get_vcf_data(aCompareFilename, i_statsDict, aCompareDict, "cmp", anIsDebug)
    #(i_radDict, i_statsDict) = get_maf_data(aRadFilename, i_statsDict, aCompareDict, "rad", anIsDebug)
    #(i_cmpDict, i_statsDict) = get_maf_data(aCompareFilename, i_statsDict, aCompareDict, "cmp", anIsDebug)
    #(i_cmpDict, i_statsDict) = get_validation_data(aCompareFilename, i_statsDict, aCompareDict, "cmp", anIsDebug)
    #(i_cmpDict, i_statsDict) = get_simulation_data(aCompareFilename, i_statsDict, aCompareDict, "cmp", anIsDebug)
    
    if (anOverlapFilename != None):
        overlapFileHandler = get_write_fileHandler(anOverlapFilename)
    if (aNonOverlapFilename != None):
        nonOverlapFileHandler = get_write_fileHandler(aNonOverlapFilename)

    # initialize some variables
    startTime = time.time()
        
    # for each cmp event
    for (cmpCoordinate, cmpLine) in i_cmpDict.iteritems():
        
        # this one is for comparing blacklist results
        #if ("SNP" in cmpLine and ("bldp" in cmpLine or "blq" in cmpLine) and cmpCoordinate not in i_radDict):
        # this one is for comparing BB, Radia, or Maf results
        if ("PASS" in cmpLine and ("SOM" in cmpLine or "EDIT" in cmpLine or "RNA_TUM_VAR" in cmpLine or "RNA_NOR_VAR" in cmpLine) and cmpCoordinate not in i_radDict):
        #if ("SNP" in cmpLine and "Somatic" in cmpLine and cmpCoordinate not in i_radDict):
        # this one is for validation data
        #if (cmpCoordinate not in i_radDict):
        #if ("PASS" in cmpLine and "SNP" in cmpLine and "SOM" in cmpLine and cmpCoordinate not in i_radDict):
        #if ("Somatic" in cmpLine and "SNP" in cmpLine and cmpCoordinate not in i_radDict):
            if (anIsDebug):
                logging.debug("no radia call %s", cmpLine)
            
            #if (aNonOverlapFilename != None):
            #    nonOverlapFileHandler.write(cmpLine + "\n")
            
            # add to maf
            #if (anOverlapFilename != None):
            #    overlapFileHandler.write(cmpLine + "\n")
            
    # for each rad event
    for (radCoordinate, radLine) in i_radDict.iteritems():
        
        #if (("bldp" in radLine or "blq" in radLine) and radCoordinate not in i_cmpDict):
        #if ("PASS" in radLine and "SNP" in radLine and radCoordinate not in i_cmpDict):
        #if (radCoordinate not in i_cmpDict):
        #if ("PASS" in radLine and "SOM" in radLine):
        #if ("SOM" in radLine and radCoordinate not in i_cmpDict):
        if ("PASS" in radLine and ("SOM" in radLine or "EDIT" in radLine or "RNA_TUM_VAR" in radLine or "RNA_NOR_VAR" in radLine) and "SNP" in radLine and radCoordinate not in i_cmpDict):
        #if ("PASS" in radLine and "SNP" in radLine and "Somatic" in radLine and radCoordinate not in i_cmpDict):
        #if ("SNP" in radLine and "Somatic" in radLine and radCoordinate not in i_cmpDict):
        #if ("PASS" in radLine and "SOM" in radLine and radCoordinate not in i_cmpDict):
        #if ("SOM" in radLine and radCoordinate not in i_cmpDict):
            if (anIsDebug):
                logging.debug("new radia call %s", radLine)
            
            if (aNonOverlapFilename != None):
                nonOverlapFileHandler.write(radLine + "\n")
            
            # add to maf
            #if (anOverlapFilename != None):
                #caller = "ucsc;"
                #if ("radia" in radLine):
                #    caller += "radia;"
                #if ("bambam" in radLine):
                #    caller += "bambam;"
                
                #caller = "rnaCall;"
                # split the line on the tab
                #splitLine = radLine.split("\t")
                #chrom = splitLine[0]
                #stopCoordinate = int(splitLine[1])
                #startCoordinate = stopCoordinate-1
                #output = ["gene", "score", caller, "score", chrom, str(startCoordinate), str(stopCoordinate), "+", "mutClass", "SNP", "Somatic"]
                #overlapFileHandler.write("\t".join(output) + "\n")
                
            
        # if the coordinates overlap, then count them
        if (radCoordinate in i_cmpDict):
                
            i_statsDict["overlap_events"] += 1
            compareLine = i_cmpDict[radCoordinate]
            
            # this one is for BB and Maf comparisons
            #if ("PASS" in radLine and "SNP" in radLine):
            # this one is for Radia to Radia comparisons
            #if ("PASS" in radLine and "PASS" in compareLine):
            # this one is for Radia and validation
            if ("PASS" in radLine and ("SOM" in radLine or "EDIT" in radLine or "RNA_TUM_VAR" in radLine or "RNA_NOR_VAR" in radLine) and 
                ("PASS" in compareLine and ("SOM" in compareLine or "EDIT" in compareLine or "RNA_TUM_VAR" in compareLine or "RNA_NOR_VAR" in compareLine))):
            #if ("PASS" in radLine and "SOM" in radLine and "SNP" in compareLine and "Somatic" in compareLine):
            #if ("Somatic" in radLine and "SNP" in radLine and "Somatic" in compareLine):
            #if ("SOM" in radLine and "Somatic" in compareLine and "SNP" in compareLine): 
            #if ("PASS" in radLine and "SOM" in radLine):
            #if ("SOM" in radLine):
            #if ("PASS" in radLine and "SNP" in compareLine and "Somatic" in compareLine):
                i_statsDict["overlap_pass_events"] += 1
            
            # for each key to compare
            # their can be multiple keys for one filter such as blq and bldp for blacklists            
            for (radKeyString, cmpKeyString) in aCompareDict.iteritems():
                # break up the strings to get the individual keys
                radKeyList = radKeyString.split(",")
                cmpKeyList = cmpKeyString.split(",")
                
                # set some booleans
                foundInRad = False
                foundInCmp = False
                # search for one of them
                for radKey in radKeyList:
                    # if we find one
                    if (radKey in radLine):
                        foundInRad = True
                        break;
                # search for one of them
                for cmpKey in cmpKeyList:
                    # if we find one
                    if (cmpKey in compareLine):
                        foundInCmp = True
                        break;
                
                # if the keys exist in both files at the same position, then count them
                if (foundInRad and foundInCmp):
                    # if these are germline or they haven't been found in dbSnp, then count them
                    #if (((radKey == "GERM") or ("DB" not in radLine and "DB" not in compareLine))):
                    #if ("SNP" in compareLine):
                    #if ("SNP" in compareLine and ((radKey == "GERM") or ("DB" not in radLine and "DB" not in compareLine))):
                    #if ("SNP" in compareLine):
                    if ("PASS" in compareLine):
                    #if ("PASS" in compareLine and "SNP" in compareLine):
                    #if (True):
                        i_statsDict["overlap_" + radKey] += 1
                        
                        splitLine = radLine.split("\t")    
                        filterString = splitLine[6]
                        filterList = filterString.split(";")
                    
                        #if ("PASS" in radLine and "SNP" in radLine):
                        #if ("PASS" in radLine and "SNP" in radLine):
                        if ("PASS" in radLine):
                        #if ("SNP" in radLine):
                        #if (True):
                            i_statsDict["overlap_pass_" + radKey] += 1
                            if (anIsDebug):
                                logging.debug("found call %s", compareLine)
                            if (anOverlapFilename != None):
                                
                                # add to maf
                                
                                #caller = ";ucsc;"
                                #if ("radia" in radLine):
                                #    caller += "radia;"
                                #if ("bambam" in radLine):
                                #    caller += "bambam;"
                                
                                #caller = ";rnaCall"
                                #splitLine = compareLine.split("\t")
                                #splitLine[2] += caller
                                #overlapFileHandler.write("\t".join(splitLine) + "\n")
                                
                                #caller = ";rnaCall"
                                #cmpSplitLine = compareLine.split("\t")
                                #callers = cmpSplitLine[2] + caller
                                #callers = callers.replace(";;", ",")
                                #callers = callers.replace(";", ",")
                                
                                #radSplitLine = radLine.split("\t")
                                #radSplitLine[7] += ";CALLER=" + callers
                                #overlapFileHandler.write("\t".join(radSplitLine) + "\n")
                                
                                overlapFileHandler.write(radLine + "\n")
                                # we only want to write the line to the overlap file once
                                # even if it matches as a SOM and an EDIT
                                break;
                                #overlapFileHandler.write(compareLine + "\n")
                        else:
                            if (anIsDebug):
                                logging.debug("found but no radia pass %s %s", radLine, compareLine)
                            #overlapFileHandler.write(compareLine + "\n")
                            splitLine = radLine.split("\t")
                            
                            filterString = splitLine[6]
                            filterList = filterString.split(";")
                            for filterKey in filterList:
                                i_filterDict[filterKey] += 1
                                
                            if (aNonOverlapFilename != None):
                                #nonOverlapFileHandler.write(compareLine + "\n")
                                nonOverlapFileHandler.write(radLine + "\n")
                        
                elif (anIsDebug and foundInRad):
                    logging.debug("overlap but not found in compare file %s %s %s", radKey, radLine, compareLine)
                    #overlapFileHandler.write(compareLine + "\n")
                elif (anIsDebug and foundInCmp):
                    logging.debug("overlap but not found in RADIA %s %s %s", cmpKey, radLine, compareLine)
                    #overlapFileHandler.write(compareLine + "\n")
                elif (anIsDebug):
                    logging.debug("overlap but not same type %s %s %s %s", radKeyList, cmpKeyList, radLine, compareLine)
                    #overlapFileHandler.write(compareLine + "\n")

    # aTCGAId, aChrom, rad_events, cmp_events, overlap_events, [rad_key, cmp_key, overlap_radKey]{n}
    #outputHeader = ["PatientId", "Chrom", "rad_events", "cmp_events", "overlap_events", "rad_pass_events", "cmp_pass_events", "overlap_pass_events"]
    outputList = [aTCGAId, aChrom]
    outputList += [str(i_statsDict["rad_events"]), str(i_statsDict["cmp_events"]), str(i_statsDict["overlap_events"])]
    outputList += [str(i_statsDict["rad_pass_events"]), str(i_statsDict["cmp_pass_events"]), str(i_statsDict["overlap_pass_events"])]
    
    # for each key to compare, get the total radias, total cmps, and overlaps
    for radKey in sorted(aCompareDict.iterkeys()):
        cmpKey = aCompareDict[radKey]
        outputList += [str(i_statsDict["rad_" + radKey]), str(i_statsDict["cmp_" + cmpKey]), str(i_statsDict["overlap_" + radKey])]
        outputList += [str(i_statsDict["rad_pass_" + radKey]), str(i_statsDict["cmp_pass_" + cmpKey]), str(i_statsDict["overlap_pass_" + radKey])]
          
    for (filterKey, count) in i_filterDict.iteritems():
        logging.debug("filter: %s\t%s", filterKey, count)
        
    # get the files
    i_statsFileHandler = None
    if (aStatsFilename != None):
        i_statsFileHandler = get_append_fileHandler(aStatsFilename)
        i_statsFileHandler.write("\t".join(outputList) + "\n")
        i_statsFileHandler.close() 
                        
    stopTime = time.time()
    logging.info("Chrom %s and Id %s: Total time=%s hrs, %s mins, %s secs", aChrom, aTCGAId, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime))
    logging.info("\t".join(outputList))
    
    if (anOverlapFilename != None):
        overlapFileHandler.close()
                      
    if (aNonOverlapFilename != None):
        nonOverlapFileHandler.close()
            
    return


def main():
    
    #python radiaCompare.py TCGA-AB-2995 12 ../data/test/TCGA-AB-2995.vcf ../data/test/TCGA-AB-2995.vcf -c "SOM=Somatic" -s ../stats/radia/cmpRadBB.tab --log=DEBUG 
        
    # create the usage statement
    usage = "usage: python %prog id chrom radFile compareFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    i_cmdLineParser.add_option("-c", "--compareList", dest="compareList", metavar="COMPARE_LIST", help="a comma separated list of key/values comparisons where the key is in RADIA and the value is in the compare file")
    i_cmdLineParser.add_option("-s", "--statsFilename", dest="statsFilename", metavar="STATS_FILE", help="the name of the stats file, sys.stdout by default")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDERR by default")
    i_cmdLineParser.add_option("-o", "--overlapFilename", dest="overlapFilename", metavar="OVERLAP_FILE", help="the name of the overlap file")
    i_cmdLineParser.add_option("-n", "--nonOverlapFilename", dest="nonOverlapFilename", metavar="NON_OVERLAP_FILE", help="the name of the non-overlap file")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(4,17,1)
    i_argLength = len(sys.argv)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)
    
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_id = str(i_cmdLineArgs[0])
    i_chr = str(i_cmdLineArgs[1])
    i_radiaFilename = str(i_cmdLineArgs[2])
    i_compareFilename = str(i_cmdLineArgs[3])
    
    # get the optional params with default values   
    i_logLevel = i_cmdLineOptions.logLevel
    
    # try to get any optional parameters with no defaults   
    # check for any errors
    writeFilenameList = []
    readFilenameList = [i_radiaFilename, i_compareFilename]
     
    i_statsFilename = None
    i_logFilename = None
    i_compareString = None
    i_overlapFilename = None
    i_nonOverlapFilename = None
    i_compareDict = collections.defaultdict(list)
    if (i_cmdLineOptions.overlapFilename != None):
        i_overlapFilename = str(i_cmdLineOptions.overlapFilename)
        writeFilenameList += [i_overlapFilename]
    if (i_cmdLineOptions.nonOverlapFilename != None):
        i_nonOverlapFilename = str(i_cmdLineOptions.nonOverlapFilename)
        writeFilenameList += [i_nonOverlapFilename]
    if (i_cmdLineOptions.statsFilename != None):
        i_statsFilename = str(i_cmdLineOptions.statsFilename)
        writeFilenameList += [i_statsFilename]
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
        writeFilenameList += [i_logFilename]
    if (i_cmdLineOptions.compareList != None):
        i_compareString = str(i_cmdLineOptions.compareList)
        i_compareList = i_compareString.split(",")
        
        for keyValue in i_compareList:
            (key, value) = keyValue.split("=")
            i_compareDict[key] = value
        
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

    # do some debugging
    if (i_debug):
        logging.debug("id=%s", i_id)
        logging.debug("chr=%s", i_chr)
        logging.debug("radiaFile=%s", i_radiaFilename)
        logging.debug("overlapFilename=%s" % i_overlapFilename)
        logging.debug("nonOverlapFilename=%s" % i_nonOverlapFilename)
        logging.debug("compareFile=%s", i_compareFilename)
        logging.debug("statsFile=%s", i_statsFilename)
        logging.debug("logLevel=%s", i_logLevel)
        logging.debug("logFile=%s", i_logFilename)
        logging.debug("compareDict=%s", i_compareDict)
        
    if (not radiaUtil.check_for_argv_errors(None, readFilenameList, writeFilenameList)):
        sys.exit(1)           
    
    compare_events(i_id, i_chr, i_radiaFilename, i_compareFilename, i_statsFilename, i_overlapFilename, i_nonOverlapFilename, i_compareDict, i_debug)
       
    return

main()
sys.exit(0)
