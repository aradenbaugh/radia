#!/usr/bin/env python

import sys
from optparse import OptionParser
import radiaUtil
import logging
import gzip
import collections
from itertools import izip


'''
'    RNA and DNA Integrated Analysis (RADIA) identifies RNA and DNA variants in NGS data.
'    Copyright (C) 2010-2015  Amie Radenbaugh
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


def get_vcf_data(aVcfFile, anIsDebug):
    
    # open the VCF file
    fileHandler = get_read_fileHandler(aVcfFile)
     
    for line in fileHandler:
          
        # here are some examples of .vcf lines:
        # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31.4;DP=53;FA=0.04;INDEL=0;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP      
        # GT:DP:INDEL:START:STOP:AD:AF:BQ:SB      0/0:2:0:0:0:2:1.0:35.5:0.0      0/0:1:0:0:0:1:1.0:39.0:1.0      0/1:50:0:2:0:48,2:0.96,0.04:31.65,17.5:0.75,0.5
        
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug):
            logging.debug("VCF: %s", line)    
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
                
        # these are header lines, so just continue    
        elif (line.startswith("#")):
            continue
        
        # split the line on the tab
        splitLine = line.split("\t")

        # the coordinate is the second element
        chrom = splitLine[0]
        stopCoordinate = int(splitLine[1])
        ids = splitLine[2]
        ref = splitLine[3]
        alt = splitLine[4]
        score = splitLine[5]
        filterSet = set(splitLine[6].split(";"))
        
        # if there are no filters so far, then clear the list
        if (len(filterSet) == 1 and "PASS" in filterSet):
            filterSet = set()
            
        infoList = splitLine[7].split(";")
        infoDict = collections.defaultdict(list)
        for info in infoList:
            keyValueList = info.split("=")
            # some keys are just singular without a value (e.g. DB, SOMATIC, etc.)
            if (len(keyValueList) == 1):
                infoDict[keyValueList[0]] = ["True"]
            else:
                # the value can be a comma separated list
                infoDict[keyValueList[0]] = keyValueList[1].split(",")
        
        # yield all the information about the current coordinate
        yield (chrom, stopCoordinate, ids, ref, alt, score, filterSet, infoDict, "\t".join(splitLine[8:]))
        
    fileHandler.close()
    return


def parse_blat_input(aBlatFile, anIsDebug):
    '''
    ' This function parses the input to BLAT.
    '
    ' aBlatFile:  An input file to BLAT
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    # open the file
    fileHandler = get_read_fileHandler(aBlatFile)
    blatDict = dict()
     
    for line in fileHandler:
          
        # we can ignore the lines that start with # for now
        if (line.isspace()):
            continue;

        # we only want the header lines        
        if (line.startswith(">")):
            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")
    
            if (anIsDebug):
                logging.debug("BLAT: %s", line)    
                
            # split the line on the tab
            splitLine = line.split()
    
            # get the coordinate data = rnaTumor_7_55196749_HS2144:2:1108:17342:164248
            blatId = splitLine[1]
            blatSplitId = blatId.split("_")
            prefix = blatSplitId[0]
            coordinateId = "_".join(blatSplitId[1:3])
            
            if coordinateId not in blatDict:
                blatDict[coordinateId] = collections.defaultdict(dict)
            if prefix not in blatDict[coordinateId]:
                blatDict[coordinateId][prefix] = list()
                
            blatDict[coordinateId][prefix].append(line)
        
    fileHandler.close()
    
    return blatDict
    
               
def main():
        
    # create the usage statement
    usage = "usage: python %prog id chrom vcfFile blatInputFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, STDOUT by default")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    
    i_cmdLineParser.add_option("-d", "--readDepthCutoff", type="int", default=int(4), dest="readDepthCutoff", metavar="READ_DP_CUTOFF", help="the minimum number of reads that are necessary before applying this filter, %default by default")
    i_cmdLineParser.add_option("-p", "--readPercentCutoff", type="float", default=float(0.95), dest="readPercentCutoff", metavar="READ_PERCENT_CUTOFF", help="the maximum percentage of reads with the alternative allele at the beginning or end of the reads, %default by default")
    
    i_cmdLineParser.add_option("-k", "--keepPreviousFilters", action="store_true", default=False, dest="keepPreviousFilters", help="by default the previous filters are overwritten with the pbias filter, include this argument if the previous filters should be kept")
    
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
    i_vcfFilename = str(i_cmdLineArgs[1])
    i_blatInputFilename = str(i_cmdLineArgs[2])
    
    # get the optional params with default values
    i_logLevel = i_cmdLineOptions.logLevel
    i_readDepthCutoff = i_cmdLineOptions.readDepthCutoff
    i_readPercentCutoff = i_cmdLineOptions.readPercentCutoff
    i_keepPreviousFiltersFlag = i_cmdLineOptions.keepPreviousFilters
    
    # try to get any optional parameters with no defaults    
    i_outputFilename = None
    i_logFilename = None
    if (i_cmdLineOptions.outputFilename != None):
        i_outputFilename = str(i_cmdLineOptions.outputFilename)
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
        logging.debug("vcfFilename=%s", i_vcfFilename)
        logging.debug("blatInputFilename=%s", i_blatInputFilename)
        logging.debug("outputFilename=%s", i_outputFilename)
        logging.debug("logFilename=%s", i_logFilename)
        logging.debug("readDepthCutoff=%s", i_readDepthCutoff)
        logging.debug("readPerentCutoff=%s", i_readPercentCutoff)
        logging.debug("keepPreviousFiltersFlag? %s", i_keepPreviousFiltersFlag)
            
    # check for any errors
    i_writeFilenameList = []
    if (i_outputFilename != None):
        i_writeFilenameList = [i_outputFilename]
    if (i_logFilename != None):
        i_writeFilenameList = [i_logFilename]
        
    i_readFilenameList = [i_vcfFilename, i_blatInputFilename]
    
    if (not radiaUtil.check_for_argv_errors(None, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)

    # open the output stream
    i_outputFileHandler = None
    if (i_outputFilename != None):
        i_outputFileHandler = get_write_fileHandler(i_outputFilename)
    
    # get the BLAT input
    i_blatCoordinateDict = parse_blat_input(i_blatInputFilename, i_debug)
    
    # get the VCF generator   
    i_vcfGenerator  = get_vcf_data(i_vcfFilename, i_debug)
    
    for (vcfChr, vcfStopCoordinate, vcfId, vcfRef, vcfAlt, vcfScore, vcfFilterSet, vcfInfoDict, restOfLine) in i_vcfGenerator:
        if (i_debug):
            logging.debug("VCF Data: %s %s %s %s %s %s %s %s %s", vcfChr, str(vcfStopCoordinate), vcfId, vcfRef, vcfAlt, vcfScore, str(vcfFilterSet), str(vcfInfoDict), restOfLine) 
            
        modTypes = vcfInfoDict["MT"]
        modTypeFilters = dict()
        atLeastOnePass = False
        blatHitsList = list()
    
        # get the coordinate for this position
        coordinate = vcfChr + "_" + str(vcfStopCoordinate)
        for modType in modTypes:
            if (modType == "NOR_EDIT"):
                if (coordinate in i_blatCoordinateDict and "rnaNormal" in i_blatCoordinateDict[coordinate]):
                    # for each coordinate, get a dict of corresponding blat hits
                    blatHitsList = i_blatCoordinateDict[coordinate]["rnaNormal"]
            elif (modType == "SOM" or modType == "TUM_EDIT"):
                if (coordinate in i_blatCoordinateDict and "rnaTumor" in i_blatCoordinateDict[coordinate]):
                    # for each coordinate, get a dict of corresponding blat hits
                    blatHitsList = i_blatCoordinateDict[coordinate]["rnaTumor"]
            
            if (i_debug):
                logging.debug("coordinate=%s", coordinate)
            
            starts = 0
            ends = 0
            middles = 0
            total = 0
            
            # for each read, investigate the blat input
            for readId in blatHitsList:
                
                readIdList = readId.split("_")
                total += 1
                position = readIdList[7]
                
                if (i_debug):
                    logging.debug("readId=%s, position=%s", readId, position)
                
                if (position == "start"):
                    starts += 1
                elif (position == "end"):
                    ends += 1
                elif (position == "middle"):
                    middles += 1    
                            
            # if we have enough reads
            if (total > i_readDepthCutoff):        
                if (round(starts/float(total),2) >= i_readPercentCutoff):
                    modTypeFilters[modType] = "pbias"
                elif (round(ends/float(total),2) >= i_readPercentCutoff):
                    modTypeFilters[modType] = "pbias"
                #elif (round(middles/float(total),2) >= i_readPercentCutoff):
                #    modTypeFilters[modType] = "pbias"
                else:
                    modTypeFilters[modType] = "PASS"
                    atLeastOnePass = True
            # if we don't have the minimum number of reads, then pass
            else:
                modTypeFilters[modType] = "PASS"
                atLeastOnePass = True
                
            if (i_debug):
                logging.debug("coordinate=%s, starts=%s, middles=%s, ends=%s, total=%s, positionalBias=%s", coordinate, starts, middles, ends, total, str(modTypeFilters))
           
        # make a copy of the list to manipulate
        modTypesTmpList = list(modTypes)
        modChanges = vcfInfoDict["MC"]
        # if at least one passed, then remove the ones that didn't
        for (modType, modChange) in izip(modTypes, modChanges):
            # if at least one passed, then remove the ones that didn't
            if (atLeastOnePass):
                if (modTypeFilters[modType] == "pbias"):
                    modTypesTmpList.remove(modType)
                    modChanges.remove(modChange)
        
        # set the modTypes and modChanges
        vcfInfoDict["MT"] = modTypesTmpList
        vcfInfoDict["MC"] = modChanges 
         
        # if at least one passed, then set pass
        if (atLeastOnePass):
            vcfFilterSet = ["PASS"]
        else:
            # if the user wants to keep the previous filters
            if (i_keepPreviousFiltersFlag):
                # if the call previous passed, then just set pbias
                if (len(vcfFilterSet) == 1 and "PASS" in vcfFilterSet):
                    vcfFilterSet = ["pbias"] 
                # otherwise, add it to the previous filters
                else:
                    vcfFilterSet.add("pbias")
            # otherwise, just set the pbias filter
            else:
                vcfFilterSet = ["pbias"]
            
            # update the mod filters    
            modTypes = vcfInfoDict["MT"]
            modChanges = vcfInfoDict["MC"]
            origins = vcfInfoDict["ORIGIN"]
            modFilters = [] if vcfInfoDict["MF"] is None else vcfInfoDict["MF"]
            modFilterTypes = [] if vcfInfoDict["MFT"] is None else vcfInfoDict["MFT"]
                    
            for origin in origins:
                for (modType, modChange) in izip(modTypes, modChanges):
                    modFilterTypes.append("_".join([origin, modType, modChange]))
                    modFilters.append("_".join(vcfFilterSet))
                
            vcfInfoDict["MF"] = modFilters
            vcfInfoDict["MFT"] = modFilterTypes
            
        output = [vcfChr, str(vcfStopCoordinate), vcfId, vcfRef, vcfAlt, vcfScore]
            
        # if there are no filters so far, then this call passes
        if (len(vcfFilterSet) == 0):
            vcfFilterSet.add("PASS")
            
        output.append(";".join(vcfFilterSet))
        
        # add the modified info dict
        infoField = ""
        for key in sorted(vcfInfoDict.iterkeys()):
            if (len(vcfInfoDict[key]) == 0):
                continue
            elif ("True" in vcfInfoDict[key]):
                infoField += key + ";"
            else:    
                infoField += key + "=" + ",".join(vcfInfoDict[key]) + ";"
        
        output.append(infoField.rstrip(";"))
        
        output.append(restOfLine)
        
        if (i_outputFilename != None):
            i_outputFileHandler.write("\t".join(output) + "\n")
        else:
            print >> sys.stdout, "\t".join(output)
                
    # close the files 
    if (i_outputFilename != None):
        i_outputFileHandler.close()
        
    return
 

main()    
sys.exit(0)
