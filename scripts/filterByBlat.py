#!/usr/bin/env python

import sys
import time
from optparse import OptionParser
from itertools import izip
import radiaUtil
import collections
import logging
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


def get_vcf_data(aVcfFile, aPassOnlyFlag, anIsDebug):
    '''
    ' This function reads from a .vcf input file and uses the python generator to yield the information
    ' one line at a time.  It ignores empty lines and strips trailing \r\n characters.  This function
    ' yields all the information from the VCF file.
    '
    ' aVcfFile:  A VCF file
    ' aPassOnlyFlag:  If all calls should be processed or only those calls that passed the filters thus far
    ' anIsDebug: A flag for outputting debug messages to STDERR
    ''' 
    
    # open the VCF file
    fileHandler = get_read_fileHandler(aVcfFile)
     
    for line in fileHandler:
          
        # here are some examples of .vcf lines:
        # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31;DP=53;FA=0.04;INDEL=0;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP      
        # GT:DP:AD:AF:INDEL:START:STOP:BQ:SB      0/0:2:2,0:1.0,0.0:0:0:0::36,0:0.0,0.0      0/0:1:1,0:1.0,0.0:0:0:0:39,0:1.0,0.0      0/1:50:48,2:0.96,0.04:0:2:0:32,18:0.75,0.5
        
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
        
        # if we are only suppose to process the passed calls
        # and this call has not passed, then skip it
        elif (aPassOnlyFlag and "PASS" not in line):    
            continue;

        # split the line on the tab
        splitLine = line.split("\t")

        # the coordinate is the second element
        chrom = splitLine[0]
        stopCoordinate = int(splitLine[1])
        idList = splitLine[2]
        refList = splitLine[3]
        altList = splitLine[4]
        score = splitLine[5]
        filterSet = set(splitLine[6].split(";"))
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
        yield (chrom, stopCoordinate, idList, refList, altList, score, filterSet, infoDict, "\t".join(splitLine[8:]))

    fileHandler.close()
    return


def parse_blat_output(aBlatFile, anOutputFormat, anIsDebug):
    '''
    ' This function parses the output from BLAT.  Two formats are supported:  BLAST NCBI-8 and PSL.  It groups 
    ' all of the information from one query sequence and uses the python generator to yield the information.  
    ' It ignores empty lines and strips trailing \r\n characters.
    '
    ' aBlatFile:  A output file from BLAT
    ' anOutputFormat:  BLAST or PSL
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    # open the file
    fileHandler = get_read_fileHandler(aBlatFile)
    blatHitsDict = collections.defaultdict(dict)
     
    for line in fileHandler:
          
        # we can ignore the lines that start with # for now
        if (line.isspace()):
            continue;
        
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        if (anIsDebug):
            logging.debug("BLAT: %s", line)    
            
        # split the line on the tab
        splitLine = line.split("\t")

        # get the coordinate data = rnaTumor_7_55196749_HS2144:2:1108:17342:164248
        if (anOutputFormat == "PSL"):
            blatId = splitLine[9]
        elif (anOutputFormat == "BLAST"):
            blatId = splitLine[0]
        
        blatSplitId = blatId.split("_")
        prefix = blatSplitId[0]
        coordinateId = "_".join(blatSplitId[1:3])
        readId = "_".join(blatSplitId[0:4])
        
        if coordinateId not in blatHitsDict:
            blatHitsDict[coordinateId] = collections.defaultdict(dict)
        if prefix not in blatHitsDict[coordinateId]:
            blatHitsDict[coordinateId][prefix] = collections.defaultdict(list)
            
        blatHitsDict[coordinateId][prefix][readId].append(line)
        
    fileHandler.close()
    return blatHitsDict
    
    
def is_valid_read_blast_format(aBlatHitsList, aVCFChrom, aVCFCoordinate, anOrderMagnitude, anIsDebug):
    '''
    ' This method determines if the read is valid.  It compares all of the BLAT hits to determine if 
    ' the read mapped to another location in the genome with a better score.
    '
    ' aBlatHitsList:        The list of blat hits from the output
    ' aVCFChrom:            The chrom
    ' aVCFCoordinate:       The coordinate
    ' anOrderMagnitude:     The order of magnitude
    ' anIsDebug:            The debug flag
    '''
    # set defaults
    validRead = ""
    minOverlappingEValue = sys.float_info.max
    otherEValues = []
    maxOverlappingIdentity = sys.float_info.min
    otherIdentities = []
    
    # for each read, investigate the blat hits to see if this read is valid
    for blatHit in aBlatHitsList:
        if (anIsDebug):
            logging.debug("blatHit: %s", blatHit)

        # split the line on the tab
        splitLine = blatHit.split("\t")
        
        blatReadIdFull = splitLine[0]
        blatReadIdSplit = blatReadIdFull.split("_")
        readFlag = int(blatReadIdSplit[9])
        readLength = int(blatReadIdSplit[10])
        readLengthHalf = readLength/2
       
        # if this read is properly paired and it is the primary alignment
        readPaired = int(readFlag) & 0x1
        properlyPaired = int(readFlag) & 0x2
        if ((readPaired and properlyPaired) or (not readPaired)):
            
            blatChrom = splitLine[1]
            blatIdentity = float(splitLine[2])
            blatAlignmentLength = int(splitLine[3])
            #blatNumMismatches = int(splitLine[4])
            #blatGapOpenings = int(splitLine[5])
            #blatQueryStart = int(splitLine[6])
            #blatQueryStop = int(splitLine[7])
            blatRefStart = int(splitLine[8])
            blatRefStop = int(splitLine[9])
            blatEValue = float(splitLine[10])
            #print "eval", blatEValue
            #blatBitScore = float(splitLine[11])
            
            # if the blat hit covers the coordinate that we're investigating
            if ((("chr" + aVCFChrom) == blatChrom or aVCFChrom == blatChrom) and aVCFCoordinate >= blatRefStart and aVCFCoordinate <= blatRefStop):
                # if the evalue is less than the min, then move the old min to the list and set the new min
                if (blatEValue < minOverlappingEValue):
                    otherEValues.append(minOverlappingEValue)
                    minOverlappingEValue = blatEValue
                    validRead = blatHit
                else:
                    otherEValues.append(blatEValue)
                    
                # if the blat alignment length is greater than half the read size and 
                # the percent identity is greater than the current max, then move the old max to the list and set the new max
                if (blatAlignmentLength > readLengthHalf and blatIdentity > maxOverlappingIdentity):
                    otherIdentities.append(maxOverlappingIdentity)
                    maxOverlappingIdentity = blatIdentity
                    validRead = blatHit
                elif (blatAlignmentLength > readLengthHalf):
                    otherIdentities.append(blatIdentity)
            else:
                #print "blat hit doesn't overlap position"
                otherEValues.append(blatEValue)
                if (blatAlignmentLength > readLengthHalf):
                    otherIdentities.append(blatIdentity)
            
            #print "min=", minOverlappingEValue, "others=", otherEValues
            
    # if none of the hits were properly paired, then just return false
    if (len(otherEValues) == 0):
        return (False, validRead)
    # if no read covered the desired location
    elif (minOverlappingEValue == sys.float_info.max):
        return (False, validRead)

    # if the overlapping read has the lowest e-value (within a certain order of magnitude) or the highest identity, then return true
    # otherwise, our read is not mapped to the best spot in the genome, so return false
    minOtherEvalue = min(otherEValues)
    if (len(otherIdentities) > 0):
        maxOtherIdentity = max(otherIdentities)
    else:
        maxOtherIdentity = 0.0
        
    if (anOrderMagnitude == 0):
        if (minOverlappingEValue < minOtherEvalue or maxOverlappingIdentity > maxOtherIdentity):
            return (True, validRead)
        else:
            return (False, validRead)
    else:
        #if (minOtherEvalue < (minOverlappingEValue/(1/(10.0**anOrderMagnitude))) and maxOtherIdentity > maxOverlappingIdentity):
        if ((minOverlappingEValue/(1/(10.0**anOrderMagnitude))) < minOtherEvalue or maxOverlappingIdentity > maxOtherIdentity):
            return (True, validRead)
        else:
            return (False, validRead)

    '''
    # if another read had a better evalue (within a certain order of magnitude), then return false
    # otherwise, our read is mapped to the best spot in the genome
    minOther = min(otherEValues)
    if (anOrderMagnitude == 0):
        if (minOther < minOverlappingEValue):
            return (False, validRead)
        else:
            return (True, validRead)
    else:
        if (minOther < (minOverlappingEValue/(1/(10.0**anOrderMagnitude)))):
            return (False, validRead)
        else:
            return (True, validRead)
    '''
    return 

    
def is_valid_read_psl_format(aBlatHitsList, aVCFChrom, aVCFCoordinate, anIsDebug):
    '''
    ' This method determines if the read is valid.  It compares all of the BLAT hits to determine if 
    ' the read mapped to another location in the genome with a better score.
    '
    ' aBlatHitsList:        The list of blat hits from the output
    ' aVCFChrom:            The chrom
    ' aVCFCoordinate:       The coordinate
    ' anIsDebug:            The debug flag
    '''
    # set defaults
    validRead = ""
    maxOverlappingScore = -sys.maxint - 1
    otherBlatScores = []
    
    # for each read, investigate the blat hits to see if this read is valid
    for blatHit in aBlatHitsList:
        if (anIsDebug):
            logging.debug("blatHit: %s", blatHit)

        # split the line on the tab
        splitLine = blatHit.split("\t")
        
        blatReadIdFull = splitLine[9]
        blatReadIdSplit = blatReadIdFull.split("_")
        readFlag = int(blatReadIdSplit[9])
       
        # if this read is properly paired and it is the primary alignment
        readPaired = int(readFlag) & 0x1
        properlyPaired = int(readFlag) & 0x2
        if ((readPaired and properlyPaired) or (not readPaired)):
            
            matches = int(splitLine[0])
            mismatches = int(splitLine[1])
            repeatMatches = int(splitLine[2])
            #numNs = int(splitLine[3])
            qNumInserts = int(splitLine[4])
            #qNumInsertBases = int(splitLine[5])
            tNumInserts = int(splitLine[6])
            #tNumInsertBases = int(splitLine[7])
            #strand = splitLine[8]
            #blatReadIdFill = splitLine[9]
            #qSize = int(splitLine[10])
            #qStart = int(splitLine[11])
            #qEnd = int(splitLine[12])
            tName = splitLine[13]
            #tSize = int(splitLine[14])
            tStart = int(splitLine[15])
            tEnd = int(splitLine[16])
            #blockCount = int(splitLine[17])
            #blockSizes = splitLine[18].split(",")
            #qStarts = splitLine[19].split(",")
            #tStarts = splitLine[20].split(",")
            
            # this is the blat score on the web interface
            blatScore = (1 * (matches + repeatMatches)) - (1 * mismatches) - qNumInserts - tNumInserts
            #print blatScore, matches, repeatMatches, mismatches, qNumInserts, tNumInserts
            
            # if the blat hit covers the coordinate that we're investigating
            if (("chr" + aVCFChrom) == tName and aVCFCoordinate >= tStart and aVCFCoordinate <= tEnd):
                # if the score is greater than the max, then move the old max to the list and set the new max
                if (blatScore > maxOverlappingScore):
                    otherBlatScores.append(maxOverlappingScore)
                    maxOverlappingScore = blatScore
                    validRead = blatHit
                else:
                    otherBlatScores.append(blatScore)
            else:
                #print "blat hit doesn't overlap position"
                otherBlatScores.append(blatScore)
            
            #print "max=", maxOverlappingScore, "others=", otherBlatScores
            
    # if none of the hits were properly paired, then just return false
    if (len(otherBlatScores) == 0):
        return (False, validRead)
    # if no read covered the desired location
    elif (maxOverlappingScore == sys.maxint):
        return (False, validRead)

    # if another read had a better score, then return false
    # otherwise, our read is mapped to the best spot in the genome
    maxOther = max(otherBlatScores)
    if (maxOther > maxOverlappingScore):
        return (False, validRead)
    else:
        return (True, validRead)
     
               
def main():
    
    # command for running this on a small test case: 
    #python filterByBlat.py TCGA-00-4454 7 ../data/test/TCGA-00-4454_EGFR.vcf ../data/test/TCGA-00-4454_EGFR.fa ../data/test/TCGA-00-4454_EGFR.psl --dnaNormalFilename=../data/test/TCGA-00-4454_EGFR.reads
    
    startTime = time.time()
    
    # create the usage statement
    usage = "usage: python %prog id chrom vcfFile blatInputFile blatOutputFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    # add the optional parameters
    i_cmdLineParser.add_option("-c", "--allVCFCalls", action="store_false", default=True, dest="passedVCFCallsOnly", help="by default only the VCF calls that have passed all filters thus far are processed, include this argument if all of the VCF calls should be processed")
    i_cmdLineParser.add_option("-k", "--keepPreviousFilters", action="store_true", default=False, dest="keepPreviousFilters", help="by default the previous filters are overwritten with the blat filter, include this argument if the previous filters should be kept")
    
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, STDOUT by default")
    i_cmdLineParser.add_option("-b", "--blatOutputFormat", dest="blatOutputFormat", metavar="OUTPUT_FORMAT", default="BLAST", help="the BLAT output format, BLAST by default")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    
    i_cmdLineParser.add_option("-n", "--blatDnaNormalReads", action="store_true", default=False, dest="blatDnaNormalReads", help="include this argument if the normal DNA reads should be processed")
    i_cmdLineParser.add_option("-x", "--blatRnaNormalReads", action="store_true", default=False, dest="blatRnaNormalReads", help="include this argument if the normal RNA reads should be processed")
    i_cmdLineParser.add_option("-t", "--blatDnaTumorReads", action="store_true", default=False, dest="blatDnaTumorReads", help="include this argument if the tumor DNA reads should be processed")
    i_cmdLineParser.add_option("-r", "--blatRnaTumorReads", action="store_true", default=False, dest="blatRnaTumorReads", help="include this argument if the tumor RNA reads should be processed")

    i_cmdLineParser.add_option("-d", "--readDepthCutoff", type="int", default=int(4), dest="readDepthCutoff", metavar="READ_DP_CUTOFF", help="the minimum number of valid reads that are necessary, %default by default")
    i_cmdLineParser.add_option("-p", "--readPercentCutoff", type="float", default=float(0.10), dest="readPercentCutoff", metavar="READ_PERCENT_CUTOFF", help="the minimum percentage of valid reads that are necessary, %default by default")
    
    #i_cmdLineParser.add_option("-e", "--eValueCutoff", type="float", default=float(10e-6), dest="eValueCutoff", metavar="EVAL_CUTOFF", help="the e-value cutoff for determining if a blat hit is significant, %default by default")
    #i_cmdLineParser.add_option("-u", "--upperIdentityCutoff", type="float", default=float(0.95), dest="upperIdentityCutoff", metavar="UPPER_CUTOFF", help="the upper cutoff for the match length adjusted identity to determine if a blat hit is significant, %default by default")
    #i_cmdLineParser.add_option("-l", "--lowerIdentityCutoff", type="float", default=float(0.5), dest="lowerIdentityCutoff", metavar="LOWER_CUTOFF", help="the lower cutoff for the match length adjusted identity to determine if a second blat hit is significant, %default by default")
          
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(5,27,1)
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
    i_blatOutputFilename = str(i_cmdLineArgs[3])
    
    # get the optional params with default values
    i_passedVCFCallsOnlyFlag = i_cmdLineOptions.passedVCFCallsOnly
    i_keepPreviousFiltersFlag = i_cmdLineOptions.keepPreviousFilters
    i_blatOutputFormat = i_cmdLineOptions.blatOutputFormat
    i_logLevel = i_cmdLineOptions.logLevel
    i_readDepthCutoff = i_cmdLineOptions.readDepthCutoff
    i_readPercentCutoff = i_cmdLineOptions.readPercentCutoff
    #i_eValueCutoff = i_cmdLineOptions.eValueCutoff
    #i_upperIdentityCutoff = i_cmdLineOptions.upperIdentityCutoff
    #i_lowerIdentityCutoff = i_cmdLineOptions.lowerIdentityCutoff
    
    i_blatDnaNormalReads = i_cmdLineOptions.blatDnaNormalReads
    i_blatDnaTumorReads = i_cmdLineOptions.blatDnaTumorReads
    i_blatRnaNormalReads = i_cmdLineOptions.blatRnaNormalReads
    i_blatRnaTumorReads = i_cmdLineOptions.blatRnaTumorReads
    
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
        logging.debug("blatOutputFilename=%s", i_blatOutputFilename)
        logging.debug("passedCallsOnly? %s", i_passedVCFCallsOnlyFlag)
        logging.debug("keepPreviousFiltersFlag? %s", i_keepPreviousFiltersFlag)
        logging.debug("blatOutputFormat=%s", i_blatOutputFormat)
        
        logging.debug("blatDnaNormal? %s", i_blatDnaNormalReads)
        logging.debug("blatDnaTumor? %s", i_blatDnaTumorReads)
        logging.debug("blatRnaNormal? %s", i_blatRnaNormalReads)
        logging.debug("blatRnaTumor? %s", i_blatRnaTumorReads)
        
        logging.debug("readDepthCutoff=%s", i_readDepthCutoff)
        logging.debug("readPerentCutoff=%s", i_readPercentCutoff)
            
    # check for any errors
    i_writeFilenameList = []
    if (i_outputFilename != None):
        i_writeFilenameList = [i_outputFilename]
    if (i_logFilename != None):
        i_writeFilenameList = [i_logFilename]
        
    i_readFilenameList = [i_vcfFilename, i_blatInputFilename, i_blatOutputFilename]
    
    if (not radiaUtil.check_for_argv_errors(None, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)

    # open the output stream
    i_outputFileHandler = None
    if (i_outputFilename != None):
        i_outputFileHandler = get_write_fileHandler(i_outputFilename)
    
    # get the BLAT results
    i_blatCoordinateDict = parse_blat_output(i_blatOutputFilename, i_blatOutputFormat, i_debug)
    
    # get the VCF generator   
    i_vcfGenerator  = get_vcf_data(i_vcfFilename, i_passedVCFCallsOnlyFlag, i_debug)
    
    for (vcfChr, vcfStopCoordinate, vcfId, vcfRef, vcfAlt, vcfScore, vcfFilterSet, vcfInfoDict, restOfLine) in i_vcfGenerator:
        if (i_debug):
            logging.debug("VCF Data: %s %s %s %s %s %s %s %s %s", vcfChr, str(vcfStopCoordinate), vcfId, vcfRef, vcfAlt, vcfScore, str(vcfFilterSet), str(vcfInfoDict), restOfLine)      
           
        modTypes = vcfInfoDict["MT"]
        modTypeFilters = dict()
        atLeastOnePass = False
        for modType in modTypes:
            
            blatHitsDict = dict()
            blatOverallReadDepth = 0
            numValidReads = 0
    
            if (modType == "NOR_EDIT" and i_blatRnaNormalReads):
                if ("rnaNormal" in i_blatCoordinateDict[vcfChr + "_" + str(vcfStopCoordinate)]):
                    # for each coordinate, get a dict of reads and corresponding blat hits
                    blatHitsDict = i_blatCoordinateDict[vcfChr + "_" + str(vcfStopCoordinate)]["rnaNormal"]
            elif ((modType == "SOM" or modType == "TUM_EDIT") and i_blatRnaTumorReads):
                if ("rnaTumor" in i_blatCoordinateDict[vcfChr + "_" + str(vcfStopCoordinate)]):
                    # for each coordinate, get a dict of reads and corresponding blat hits
                    blatHitsDict = i_blatCoordinateDict[vcfChr + "_" + str(vcfStopCoordinate)]["rnaTumor"]
                
            # for each read, investigate the blat hits to see if this read is valid
            for (readId, blatHitList) in blatHitsDict.iteritems():
                if (i_debug):
                    logging.debug("num of blat hits for read %s=%s", readId, len(blatHitList))
                
                blatOverallReadDepth +=1
    
                # find out if the read is valid or if it maps to other places in the genome
                if (i_blatOutputFormat == "PSL"):
                    (isValidRead, validRead) = is_valid_read_psl_format(blatHitList, vcfChr, vcfStopCoordinate, i_debug)
                elif (i_blatOutputFormat == "BLAST"):
                    (isValidRead, validRead) = is_valid_read_blast_format(blatHitList, vcfChr, vcfStopCoordinate, 0, i_debug)
                    #(isValidRead, validRead) = is_valid_read_blast_format(blatHitList, vcfChr, vcfStopCoordinate, 1, i_debug)
                    #(isValidRead, validRead) = is_valid_read_blast_format(blatHitList, vcfChr, vcfStopCoordinate, 2, i_debug)
                
                # if we have only one valid blat hit, then the read doesn't map to other places in the genome very well, so let's use it
                if (isValidRead):
                    numValidReads += 1
    
                    if (i_debug):
                        logging.debug("ValidRead: %s", validRead) 
            
            if (blatOverallReadDepth > 0):
                altPercent = round(numValidReads/float(blatOverallReadDepth),2)
            else:
                altPercent = 0.0
    
            if (numValidReads < i_readDepthCutoff or altPercent < i_readPercentCutoff):
                modTypeFilters[modType] = "blat"
            else:
                modTypeFilters[modType] = "PASS"
                atLeastOnePass = True
                
            if (i_debug):
                logging.debug("blatOverallReadDepth=%s, numValidReads=%s, altPercent=%s", str(blatOverallReadDepth), str(numValidReads), str(altPercent))
                logging.debug("modType=%s, passed? %s", modType, modTypeFilters[modType])
                logging.debug("blatFilter originalDepth=%s, afterBlatDepth=%s", str(blatOverallReadDepth), str(numValidReads))
            
        # make a copy of the list to manipulate
        modTypesTmpList = list(modTypes)
        modChanges = vcfInfoDict["MC"]
        # if at least one passed, then remove the ones that didn't
        for (modType, modChange) in izip(modTypes, modChanges):
            # if at least one passed, then remove the ones that didn't
            if (atLeastOnePass):
                if (modTypeFilters[modType] == "blat"):
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
                # if the call previous passed, then just set blat
                if (len(vcfFilterSet) == 1 and "PASS" in vcfFilterSet):
                    vcfFilterSet = ["blat"] 
                # otherwise, add it to the previous filters
                else:
                    vcfFilterSet.add("blat")
            # otherwise, just set the blat filter
            else:
                vcfFilterSet = ["blat"]
        
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
            
        output = [vcfChr, str(vcfStopCoordinate), vcfId, vcfRef, vcfAlt, vcfScore, ";".join(vcfFilterSet)]
        
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
            
    stopTime = time.time()  
    logging.info("Total time for Id %s: Total time=%s hrs, %s mins, %s secs", i_id, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime)) 
        
    # close the files 
    if (i_outputFilename != None):
        i_outputFileHandler.close()
        
    return
 

main()    
sys.exit(0)
