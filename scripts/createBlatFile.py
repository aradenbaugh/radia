#!/usr/bin/env python

import sys
import time
import re
import subprocess
from optparse import OptionParser
import radiaUtil
import collections
import logging
from itertools import izip
import os
import gzip

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

i_cigarRegEx = re.compile("[0-9]+[MIDNSHPX=]")


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


def get_vcf_data(aVcfFile, aHeaderFile, aPassOnlyFlag, anAltOnlyFlag, anIsDebug):
    '''
    ' This function reads from a .vcf input file and uses the python generator to yield the information
    ' one line at a time.  It ignores empty lines and strips trailing \r\n characters.  This function
    ' yields all the information from the VCF file.
    '
    ' aVcfFile:  A VCF file
    ' aPassOnlyFlag:  If all calls should be processed or only those calls that passed the filters thus far
    ' anAltOnlyFlag:  If all reads should be processed or only those that have the alternate allele
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    # open the header file
    fileHandler = get_read_fileHandler(aHeaderFile)
     
    for line in fileHandler:
          
        # here are some examples of .vcf lines:
        # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31.4;DP=53;FA=0.04;INDEL=0;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP      
        # GT:DP:INDEL:START:STOP:AD:AF:BQ:SB      0/0:2:0:0:0:2:1.0:35.5:0.0      0/0:1:0:0:0:1:1.0:39.0:1.0      0/1:50:0:2:0:48,2:0.96,0.04:31.65,17.5:0.75,0.5
        
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug):
            logging.debug("VCF Header: %s", line)    
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # if we find the column headers
        elif ("#CHROM" in line):
            columnsLine = line.lstrip("#")
            columnsList = columnsLine.split("\t")
            columnsList = columnsList[9:len(columnsList)]
        
        # if we find the vcfGenerator line, then create the dict of params
        elif ("vcfGenerator" in line):
            #generatorLine = line.rstrip(">")
            #generatorLine = generatorLine.lstrip("##vcfGenerator=<")
            generatorLine = line[0:(len(line)-1)]
            #print "generatorLine: %s", generatorLine
            generatorLine = generatorLine[16:len(generatorLine)]
            #print "generatorLine: %s", generatorLine
            generatorParamsList = generatorLine.split(",")
            generatorParamsDict = {}
            
            # create a dictionary of existing params
            for param in generatorParamsList:
                (key, value) = param.split("=")
                value = value.rstrip(">")
                value = value.lstrip("<")
                generatorParamsDict[key] = value
                
        # if we are done with the header, then stop    
        elif (not line.startswith("#")):
            break
        
    fileHandler.close()
    
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
        
        # if we find the column headers
        elif ("#CHROM" in line):
            columnsLine = line.lstrip("#")
            columnsList = columnsLine.split("\t")
            columnsList = columnsList[9:len(columnsList)]
        
        # if we find the vcfGenerator line, then create the dict of params
        elif ("vcfGenerator" in line):
            #generatorLine = line.rstrip(">")
            #generatorLine = generatorLine.lstrip("##vcfGenerator=<")
            generatorLine = line[0:(len(line)-1)]
            #print "generatorLine: %s", generatorLine
            generatorLine = generatorLine[16:len(generatorLine)]
            #print "generatorLine: %s", generatorLine
            generatorParamsList = generatorLine.split(",")
            generatorParamsDict = {}
            
            # create a dictionary of existing params
            for param in generatorParamsList:
                (key, value) = param.split("=")
                value = value.rstrip(">")
                value = value.lstrip("<")
                generatorParamsDict[key] = value
                
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
        idList = splitLine[2].split(";")
        refList = splitLine[3].split(",")
        altList = splitLine[4].split(",")
        score = float(splitLine[5])
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
        yield (chrom, stopCoordinate, idList, refList, altList, score, filterSet, infoDict, "\t".join(splitLine[8:]), generatorParamsDict)
    fileHandler.close()
    return


def get_read_data(aReadFile, anIsDebug):
    '''
    ' This function returns reads from a test file.
    '
    ' aReadFile:  An test file that has the output from a samtools view command
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    # open the file
    fileHandler = get_read_fileHandler(aReadFile)
    reads = fileHandler.readlines()
    fileHandler.close()
    
    return reads


def execute_samtools_cmd(aBamFile, aFastaFile, aMappingQuality, aChrom, aCoordinate, aUseChrPrefix, anIsDebug):
    '''
    ' This function executes an external command.  The command is the "samtools view" command which returns all 
    ' the information about the sequencing reads that overlap a specific coordinate.  Some .bam files use the 'chr' 
    ' prefix when specifying the region.  If the 'chr' prefix is required, then specify the --useChrPrefix argument.
    ' Here are some examples of the commands that can be copy/pasted to the command line to view the output:
    '
    ' samtools view -q 10 myBamfile.bam 10:8100500-8100500
    ' samtools view -q 10 myBamfile.bam chr10:8100500-8100500
    '
    ' aBamFile:            A .bam file to be read from
    ' aFastaFile:          A fasta file that can be used in the samtools view command
    ' aMappingQuality:     The mapping quality score for the samtools command
    ' aChrom:              The chromosome that we are selecting from
    ' aCoordinate:         The coordinate of the selection
    ' aUseChrPrefix:       Whether the 'chr' should be used in the samtools command
    ' anIsDebug:           A flag for outputting debug messages to STDERR
    '''
  
    if (not os.path.isfile(aBamFile)):
        logging.critical("The BAM file specified in the VCF header does not exist: %s", aBamFile)
        sys.exit(1)

    if (aUseChrPrefix):
        # create the samtools command
        samtoolsSelectStatement = "samtools view -q " + str(aMappingQuality) + " " + aBamFile + " chr" + aChrom + ":" + str(aCoordinate) + "-" + str(aCoordinate)
    else:
        # create the samtools command
        samtoolsSelectStatement = "samtools view -q " + str(aMappingQuality) + " " + aBamFile + " " + aChrom + ":" + str(aCoordinate) + "-" + str(aCoordinate)
    
    # keep track of how long it takes to run the samtools command
    if (anIsDebug):
        timeSamtoolsStart = time.time()
        logging.debug(samtoolsSelectStatement)
    
    # execute the samtools command
    samtoolsCall = subprocess.Popen(samtoolsSelectStatement, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    # communicate() waits for the process to finish
    (pileups, samtoolsStdErr) = samtoolsCall.communicate()
    # wait for the child process to terminate and set the returncode
    #samtoolsStdErr = samtoolsCall.wait()
    #samtoolsCall.poll()
    
    if (anIsDebug):
        timeSamtoolsEnd = time.time()
        logging.debug("Time spent executing samtools command: %s", str(timeSamtoolsEnd-timeSamtoolsStart)) 
    
    # get the output from executing the samtools command
    #pileups = samtoolsCall.stdout.readlines()
    
    # for some reason, we have to get the stdout.readlines() before the stderr.readlines(), otherwise the process hangs
    if (samtoolsCall.returncode != 0):
        logging.error("Error from %s\n %s", samtoolsSelectStatement, samtoolsStdErr)
        #logging.error("Error from %s\n %s", samtoolsSelectStatement, samtoolsCall.stderr.readlines())
        sys.exit(1)
        
    # kill the command
    #samtoolsCall.kill()
    
    #return pileups
    return pileups.split("\n")
   
   
def write_to_blat_file(aBlatFileHandler, aChr, aStopCoordinate, aParamsDict, anInfoDict, aPrefix, anIsDebug):
    '''
    ' This function gets all of the reads at a specific coordinate and creates a BLAT input file.
    '
    ' aBlatFileHandler:    A file handler where all the BLAT query data is written
    ' aChr:                The chromosome that we are selecting from
    ' aStopCoordinate:     The coordinate of the selection
    ' aParamsDict:         The parameters from the header
    ' anInfoDict:          A dictionary of the INFO column
    ' aParamPrefix:        The prefix for the parameters dictionary
    ' aNamePrefix:         A prefix to be added to the name (usually specifying whether it came from DNA or RNA)
    ' anIsDebug:           A flag for outputting debug messages to STDERR
    '''
    
    # get the info for executing the samtools command
    bamFile = aParamsDict[aPrefix + "Filename"]
    fastaFile = aParamsDict[aPrefix + "FastaFilename"]
    mappingQuality = aParamsDict[aPrefix + "MinMappingQuality"]
    useChrPrefixString = aParamsDict[aPrefix + "UseChrPrefix"]
    
    if (useChrPrefixString == "True"):
        useChrPrefix = True
    else:
        useChrPrefix = False

    if (not os.path.isfile(bamFile)):
        logging.critical("The BAM file specified in the VCF header does not exist: %s", bamFile)
        sys.exit(1)

    if (not os.path.isfile(fastaFile)):
        logging.critical("The FASTA file specified in the VCF header does not exist: %s", fastaFile)
        sys.exit(1)

    # execute the samtools command
    reads = execute_samtools_cmd(bamFile, fastaFile, mappingQuality, aChr, aStopCoordinate, useChrPrefix, anIsDebug)
    #reads = get_read_data(aBamFile, anIsDebug)
    
    if (anIsDebug):
        logging.debug("samtools number of reads selected from chr%s:%s=%s", aChr, aStopCoordinate, len(reads))

    modChanges = anInfoDict["MC"]
    modTypes = anInfoDict["MT"]
    altSet = []
    refSet = []
    for (modChange, modType) in izip(modChanges, modTypes):
        (ref, alt) = modChange.split(">")
        altSet.append(alt)
        refSet.append(ref)
    
    baseCounter = collections.defaultdict(int)
    numMatchesFound = 0
    numAltsFound = 0
    numSkipsFound = 0
    numNonRefAltBasesFound = 0
    numMsFound = 0
    numIsFound = 0
    numDsFound = 0
    numNsFound = 0
    numHsFound = 0
    numPsFound = 0
    numEqualsFound = 0
    numXsFound = 0
    numSsFound = 0
    
    # if there was data for this coordinate
    if (len(reads) > 0):
        # for each line representing one coordinate
        for line in reads:
            # if the samtools select statement returns no reads which can happen when the batch size is
            # small and the selection is done in an area with no reads, then a warning message will be
            # returned that starts with [main_samview], [sam_header_read2], or [fai_build_core].  We can 
            # ignore the message and move on to the next select statement.
            if (line == "" or line.isspace() or line.startswith("[")):
                continue;

            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")
            
            if (anIsDebug):
                logging.debug("samtools read: %s", line)

            # split the line on the tab
            splitLine = line.split("\t")

            readName = splitLine[0]
            readFlag = splitLine[1]
            
            # if this read is properly paired
            readPaired = int(readFlag) & 0x1
            properlyPaired = int(readFlag) & 0x2
            #unmapped = int(readFlag) & 0x4
            #mateUnmapped = int(readFlag) & 0x8
            #secondaryAlignment = int(readFlag) & 0x100
            #notPassingQC = int(readFlag) & 0x200
            #pcrDup = int(readFlag) & 0x400
            #if (readPaired and properlyPaired and not unmapped and not mateUnmapped and not secondaryAlignment and not notPassingQC and not pcrDup):
            
            if ((readPaired and properlyPaired) or (not readPaired)):
                
                if (anIsDebug):
                    logging.debug("read properly paired %s", readFlag)
                
                #readChr = splitLine[2]
                readStartCoordinate = int(splitLine[3])
                readMapQuality = int(splitLine[4])
                
                if (anIsDebug):
                    logging.debug("readMapQual=%s", readMapQuality)
                
                if (readMapQuality >= 10):
                    readCigar = splitLine[5]
                    #readMateChr = splitLine[6]
                    #readMateStartCoordinate = int(splitLine[7])
                    readInsertSize = int(splitLine[8])
                    readSequence = splitLine[9]
                    readQualScores = splitLine[10]
                    
                    cigarIter = i_cigarRegEx.finditer(readCigar)
                    cigarList = []
                    skipsBeforeFirstMatch = 0
                    foundFirstMatch = False
                    for regEx in cigarIter:
                        match = regEx.group()
                        number = int(match[0:len(match)-1])
                        #print match, number, match[len(match)-1:]
                        
                        # M alignment match (can be match or mismatch)
                        # I insertion to the ref
                        # D deletion from the ref
                        # N skipped region from the ref
                        # S soft-clipping
                        # H hard-clipping
                        # P padding
                        # = sequence match
                        # X sequence mismatch
                        if ("M" in match):
                            # alignment match (can be sequence match or mismatch)
                            cigarList += ["M"] * number
                            numMsFound += 1
                            foundFirstMatch = True
                        elif ("I" in match):
                            # insertion to the reference
                            cigarList += ["S"] * number
                            numIsFound += 1
                        elif ("D" in match):
                            # deletion from the reference
                            cigarList += ["S"] * number
                            numDsFound += 1
                        elif ("N" in match):
                            # for mRNA, N represents intronic sequences
                            cigarList += ["S"] * number
                            numNsFound += 1
                        elif ("H" in match):
                            # hard-clipping where the clipped sequence is not a part of our sequence,
                            # so no need to skip our sequence
                            #cigarList += ["S"] * number
                            numHsFound += 1
                        elif ("P" in match):
                            # silent deletion from padded reference
                            cigarList += ["S"] * number
                            numPsFound += 1
                        elif ("=" in match):
                            # sequence match
                            cigarList += ["M"] * number
                            numEqualsFound += 1
                        elif ("X" in match):
                            # sequence mismatch
                            cigarList += ["M"] * number
                            numXsFound += 1
                        elif ("S" in match):
                            # soft-clipping where the clipped sequence is part of our sequence
                            # this happens sometimes where the first few bases of our sequence are
                            # clipped off and then the rest is aligned
                            cigarList += ["S"] * number
                            numSsFound += 1
                        
                        # if we haven't found the first "M" yet, then count the number of skips before the first "M"
                        if (not foundFirstMatch):
                            skipsBeforeFirstMatch += number
                            
                    #print readCigar, cigarList
            
                    # the VCF file is 1-based
                    # the samtools view returns reads that are 1-based.
                    # The SAM, GFF and Wiggle formats are using the 1-based coordinate system.
                    # The BAM, BED, and PSL formats are using the 0-based coordinate system.
                    # when using pysam:  the readSequence and readQualScores sequences are 0-based, so we need to subtract 1 to get the index
                    
                    # the readStartCoordinate corresponds to the first M, so if the cigar starts with skips, we need to add them here
                    cigarListIndex = aStopCoordinate-readStartCoordinate+skipsBeforeFirstMatch
                    
                    if (anIsDebug):
                        logging.debug("cigarListIndex=%s, cigar=%s, skipsBeforeFirstMatch=%s, lenCigarList=%s", str(cigarListIndex), readCigar, str(skipsBeforeFirstMatch), str(len(cigarList)))
                        logging.debug("cigarListIndex=%s, cigar=%s, cigarAtIndex=%s", str(cigarListIndex), readCigar, cigarList[cigarListIndex])
                                        
                    if (cigarList[cigarListIndex] == "M" or cigarList[cigarListIndex] == "X" or cigarList[cigarListIndex] == "="):
                            
                        # we only want to count the number of skips after the skipsBeforeFirstMatch to the cigarListIndex
                        numOfSkips = cigarList[skipsBeforeFirstMatch:cigarListIndex].count("S")
                        # the sequence index is the index of the variant position in the sequence, 
                        sequenceIndex = aStopCoordinate-readStartCoordinate-numOfSkips+skipsBeforeFirstMatch    
                        
                        if (anIsDebug):
                            logging.debug("cigarListIndex=%s, numOfSkips=%s, skipsBeforeFirstMatch=%s, sequenceIndex=%s, lenReadSequence=%s", str(sequenceIndex), str(numOfSkips), str(skipsBeforeFirstMatch), str(sequenceIndex), len(readSequence))
                            logging.debug("cigarListIndex=%s, numOfSkips=%s, skipsBeforeFirstMatch=%s, sequenceIndex=%s, baseAtIndex=%s", str(sequenceIndex), str(numOfSkips), str(skipsBeforeFirstMatch), str(sequenceIndex), str(readSequence[sequenceIndex]))
                            
                        #print "bases", readSequence[index], readSequence[otherIndex] 
                        base = readSequence[sequenceIndex]
                        numMatchesFound += 1
                        
                        # count the number of time a base is not the ref nor the alt
                        if (base not in refSet and base not in altSet):
                            numNonRefAltBasesFound += 1
                        
                        if (base in altSet):
                            baseQuality = readQualScores[sequenceIndex]
                            baseQualityConverted = ord(baseQuality)-33
                            if (anIsDebug):
                                logging.debug("baseQualConverted=%s", str(baseQualityConverted))
                            
                            if (baseQualityConverted >= 10):
                                numAltsFound += 1
                                if (anIsDebug):
                                    logging.debug("found an alt for stopCoordinate=%s, base=%s, baseQual=%s, readMapQual=%s", str(aStopCoordinate), base, str(ord(baseQuality)-33), str(readMapQuality))
                                
                                baseCounter[base] += 1
                                    
                                # position in read
                                readLength = len(readSequence)
                                if (sequenceIndex/float(readLength) <= 0.33):
                                    position = "start"
                                elif (sequenceIndex/float(readLength) <= 0.66):
                                    position = "middle"
                                else:
                                    position = "end"
                            
                                outputList = [aPrefix, aChr, str(aStopCoordinate), readName.replace("_", ""), base, str(ord(baseQuality)-33), str(readMapQuality), position, readFlag, str(readInsertSize), str(readLength)]
                                if (aBlatFileHandler != None):
                                    aBlatFileHandler.write("> " + "_".join(outputList) + "\n")
                                    aBlatFileHandler.write(readSequence + "\n")
                                else:
                                    print >> sys.stdout, "> " + "_".join(outputList)
                                    print >> sys.stdout, readSequence     
                    else:
                        numSkipsFound += 1   
                        
    if (anIsDebug):  
        logging.debug("For stopCoordinate=%s, numSkipsFound=%s, numSkipsBeforeFirstMatch=%s, numNonRefAltBasesFound=%s, numAltsFound=%s, numTotalMatches=%s", aStopCoordinate, numSkipsFound, skipsBeforeFirstMatch, numNonRefAltBasesFound, numAltsFound, numMatchesFound)
        logging.debug("numMs=%s, numSs=%s, numIs=%s, numDs=%s, numHs=%s, numPs=%s, numNs=%s, numXs=%s, numEquals=%s", numMsFound, numSsFound, numIsFound, numDsFound, numHsFound, numPsFound, numNsFound, numXsFound, numEqualsFound)
        
    return 
    
               
def main():
    
    # command for running this on a small test case: 
    #python createBlatFile.py TCGA-00-4454 7 ../data/test/TCGA-00-4454_EGFR.vcf ../data/test/tmp/ --dnaNormalFilename=../data/test/TCGA-00-4454_EGFR.reads
 
    startTime = time.time()
    
    # create the usage statement
    usage = "usage: python %prog id vcfFile headerFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    # add the optional parameters
    i_cmdLineParser.add_option("-c", "--allVCFCalls", action="store_false", default=True, dest="passedVCFCallsOnly", help="by default only the VCF calls that have passed all filters thus far are processed, include this argument if all of the VCF calls should be processed")
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, STDOUT by default")
    i_cmdLineParser.add_option("-b", "--allReadBases", action="store_false", default=True, dest="altBasesOnly", help="by default only the reads with the alternate base are processed, include this argument if all of the reads should be processed")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    
    i_cmdLineParser.add_option("-n", "--blatDnaNormalReads", action="store_true", default=False, dest="blatDnaNormalReads", help="include this argument if the normal DNA reads should be processed")
    i_cmdLineParser.add_option("-x", "--blatRnaNormalReads", action="store_true", default=False, dest="blatRnaNormalReads", help="include this argument if the normal RNA reads should be processed")
    i_cmdLineParser.add_option("-t", "--blatDnaTumorReads", action="store_true", default=False, dest="blatDnaTumorReads", help="include this argument if the tumor DNA reads should be processed")
    i_cmdLineParser.add_option("-r", "--blatRnaTumorReads", action="store_true", default=False, dest="blatRnaTumorReads", help="include this argument if the tumor RNA reads should be processed")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(3,22,1)
    i_argLength = len(sys.argv)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)
    
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_id = i_cmdLineArgs[0]
    i_vcfFilename = i_cmdLineArgs[1]
    i_headerFilename = i_cmdLineArgs[2]
    
    # get the optional params with default values
    i_passedVCFCallsOnlyFlag = i_cmdLineOptions.passedVCFCallsOnly
    i_altBasesOnlyFlag = i_cmdLineOptions.altBasesOnly
    i_logLevel = i_cmdLineOptions.logLevel
    
    i_blatDnaNormalReads = i_cmdLineOptions.blatDnaNormalReads
    i_blatDnaTumorReads = i_cmdLineOptions.blatDnaTumorReads
    i_blatRnaNormalReads = i_cmdLineOptions.blatRnaNormalReads
    i_blatRnaTumorReads = i_cmdLineOptions.blatRnaTumorReads
    
    # try to get any optional parameters with no defaults    
    i_readFilenameList = [i_vcfFilename, i_headerFilename]
    i_writeFilenameList = []
    
    i_logFilename = None
    i_outputFilename = None
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
        i_writeFilenameList += [i_logFilename]
    if (i_cmdLineOptions.outputFilename != None):
        i_outputFilename = str(i_cmdLineOptions.outputFilename)
        i_writeFilenameList += [i_outputFilename]
        
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
        
    # set the debug flag    
    i_debug = (i_numericLogLevel == logging.DEBUG)

    # output some debug info
    if (i_debug):
        logging.debug("id=%s", i_id)
        logging.debug("vcfFilename=%s", i_vcfFilename)
        logging.debug("headerFilename=%s", i_headerFilename)
        logging.debug("outputFilename=%s", i_outputFilename)
        
        logging.debug("passedCallsOnly? %s", i_passedVCFCallsOnlyFlag)
        logging.debug("altBasesOnlyFlag? %s", i_altBasesOnlyFlag)
        
        logging.debug("blatDnaNormal? %s", i_blatDnaNormalReads)
        logging.debug("blatDnaTumor? %s", i_blatDnaTumorReads)
        logging.debug("blatRnaNormal? %s", i_blatRnaNormalReads)
        logging.debug("blatRnaTumor? %s", i_blatRnaTumorReads)
                    
    if (not radiaUtil.check_for_argv_errors(None, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)
        
    # open the output stream
    i_outputFileHandler = None
    if (i_outputFilename != None):
        i_outputFileHandler = get_write_fileHandler(i_outputFilename)
                
    # get the VCF generator
    i_vcfGenerator  = get_vcf_data(i_vcfFilename, i_headerFilename, i_passedVCFCallsOnlyFlag, i_altBasesOnlyFlag, i_debug)    
   
    # for each VCF call that should be investigated   
    for (vcfChr, vcfStopCoordinate, vcfId, vcfRef, vcfAlt, vcfScore, vcfFilterSet, vcfInfoDict, restOfLine, vcfParamsDict) in i_vcfGenerator:
        if (i_debug):
            logging.debug("VCF Data: %s %s %s %s %s %s %s %s %s", vcfChr, str(vcfStopCoordinate), vcfId, vcfRef, vcfAlt, vcfScore, str(vcfFilterSet), str(vcfInfoDict), restOfLine) 
        
        modTypes = vcfInfoDict["MT"]
        for modType in modTypes:
            # get the reads contributing to a call and put them in a blat query file
            if (i_blatDnaNormalReads):
                write_to_blat_file(i_outputFileHandler, vcfChr, vcfStopCoordinate, vcfParamsDict, vcfInfoDict, "dnaNormal", i_debug)                      
                
            if (modType == "NOR_EDIT" and i_blatRnaNormalReads):
                write_to_blat_file(i_outputFileHandler, vcfChr, vcfStopCoordinate, vcfParamsDict, vcfInfoDict, "rnaNormal", i_debug)
            
            if (i_blatDnaTumorReads):
                write_to_blat_file(i_outputFileHandler, vcfChr, vcfStopCoordinate, vcfParamsDict, vcfInfoDict, "dnaTumor", i_debug)    
                
            if ((modType == "SOM" or modType == "TUM_EDIT") and i_blatRnaTumorReads):
                write_to_blat_file(i_outputFileHandler, vcfChr, vcfStopCoordinate, vcfParamsDict, vcfInfoDict, "rnaTumor", i_debug)
            
    stopTime = time.time()       
    logging.info("Id %s: Total time=%s hrs, %s mins, %s secs", i_id, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime))         
        
    # close the files 
    if (i_outputFilename != None):
        i_outputFileHandler.close()
        
    return
 

main()    
sys.exit(0)
