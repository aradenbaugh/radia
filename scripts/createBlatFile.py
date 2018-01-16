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

i_reverseCompDict = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

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
    
    
def reverse_complement_nucleotide(aNucleotide):
    '''
    ' This function returns the reverse complement of the parameter aNucleotide
    '
    ' aNucleotide:    A nucleotide to reverse complement
    '''
    if aNucleotide in i_reverseCompDict:
        return i_reverseCompDict[aNucleotide]
    else:
        logging.error("Trying to reverse complement an unknown nucleotide: %s", aNucleotide)
        sys.exit(1)
        
    return None


def get_vcf_data(aVcfFile, aHeaderFile, aPassOnlyFlag, anIsDebug):
    '''
    ' This function reads from a .vcf input file and uses the python generator to yield the information
    ' one line at a time.  It ignores empty lines and strips trailing \r\n characters.  This function
    ' yields all the information from the VCF file.
    '
    ' aVcfFile:  A VCF file
    ' aPassOnlyFlag:  If all calls should be processed or only those calls that passed the filters thus far
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
        
        #if (anIsDebug):
        #    logging.debug("VCF Header: %s", line)    
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # if we find the column headers
        elif ("#CHROM" in line):
            columnsLine = line.lstrip("#")
            columnsList = columnsLine.split("\t")
            columnsList = columnsList[9:len(columnsList)]
            continue
        
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
            continue
        
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
        
        #if (anIsDebug):
        #    logging.debug("VCF: %s", line)    
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # if we find the column headers
        elif ("#CHROM" in line):
            columnsLine = line.lstrip("#")
            columnsList = columnsLine.split("\t")
            columnsList = columnsList[9:len(columnsList)]
            continue
        
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
            continue
                
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


def execute_samtools_cmd(aBamFile, aMappingQuality, aChrom, aCoordinate, aUseChrPrefix, anRnaIncludeSecondaryAlignmentsFlag, anIsDebug):
    '''
    ' This function executes an external command.  The command is the "samtools view" command which returns all 
    ' the information about the sequencing reads that overlap a specific coordinate.  Some .bam files use the 'chr' 
    ' prefix when specifying the region.  If the 'chr' prefix is required, then specify the --useChrPrefix argument.
    ' Here are some examples of the commands that can be copy/pasted to the command line to view the output:
    '
    ' samtools view -q 10 myBamfile.bam 10:8100500-8100500
    ' samtools view -q 10 myBamfile.bam chr10:8100500-8100500
    '
    ' aBamFile:                              A .bam file to be read from
    ' aMappingQuality:                       The mapping quality score for the samtools command
    ' aChrom:                                The chromosome or transcript name that we are selecting from
    ' aCoordinate:                           The coordinate of the selection
    ' aUseChrPrefix:                         Whether the 'chr' should be used in the samtools command
    ' anRnaIncludeSecondayAlignmentsFlag:    If you align the RNA to transcript isoforms, then you may want to include RNA secondary alignments in the samtools command
    ' anIsDebug:                             A flag for outputting debug messages to STDERR
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
    
    # for the samtools mpileup command, the flags are:
    # --ff (exclude flags):  unmapped reads, reads that fail quality checks, pcr duplicates
    # --rf (include flags):  everything else (including secondary alignments)
    # --ff, --excl-flags STR|INT  filter flags: skip reads with mask bits set [UNMAP,SECONDARY,QCFAIL,DUP]
    # --rf, --incl-flags STR|INT  required flags: skip reads with mask bits unset []
    
    # for the samtools view command, the flags are:
    # -F (exclude flags):  unmapped reads, reads that fail quality checks, pcr duplicates
    # -f (include flags):  everything else (including secondary alignments)
    # -F INT   only include reads with none of the bits set in INT set in FLAG [0]
    # -f INT   only include reads with all bits set in INT set in FLAG [0]
    if (anRnaIncludeSecondaryAlignmentsFlag):
        # for the samtools mpileup command:
        #samtoolsSelectStatement += " --ff 1540 --rf 2555"
        # for the samtools view command, we only need the -F flag
        samtoolsSelectStatement += " -F 1540"
    
    # keep track of how long it takes to run the samtools command
    if (anIsDebug):
        timeSamtoolsStart = time.time()
        logging.debug(samtoolsSelectStatement)
    
    # execute the samtools command
    samtoolsCall = subprocess.Popen(samtoolsSelectStatement, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    # communicate() waits for the process to finish
    (pileups, samtoolsStdErr) = samtoolsCall.communicate()
    
    if (anIsDebug):
        timeSamtoolsEnd = time.time()
        logging.debug("Time spent executing samtools command: %s", str(timeSamtoolsEnd-timeSamtoolsStart)) 
    
    # for some reason, we have to get the stdout.readlines() before the stderr.readlines(), otherwise the process hangs
    if (samtoolsCall.returncode != 0):
        logging.error("Error from '%s':\n %s", samtoolsSelectStatement, samtoolsStdErr)
        sys.exit(1)
    # if the stdErr is not empty, then just warn the user
    elif samtoolsStdErr:
        logging.warning("Message from '%s':\n %s", samtoolsSelectStatement, samtoolsStdErr)
        
    return pileups.split("\n")


def group_reads_by_name(aReadsList, aMinMapQual, aGenomicCoordinate, aTranscriptName, aTranscriptCoordinate, aTranscriptStrand, anIsDebug):
    '''
    ' This function loops through the reads that were returned from the "samtools view" command and
    ' creates a dictionary where the reads are grouped by their read name.  Reads with the following 
    ' conditions are ignored:
    '    - unmapped
    '    - not passing QC
    '    - PCR duplicate
    '    - paired but not properly paired
    '    - mapping quality < aMinMapQual
    '
    ' aReadsList:               The list of reads returned from the samtools view command
    ' aMinMapQuality:           A minimum mapping quality score for the reads
    ' aGenomicCoordinate:       The genomic coordinate
    ' aTranscriptName:          The transcript name
    ' aTranscriptCoordinate:    The transcript coordinate
    ' aTranscriptStrand:        The transcript strand
    ' anIsDebug:                A flag for outputting debug messages to STDERR
    '''
    
    readsDict = collections.defaultdict(list)
    numReadsSkipped = 0
    samtoolsViewBases = ""
    samtoolsViewQuals = ""
    
    # if there was data for this coordinate
    if (len(aReadsList) > 0):
        # for each line representing one coordinate
        for line in aReadsList:
            # if the samtools view statement returns no reads which can happen when the batch size is
            # small and the selection is done in an area with no reads, then a warning message will be
            # returned that starts with [main_samview], [sam_header_read2], or [fai_build_core].  We can 
            # ignore the message and move on to the next view statement.
            if (line == "" or line.isspace() or line.startswith("[")):
                continue;

            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")
            
            if (anIsDebug):
                logging.debug("samtools read: %s", line)

            # split the line on the tab
            splitLine = line.split("\t")

            readName = splitLine[0]
            readFlag = int(splitLine[1])
            #readRefName = splitLine[2]
            readStart = int(splitLine[3])
            readMapQual = int(splitLine[4])
            readCigar = splitLine[5]
            #readMateRef = splitLine[6]
            readMateStart = int(splitLine[7])
            readInsertSize = int(splitLine[8])
            readSequence = splitLine[9]
            readQualScores = splitLine[10]
                        
            # if this read is properly paired
            readPaired = readFlag & 0x1
            properlyPaired = readFlag & 0x2
            unmapped = readFlag & 0x4
            mateUnmapped = readFlag & 0x8
            seqReversed = readFlag & 0x10
            mateSeqReversed = readFlag & 0x20
            firstSegment = readFlag & 0x40
            lastSegment = readFlag & 0x80
            secondaryAlignment = readFlag & 0x100
            notPassingQC = readFlag & 0x200
            pcrDup = readFlag & 0x400
            supplementaryAlignment = readFlag & 0x800
            #if (readPaired and properlyPaired and not unmapped and not mateUnmapped and not secondaryAlignment and not notPassingQC and not pcrDup):
            
            if (anIsDebug):
                logging.debug("readPaired=%s, properlyPaired=%s, unmapped=%s, mateUnmapped=%s, seqReversed=%s, mateSeqReversed=%s, firstSeg=%s, lastSeg=%s, secondaryAlignment=%s, notPassingQC=%s, pcrDup=%s, suppAlign=%s", readPaired, properlyPaired, unmapped, mateUnmapped, seqReversed, mateSeqReversed, firstSegment, lastSegment, secondaryAlignment, notPassingQC, pcrDup, supplementaryAlignment)
            
            # don't process these reads
            if (unmapped or notPassingQC or pcrDup):
                numReadsSkipped += 1
                continue
            
            # MapSplice doesn't set the proper pair flag for RNA-Seq reads, so only do this for DNA reads
            if ((readPaired and properlyPaired) or (not readPaired)):
                if (anIsDebug):
                    logging.debug("read properly paired or not paired at all %s", readFlag)
                    logging.debug("readMapQual=%s, minMapQual=%s", readMapQual, aMinMapQual)
                
                # if the mapping quality is good
                if (readMapQual >= int(aMinMapQual)):
                    # create an iterator over the cigar string for this read
                    cigarIter = i_cigarRegEx.finditer(readCigar)
                    cigarList = []
                    skipsBeforeFirstMatch = 0
                    foundFirstMatch = False
                    alignedLength = 0
                   
                    numMsFound = 0
                    numIsFound = 0
                    numDsFound = 0
                    numNsFound = 0
                    numHsFound = 0
                    numPsFound = 0
                    numEqualsFound = 0
                    numXsFound = 0
                    numSsFound = 0
                                    
                    # for each item in the cigar string
                    for regEx in cigarIter:
                        
                        cigar = regEx.group()
                        number = int(cigar[0:len(cigar)-1])
                        alignedLength += number
                        
                        if (anIsDebug):
                            logging.debug("cigar=%s, number=%s, cigarAtNumber=%s", cigar, number, cigar[len(cigar)-1:])
                        
                        # M alignment match (can be match or mismatch)
                        # I insertion to the ref
                        # D deletion from the ref
                        # N skipped region from the ref
                        # S soft-clipping
                        # H hard-clipping
                        # P padding
                        # = sequence match
                        # X sequence mismatch
                        if ("M" in cigar):
                            # alignment match (can be sequence match or mismatch)
                            cigarList += ["M"] * number
                            numMsFound += 1
                            foundFirstMatch = True
                        elif ("I" in cigar):
                            # insertion to the reference
                            cigarList += ["S"] * number
                            numIsFound += 1
                        elif ("D" in cigar):
                            # deletion from the reference
                            cigarList += ["S"] * number
                            numDsFound += 1
                        elif ("N" in cigar):
                            # for mRNA, N represents intronic sequences
                            cigarList += ["S"] * number
                            numNsFound += 1
                        elif ("H" in cigar):
                            # hard-clipping where the clipped sequence is not a part of our sequence,
                            # so no need to skip our sequence
                            #cigarList += ["S"] * number
                            numHsFound += 1
                        elif ("P" in cigar):
                            # silent deletion from padded reference
                            cigarList += ["S"] * number
                            numPsFound += 1
                        elif ("=" in cigar):
                            # sequence match
                            cigarList += ["M"] * number
                            numEqualsFound += 1
                        elif ("X" in cigar):
                            # sequence mismatch
                            cigarList += ["M"] * number
                            numXsFound += 1
                        elif ("S" in cigar):
                            # soft-clipping where the clipped sequence is part of our sequence
                            # this happens sometimes where the first few bases of our sequence are
                            # clipped off and then the rest is aligned
                            cigarList += ["S"] * number
                            numSsFound += 1
                        
                        # if we haven't found the first "M" yet, then count the number of skips before the first "M"
                        if (not foundFirstMatch) and ("H" not in cigar):
                            skipsBeforeFirstMatch += number
                     
                    if (anIsDebug):       
                        logging.debug("readCigar=%s, cigarList=%s", readCigar, cigarList)
            
                    # the VCF file is 1-based
                    # the samtools view returns reads that are 1-based.
                    # pysam is 0-based.
                    # The SAM, GFF and Wiggle formats are using the 1-based coordinate system.
                    # The BAM, BED, and PSL formats are using the 0-based coordinate system.
                    
                    # if we are processing the RNA and
                    # the transcript name and coordinate should be used instead of the genomic chrom and coordinate
                    if (aTranscriptName != None and aTranscriptCoordinate != None):
                        # the readStart corresponds to the first M, so if the cigar starts with skips, we need to add them here
                        cigarListIndex = aTranscriptCoordinate - readStart + skipsBeforeFirstMatch
                    else:
                        # the readStart corresponds to the first M, so if the cigar starts with skips, we need to add them here
                        cigarListIndex = aGenomicCoordinate - readStart + skipsBeforeFirstMatch
                 
                    if (anIsDebug):
                        logging.debug("numMs=%s, numSs=%s, numIs=%s, numDs=%s, numHs=%s, numPs=%s, numNs=%s, numXs=%s, numEquals=%s", numMsFound, numSsFound, numIsFound, numDsFound, numHsFound, numPsFound, numNsFound, numXsFound, numEqualsFound)
                        logging.debug("cigar=%s, cigarListIndex=%s, cigarAtIndex=%s, skipsBeforeFirstMatch=%s, lenCigarList=%s", readCigar, str(cigarListIndex), cigarList[cigarListIndex], str(skipsBeforeFirstMatch), str(len(cigarList)))
                        
                    if (cigarList[cigarListIndex] == "M" or cigarList[cigarListIndex] == "X" or cigarList[cigarListIndex] == "="):
                            
                        # we only want to count the number of skips after the skipsBeforeFirstMatch to the cigarListIndex
                        numOfSkips = cigarList[skipsBeforeFirstMatch:cigarListIndex].count("S")
                        # if we are processing the RNA and
                        # the transcript name and coordinate should be used instead of the genomic chrom and coordinate
                        if (aTranscriptName != None and aTranscriptCoordinate != None):
                            # the sequence index is the index of the variant position in the sequence, 
                            sequenceIndex = aTranscriptCoordinate - readStart - numOfSkips + skipsBeforeFirstMatch    
                        else:
                            # the sequence index is the index of the variant position in the sequence, 
                            sequenceIndex = aGenomicCoordinate - readStart - numOfSkips + skipsBeforeFirstMatch
                        
                        base = readSequence[sequenceIndex]
                        baseQual = readQualScores[sequenceIndex]
                        
                        # is the transcript fasta strand was on the reverse, then reverse comp the base
                        if (aTranscriptStrand != None and aTranscriptStrand == "-"):
                            base = reverse_complement_nucleotide(base)
                            
                        if (anIsDebug):
                            logging.debug("cigarListIndex=%s, alignedBaseAtIndex=%s, numOfSkips=%s, skipsBeforeFirstMatch=%s, alignedLength=%s, lenReadSequence=%s, ", str(sequenceIndex), base, str(numOfSkips), str(skipsBeforeFirstMatch), str(alignedLength), len(readSequence))
                        
                        samtoolsViewBases += base
                        samtoolsViewQuals += baseQual
                        
                        # keep a dictionary of all reads, using the readName as the key
                        # due to the inclusion of secondary alignments for RNA-Seq, there could be more than 2 reads that are paired
                        oneReadDict = {}
                        oneReadDict["name"] = readName
                        oneReadDict["read"] = line
                        oneReadDict["base"] = base
                        oneReadDict["baseQual"] = baseQual
                        oneReadDict["mapQual"] = readMapQual
                        oneReadDict["flag"] = readFlag
                        oneReadDict["insertSize"] = readInsertSize
                        oneReadDict["sequenceIndex"] = sequenceIndex
                        oneReadDict["sequence"] = readSequence
                        oneReadDict["start"] = readStart
                        oneReadDict["mateStart"] = readMateStart
                        oneReadDict["alignedLength"] = alignedLength
                        # add it to the dictionary of reads
                        readsDict[readName].append(oneReadDict)
                    # base at this coordinate for this read doesn't align
                    else:
                        numReadsSkipped += 1      
                # read low mapping quality
                else:
                    numReadsSkipped += 1
            # read not properly paired
            else:
                numReadsSkipped += 1
    
    if (anIsDebug):
        logging.debug("numReadsSkipped(unmap, QC, dupe, pairing, MQ, non-align)=%s", numReadsSkipped)
        logging.debug("samtools view aligned bases: %s", samtoolsViewBases)
        logging.debug("samtools view aligned quals: %s", samtoolsViewQuals)
    
    return readsDict


def find_non_overlapping_reads(aReadsDict, aMinBaseQual, aGenomicCoordinate, aTranscriptCoordinate, anIsDebug):
    '''
    ' This function loops through the reads with the same name.  Due to the inclusion of secondary
    ' alignments for RNA-Seq, there could be more than 2 reads that are paired.  If the read pairs overlap and
    ' the bases agree, then we keep the read with the highest base quality.  If the read pairs overlap and
    ' the bases disagree, then we only keep the read with the highest base quality if the other reads have a 
    ' low quality.  If the read pairs don't overlap, then keep them. 
    '
    ' aReadsDict:               The dictionary of reads grouped by their read name
    ' aMinBaseQual:             A minimum base quality score for the base
    ' aGenomicCoordinate:       The genomic coordinate
    ' aTranscriptCoordinate:    The transcript coordinate
    ' anIsDebug:                A flag for outputting debug messages to STDERR
    '''
    # loop through read pairs
    #     if read pairs overlap:
    #         if the bases agree:
    #             keep the one with the highest base quality
    #         else bases don't agree:
    #            keep the read with the high quality only if the other reads have a low quality
    #     else:
    #        keep all non-overlaps
     
    nonOverlappingReadsList = []
    for (readName, readPairsList) in aReadsDict.iteritems():
        if (anIsDebug):
            logging.debug("readList=%s", readPairsList)
        
        # if there is only one read, then no overlapping reads exist
        if (len(readPairsList) == 1):
            if (anIsDebug):
                logging.debug("only one read %s", readName)
            nonOverlappingReadsList.append(readPairsList[0])
        # check for overlapping reads
        else:
            maxBaseQual = -1
            maxBaseQualIndex = -1
            nextMaxBaseQual = -1
            basesSet = set()
            
            # for each read
            for index in range(0, len(readPairsList)):
                readStart = readPairsList[index]["start"]
                #readAlignedLength = readPairsList[index]["alignedLength"]
                readLength = len(readPairsList[index]["sequence"])
                readSequenceIndex = readPairsList[index]["sequenceIndex"]
                readMateStart = readPairsList[index]["mateStart"]
                # abs(p->b->core.pos + p->qpos - p->b->core.mpos) < p->b->core.l_qseq)
                # samtools view pos = reference_start = 1-based leftmost coordinate
                # pysam pos = reference_start = 0-based leftmost coordinate
                # samtools view pnext (1-based) = pysam mpos (0-based) = next_reference_start = the position of the mate/next read.
                # sequenceIndex pysam qpos = query_position = position of the read base at the pileup site, 0-based. None if is_del or is_refskip is set.
                # len(sequence) or pysam qlen = query_length length of the query sequence
                
                if (anIsDebug):
                    if (aTranscriptCoordinate):
                        logging.debug("aTranscriptCoordinate=%s, readStart(pos)=%s, readSequenceIndex(qpos)=%s, readMateStart(mpos)=%s, readLength(l_qseq)=%s", aTranscriptCoordinate, readStart, readSequenceIndex, readMateStart, readLength)
                    else:
                        logging.debug("aGenomicCoordinate=%s, readStart(pos)=%s, readSequenceIndex(qpos)=%s, readMateStart(mpos)=%s, readLength(l_qseq)=%s", aGenomicCoordinate, readStart, readSequenceIndex, readMateStart, readLength)
                    logging.debug("%s (readStart + readSequenceIndex): abs(readStart + readSequenceIndex - readMateStart) %s <? (readLength) %s", (readStart + readSequenceIndex), abs(readStart + readSequenceIndex - readMateStart), readLength)
                    
                # if the reads overlap
                # if abs(p->b->core.pos + p->qpos - p->b->core.mpos) < p->b->core.l_qseq)           
                # if (abs(readStart + readSequenceIndex - readMateStart) < readAlignedLength):
                if (abs(readStart + readSequenceIndex - readMateStart) < readLength):            
                    
                    if (anIsDebug):
                        logging.debug("reads overlap readStart=%s, sequenceIndex=%s, readMateStart=%s, abs=%s, readLength=%s", str(readStart), str(readSequenceIndex), str(readMateStart), str(abs(readStart + readSequenceIndex - readMateStart)), str(readLength))
                    base = readPairsList[index]["base"]
                    qual = ord(readPairsList[index]["baseQual"])-33
                    basesSet.add(base)
                    if (qual > maxBaseQual):
                        # if the current max is higher than the next max,
                        # then propagate the current max down
                        if (maxBaseQual > nextMaxBaseQual):
                            nextMaxBaseQual= maxBaseQual
                        maxBaseQual = qual
                        maxBaseQualIndex = index
                    elif (qual > nextMaxBaseQual):
                        nextMaxBaseQual = qual
                    if (anIsDebug):
                        logging.debug("base=%s, qual=%s, maxBaseQual=%s, maxBaseQualIndex=%s, basesSet=%s", base, str(qual), str(maxBaseQual), str(maxBaseQualIndex), basesSet)    
                # if they don't overlap, then just add the reads
                else:
                    if (anIsDebug):
                        logging.debug("reads don't overlap readStart=%s, sequenceIndex=%s, readMateStart=%s, abs=%s, readLength=%s", str(readStart), str(readSequenceIndex), str(readMateStart), str(abs(readStart + readSequenceIndex - readMateStart)), str(readLength))
                    nonOverlappingReadsList.append(readPairsList[index])
                    
            # if all the bases agree, then keep the one with the highest base quality and avoid double counting
            if (len(basesSet) == 1):
                if (anIsDebug):
                    logging.debug("all bases in overlapping reads agree, keeping=%s, %s %s", str(maxBaseQualIndex), readName, readPairsList)
                nonOverlappingReadsList.append(readPairsList[maxBaseQualIndex])
            # if the bases disagree, then it's likely a sequencing error. only keep the read with a high qual if the other reads have a low qual
            elif(maxBaseQual >= aMinBaseQual and nextMaxBaseQual < aMinBaseQual):
                if (anIsDebug):
                    logging.debug("not all bases in overlapping reads agree, but keeping a high qual one, maxBaseQual=%s, nextMaxBaseQual=%s, minBaseQual=%s, keeping=%s, %s %s", str(maxBaseQual), str(nextMaxBaseQual), str(aMinBaseQual), str(maxBaseQualIndex), readName, readPairsList)
                nonOverlappingReadsList.append(readPairsList[maxBaseQualIndex])
            elif anIsDebug:
                logging.debug("bases in overlapping reads don't agree and no obvious read to select based on qual scores, maxBaseQual=%s, nextMaxBaseQual=%s, minBaseQual=%s, %s %s", str(maxBaseQual), str(nextMaxBaseQual), str(aMinBaseQual), readName, readPairsList)
    
    return nonOverlappingReadsList
   
   
def write_to_blat_file(aBlatFileHandler, aChr, aGenomicCoordinate, aTranscriptNameTag, aTranscriptCoordinateTag, aTranscriptStrandTag, aParamsDict, anInfoDict, aPrefix, anAltOnlyFlag, anRnaIncludeSecondaryAlignmentsFlag, anIsDebug):
    '''
    ' This function gets all of the reads at a specific coordinate and creates a BLAT input file.
    '
    ' aBlatFileHandler:                        A file handler where all the BLAT query data is written
    ' aChr:                                    The chromosome that we are selecting from
    ' aGenomicCoordinate:                      The genomic stop coordinate of the selection
    ' aTranscriptNameTag:                      The transcript name
    ' aTranscriptCoordinateTag:                The transcript coordinate
    ' aTranscriptStrandTag:                    The transcript strand
    ' aParamsDict:                             The parameters from the header
    ' anInfoDict:                              A dictionary of the INFO column
    ' aPrefix:                                 The prefix for the parameters dictionary
    ' anAltOnlyFlag:                           If all reads should be processed or only those that have the alternate allele
    ' anRnaIncludeSecondaryAlignmentsFlag:     A flag to include the RNA secondary alignments
    ' anIsDebug:                               A flag for outputting debug messages to STDERR
    '''
    
    # get the info for executing the samtools command
    bamFile = aParamsDict[aPrefix + "Filename"]
    minMapQual = aParamsDict[aPrefix + "MinMappingQuality"]
    minBaseQual = aParamsDict[aPrefix + "MinBaseQuality"]
    useChrPrefixString = aParamsDict[aPrefix + "UseChrPrefix"]
    
    if (useChrPrefixString == "True"):
        useChrPrefix = True
    else:
        useChrPrefix = False

    if (not os.path.isfile(bamFile)):
        logging.critical("The BAM file specified in the VCF header does not exist: %s", bamFile)
        sys.exit(1)

    # if we are processing the RNA and
    # the transcript name and coordinate should be used instead of the genomic chrom and coordinate
    transcriptName = None
    transcriptCoordinate = None
    transcriptStrand = None
    if ((aPrefix == "rnaNormal" or aPrefix == "rnaTumor") and (aTranscriptNameTag != None)):
        transcriptName = anInfoDict[aTranscriptNameTag][0]
        transcriptCoordinate = int(anInfoDict[aTranscriptCoordinateTag][0])
        transcriptStrand = anInfoDict[aTranscriptStrandTag][0]
        
        # execute the samtools command
        reads = execute_samtools_cmd(bamFile, minMapQual, transcriptName, transcriptCoordinate, useChrPrefix, anRnaIncludeSecondaryAlignmentsFlag, anIsDebug)
        #reads = get_read_data(aBamFile, anIsDebug)
        
        if (anIsDebug):
            logging.debug("samtools number of reads selected from %s:%s=%s", transcriptName, transcriptCoordinate, len(reads))
    else:            
        # execute the samtools command
        reads = execute_samtools_cmd(bamFile, minMapQual, aChr, aGenomicCoordinate, useChrPrefix, False, anIsDebug)
        #reads = get_read_data(aBamFile, anIsDebug)
        
        if (anIsDebug):
            logging.debug("samtools number of reads selected from %s:%s=%s", aChr, aGenomicCoordinate, len(reads))
    
    # group all the reads by readName and add the base and baseQual
    readsDict = group_reads_by_name(reads, minMapQual, aGenomicCoordinate, transcriptName, transcriptCoordinate, transcriptStrand, anIsDebug)
    
    if (anIsDebug):     
        logging.debug("readsDictLen=%s", len(readsDict.keys()))
    
    # get all of the non-overlapping reads
    nonOverlappingReadsList = find_non_overlapping_reads(readsDict, minBaseQual, aGenomicCoordinate, transcriptCoordinate, anIsDebug)
    
    if (anIsDebug):
        logging.debug("nonOverlappingReadsListLen=%s", len(nonOverlappingReadsList))

    modChanges = anInfoDict["MC"]
    modTypes = anInfoDict["MT"]
    altSet = []
    refSet = []
    for (modChange, modType) in izip(modChanges, modTypes):
        (ref, alt) = modChange.split(">")
        altSet.append(alt)
        refSet.append(ref)
    
    numRefsFound = 0
    numAltsFound = 0
    numNonRefAltBasesFound = 0
    samtoolsViewBases = ""
    samtoolsViewQuals = ""
                        
    for readDict in nonOverlappingReadsList:
        readBase = readDict["base"]
        readBaseQual = readDict["baseQual"] 
        readSequence = readDict["sequence"]
        readSequenceIndex = readDict["sequenceIndex"]
        readName = readDict["name"]
        readMapQual = readDict["mapQual"]
        readInsertSize = readDict["insertSize"]
        readFlag = readDict["flag"]
        
        samtoolsViewBases += readBase
        samtoolsViewQuals += readBaseQual
        
        # count the number of times a base is not the ref nor the alt
        if (readBase not in refSet and readBase not in altSet):
            if (anIsDebug):
                logging.debug("base=%s is not a ref=%s, nor an alt=%s", readBase, refSet, altSet)
            numNonRefAltBasesFound += 1
            
        # count the refs
        if (readBase in refSet):
            numRefsFound += 1
        # if only the reads with the alts should be processed and the base is an alt, or
        # if all reads should be processed
        elif ((anAltOnlyFlag and readBase in altSet) or (not anAltOnlyFlag)):
            baseQualityConverted = ord(readBaseQual)-33
            if (anIsDebug):
                logging.debug("baseQualOrg=%s, baseQualityOrgOrd=%s, baseQualConverted=%s, minBaseQual=%s, greaterThan? %s", readBaseQual, str(ord(readBaseQual)), str(baseQualityConverted), str(minBaseQual), str(baseQualityConverted >= int(minBaseQual)))
            
            if (baseQualityConverted >= int(minBaseQual)):
                numAltsFound += 1
                if (anIsDebug):
                    if (transcriptName != None and transcriptCoordinate != None):
                        logging.debug("found an alt for %s:%s, base=%s, baseQual=%s", transcriptName, str(transcriptCoordinate), readBase, str(baseQualityConverted))
                    else:
                        logging.debug("found an alt for %s:%s, base=%s, baseQual=%s", aChr, str(aGenomicCoordinate), readBase, str(baseQualityConverted))
                        
                # position in read
                readLength = len(readSequence)
                if (readSequenceIndex/float(readLength) <= 0.33):
                    position = "start"
                elif (readSequenceIndex/float(readLength) <= 0.66):
                    position = "middle"
                else:
                    position = "end"
            
                # we want to use the genomic chr and stop coordinate in the output instead of the transcript coordinate
                outputList = [aPrefix, aChr, str(aGenomicCoordinate), readName.replace("_", ""), readBase, str(baseQualityConverted), str(readMapQual), position, str(readFlag), str(readInsertSize), str(readLength)]
                if (aBlatFileHandler != None):
                    aBlatFileHandler.write("> " + "_".join(outputList) + "\n")
                    aBlatFileHandler.write(readSequence + "\n")
                else:
                    print >> sys.stdout, "> " + "_".join(outputList)
                    print >> sys.stdout, readSequence
                           
    if (anIsDebug): 
        if (transcriptCoordinate): 
            logging.debug("For transcriptCoordinate=%s, numNonRefAltBasesFound=%s, numAltsFound=%s, numRefsFound=%s, numTotal=%s", transcriptCoordinate, str(numNonRefAltBasesFound), str(numAltsFound), str(numRefsFound), str(numRefsFound + numAltsFound + numNonRefAltBasesFound))
        else:
            logging.debug("For stopCoordinate=%s, numNonRefAltBasesFound=%s, numAltsFound=%s, numRefsFound=%s, numTotal=%s", aGenomicCoordinate, str(numNonRefAltBasesFound), str(numAltsFound), str(numRefsFound), str(numRefsFound + numAltsFound + numNonRefAltBasesFound))
        logging.debug("samtools view bases=%s", samtoolsViewBases)
        logging.debug("samtools view quals=%s", samtoolsViewQuals)
                                    
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
    i_cmdLineParser.add_option("-b", "--allReadBases", action="store_false", default=True, dest="altBasesOnly", help="by default only the reads with the alternate base are processed, include this argument if all of the reads should be processed")
    
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, STDOUT by default")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    
    i_cmdLineParser.add_option("", "--transcriptNameTag", dest="transcriptNameTag", help="the INFO key where the original transcript name can be found")
    i_cmdLineParser.add_option("", "--transcriptCoordinateTag", dest="transcriptCoordinateTag", help="the INFO key where the original transcript coordinate can be found")
    i_cmdLineParser.add_option("", "--transcriptStrandTag", dest="transcriptStrandTag", help="the INFO key where the original transcript strand can be found")
    i_cmdLineParser.add_option("", "--rnaIncludeSecondaryAlignments", action="store_true", default=False, dest="rnaIncludeSecondaryAlignments", help="if you align the RNA to transcript isoforms, then you may want to include RNA secondary alignments in the samtools mpileups")
    
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
    i_rnaIncludeSecondaryAlignments = i_cmdLineOptions.rnaIncludeSecondaryAlignments
    
    i_blatDnaNormalReads = i_cmdLineOptions.blatDnaNormalReads
    i_blatDnaTumorReads = i_cmdLineOptions.blatDnaTumorReads
    i_blatRnaNormalReads = i_cmdLineOptions.blatRnaNormalReads
    i_blatRnaTumorReads = i_cmdLineOptions.blatRnaTumorReads
    
    # try to get any optional parameters with no defaults    
    i_readFilenameList = [i_vcfFilename, i_headerFilename]
    i_writeFilenameList = []
    
    i_logFilename = None
    i_outputFilename = None
    i_transcriptNameTag = None
    i_transcriptCoordinateTag = None
    i_transcriptStrandTag = None
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
        i_writeFilenameList += [i_logFilename]
    if (i_cmdLineOptions.outputFilename != None):
        i_outputFilename = str(i_cmdLineOptions.outputFilename)
        i_writeFilenameList += [i_outputFilename]
    if (i_cmdLineOptions.transcriptNameTag != None):
        i_transcriptNameTag = i_cmdLineOptions.transcriptNameTag
    if (i_cmdLineOptions.transcriptCoordinateTag != None):
        i_transcriptCoordinateTag = i_cmdLineOptions.transcriptCoordinateTag
    if (i_cmdLineOptions.transcriptStrandTag != None):
        i_transcriptStrandTag = i_cmdLineOptions.transcriptStrandTag
           
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
        
        logging.debug("transcriptNameTag %s", i_transcriptNameTag)
        logging.debug("transcriptCoordinateTag %s", i_transcriptCoordinateTag)
        logging.debug("transcriptStrandTag %s", i_transcriptStrandTag)
        logging.debug("rnaIncludeSecondaryAlignments=%s" % i_rnaIncludeSecondaryAlignments)
        
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
    i_vcfGenerator  = get_vcf_data(i_vcfFilename, i_headerFilename, i_passedVCFCallsOnlyFlag, i_debug)    
   
    # for each VCF call that should be investigated   
    for (vcfChr, vcfStopCoordinate, vcfId, vcfRef, vcfAlt, vcfScore, vcfFilterSet, vcfInfoDict, restOfLine, vcfParamsDict) in i_vcfGenerator:
        if (i_debug):
            logging.debug("VCF Data: %s %s %s %s %s %s %s %s %s", vcfChr, str(vcfStopCoordinate), vcfId, vcfRef, vcfAlt, vcfScore, str(vcfFilterSet), str(vcfInfoDict), restOfLine) 
        
        modTypes = vcfInfoDict["MT"]
        for modType in modTypes:
            # get the reads contributing to a call and put them in a blat query file
            if (i_blatDnaNormalReads):
                write_to_blat_file(i_outputFileHandler, vcfChr, vcfStopCoordinate, i_transcriptNameTag, i_transcriptCoordinateTag, i_transcriptStrandTag, vcfParamsDict, vcfInfoDict, "dnaNormal", i_altBasesOnlyFlag, False, i_debug)                      
                
            if (modType == "NOR_EDIT" and i_blatRnaNormalReads):
                write_to_blat_file(i_outputFileHandler, vcfChr, vcfStopCoordinate, i_transcriptNameTag, i_transcriptCoordinateTag, i_transcriptStrandTag, vcfParamsDict, vcfInfoDict, "rnaNormal", i_altBasesOnlyFlag, i_rnaIncludeSecondaryAlignments, i_debug)
            
            if (i_blatDnaTumorReads):
                write_to_blat_file(i_outputFileHandler, vcfChr, vcfStopCoordinate, i_transcriptNameTag, i_transcriptCoordinateTag, i_transcriptStrandTag, vcfParamsDict, vcfInfoDict, "dnaTumor", i_altBasesOnlyFlag, False, i_debug)    
                
            if ((modType == "SOM" or modType == "TUM_EDIT") and i_blatRnaTumorReads):
                write_to_blat_file(i_outputFileHandler, vcfChr, vcfStopCoordinate, i_transcriptNameTag, i_transcriptCoordinateTag, i_transcriptStrandTag, vcfParamsDict, vcfInfoDict, "rnaTumor", i_altBasesOnlyFlag, i_rnaIncludeSecondaryAlignments, i_debug)
            
    stopTime = time.time()       
    logging.info("createBlatFile.py Id %s: Total time=%s hrs, %s mins, %s secs", i_id, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime))         
        
    # close the files 
    if (i_outputFilename != None):
        i_outputFileHandler.close()
        
    return
 

main()    
sys.exit(0)
