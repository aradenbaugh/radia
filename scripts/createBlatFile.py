#!/usr/bin/env python
__requires__=['pysam>=0.8.1']
import pkg_resources
import pysam
import sys
import time
from optparse import OptionParser
import radiaUtil
import collections
import logging
from itertools import izip
import os
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


def group_reads_by_name(aChromList, aPosList, aTranscriptStrandList, aBamFile, aFastaFile, aBamOrigin, anRnaIncludeSecondaryAlignmentsFlag, aMaxDepth, anIsDebug):
        
        readsDict = collections.defaultdict(list)
        # loop through all of the transcripts
        for (chrom, pos, strand) in izip(aChromList, list(map(int, aPosList)), aTranscriptStrandList):
            
            # vcfs are 1-based and pysam requires a 0-based coordinate
            pos = pos - 1
            
            # get the reference base from the fasta
            refBase = aFastaFile.fetch(chrom, pos, pos+1).upper()
            
            if (anIsDebug):
                logging.debug("getting pileups for chrom=%s, pos=%s", chrom, pos)
            
            # get the pileups
            for pileupColumn in aBamFile.pileup(chrom, pos, pos+1, stepper="nofilter", max_depth=aMaxDepth):
                
                # move through pileup until at the correct position
                if pileupColumn.pos < pos:
                    #if (anIsDebug):
                    #    logging.debug("continue pileupColumn.pos=%s, pos=%s", pileupColumn.pos, pos)
                    continue
                if pileupColumn.pos > pos:
                    #if (anIsDebug):
                    #    logging.debug("break out pileupColumn.pos=%s, pos=%s", pileupColumn.pos, pos)
                    break
                
                totalReads = 0
                keptReads = 0
                
                # loop through the reads and create a dictionary
                for pileupRead in pileupColumn.pileups:
                    
                    totalReads += 1
                    alignedRead = pileupRead.alignment
                    
                    # skip over reads with problems
                    if (alignedRead.is_qcfail or alignedRead.is_unmapped or alignedRead.is_duplicate or pileupRead.is_del or pileupRead.is_refskip):
                        #if (anIsDebug):
                        #    logging.debug("read is unmapped, duplicate, qcfail, del, or refskip at %s:%s", chrom, pos)
                        continue;
                    # no longer needed?? indel length for the position following the current pileup site
                    if pileupRead.indel != 0:
                    #    #if (anIsDebug):
                    #    #    logging.debug("base is an indel at %s:%s", chrom, pos)
                        continue;
                    # skip over secondary mappings for RNA if the param is not set to include them
                    if (aBamOrigin == "RNA" and not anRnaIncludeSecondaryAlignmentsFlag and pileupRead.alignment.is_secondary):
                        #if (anIsDebug):
                        #    logging.debug("read is secondary alignment but flag to include secondary alignments for RNA is not set %s:%s", chrom, pos)
                        continue;
                    
                    # keep a dictionary of all reads, using the readName as the key
                    # due to the inclusion of secondary alignments for RNA-Seq, there could be more than 2 reads that are paired
                    oneReadDict = {}
                    keptReads += 1
                    oneReadDict["alignedRead"] = alignedRead
                    oneReadDict["pileupRead"] = pileupRead
                    #oneReadDict["qname"] = alignedRead.query_name                          # qname
                    oneReadDict["flag"] = alignedRead.flag                                  # flag
                    oneReadDict["name"] = alignedRead.query_name                            # qname
                    oneReadDict["start"] = alignedRead.reference_start                      # pos
                    oneReadDict["mapQual"] = alignedRead.mapping_quality                    # mapq
                    #oneReadDict["cigar"] = alignedRead.cigar                               # cigar
                    #oneReadDict["mateName"] = alignedRead.next_reference_name              # rnext
                    oneReadDict["mateStart"] = alignedRead.next_reference_start             # pnext or mpos
                    oneReadDict["insertSize"] = alignedRead.template_length                 # isize or tlen
                    oneReadDict["sequence"] = alignedRead.seq                               # seq
                    #oneReadDict["qualities"] = alignedRead.qual                            # qual
                    oneReadDict["qlen"] = alignedRead.query_length                          # qlen
                    oneReadDict["base"] = alignedRead.seq[pileupRead.query_position]
                    oneReadDict["baseQual"] = alignedRead.qual[pileupRead.query_position]
                    oneReadDict["sequenceIndex"] = pileupRead.query_position
                    oneReadDict["chrom"] = chrom
                    oneReadDict["pos"] = pos
                    oneReadDict["strand"] = strand
                    oneReadDict["refBase"] = refBase
                    
                    #if (anIsDebug):
                    #    logging.debug("name=%s, oneReadDict=%s", alignedRead.query_name, oneReadDict)
                    
                    # add it to the dictionary of reads
                    readsDict[alignedRead.query_name].append(oneReadDict)
                    
                if (anIsDebug):
                    logging.debug("group_reads_by_name(): %s:%s added %s out of %s reads to initial reads dict, readsDictLen=%s", chrom, pos, keptReads, totalReads, len(readsDict.keys()))
        
        return readsDict


def find_non_overlapping_reads(aReadsDict, aMinBaseQual, anIsDebug):
    '''
    ' This function loops through the reads with the same name.  Due to the inclusion of secondary
    ' alignments for RNA-Seq, there could be more than 2 reads that are paired.  If the read pairs overlap and
    ' the bases agree, then we keep the read with the highest base quality.  If the read pairs overlap and
    ' the bases disagree, then we only keep the read with the highest base quality if the other reads have a 
    ' low quality.  If the read pairs don't overlap, then keep them. 
    '
    ' aReadsDict:               The dictionary of reads grouped by their read name
    ' aMinBaseQual:             A minimum base quality score for the base
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
        #if (anIsDebug):
        #    logging.debug("readName=%s, len=%s, readList=%s", readName, len(readPairsList), readPairsList)
        
        # if there is only one read, then no overlapping reads exist
        if (len(readPairsList) == 1):
            #if (anIsDebug):
            #    logging.debug("only one read %s", readName)
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
                readSequenceIndex = readPairsList[index]["sequenceIndex"]
                readMateStart = readPairsList[index]["mateStart"]
                readLength = readPairsList[index]["qlen"]
                # abs(p->b->core.pos + p->qpos - p->b->core.mpos) < p->b->core.l_qseq)
                # samtools view pos = reference_start = 1-based leftmost coordinate
                # pysam pos = reference_start = 0-based leftmost coordinate
                # samtools view pnext (1-based) = pysam mpos (0-based) = next_reference_start = the position of the mate/next read.
                # sequenceIndex pysam qpos = query_position = position of the read base at the pileup site, 0-based. None if is_del or is_refskip is set.
                # len(sequence) or pysam qlen = query_length length of the query sequence
                    
                #if (anIsDebug):
                #    logging.debug("readStart(pos)=%s, readSequenceIndex(qpos)=%s, readMateStart(mpos)=%s, readLength(l_qseq)=%s", readStart, readSequenceIndex, readMateStart, readLength)
                #    logging.debug("%s (readStart + readSequenceIndex): abs(readStart + readSequenceIndex - readMateStart) %s <? (readAlignedLength) %s", (readStart + readSequenceIndex), abs(readStart + readSequenceIndex - readMateStart), readLength)
                
                # if the reads overlap
                # if abs(p->b->core.pos + p->qpos - p->b->core.mpos) < p->b->core.l_qseq)
                if (abs(readStart + readSequenceIndex - readMateStart) < readLength):            
                    
                    #if (anIsDebug):
                    #    logging.debug("reads overlap readStart=%s, sequenceIndex=%s, readMateStart=%s, abs=%s, alignedLen=%s", str(readStart), str(readSequenceIndex), str(readMateStart), str(abs(readStart + readSequenceIndex - readMateStart)), str(readLength))
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
                    #if (anIsDebug):
                    #    logging.debug("base=%s, qual=%s, maxBaseQual=%s, maxBaseQualIndex=%s, basesSet=%s", base, str(qual), str(maxBaseQual), str(maxBaseQualIndex), basesSet)
                # if they don't overlap, then just add the reads
                else:
                    #if (anIsDebug):
                    #    logging.debug("reads don't overlap readStart=%s, sequenceIndex=%s, readMateStart=%s, abs=%s, alignedLen=%s", str(readStart), str(readSequenceIndex), str(readMateStart), str(abs(readStart + readSequenceIndex - readMateStart)), str(readLength))
                    nonOverlappingReadsList.append(readPairsList[index])
            
            # if all the bases agree, then keep the one with the highest base quality and avoid double counting
            if (len(basesSet) == 1):
                #if (anIsDebug):
                #    logging.debug("all bases in overlapping reads agree, keeping=%s, %s %s", str(maxBaseQualIndex), readName, readPairsList)
                nonOverlappingReadsList.append(readPairsList[maxBaseQualIndex])
            # if the bases disagree, then it's likely a sequencing error. only keep the read with a high qual if the other reads have a low qual
            elif(maxBaseQual >= aMinBaseQual and nextMaxBaseQual < aMinBaseQual):
                #if (anIsDebug):
                #    logging.debug("not all bases in overlapping reads agree, but keeping a high qual one, maxBaseQual=%s, nextMaxBaseQual=%s, minBaseQual=%s, keeping=%s, %s %s", str(maxBaseQual), str(nextMaxBaseQual), str(aMinBaseQual), str(maxBaseQualIndex), readName, readPairsList)
                nonOverlappingReadsList.append(readPairsList[maxBaseQualIndex])
            #elif anIsDebug:
            #    logging.debug("bases in overlapping reads don't agree and no obvious read to select based on qual scores, maxBaseQual=%s, nextMaxBaseQual=%s, minBaseQual=%s, %s %s", str(maxBaseQual), str(nextMaxBaseQual), str(aMinBaseQual), readName, readPairsList)
    
    if (anIsDebug):
        logging.debug("find_non_overlapping_reads(): nonOverlappingReadsListLen=%s", len(nonOverlappingReadsList))
    
    return nonOverlappingReadsList

   
def write_to_blat_file(aBlatFileHandler, aGenomicChr, aGenomicCoordinate, aChromList, aPosList, aTranscriptStrandList, aParamsDict, anInfoDict, aPrefix, anAltOnlyFlag, anRnaIncludeSecondaryAlignmentsFlag, aMaxDepth, anIsDebug):
    '''
    ' This function gets all of the reads at a specific coordinate and creates a BLAT input file.
    '
    ' aBlatFileHandler:                        A file handler where all the BLAT query data is written
    ' aGenomicChr:                             The genomic chromosome
    ' aGenomicCoordinate:                      The genomic stop coordinate
    ' aChromList:                              A list of chromosomes or transcript names to process
    ' aPosList:                                A list of coordinates to process
    ' aTranscriptStrandList:                   A list of transcript strands
    ' aParamsDict:                             The parameters from the header
    ' anInfoDict:                              A dictionary of the INFO column
    ' aPrefix:                                 The prefix for the parameters dictionary
    ' anAltOnlyFlag:                           If all reads should be processed or only those that have the alternate allele
    ' anRnaIncludeSecondaryAlignmentsFlag:     A flag to include the RNA secondary alignments
    ' aMaxDepth                                The max depth for the pysam pileup command
    ' anIsDebug:                               A flag for outputting debug messages to STDERR
    '''
    
    # get the info for executing the samtools command
    bamFilename = aParamsDict[aPrefix + "Filename"]
    fastaFilename = aParamsDict[aPrefix + "FastaFilename"]
    minBaseQual = aParamsDict[aPrefix + "MinBaseQuality"]
    bamOrigin = anInfoDict["ORIGIN"]
    
    if (not os.path.isfile(bamFilename)):
        logging.critical("The BAM file specified in the VCF header does not exist: %s", bamFilename)
        sys.exit(1)
    
    if (not os.path.isfile(fastaFilename)):
        logging.critical("The FASTA file specified in the VCF header does not exist: %s", fastaFilename)
        sys.exit(1)
    
    bamFile = pysam.Samfile(bamFilename, 'rb')
    fastaFile = pysam.Fastafile(fastaFilename)
       
    # group all of the reads by name
    readsDict = group_reads_by_name(aChromList, aPosList, aTranscriptStrandList, bamFile, fastaFile, bamOrigin, anRnaIncludeSecondaryAlignmentsFlag, aMaxDepth, anIsDebug)
        
    # get all of the non-overlapping reads
    nonOverlappingReadsList = find_non_overlapping_reads(readsDict, minBaseQual, anIsDebug)
    
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
            baseQualityConverted = ord(readBaseQual)-33
            if (baseQualityConverted >= int(minBaseQual)):
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
                    logging.debug("found an alt for %s:%s, base=%s, baseQual=%s", aChromList, aPosList, readBase, str(baseQualityConverted))
            
                # we want to use the genomic chr and stop coordinate in the output instead of the transcript coordinate
                outputList = [aPrefix, aGenomicChr, str(aGenomicCoordinate), readName.replace("_", ""), readBase, str(baseQualityConverted), str(readMapQual), str(readFlag), str(readInsertSize)]
                if (aBlatFileHandler != None):
                    aBlatFileHandler.write("> " + "_".join(outputList) + "\n")
                    aBlatFileHandler.write(readSequence + "\n")
                else:
                    print >> sys.stdout, "> " + "_".join(outputList)
                    print >> sys.stdout, readSequence

    if (anIsDebug): 
        logging.debug("For %s:%s, numNonRefAltBasesFound=%s, numAltsFound=%s, numRefsFound=%s, numTotal=%s", aChromList, aPosList, str(numNonRefAltBasesFound), str(numAltsFound), str(numRefsFound), str(numRefsFound + numAltsFound + numNonRefAltBasesFound))
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
    i_cmdLineParser.add_option("-d", "--maxReadDepth", type="int", default=int(8000), dest="maxReadDepth", metavar="MAX_READ_DEPTH", help="the maximum read depth to process from the samtools view command, %default by default")
    
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
    i_maxReadDepth = i_cmdLineOptions.maxReadDepth
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
        logging.debug("maxReadDepth %s", i_maxReadDepth)
        
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
                write_to_blat_file(i_outputFileHandler, 
                                   vcfChr, 
                                   vcfStopCoordinate, 
                                   [vcfChr], 
                                   [vcfStopCoordinate], 
                                   [None], 
                                   vcfParamsDict, 
                                   vcfInfoDict, 
                                   "dnaNormal", 
                                   i_altBasesOnlyFlag, 
                                   False,
                                   i_maxReadDepth, 
                                   i_debug)
                
            if (modType == "NOR_EDIT" and i_blatRnaNormalReads):
                # if we should process the transcripts
                if ((i_transcriptNameTag != None) and (i_transcriptNameTag in vcfInfoDict)):
                    write_to_blat_file(i_outputFileHandler, 
                                       vcfChr, 
                                       vcfStopCoordinate, 
                                       vcfInfoDict[i_transcriptNameTag], 
                                       vcfInfoDict[i_transcriptCoordinateTag], 
                                       vcfInfoDict[i_transcriptStrandTag], 
                                       vcfParamsDict, 
                                       vcfInfoDict, 
                                       "rnaNormal", 
                                       i_altBasesOnlyFlag, 
                                       i_rnaIncludeSecondaryAlignments,
                                       i_maxReadDepth, 
                                       i_debug)
                else:
                    write_to_blat_file(i_outputFileHandler, 
                                       vcfChr, 
                                       vcfStopCoordinate, 
                                       [vcfChr], 
                                       [vcfStopCoordinate], 
                                       [None], 
                                       vcfParamsDict, 
                                       vcfInfoDict, 
                                       "rnaNormal", 
                                       i_altBasesOnlyFlag, 
                                       i_rnaIncludeSecondaryAlignments, 
                                       i_maxReadDepth,
                                       i_debug)
            
            if (i_blatDnaTumorReads):
                write_to_blat_file(i_outputFileHandler, 
                                   vcfChr, 
                                   vcfStopCoordinate, 
                                   [vcfChr], 
                                   [vcfStopCoordinate], 
                                   [None], 
                                   vcfParamsDict, 
                                   vcfInfoDict, 
                                   "dnaTumor", 
                                   i_altBasesOnlyFlag, 
                                   False, 
                                   i_maxReadDepth,
                                   i_debug)
                
            if ((modType == "SOM" or modType == "TUM_EDIT") and i_blatRnaTumorReads):
                # if we should process the transcripts
                if ((i_transcriptNameTag != None) and (i_transcriptNameTag in vcfInfoDict)):
                    write_to_blat_file(i_outputFileHandler, 
                                       vcfChr, 
                                       vcfStopCoordinate,
                                       list(vcfInfoDict[i_transcriptNameTag]), 
                                       vcfInfoDict[i_transcriptCoordinateTag], 
                                       vcfInfoDict[i_transcriptStrandTag], 
                                       vcfParamsDict, 
                                       vcfInfoDict, 
                                       "rnaTumor", 
                                       i_altBasesOnlyFlag, 
                                       i_rnaIncludeSecondaryAlignments,
                                       i_maxReadDepth, 
                                       i_debug)
                else:
                    write_to_blat_file(i_outputFileHandler, 
                                       vcfChr, 
                                       vcfStopCoordinate, 
                                       [vcfChr], 
                                       [vcfStopCoordinate], 
                                       [None], 
                                       vcfParamsDict, 
                                       vcfInfoDict, 
                                       "rnaTumor", 
                                       i_altBasesOnlyFlag, 
                                       i_rnaIncludeSecondaryAlignments,
                                       i_maxReadDepth, 
                                       i_debug)
            
    stopTime = time.time()       
    logging.info("createBlatFile.py Id %s: Total time=%s hrs, %s mins, %s secs", i_id, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime))         
        
    # close the files 
    if (i_outputFilename != None):
        i_outputFileHandler.close()
        
    return
 

main()    
sys.exit(0)
