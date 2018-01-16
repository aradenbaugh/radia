#!/usr/bin/env python

import sys                          # system module
from optparse import OptionParser   # used for parsing command line arguments
import radiaUtil                    # utility functions for rna editing
import logging
import time
from pybed import pybed
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


def get_vcf_data(anInputFileHandler, anOutputFileHandler, aHeaderLine, anIsDebug):
    '''
    ' The .vcf files must have at least 10 fields:  chromosome, coordinate, id
    ' references, alts, quality score, filters, infos, format, and summary info
    ' for at least one .bam file.
    '
    ' Note: vcf files are 1-based, so create a fake startCoordinate to compare to 0-based beds
    '
    ' anInputFileHandler: The input stream for the file
    ' anOutputFileHandler: The output stream for the file
    ' aHeaderLine: The header line that should be added to the output file
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''

    hasAddedHeader = False
         
    for line in anInputFileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug):
            logging.debug("VCF Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        # if we find the INFO or FILTER section, then add the filters from here
        elif (((aHeaderLine != None) and (not hasAddedHeader)) and
              ((aHeaderLine.startswith("##INFO") and line.startswith("##INFO")) or
              (aHeaderLine.startswith("##FILTER") and line.startswith("##FILTER")))):
            hasAddedHeader = True
            if (anOutputFileHandler != None):        
                anOutputFileHandler.write(aHeaderLine + "\n")
                anOutputFileHandler.write(line + "\n")
            else:
                print >> sys.stdout, aHeaderLine
                print >> sys.stdout, line
                
        # these lines are from previous scripts in the pipeline, so output them    
        elif (line.startswith("#")):
            if (anOutputFileHandler != None):
                anOutputFileHandler.write(line + "\n")
            else:
                print >> sys.stdout, line
        # now we are to the data
        else:    
            
            # split the line on the tab
            splitLine = line.split("\t")
            
            # get the fields to yield
            #columnHeaders = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            chrom = splitLine[0]
            stopCoordinate = int(splitLine[1])
            # since vcf files are 1-based and bed files are 0-based,
            # make a fake starting coordinate so that we can compare 
            # on the 0-based level
            startCoordinate = (stopCoordinate-1)
            ids = splitLine[2]
            ref = splitLine[3]
            alt = splitLine[4]
            qual = splitLine[5]
            filters = splitLine[6]
            info = splitLine[7]
            restLine = splitLine[8:len(splitLine)]
            
            yield chrom, startCoordinate, stopCoordinate, ids, ref, alt, qual, filters, info, restLine, line
        
    return


def add_filter(aVCFFilter, aVCFInfo, aFilterName, aFilterField, anIncludeCount, aCount):
    '''
    ' Add the filter name to the filter or info column.
    '
    ' aVCFFilter: The filter column from the VCF
    ' aVCFInfo:  The info column from the VCF
    ' aFilterName: The name of the filter or info tag that should be added
    ' aFilterField: The column that should be altered (filter or info)
    ' anIncludeCount: A flag for whether the count should be included or not
    ' aCount: The number of overlaps (e.g. self-chain hits)
    '''
    # if we should add the filter to the FILTER field
    if (aFilterField == "FILTER"):
        # if the call has passed so far, then ignore the pass and return the filter
        if (aVCFFilter == "PASS"):
            aVCFFilter = aFilterName
        # if there are already filters for this call, then add this one
        else:
            vcfFilterList = aVCFFilter.split(";")
            vcfFilterList.append(aFilterName)
            aVCFFilter = ";".join(vcfFilterList)
    # if we should add the filter to the INFO field
    else:
        # if the info field is empty so far, then return the filter name
        if (aVCFInfo == "."):
            aVCFInfo = aFilterName
        # if an id already exists, then add this one
        else:
            vcfInfoList = aVCFInfo.split(";")
            if (anIncludeCount):
                filterName = aFilterName + "=" + str(aCount)
                vcfInfoList.append(filterName)
            else:
                vcfInfoList.append(aFilterName)
            aVCFInfo = ";".join(sorted(vcfInfoList))
            
    return (aVCFFilter, aVCFInfo)    
    
    
def add_id(aVCFId, anIdName):
    '''
    ' Add the id (e.g. dbSNP rs Ids) to the id column
    '
    ' aVCFId: The id column from the VCF
    ' anIdName: The id that should be added to the id column
    '''
    # if the id is empty so far, then return the id name
    if (aVCFId == "."):
        return anIdName
    # if an id already exists, then add this one
    else:
        vcfIdList = aVCFId.split(";")
        vcfIdList.append(anIdName)        
        return ";".join(vcfIdList)
        
        
def filter_events(aTCGAId, aChrom, aBedFilename, aVCFFilename, anOutputFilename, aFilterName, aFilterField, anIncludeOverlapInfo, anIncludeFilterName, anIncludeIdName, anIncludeCount, aFilterHeaderLine, aBinSize, anIsDebug):
    '''
    ' The function reads from a .bed file and a .vcf file line by line and looks for variants that should be
    ' filtered out.  The .bed file specifies coordinates for areas where variants should either be included
    ' or excluded.  For example, a .bed file specifying transcription or exon start and stop coordinates can 
    ' be provided along with the --includeOverlaps flag to indicate that the variants in these regions 
    ' should be kept, and all the others should be filtered out.  Conversely, a bed file specifying areas of 
    ' the genome that are accessible (as defined by the 1000 Genomes project) can be given without the
    ' --includeOverlaps flat to indicate that the variants outside of the accessible genome should be filtered 
    ' out, and all the other ones should be kept. 
    '
    ' aTCGAId: The TCGA Id for this sample
    ' aChrom: The chromosome being filtered
    ' aBedFilename: A .bed file with at least 3 columns specifying the chrom, start, and stop coordinates and possibly a 4th column with an id
    ' aVCFFilename: A .vcf file with variants that will be either included or excluded
    ' anOutputFilename: An output .vcf file where the filtered variants should be output
    ' anIncludeOverlapsFlag: A flag specifying whether the variants should be included or excluded when they overlap
    ' anIncludeFilterName: A flag specifying whether the filtering name should be included in the output or not
    ' anIncludeIdName: A flag specifying whether the id name should be included in the output or not
    ' anIncludeCount: A flag specifying whether the number of overlaps should be included in the output or not
    ' aFilterHeaderLine: A filter header line that should be added to the VCF header describing this filter
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    # initialize pybed with the filtering file
    filterPybed = pybed(binsize=aBinSize)
    filterPybed.loadfromfile(aBedFilename)
    
    # get the vcf file
    i_vcfFileHandler = get_read_fileHandler(aVCFFilename)
    
    # get the output file
    i_outputFileHandler = None
    if (anOutputFilename != None):
        i_outputFileHandler = get_write_fileHandler(anOutputFilename)  
    
    # create the generator for the vcf file
    vcfGenerator = get_vcf_data(i_vcfFileHandler, i_outputFileHandler, aFilterHeaderLine, anIsDebug)
    
    # initialize some variables
    overlappingEvents = 0
    nonOverlappingEvents = 0
    totalEvents = 0
    startTime = time.time()
    
    # for each vcf line
    for (vcf_chr, vcf_startCoordinate, vcf_stopCoordinate, vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_restLine, vcf_line) in (vcfGenerator):
    
        totalEvents += 1
        
        if (anIsDebug):
            logging.debug("VCF: %s", vcf_line)
        
        # check if this vcf coordinate overlaps with the filter coordinates    
        (isOverlapping, filter_id, count) = filterPybed.overlapswith((vcf_chr, vcf_startCoordinate, vcf_stopCoordinate), anIncludeCount)
        #print vcf_chr, vcf_startCoordinate, vcf_stopCoordinate, isOverlapping, filter_id

        # if an event overlaps with the filters
        if (isOverlapping):
            # count the overlap
            overlappingEvents += 1
            
            # if we want to add info about overlaps
            if (anIncludeOverlapInfo):
            
                # alter the filter and id name if appropriate
                if (anIncludeFilterName):
                    (vcf_filter, vcf_info) = add_filter(vcf_filter, vcf_info, aFilterName, aFilterField, anIncludeCount, count)
                        
                if (anIncludeIdName):
                    vcf_id = add_id(vcf_id, filter_id)
                
                # output the event
                outputList = (vcf_chr, str(vcf_stopCoordinate), vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info)
                if (anOutputFilename != None):
                    i_outputFileHandler.write("\t".join(outputList) + "\t" + "\t".join(vcf_restLine) + "\n")
                else:
                    print >> sys.stdout, "\t".join(outputList) + "\t" + "\t".join(vcf_restLine)
            # we don't want to add info about overlaps, just output them
            else:
                # output the event
                if (anOutputFilename != None):
                    i_outputFileHandler.write(vcf_line + "\n")
                else:
                    print >> sys.stdout, vcf_line
        # these events don't overlap with the filters
        else:
            # count the non overlap
            nonOverlappingEvents += 1
                        
            # if we don't want to add info about overlaps, then we do want to add info about non-overlaps
            if (not anIncludeOverlapInfo):
            
                # alter the filter and id name if appropriate
                if (anIncludeFilterName):
                    (vcf_filter, vcf_info) = add_filter(vcf_filter, vcf_info, aFilterName, aFilterField, anIncludeCount, count)
                        
                if (anIncludeIdName):
                    vcf_id = add_id(vcf_id, filter_id)
                
                # output the event
                outputList = (vcf_chr, str(vcf_stopCoordinate), vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info)
                if (anOutputFilename != None):
                    i_outputFileHandler.write("\t".join(outputList) + "\t" + "\t".join(vcf_restLine) + "\n")
                else:
                    print >> sys.stdout, "\t".join(outputList) + "\t" + "\t".join(vcf_restLine)
            # we do want to add info about overlaps, so just output non-overlaps
            else:
                # output the event
                if (anOutputFilename != None):
                    i_outputFileHandler.write(vcf_line + "\n")
                else:
                    print >> sys.stdout, vcf_line
                    
    stopTime = time.time()
    logging.info("Chrom %s and Id %s: Total time=%s hrs, %s mins, %s secs", aChrom, aTCGAId, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime)) 
           
    if (overlappingEvents + nonOverlappingEvents == totalEvents):
        logging.info("For chrom %s and Id %s: %s (overlapping events) + %s (non-overlapping events) = %s", aChrom, aTCGAId, overlappingEvents, nonOverlappingEvents, totalEvents) 
    else:
        logging.info("FilterByCoordinateBruteForce Warning: For chrom %s and Id %s: %s (overlapping events) + %s (non-overlapping events) = %s", aChrom, aTCGAId, overlappingEvents, nonOverlappingEvents, totalEvents)
        
    # close the files
    i_vcfFileHandler.close()   
    if (anOutputFilename != None):
        i_outputFileHandler.close() 
    return


def main():
    
    #python filterByPybed.py TCGA-AB-2995 12 ../data/test/filterBlacklist.bed ../data/test/TCGA-AB-2995.vcf blck -d FILTER --includeOverlaps --includeFilterName --log=DEBUG -f "##FILTER=<ID=blck,Description=\"Position overlaps 1000 Genomes Project blacklist\">"
     
    # create the usage statement
    usage = "usage: python %prog id chrom filterFile vcfFile filterName [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    i_cmdLineParser.add_option("-f", "--filterHeader", dest="filterHeader", metavar="FILTER_HEADER", help="the INFO or FORMAT line that should be included in the VCF header")
    i_cmdLineParser.add_option("-p", "--includeOverlaps", action="store_true", default=False, dest="includeOverlaps", help="whether the events that overlap should be considered to PASS (True) or FILTER (False), %default by default")
    i_cmdLineParser.add_option("-n", "--includeFilterName", action="store_true", default=False, dest="includeFilterName", help="whether the filter name should be included in the INFO or FILTER fields of the VCF output, %default by default")
    i_cmdLineParser.add_option("-c", "--includeFilterCount", action="store_true", default=False, dest="includeFilterCount", help="whether the number of overlaps with the filters should be included, %default by default")
    i_cmdLineParser.add_option("-d", "--filterField", default="FILTER", dest="filterField", help="the column where the filter name should be included (INFO or FILTER) field of the VCF output, %default by default")
    i_cmdLineParser.add_option("-i", "--includeIdName", action="store_true", default=False, dest="includeIdName", help="whether the name found in the filtering file should be included in the ID field of the VCF output, %default by default")
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, sys.stdout by default")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    i_cmdLineParser.add_option("-b", "--binSize", dest="binSize", default=int(10000), metavar="BIN_SIZE", help="the size of the interval between each bin, %default by default")
    
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
    i_chr = str(i_cmdLineArgs[1])
    i_filterFilename = str(i_cmdLineArgs[2])
    i_vcfFilename = str(i_cmdLineArgs[3])
    i_filterName = str(i_cmdLineArgs[4])
    
    # get the optional params with default values   
    i_includeOverlapsFlag = i_cmdLineOptions.includeOverlaps
    i_includeFilterName = i_cmdLineOptions.includeFilterName
    i_includeFilterCount = i_cmdLineOptions.includeFilterCount
    i_filterField = i_cmdLineOptions.filterField
    i_includeIdName = i_cmdLineOptions.includeIdName
    i_logLevel = i_cmdLineOptions.logLevel
    i_binSize = i_cmdLineOptions.binSize
    
    # try to get any optional parameters with no defaults    
    i_outputFilename = None
    i_logFilename = None
    i_filterHeader = None
    if (i_cmdLineOptions.outputFilename != None):
        i_outputFilename = str(i_cmdLineOptions.outputFilename)
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
    if (i_cmdLineOptions.filterHeader != None):
        i_filterHeader = str(i_cmdLineOptions.filterHeader)
        
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
        logging.debug("filterFile=%s", i_filterFilename)
        logging.debug("vcfFile=%s", i_vcfFilename)
        logging.debug("output=%s", i_outputFilename)
        logging.debug("logLevel=%s", i_logLevel)
        logging.debug("logFile=%s", i_logFilename)
        logging.debug("filterName=%s", i_filterName)
        logging.debug("filterHeader=%s", i_filterHeader)
        logging.debug("includeOverlapsFlag=%s", i_includeOverlapsFlag)
        logging.debug("includeFilterName=%s", i_includeFilterName)
        logging.debug("includeFilterCount=%s", i_includeFilterCount)
        logging.debug("filterField=%s", i_filterField)
        logging.debug("includeIdName=%s", i_includeIdName)
        logging.debug("binSize=%s", i_binSize)
    
    # check for any errors
    writeFilenameList = []
    if (i_outputFilename != None):
        writeFilenameList += [i_outputFilename]
    if (i_logFilename != None):
        writeFilenameList += [i_logFilename]
        
    readFilenameList = [i_filterFilename, i_vcfFilename]        
    if (not radiaUtil.check_for_argv_errors(None, readFilenameList, writeFilenameList)):
        sys.exit(1)           
    
    filter_events(i_id, i_chr, i_filterFilename, i_vcfFilename, i_outputFilename, i_filterName, i_filterField, i_includeOverlapsFlag, i_includeFilterName, i_includeIdName, i_includeFilterCount, i_filterHeader, i_binSize, i_debug)
       
    return

main()
sys.exit(0)
