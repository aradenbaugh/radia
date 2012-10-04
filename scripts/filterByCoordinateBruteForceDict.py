#!/usr/bin/env python2.7

import sys                          # system module
from optparse import OptionParser   # used for parsing command line arguments
import rnaEditingUtil               # utility functions for rna editing
import logging
import time
import collections

'''
'   Amie Radenbaugh - 03/03/2011
'   UCSC - RNA Editing  
'   Program name: "filterByCoordinate.py"
'''

def get_bed_data(anInputStream, anIsDebug):
    '''
    ' The bed files must have at least 3 fields:  chromosome, 
    ' start coordinate, stop coordinate.  Any other column after 
    ' that will be only be used as extra information if appropriate.
    '
    ' anInputStream: The input stream for the file
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
     
    dbSnpDict = collections.defaultdict(list)
    for line in anInputStream:
        # if it is an empty line, then just continue
        if (line.isspace() or line.startswith("#")):
            continue;

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug):
            logging.debug("Filter Line: %s", line)
        
        # split the line on the tab
        splitLine = line.split("\t")
        
        # get the fields to yield
        chr = splitLine[0]
        startCoordinate = int(splitLine[1])
        # a bunch of the entries in dbSnp have the same start and stop coordinate
        # they will never be matched to the vcf entries, one way to deal with this
        # is to just take the startCoordinate and add one for the stop coordinate
        # bed files are 0-based, and vcf files are 1-based 
        #stopCoordinate = (startCoordinate + 1)
        stopCoordinate = int(splitLine[2])
        
        # the blacklist only has 3 columns
        # the dbSNP has the name of the SNP in the 4th column
        # the transcript and exon filters have the name of the gene and strand at the 4th and 5th columns
        if (len(splitLine) == 3):
            name = ""
        elif (len(splitLine) == 4):
            name = splitLine[3]
        elif (len(splitLine) == 5):
            name = splitLine[3]
        
        dbSnpDict[str(startCoordinate) + "_" + str(stopCoordinate)].append(name)
        
    return(dbSnpDict)


def get_vcf_data(anInputFileHandler, anOutputFileHandler, aHeaderLine, anIsDebug):
    '''
    ' The .vcf files must have at least 10 fields:  chromosome, coordinate, id
    ' references, alts, quality score, filters, infos, format, and summary info
    ' for at least one .bam file.
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
            ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 130">
            ##FILTER=<ID=bldp,Description="Position overlap 1000 Genomes Project depth blacklist">
            ##FILTER=<ID=blq,Description="Position overlaps 1000 Genomes Project mapping quality blacklist">
            ##FILTER=<ID=blck,Description="Position overlaps 1000 Genomes Project blacklist">
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
            chr = splitLine[0]
            stopCoordinate = int(splitLine[1])
            id = splitLine[2]
            ref = splitLine[3]
            alt = splitLine[4]
            qual = splitLine[5]
            filter = splitLine[6]
            info = splitLine[7]
            restLine = splitLine[8:len(splitLine)]
            
            yield chr, (stopCoordinate-1), stopCoordinate, id, ref, alt, qual, filter, info, restLine, line
        
    return


def add_filter(aVCFFilter, aVCFInfo, aFilterName, aFilterField):
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
            vcfInfoList.append(aFilterName)
            aVCFInfo = ";".join(vcfInfoList)
            
    return (aVCFFilter, aVCFInfo)    
    
    
def add_id(aVCFId, aNamesList):
    # if the id is empty so far, then return the id name
    if (aVCFId == "."):
        return ";".join(aNamesList)
    # if an id already exists, then add this one
    else:
        vcfIdList = aVCFId.split(";")
        vcfIdList += aNamesList        
        return ";".join(vcfIdList)
        
        
def filter_events(aTCGAId, aChrom, aBedFilename, aVCFFilename, anOutputFilename, aFilterName, aFilterField, anIncludeOverlapInfo, anIncludeFilterName, anIncludeIdName, aFilterHeaderLine, anIsDebug):
    '''
    ' The function reads from a bed file and a .vcf file line by line and looks for editing events that should be
    ' filtered out.  The bed file specifies coordinates for areas where editing events should either be included
    ' or excluded.  For example, a bed file specifying transcription or exon start and stop coordinates can be provided
    ' along with the anIncludeOverlapsFlag=True to indicate that the editing events in these regions should be kept, 
    ' and all the others should be filtered out.  Conversely, a bed file specifying areas of the genome that are blacklisted 
    ' can be given along with the anIncludeOverlapsFlag=False to indicate that the editing events in these regions should
    ' be filtered out, and all the other ones should be kept. 
    '
    ' aTCGAId: The TCGA Id for this sample
    ' aChrom: The chromosome being filtered
    ' aBedFilename: A .bed file with at least 3 columns specifying the chrom, start, and stop coordinates
    ' aVCFFilename: A .vcf file with RNA-editing events that will be either included or excluded
    ' anOutputFilename: A output file where the filtered RNA-editing events should be output
    ' anIncludeOverlapsFlag: A flag specifying whether the events should be included or excluded when they overlap
    ' anIncludeFilterName: A flag specifying whether the filtering name should be included in the output or not
    ' anIncludeIdName: A flag specifying whether the id name should be included in the output or not
    ' aFilterHeaderLine: A filter header line that should be added to the VCF header describing this filter
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    # get the files
    i_vcfFileHandler = open(aVCFFilename, "r")
    i_filterFileHandler = open(aBedFilename, "r")
    
    i_outputFileHandler = None
    if (anOutputFilename != None):
        i_outputFileHandler = open(anOutputFilename, "w")  
    
    # create the generators for the filter and vcf files
    (i_dbSnpDict) = get_bed_data(i_filterFileHandler, anIsDebug)
    vcfGenerator = get_vcf_data(i_vcfFileHandler, i_outputFileHandler, aFilterHeaderLine, anIsDebug)
    
    # initialize some variables
    overlappingEvents = 0
    nonOverlappingEvents = 0
    totalEvents = 0
    startTime = time.time()
    
    # for each vcf line
    for (vcf_chr, vcf_startCoordinate, vcf_stopCoordinate, vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_restLine, vcf_line) in (vcfGenerator):
    
        isOverlapping = False
        totalEvents += 1
        
        if (anIsDebug):
            logging.debug("VCF: %s", vcf_line)
            
        vcf_coordinateKey = str(vcf_startCoordinate) + "_" + str(vcf_stopCoordinate)
        if (vcf_coordinateKey in i_dbSnpDict):
            isOverlapping = True
            filterNamesList = i_dbSnpDict[vcf_coordinateKey]
            
        # if an event overlaps with the filters
        if (isOverlapping):
            # count the overlap
            overlappingEvents += 1
            
            # if we want to add info about overlaps
            if (anIncludeOverlapInfo):
            
                # alter the filter and id name if appropriate
                if (anIncludeFilterName):
                    (vcf_filter, vcf_info) = add_filter(vcf_filter, vcf_info, aFilterName, aFilterField)
                        
                if (anIncludeIdName):
                    vcf_id = add_id(vcf_id, filterNamesList)
                
                # output the event
                outputList = (vcf_chr, str(vcf_stopCoordinate), vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info)
                if (anOutputFilename != None):
                    i_outputFileHandler.write("\t".join(outputList) + "\t" + "\t".join(vcf_restLine) + "\n")
                else:
                    print >> sys.stdout, "\t".join(outputList) + "\t" + "\t".join(vcf_restLine)
            # we don't want to add info about overlaps, just output them
            else:
                # output the editing event
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
                    (vcf_filter, vcf_info) = add_filter(vcf_filter, vcf_info, aFilterName, aFilterField)
                        
                if (anIncludeIdName):
                    vcf_id = add_id(vcf_id, filterNamesList)
                
                # output the event
                outputList = (vcf_chr, str(vcf_stopCoordinate), vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info)
                if (anOutputFilename != None):
                    i_outputFileHandler.write("\t".join(outputList) + "\t" + "\t".join(vcf_restLine) + "\n")
                else:
                    print >> sys.stdout, "\t".join(outputList) + "\t" + "\t".join(vcf_restLine)
            # we do want to add info about overlaps, so just output non-overlaps
            else:
                # output the editing event
                if (anOutputFilename != None):
                    i_outputFileHandler.write(vcf_line + "\n")
                else:
                    print >> sys.stdout, vcf_line
                    
    stopTime = time.time()
    logging.warning("Chrom %s and Id %s: Total time=%s hrs, %s mins, %s secs", aChrom, aTCGAId, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime)) 
           
    if (overlappingEvents + nonOverlappingEvents == totalEvents):
        logging.warning("For chrom %s and Id %s: %s (overlapping events) + %s (non-overlapping events) = %s", aChrom, aTCGAId, overlappingEvents, nonOverlappingEvents, totalEvents) 
    else:
        logging.warning("FilterByCoordinateBruteForce Warning: For chrom %s and Id %s: %s (overlapping events) + %s (non-overlapping events) = %s", aChrom, aTCGAId, overlappingEvents, nonOverlappingEvents, totalEvents)
        
    # close the files 
    i_filterFileHandler.close()
    i_vcfFileHandler.close()   
    if (anOutputFilename != None):
        i_outputFileHandler.close() 
    return


def main():
    
    #python2.7 filterByCoordinateBruteForce.py TCGA-AB-2995 12 ../data/test/filterBlacklist.bed ../data/test/TCGA-AB-2995.vcf blck -d FILTER --includeOverlaps --includeFilterName --log=DEBUG -f "##FILTER=<ID=blck,Description=\"Position overlaps 1000 Genomes Project blacklist\">"
    #python2.7 filterByCoordinateBruteForce.py TCGA-AB-2995 12 ../data/test/filterSNP.bed ../data/test/TCGA-AB-2995.vcf DB -d INFO --includeOverlaps --includeFilterName --includeIdName --log=DEBUG -f "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 130\">"
    
    #python2.7 filterByCoordinateBruteForce.py TCGA-AB-2897 1 ../data/hg18_chr1_tx_sorted_uniq.bed ../data/TCGA-AB-2897_chr1.rad True -i True -s True
    #python2.7 filterByCoordinateBruteForce.py TCGA-XXXX 1 ../data/filter.bed ../data/filter.rad False -i True -s False
    #python2.7 filterByCoordinateBruteForce.py TCGA-AB-2995 12 ../data/blacklists/chr12_sorted.bed /inside/grotto/users/aradenba/data/hg19/aml/rnaEditing/TCGA-AB-2995_chr12.rad False -i True -s False
    
    # create the usage statement
    usage = "usage: python2.7 %prog id chrom filterFile vcfFile filterName [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    i_cmdLineParser.add_option("-f", "--filterHeader", dest="filterHeader", metavar="FILTER_HEADER", help="the INFO or FORMAT line that should be included in the VCF header")
    i_cmdLineParser.add_option("-p", "--includeOverlaps", action="store_true", default=False, dest="includeOverlaps", help="whether the events that overlap should be considered to PASS (True) or FILTER (False), %default by default")
    i_cmdLineParser.add_option("-n", "--includeFilterName", action="store_true", default=False, dest="includeFilterName", help="whether the filter name should be included in the INFO or FILTER fields of the VCF output, %default by default")
    i_cmdLineParser.add_option("-d", "--filterField", default="FILTER", dest="filterField", help="the column where the filter name should be included (INFO or FILTER) field of the VCF output, %default by default")
    i_cmdLineParser.add_option("-i", "--includeIdName", action="store_true", default=False, dest="includeIdName", help="whether the name found in the filtering file should be included in the ID field of the VCF output, %default by default")
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, sys.stdout by default")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    #i_cmdLineParser.add_option("-s", "--includeStrand", help="whether the strand found in the filtering file should be used to convert the events and added to the .rad output, default=False", dest="includeStrand", default=False)
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(6,21,1)
    i_argLength = len(sys.argv)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        #print "usage: python2.7 filterByCoordinate.py id chr filterFile vcfFile filterName --includeOverlaps --includeFilterName --includeIdName -f filterHeader -o outputFile"
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
    i_filterField = i_cmdLineOptions.filterField
    i_includeIdName = i_cmdLineOptions.includeIdName
    i_logLevel = i_cmdLineOptions.logLevel
    #i_includeStrand = i_cmdLineOptions.includeStrand
    
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
    i_debug = (i_numericLogLevel < logging.WARNING)

    # do some debugging
    if (i_debug):
        logging.debug("id=%s", i_id)
        logging.debug("chr=%s", i_chr)
        logging.debug("filterFile=%s", i_filterFilename)
        logging.debug("vcfFile=%s", i_vcfFilename)
        logging.debug("output=%s", i_outputFilename)
        logging.debug("i_logLevel=%s", i_logLevel)
        logging.debug("logFile=%s", i_logFilename)
        logging.debug("filterName=%s", i_filterName)
        logging.debug("i_filterHeader=%s", i_filterHeader)
        logging.debug("i_includeOverlapsFlag=%s", i_includeOverlapsFlag)
        logging.debug("i_includeFilterName=%s", i_includeFilterName)
        logging.debug("i_filterField=%s", i_filterField)
        logging.debug("i_includeIdName=%s", i_includeIdName)
    
    # check for any errors
    writeFilenameList = []
    if (i_outputFilename != None):
        writeFilenameList += [i_outputFilename]
    if (i_logFilename != None):
        writeFilenameList += [i_logFilename]
        
    readFilenameList = [i_filterFilename, i_vcfFilename]        
    if (not rnaEditingUtil.check_for_argv_errors(None, readFilenameList, writeFilenameList)):
        sys.exit(1)           
    
    filter_events(i_id, i_chr, i_filterFilename, i_vcfFilename, i_outputFilename, i_filterName, i_filterField, i_includeOverlapsFlag, i_includeFilterName, i_includeIdName, i_filterHeader, i_debug)
       
    return

main()
sys.exit(0)
