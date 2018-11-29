#!/usr/bin/env python

import sys
import time
from optparse import OptionParser
import radiaUtil
import logging
import gzip


'''
'    RNA and DNA Integrated Analysis (RADIA):
'    Identifies RNA and DNA variants in NGS data.
'
'    Copyright (C) 2010  Amie J. Radenbaugh, Ph.D.
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


def get_vcf_data(aVCFFile, anIsDebug):
    
    headerList = list()
    chromLine = None 
    infoList = list() 
    filterList = list() 
    coordinateDict = dict()
    
    vcfFileHandler = get_read_fileHandler(aVCFFile)
    
    for line in vcfFileHandler:
        
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        if (anIsDebug):
            logging.debug("vcfLine: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # if we find the FILTER section, then record the filters
        elif (line.startswith("##FILTER")):
            filterList.append(line)
            
        # if we find the INFO section, then record the info
        elif (line.startswith("##INFO")):
            infoList.append(line)
              
        # if we find the header line section    
        elif (line.startswith("#CHROM")):
            chromLine = line
        
        # if we find the header line section    
        elif (line.startswith("#")):
            headerList.append(line)
           
        # now we are to the data
        else:
            # split the line on the tab
            splitLine = line.split("\t")
            
            # the coordinate is the second element
            #chrom = splitLine[0]
            stopCoordinate = splitLine[1]
            coordinateDict[stopCoordinate] = line + "\n"
                
    return (headerList, chromLine, infoList, filterList, coordinateDict)

               
def main():
           
    startTime = time.time()
    
    # create the usage statement
    usage = "usage: python %prog passingFile originalFile outputFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(3,10,1)
    i_argLength = len(sys.argv)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)
    
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_passingFilename = i_cmdLineArgs[0]
    i_originalFilename = i_cmdLineArgs[1]
    i_outputFilename = i_cmdLineArgs[2]
    
    # get the optional params with default values
    i_logLevel = i_cmdLineOptions.logLevel
    
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
        logging.debug("passingFile=%s", i_passingFilename)
        logging.debug("originalFile=%s", i_originalFilename)
        logging.debug("outputFilename=%s", i_outputFilename)
                    
    # check for any errors
    i_readFilenameList = [i_passingFilename, i_originalFilename]
    i_writeFilenameList = [i_outputFilename]
    i_dirList = None
    
    if (not radiaUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)
                
    # get the VCF generator
    (passHeaderList, chromLine, passInfoList, passFilterList, passCoordinateDict) = get_vcf_data(i_passingFilename, i_debug)
    (orgHeaderList, chromLine, orgInfoList, orgFilterList, orgCoordinateDict) = get_vcf_data(i_originalFilename, i_debug)    
    
    outputFileHandler = get_write_fileHandler(i_outputFilename)
    
    for headerLine in orgHeaderList:
        outputFileHandler.write(headerLine + "\n")
    
    for headerLine in orgInfoList:
        outputFileHandler.write(headerLine + "\n")
        
    for headerLine in passInfoList:
        if (headerLine not in orgInfoList):
            outputFileHandler.write(headerLine + "\n")
            
    for headerLine in orgFilterList:
        outputFileHandler.write(headerLine + "\n")
        
    for headerLine in passFilterList:
        if (headerLine not in orgFilterList):
            outputFileHandler.write(headerLine + "\n")
        
    outputFileHandler.write(chromLine + "\n")
    
    numericKeys = orgCoordinateDict.keys()
    numericKeys.sort(key=int)
    for coordinate in numericKeys:
        if (coordinate in passCoordinateDict):
            line = passCoordinateDict[coordinate]
        else:
            line = orgCoordinateDict[coordinate]
        outputFileHandler.write(line)
        
    stopTime = time.time() 
    logging.info("Total time=%s hrs, %s mins, %s secs", ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime))    
    # close the files 
    outputFileHandler.close()
    
    return
 

main()    
sys.exit(0)
