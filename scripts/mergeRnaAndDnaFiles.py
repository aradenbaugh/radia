#!/usr/bin/env python

import sys
import time
from optparse import OptionParser
import radiaUtil
import os
import logging
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
    

def get_vcf_data(aDnaFile, anRnaOverlapsFile, anRnaNonOverlapsFile, aDnaHeaderOnlyFlag, anIsDebug):
    
    # open the header file
    dnaFileHandler = get_read_fileHandler(aDnaFile)
    rnaOverlapsFileHandler = get_read_fileHandler(anRnaOverlapsFile)
    if (os.path.isfile(anRnaNonOverlapsFile)):
        rnaNonOverlapsFileHandler = get_read_fileHandler(anRnaNonOverlapsFile)
    
    headerList = list()
    coordinateDict = dict()
    
    for line in dnaFileHandler:
          
        # here are some examples of .vcf lines:
        # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31.4;DP=53;FA=0.04;INDEL=0;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP      
        # GT:DP:INDEL:START:STOP:AD:AF:BQ:SB      0/0:2:0:0:0:2:1.0:35.5:0.0      0/0:1:0:0:0:1:1.0:39.0:1.0      0/1:50:0:2:0:48,2:0.96,0.04:31.65,17.5:0.75,0.5
        
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug):
            print >> sys.stderr, "VCF Line: ", line    
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # if we find the vcfGenerator line, then create the dict of params
        elif (line.startswith("#")): 
            
            # keep all the header lines
            headerList.append(line+"\n")
                
        # now we are to the data
        # if we only want the header, then break     
        elif (aDnaHeaderOnlyFlag):
            break
        # if we want the DNA data then process it
        else:
            # split the line on the tab
            splitLine = line.split("\t")
    
            # the coordinate is the second element
            #chr = splitLine[0]
            stopCoordinate = splitLine[1]
            
            coordinateDict[stopCoordinate] = line + "\n"
            
                
    for line in rnaOverlapsFileHandler:
          
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug):
            print >> sys.stderr, "VCF Line: ", line    
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # if we find the vcfGenerator line, then create the dict of params
        elif (line.startswith("#")): 
            continue;
        
        # now we are to the data    
        else:
            # split the line on the tab
            splitLine = line.split("\t")
    
            # the coordinate is the second element
            #chr = splitLine[0]
            stopCoordinate = splitLine[1]
            #chrCoordinate = chr + "_" + stopCoordinate
            
            # if this was found by both the RNA and DNA, then adjust the origin
            if (stopCoordinate in coordinateDict):
                dnaLine = coordinateDict[stopCoordinate]
                dnaLine = dnaLine.replace("ORIGIN=DNA", "ORIGIN=DNA,RNA")
                coordinateDict[stopCoordinate] = dnaLine
            else:
                coordinateDict[stopCoordinate] = line + "\n"
        
    if (os.path.isfile(anRnaNonOverlapsFile)):
        for line in rnaNonOverlapsFileHandler:
          
            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")
        
            if (anIsDebug):
                print >> sys.stderr, "VCF Line: ", line    
            
            # if it is an empty line, then just continue
            if (line.isspace()):
                continue;
        
            # if we find the vcfGenerator line, then create the dict of params
            elif (line.startswith("#")): 
                continue;
        
            # now we are to the data    
            else:
                # split the line on the tab
                splitLine = line.split("\t")
    
                # the coordinate is the second element
                #chr = splitLine[0]
                stopCoordinate = splitLine[1]
                #chrCoordinate = chr + "_" + stopCoordinate
            
                # if this was found only by the RNA, then add it
                coordinateDict[stopCoordinate] = line + "\n"
    
    dnaFileHandler.close()
    rnaOverlapsFileHandler.close()
    if (os.path.isfile(anRnaNonOverlapsFile)):
        rnaNonOverlapsFileHandler.close()
    
    return (headerList, coordinateDict)

               
def main():
    
    #python mergeRnaAndDnaFiles.py TCGA-AB-2995 5 ../data/test/TCGA-AB-2995_dnaFile.vcf ../data/test/TCGA-AB-2995_rnaFile.vcf ../data/test/TCGA-AB-2995_rnaFile.vcf ../data/test/
       
    startTime = time.time()
    
    # create the usage statement
    usage = "usage: python %prog id chrom dnaFile rnaOverlapsFile rnaNonOverlapsFile outputFile [Options]"
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
    i_rnaOverlapsFilename = i_cmdLineArgs[3]
    i_rnaNonOverlapsFilename = i_cmdLineArgs[4]
    i_outputFilename = i_cmdLineArgs[5]
    
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
    i_debug = (i_numericLogLevel < logging.WARNING)
    
    # output some debug info
    if (i_debug):
        logging.debug("id=%s", i_id)
        logging.debug("chrom=%s", i_chrom)
        logging.debug("dnaFilename=%s", i_dnaFilename)
        logging.debug("rnaOverlapsFilename=%s", i_rnaOverlapsFilename)
        logging.debug("rnaNonOverlapsFilename=%s", i_rnaNonOverlapsFilename)
        logging.debug("outputFilename=%s", i_outputFilename)
                    
    # check for any errors
    i_readFilenameList = [i_dnaFilename, i_rnaOverlapsFilename]
    i_writeFilenameList = [i_outputFilename]
    i_dirList = None
    
    if (not radiaUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)
                
    # get the VCF generator
    (headerList, coordinateDict) = get_vcf_data(i_dnaFilename, i_rnaOverlapsFilename, i_rnaNonOverlapsFilename, i_dnaHeaderOnly, i_debug)    
    
    outputFileHandler = get_write_fileHandler(i_outputFilename)
    
    for headerLine in headerList:
        outputFileHandler.write(headerLine)

    numericKeys = coordinateDict.keys()
    numericKeys.sort(key=int)
    for coordinate in numericKeys:
        outputFileHandler.write(coordinateDict[coordinate])
            
    stopTime = time.time() 
    logging.info("Total time for Id %s: Total time=%s hrs, %s mins, %s secs", i_id, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime))    
    # close the files 
    outputFileHandler.close()
    
    return
 

main()    
sys.exit(0)
