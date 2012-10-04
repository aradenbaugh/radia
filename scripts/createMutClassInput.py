#!/usr/bin/env python2.7

from optparse import OptionParser   # used for parsing command line arguments
import rnaEditingUtil               # utility functions for rna editing
import sys                          # system module
import re
import logging

'''
'   Amie Radenbaugh - 02/29/2012
'   UCSC - RNA Editing  
'   Program name: "createMutClassInput.py"
'''

def main():
    
    # python2.7 createMutClassInput.py TCGA-AB-2995 12 ../data/test/TCGA-AB-2995.vcf ../data/test/TCGA-AB-2995_mutClass_input.bed --log=DEBUG
    
    # create the usage statement
    usage = "usage: python2.7 %prog id chrom vcfFile mutClassInputBedFile"
    i_cmdLineParser = OptionParser(usage=usage)
    
    # add the optional params    
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(4,8,1)
    i_argLength = len(sys.argv)
    #print >> sys.stderr, "argLen", i_argLength
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)
        
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_id = str(i_cmdLineArgs[0])
    i_chrom = str(i_cmdLineArgs[1])
    i_vcfFilename = str(i_cmdLineArgs[2])
    i_mutClassInputFilename = str(i_cmdLineArgs[3])
    
    # get the optional params with default values
    i_logLevel = i_cmdLineOptions.logLevel
    
    # try to get any optional parameters with no defaults    
    i_readFilenameList = [i_vcfFilename]
    i_writeFilenameList = [i_mutClassInputFilename]
    i_dirList = []
    
    i_logFilename = None
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
        i_writeFilenameList += [i_logFilename]
            
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
        logging.debug("id=%s" % i_id)
        logging.debug("chrom=%s" % i_chrom)
        logging.debug("i_vcfFilename=%s" % i_vcfFilename)
        logging.debug("i_mutClassInputFilename=%s" % i_mutClassInputFilename)
        logging.debug("logLevel=%s" % i_logLevel)
        logging.debug("logFilename=%s" % i_logFilename)
                        
    # check for any errors
    if (not rnaEditingUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)

    
    # get the files
    i_vcfFileHandler = open(i_vcfFilename, "r")
    i_outputFileHandler = open(i_mutClassInputFilename, "w")
            
    # for each event in the vcf file 
    for line in i_vcfFileHandler:
        # here are some examples of .vcf lines:
        # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31.4;DP=53;FA=0.04;INDEL=0;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP      
        # GT:DP:INDEL:START:STOP:AD:AF:BQ:SB      0/0:2:0:0:0:2:1.0:35.5:0.0      0/0:1:0:0:0:1:1.0:39.0:1.0      0/1:50:0:2:0:48,2:0.96,0.04:31.65,17.5:0.75,0.5
        
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (i_debug):
            logging.debug("VCF Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # these lines are not necessary for the mut class info    
        elif (line.startswith("#")):
            continue;
        
        # since mutClass takes soooo long, only do the passing calls    
        elif ("PASS" not in line):
            continue;
        
        # now we are to the data
        else:    
            
            # split the line on the tab
            splitLine = re.split("\t", line)
    
            # the coordinate is the second element
            event_chr = splitLine[0]
            event_stopCoordinate = int(splitLine[1])
            event_altList = re.split(",", splitLine[4])
            
            for alt in event_altList:
                # create the output list                        
                mutClassOutputList = [("chr" + event_chr), str(event_stopCoordinate-1), str(event_stopCoordinate-1), alt]
                
                # output the list        
                i_outputFileHandler.write("\t".join(mutClassOutputList) + "\n")
    
    i_vcfFileHandler.close()
    i_outputFileHandler.close()
            
    return

main()
sys.exit(0)