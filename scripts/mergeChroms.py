#!/usr/bin/env python

import sys
import time
from optparse import OptionParser
import radiaUtil
import os
import glob
import logging


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


def get_vcf_data(anId, anInputDir, anIsDebug):

    # for each file that starts with this id
        # load the first file to get the header
        # get the coordinates for all

    processedHeader = False
    headerDict = dict()
    headerDict["metadata"] = list()
    headerDict["format"] = list()
    headerDict["info"] = list()
    headerDict["filter"] = list()
    headerDict["chrom"] = list()
    coordinateDict = dict()
    coordinateDict["numbers"] = dict()
    coordinateDict["letters"] = dict()

    # if the input directory doesn't end with a forward slash,
    # then add one so that glob.glob will work
    if (not anInputDir.endswith("/")):
        anInputDir = anInputDir + "/"

    # for each vcf file
    # they might be gzipped, they might not
    for vcfFile in (glob.glob(anInputDir + anId + "_chr*.vcf*")):

        # open the file
        vcfFileHandler = radiaUtil.get_read_fileHandler(vcfFile)

        for line in vcfFileHandler:

            # if it is an empty line, then just continue
            if (line.isspace()):
                continue

            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")

            if (anIsDebug):
                logging.debug("vcfLine: %s", line)

            # if we haven't processed the header yet, then do it here
            if (not processedHeader):
                # extract the metadata
                if (line.startswith("##FORMAT")):
                    headerDict["format"].append(line)
                elif (line.startswith("##INFO")):
                    headerDict["info"].append(line)
                elif (line.startswith("##FILTER")):
                    headerDict["filter"].append(line)
                elif (line.startswith("##")):
                    headerDict["metadata"].append(line)
                elif (line.startswith("#CHROM")):
                    headerDict["chrom"].append(line)
                    # now we've processed the header
                    processedHeader = True

            if (line.startswith("#")):
                continue
            else:
                # split the line on the tab
                splitLine = line.split("\t")

                # the coordinate is the second element
                chrom = splitLine[0]

                # we want to sort everything at the end, so keep track
                # of the chroms that are numbers and letters separately
                if (is_number(chrom)):
                    if chrom not in coordinateDict["numbers"]:
                        coordinateDict["numbers"][chrom] = list()
                    coordinateDict["numbers"][chrom].append(line)
                else:
                    if chrom not in coordinateDict["letters"]:
                        coordinateDict["letters"][chrom] = list()
                    coordinateDict["letters"][chrom].append(line)

        # close the file and move onto the next one
        vcfFileHandler.close()

    return (headerDict, coordinateDict)


def is_number(aChrom):
    try:
        int(aChrom)
        return True
    except ValueError:
        return False


def main():

    # command for running this on a small test case:
    # python mergeChroms.py TCGA-BH-A18P
    # ../data/test/ ../data/test/ --log=DEBUG

    startTime = time.time()

    # create the usage statement
    usage = "usage: python %prog id inputDir outputDir [Options]"
    i_cmdLineParser = OptionParser(usage=usage)

    i_cmdLineParser.add_option(
        "-l", "--log",
        dest="logLevel", default="WARNING", metavar="LOG",
        help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), " +
             "%default by default")
    i_cmdLineParser.add_option(
        "-g", "--logFilename",
        dest="logFilename", metavar="LOG_FILE",
        help="the name of the log file, STDERR by default")
    i_cmdLineParser.add_option(
        "", "--gzip",
        dest="gzip", action="store_true", default=False,
        help="include this argument if the final VCF should be " +
             "compressed with gzip")

    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(3, 10, 1)
    i_argLength = len(sys.argv)

    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)

    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_id = i_cmdLineArgs[0]
    i_inputDir = i_cmdLineArgs[1]
    i_outputDir = i_cmdLineArgs[2]

    # get the optional params with default values
    i_logLevel = i_cmdLineOptions.logLevel
    i_gzip = i_cmdLineOptions.gzip

    i_logFilename = None
    if (i_cmdLineOptions.logFilename is not None):
        i_logFilename = str(i_cmdLineOptions.logFilename)

    # assuming loglevel is bound to the string value obtained from the
    # command line argument. Convert to upper case to allow the user to
    # specify --log=DEBUG or --log=debug
    i_numericLogLevel = getattr(logging, i_logLevel.upper(), None)
    if not isinstance(i_numericLogLevel, int):
        raise ValueError("Invalid log level: '%s' must be one of the " +
                         "following:  DEBUG, INFO, WARNING, ERROR, CRITICAL",
                         i_logLevel)

    # set up the logging
    if (i_logFilename is not None):
        logging.basicConfig(
            level=i_numericLogLevel,
            filename=i_logFilename,
            filemode='w',
            format='%(asctime)s\t%(levelname)s\t%(message)s',
            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(
            level=i_numericLogLevel,
            format='%(asctime)s\t%(levelname)s\t%(message)s',
            datefmt='%m/%d/%Y %I:%M:%S %p')

    # set the debug
    i_debug = (i_numericLogLevel == logging.DEBUG)

    # output some debug info
    if (i_debug):
        logging.debug("id=%s", i_id)
        logging.debug("inputDir=%s", i_inputDir)
        logging.debug("outputDir=%s", i_outputDir)
        logging.debug("logLevel=%s", i_logLevel)
        logging.debug("logFile=%s", i_logFilename)
        logging.debug("gzip=%s", i_gzip)

    # check for any errors
    i_readFilenameList = None
    if (i_logFilename is not None):
        i_writeFilenameList = [i_logFilename]
    else:
        i_writeFilenameList = None
    i_dirList = [i_inputDir, i_outputDir]

    if (not radiaUtil.check_for_argv_errors(i_dirList,
                                            i_readFilenameList,
                                            i_writeFilenameList)):
        sys.exit(1)

    # get the VCF generator
    (headerDict, coordDict) = get_vcf_data(i_id, i_inputDir, i_debug)

    if (i_gzip):
        i_outputFilename = os.path.join(i_outputDir, i_id + ".vcf.gz")
    else:
        i_outputFilename = os.path.join(i_outputDir, i_id + ".vcf")

    outputFileHandler = radiaUtil.get_write_fileHandler(i_outputFilename)

    # if we have header info to output
    if (len(headerDict["metadata"]) > 0):
        # output the header information
        outputFileHandler.write("\n".join(headerDict["metadata"]) + "\n")
        outputFileHandler.write("\n".join(headerDict["filter"]) + "\n")
        outputFileHandler.write("\n".join(headerDict["info"]) + "\n")
        outputFileHandler.write("\n".join(headerDict["format"]) + "\n")
        outputFileHandler.write("".join(headerDict["chrom"]) + "\n")

    # first output the numerical chroms in order
    numericChromKeys = coordDict["numbers"].keys()
    numericChromKeys.sort(key=int)
    for chrom in numericChromKeys:
        outputFileHandler.write("\n".join(coordDict["numbers"][chrom]) + "\n")

    # then output the alphabetical chroms in order
    letterChromKeys = coordDict["letters"].keys()
    letterChromKeys.sort(key=str)
    for chrom in letterChromKeys:
        outputFileHandler.write("\n".join(coordDict["letters"][chrom]) + "\n")

    stopTime = time.time()
    logging.info("Total time for Id %s: Total time=%s hrs, %s mins, %s secs",
                 i_id, ((stopTime-startTime)/(3600)),
                 ((stopTime-startTime)/60), (stopTime-startTime))

    # close the files
    outputFileHandler.close()

    return


main()
sys.exit(0)
