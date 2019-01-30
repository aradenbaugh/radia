#!/usr/bin/env python

import sys
from optparse import OptionParser
import radiaUtil
import logging
import collections
import re


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


def get_rna_genes(anRnaGeneFile, anRnaGeneFamilyFile, anIsDebug):
    '''
    ' This function parses the RNA gene and RNA gene family blacklist files.
    '
    ' anRnaGeneFile:  An RNA gene file
    ' anRnaGeneFamilyFile:  An RNA gene family file
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''

    # open the file
    geneFileHandler = radiaUtil.get_read_fileHandler(anRnaGeneFile)
    geneFamilyFileHandler = radiaUtil.get_read_fileHandler(anRnaGeneFamilyFile)
    rnaGeneList = list()
    rnaGeneFamilyList = list()

    for line in geneFileHandler:

        # we can ignore the lines that start with # for now
        if (line.startswith("#") or line.isspace()):
            continue

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        if (anIsDebug):
            logging.debug("RNA Blacklist: %s", line)

        rnaGeneList.append(line)

    for line in geneFamilyFileHandler:

        # we can ignore the lines that start with # for now
        if (line.startswith("#") or line.isspace()):
            continue

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        if (anIsDebug):
            logging.debug("RNA Blacklist: %s", line)

        rnaGeneFamilyList.append(line)

    geneFileHandler.close()
    geneFamilyFileHandler.close()

    return rnaGeneList, rnaGeneFamilyList


def main():

    # create the usage statement
    usage = "usage: python %prog vcfFile rnaGeneFile rnaGeneFamilyFile [Opts]"
    i_cmdLineParser = OptionParser(usage=usage)

    i_cmdLineParser.add_option(
        "-o", "--outputFilename", default=sys.stdout,
        dest="outputFilename", metavar="OUTPUT_FILE",
        help="the name of the output file, STDOUT by default")
    i_cmdLineParser.add_option(
        "-l", "--log",
        dest="logLevel", default="WARNING", metavar="LOG",
        help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), " +
             "%default by default")
    i_cmdLineParser.add_option(
        "-g", "--logFilename",
        dest="logFilename", metavar="LOG_FILE",
        help="the name of the log file, STDOUT by default")
    i_cmdLineParser.add_option(
        "-c", "--allVCFCalls", action="store_false", default=True,
        dest="passedVCFCallsOnly",
        help="by default only the VCF calls that have passed all filters " +
             "thus far are processed, include this argument if all of the " +
             "VCF calls should be processed")

    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(3, 14, 1)
    i_argLength = len(sys.argv)

    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)

    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_vcfFilename = str(i_cmdLineArgs[0])
    i_rnaGeneFilename = str(i_cmdLineArgs[1])
    i_rnaGeneFamilyFilename = str(i_cmdLineArgs[2])

    # get the optional params with default values
    i_logLevel = i_cmdLineOptions.logLevel
    i_passedVCFCallsOnlyFlag = i_cmdLineOptions.passedVCFCallsOnly

    # try to get any optional parameters with no defaults
    i_outputFilename = None
    i_logFilename = None
    if (i_cmdLineOptions.outputFilename is not None):
        i_outputFilename = str(i_cmdLineOptions.outputFilename)
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
    if (i_logFilename is not sys.stdout):
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
        logging.debug("vcfFilename=%s", i_vcfFilename)
        logging.debug("rnaGeneFilename=%s", i_rnaGeneFilename)
        logging.debug("rnaGeneFamilyFilename=%s", i_rnaGeneFamilyFilename)
        logging.debug("outputFilename=%s", i_outputFilename)
        logging.debug("logFilename=%s", i_logFilename)
        logging.debug("passedOnly?=%s", i_passedVCFCallsOnlyFlag)

    # check for any errors
    i_writeFilenameList = []
    if (i_outputFilename is not sys.stdout):
        i_writeFilenameList = [i_outputFilename]
    if (i_logFilename is not None):
        i_writeFilenameList = [i_logFilename]

    i_readFilenameList = [i_vcfFilename,
                          i_rnaGeneFilename,
                          i_rnaGeneFamilyFilename]

    if (not radiaUtil.check_for_argv_errors(None,
                                            i_readFilenameList,
                                            i_writeFilenameList)):
        sys.exit(1)

    # open the input stream
    i_vcfFileHandler = radiaUtil.get_read_fileHandler(i_vcfFilename)

    # open the output stream
    if i_outputFilename is not sys.stdout:
        i_outputFileHandler = radiaUtil.get_write_fileHandler(i_outputFilename)
    else:
        i_outputFileHandler = i_outputFilename

    # get the RNA gene blacklists
    (i_rnaGeneList,
     i_rnaGeneFamilyList) = get_rna_genes(i_rnaGeneFilename,
                                          i_rnaGeneFamilyFilename,
                                          i_debug)

    hasAddedFilterHeader = False

    for line in i_vcfFileHandler:

        if (i_debug):
            logging.debug("vcfLine: %s", line)

        # if it is an empty line, then just continue
        if (line.isspace()):
            continue

        # if we find the FILTER section, then add the filters from here
        elif ((not hasAddedFilterHeader) and (line.startswith("##FILTER"))):
            hasAddedFilterHeader = True
            i_outputFileHandler.write(
                "##FILTER=<ID=rgene,Description=\"This gene is on the " +
                "RNA gene blacklist\">\n")
            i_outputFileHandler.write(
                "##FILTER=<ID=rgfam,Description=\"This gene family is on " +
                "the RNA gene family blacklist\">\n")
            i_outputFileHandler.write(line)

        # these lines are from previous scripts in the pipeline, so output them
        elif (line.startswith("#")):
            i_outputFileHandler.write(line)

        # if we are only suppose to process the passed calls
        # and this call has not passed, then skip it
        elif (i_passedVCFCallsOnlyFlag and "PASS" not in line):
            i_outputFileHandler.write(line)

        # now we are to the data
        else:

            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")

            # split the line on the tab
            splitLine = line.split("\t")

            filterSet = set(splitLine[6].split(";"))

            # if there are no filters so far, then clear the list
            if (len(filterSet) == 1 and "PASS" in filterSet):
                filterSet = set()

            # parse the info column and create a dict
            infoList = splitLine[7].split(";")
            infoDict = collections.defaultdict(list)
            for info in infoList:
                keyValueList = info.split("=")
                # some keys are just singular without a value (e.g. DB, etc.)
                if (len(keyValueList) == 1):
                    infoDict[keyValueList[0]] = ["True"]
                else:
                    # the value can be a comma separated list
                    infoDict[keyValueList[0]] = keyValueList[1].split(",")

            effectList = infoDict["EFF"]
            effectRegEx = re.compile("(\\w).*\\({1}")
            ignoreEffectsList = ["UPSTREAM", "DOWNSTREAM"]

            isRnaBlacklistGene = False
            isRnaBlacklistGeneFamily = False

            for rawEffect in effectList:
                rawEffect = rawEffect.rstrip(")")
                iterator = effectRegEx.finditer(rawEffect)

                # for each match object in the iterator
                for match in iterator:
                    effect = match.group()
                    rawEffect = rawEffect.replace(effect, "")
                    effect = effect.rstrip("(")

                if (effect in ignoreEffectsList):
                    continue

                effectParts = rawEffect.split("|")
                # effectImpact = effectParts[0]
                # functionalClass = effectParts[1]
                # codonChange = effectParts[2]
                # aaChange = effectParts[3]
                # aaLength = effectParts[4]
                geneName = effectParts[5]
                transcriptBiotype = effectParts[6]
                # geneCoding = effectParts[7]
                # ensembleId = effectParts[8]
                # exonNumber = effectParts[9]
                # genotypeNumber = effectParts[10]

                # the RNA gene list can have "RP11" and that
                # should filter out any gene with RP11 in it
                for rnaGene in i_rnaGeneList:
                    if (rnaGene in geneName):
                        isRnaBlacklistGene = True
                        break

                if (transcriptBiotype in i_rnaGeneFamilyList):
                    isRnaBlacklistGeneFamily = True

            output = ["\t".join(splitLine[0:6])]

            # if the filter should be applied
            if (isRnaBlacklistGene):
                filterSet.add("rgene")
            # if the filter should be applied
            if (isRnaBlacklistGeneFamily):
                filterSet.add("rgfam")

            # if there are no filters so far, then this call passes
            if (len(filterSet) == 0):
                filterSet.add("PASS")

            output.append(";".join(filterSet))

            output.append("\t".join(splitLine[7:]))

            if (i_outputFilename is not sys.stdout):
                i_outputFileHandler.write("\t".join(output) + "\n")
            else:
                print >> sys.stdout, "\t".join(output)

    # close the files
    i_vcfFileHandler.close()
    if (i_outputFilename is not sys.stdout):
        i_outputFileHandler.close()

    return


main()
sys.exit(0)
