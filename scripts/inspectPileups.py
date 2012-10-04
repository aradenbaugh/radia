#!/usr/bin/env python2.6

import sys
import os
import re
import time
import logging
from itertools import izip
import collections
import subprocess
import rnaEditingUtil
from optparse import OptionParser


# this regular expression is used to remove insertions and deletions from raw reads
# a read could look like:  "T$TT+3AGGGT+2AG+2AG.-2AGGG..-1A"
# insertions start with a "+", deletions with a "-"
# in theory, there could be multiple digits
i_numOfIndelsRegEx = re.compile("[+-]{1}(\\d+)")

# this regular expression will match any number of valid cDNA strings
i_cDNARegEx = re.compile("[ACGTNacgtn]+")

# this regular expression will match full TCGA sample Ids, e.g. TCGA-AG-A016-01A-01R or TCGA-37-4133-10A-01D
i_tcgaNameRegEx = re.compile("TCGA-(\\w){2}-(\\w){4}-(\\w){3}-(\\w){3}")


def execute_samtools_cmd(aBamFile, aFastaFile, aBaseQuality, aMappingQuality, aUseBAQFlag, aChrom, aStartCoordinate, aStopCoordinate, aUseChrPrefix, anIsDebug):
    '''
    ' This function executes an external command.  The command is the "samtools mpileup" command which returns all 
    ' the information about the sequencing reads for specific coordinates.  There are two things to be careful about
    ' when using the samtools mpileup command.  Some .bam files use the 'chr' prefix when specifying the region to 
    ' select with the -r argument.  If the 'chr' prefix is required, then specify the --useChrPrefix argument and also
    ' make sure that the fasta file that is specified also has the 'chr' prefix.  Here are some examples of the commands
    ' that can be copy/pasted to the command line to view the output:
    '
    ' samtools mpileup -f /inside/home/dzerbino/data/hg19.fa -Q 20 -q 10 -r chr1:10000-20000 /inside/depot/wlee/bccagsc_070610/TCGA-AB-2995*/all_reads.bam
    ' samtools mpileup -f /inside/depot/fa/all_sequences.fasta -Q 20 -q 10 -r chr1:855155-1009900 /inside/depot/lusc/exome/TCGA-21-1070-01A-01W-0782-08.bam
    '
    ' aBamFile:            A .bam file to be read from
    ' aFastaFile:          The FASTA file which is needed for the reference base.
    ' aBaseQuality:        The base quality score for the samtools command
    ' aMappingQuality:     The mapping quality score for the samtools command
    ' aChrom:              The chromosome that we are selecting from
    ' aStartCoordinate:    The start coordinate of the selection
    ' aStopCoordinate:     The stop coordinate of the selection
    ' aUseChrPrefix:       Whether the 'chr' should be used in the samtools command
    ' anIsDebug:           A flag for outputting debug messages to STDERR
    '''
    
    # create the samtools command
    if (aUseChrPrefix):
        if (aUseBAQFlag):
            samtoolsSelectStatement = "samtools mpileup -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r chr" + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
        else:
            samtoolsSelectStatement = "samtools mpileup -B -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r chr" + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
            
    else:
        if (aUseBAQFlag):
            samtoolsSelectStatement = "samtools mpileup -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r " + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
        else:
            samtoolsSelectStatement = "samtools mpileup -B -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r " + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
    
    # output the samtools command
    if (anIsDebug):
        logging.debug(samtoolsSelectStatement)
        # keep track of how long it takes to run the samtools command
        timeSamtoolsStart = time.time()
        
    # execute the samtools command
    samtoolsCall = subprocess.Popen(samtoolsSelectStatement, shell=True, bufsize=-1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    #samtoolsCall = subprocess.Popen(samtoolsSelectStatement, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    samtoolsCall.poll()
    
    if (anIsDebug):
        timeSamtoolsEnd = time.time()
        timeSpent = timeSamtoolsEnd-timeSamtoolsStart
        logging.debug("Time spent executing samtools command: %s hrs, %s mins, %s secs", (timeSpent/3600), (timeSpent/60), (timeSpent))
    
    # get the output from executing the samtools command
    pileups = samtoolsCall.stdout.readlines()
    
    # for some reason, we have to get the stdout.readlines() before the stderr.readlines(), otherwise the process hangs
    if (anIsDebug and samtoolsCall.returncode != 0):
        print >> sys.stderr, "Warning from", samtoolsSelectStatement, "\n", samtoolsCall.stderr.readlines()
    
    # kill the command
    samtoolsCall.kill()
    
    return pileups

            
def get_pileup(aBamFile, aFastaFile, aBaseQuality, aMappingQuality, aUseBAQFlag, aChrom, aStartCoordinate, aStopCoordinate, aUseChrPrefix, anIsDebug):
    
    
    pileups = execute_samtools_cmd(aBamFile, aFastaFile, aBaseQuality, aMappingQuality, aUseBAQFlag, aChrom, aStartCoordinate, aStopCoordinate, aUseChrPrefix, anIsDebug)
    
    if (len(pileups) > 1):
        logging.warning("Input file is trying to query more than one coordinate at a time (%s to %s), only the first pileup is returned", aStartCoordinate, aStopCoordinate)
    
    # if there was data for these coordinates
    if (len(pileups) > 0):
        for line in pileups:
            # if the samtools select statement returns no reads which can happen when the selection is done 
            # in an area with no reads, then a warning message will be returned that starts with "[mpileup]".  
            # Let's just ignore it
            if (line.isspace() or line.startswith("[mpileup]") or line.startswith("<mpileup>")):
                continue;

            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")

            # split the line on the tab
            splitLine = re.split("\t", line)
        
            if (anIsDebug):    
                logging.debug("Pileup: %s", line)

            # the coordinate is the second element
            #chr = splitLine[0]
            #coordinate = int(splitLine[1])
            reference = splitLine[2].upper()
            #numOfReads = int(splitLine[3])
            reads = splitLine[4]
            qualScores = splitLine[5]
                
            # return all the information about the current coordinate
            return (reference, reads, qualScores)
    return


def convert_raw_reads(aStringOfRawReads, aStringOfQualScores, aReferenceBase, anIsDebug):
    '''
    ' This function returns all of the valid RNA (cDNA) or DNA bases from the given pile-up of read bases.
    ' It converts all of the samtools specific characters into human-readable bases and filters out any non 
    ' RNA/DNA characters. 
    '
    ' This is from the samtools documentation:
    '
    ' In the pileup format, each line represents a genomic position, consisting of chromosome name, 
    ' 1-based coordinate, reference base, read bases, read qualities and alignment mapping qualities. 
    ' Information on match, mismatch, indel, strand, mapping quality and start and end of a read are all 
    ' encoded at the read base column. At this column, a dot stands for a match to the reference base on 
    ' the forward strand, a comma for a match on the reverse strand, a ">" or "<" for a reference skip, 
    ' "ACGTN" for a mismatch on the forward strand and "acgtn" for a mismatch on the reverse strand. A 
    ' pattern "\+[0-9]+[ACGTNacgtn]+" indicates there is an insertion between this reference position and 
    ' the next reference position. The length of the insertion is given by the integer in the pattern, 
    ' followed by the inserted sequence. Similarly, a pattern "-[0-9]+[ACGTNacgtn]+" represents a deletion 
    ' from the reference. The deleted bases will be presented as "*" in the following lines. Also at the 
    ' read base column, a symbol "^" marks the start of a read. The ASCII of the character following "^" 
    ' minus 33 gives the mapping quality. A symbol "$" marks the end of a read segment.
    '
    ' We are converting all dots and commas to the upper case reference base.  Even though the comma represents 
    ' a match on the reverse strand, there is no need to take the complement of it, since samtools does
    ' that for us.  We are converting all mismatches on the reverse strand to upper case as well, and again 
    ' no complement is needed.
    '
    ' We are ignoring the following for now:
    ' 1) Reference skips (">" and "<") 
    ' 2) "N" in the reads
    '
    ' aStringOfRawReads: A string representing the pile-up of read bases from a samtools mpileup command 
    ' aStringOfQualScores: A string representing the raw quality scores for the read bases from the mpileup command
    ' aReferenceBase: Used to convert "." and "," from the samtools mpileup command
    '''
    
    # Note:  Reverse strand mismatches have been reverse-complemented by samtools
       
    # remove insertions and deletions
    # a read could look like:  "T$TT+3AGGGT+2AG+2AG.-2AGGG..-1A"
    # insertions start with a "+", deletions with a "-"
    # in theory, there could be multiple digits
    # i_numOfIndelsRegEx = re.compile("[+-]{1}(\\d+)")

    # if we have an indel
    if ("+" in aStringOfRawReads or "-" in aStringOfRawReads):
        # get an iterator of match objects for all indels
        iterator = i_numOfIndelsRegEx.finditer(aStringOfRawReads)
        
        # for each match object in the iterator
        for match in iterator:
            # get the pattern that matched the reg ex, i.e. +3 or -2
            indel = match.group()
            # the length of the indel is the number following the + or - sign
            lengthOfIndel = indel[1:len(indel)]
            # this reg ex will specifically match indels with the specified length, i.e. +3AGG, -2AG
            indelRegEx = re.compile("\\" + indel + "[ACGTNacgtn=]{" + lengthOfIndel + "}")
            
            # we can simply remove the indels and replace them with empty strings for now
            # there are no base quality scores for indels that need to be removed
            aStringOfRawReads = indelRegEx.sub("", aStringOfRawReads) 
                               
    # remove all start of read symbols "^" (plus the following quality score)
    # there are no base quality scores for start symbols that need to be removed
    while ("^" in aStringOfRawReads):
        start = aStringOfRawReads.find("^")
        end = start+2
        # replace will replace all unless a max is set, but we don't care, 
        # we want to get rid of all of them
        aStringOfRawReads = aStringOfRawReads.replace(aStringOfRawReads[start:end], "")
    
    # remove all end of read symbols "$"
    # there are no base quality scores for stop symbols that need to be removed
    aStringOfRawReads = aStringOfRawReads.replace("$", "")
    
    # replace all the periods for uppercase references representing the plus strand
    # replace all the commas for lowercase references representing the minus strand
    aStringOfRawReads = aStringOfRawReads.replace(".", aReferenceBase.upper())
    aStringOfRawReads = aStringOfRawReads.replace(",", aReferenceBase.lower())
        
    # get an iterator of match objects for all valid cDNA
    # this regular expression will match any number of valid cDNA strings
    # i_cDNARegEx = re.compile("[ACGTacgt]+")
    iterator = i_cDNARegEx.finditer(aStringOfRawReads)
    
    # create final strings consisting of just the valid cDNA and corresponding qual scores
    finalPileups = ""
    finalQuals = ""
    
    # only extract the valid cDNA and corresponding qual scores
    # ignore >", "<", etc.
    for match in iterator:
        start = match.start()
        end = match.end()
        finalPileups += aStringOfRawReads[start:end]
        finalQuals += aStringOfQualScores[start:end]
              
    # get the lengths
    lenFinalPileups = len(finalPileups)
    lenFinalQuals = len(finalQuals) 
    
    # at this point, the length of the pileups string should be equal to the length of the quality scores
    if (lenFinalPileups != lenFinalQuals):
        logging.error("convert_raw_reads() Error:  The length %s of the final pileup of reads is != the length %s of the final quality scores.  Original Pileup=%s, Final Pileup=%s, Original QualScores=%s, Final QualScores=%s", lenFinalPileups, lenFinalQuals, aStringOfRawReads, finalPileups, aStringOfQualScores, finalQuals)
     
    return (finalPileups, finalQuals)     


def filter_by_base_quality(aStringOfReads, aStringOfQualScores, aMinBaseQualityScore, anIsDebug):
    '''
    ' This function filters out all the bases that don't meet the user-specified minimum 
    ' base quality score which is specified here with the "aMinBaseQualityScore" parameter.
    '
    ' aStringOfReads: A string representing the pile-up of reads from a samtools mpileup command
    ' aStringOfQualScores: A string representing the raw quality scores for the read bases from the mpileup command 
    ' aMinBaseQualityScore: An integer with the user-specified minimum base quality score (also used as -Q parameter to samtools mpileup command)
    '''
    
    # create strings consisting of just the reads that are greater than or equal to the minimum base quality score 
    pileups = ""
    qualScores = ""
                
    # loop through the reads and the corresponding quality scores
    for (base, rawScore) in izip(aStringOfReads, aStringOfQualScores):
        convertedScore = ord(rawScore)-33
        # the scores are in ascii, so convert them to integers
        if (convertedScore >= aMinBaseQualityScore):
            
            # convert all to plus strand after counting
            base = base.upper()
            
            # keep the base and quality
            pileups += base
            qualScores += str(convertedScore) + ","
                
    return (pileups, qualScores)               


def get_input(anInputFile, anIsDebug):
    
    # open the input file
    inputFileHandler = open(anInputFile, "r")
     
    for line in inputFileHandler:
          
        if (line.isspace() or line.startswith("#")):
            continue;

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        if (anIsDebug):
            logging.debug("Input line: %s", line)    
        
        # split the .sam line on the tab
        splitLine = re.split("\t", line)

        # the coordinate is the second element
        id = splitLine[0]
        chr = splitLine[1]
        startCoordinate = splitLine[2]
        stopCoordinate = splitLine[3]
        bamFile = splitLine[4]
        refFastaFile = splitLine[5]
        baseQuality = splitLine[6]
        mappingQuality = splitLine[7]
        useChrPrefix = splitLine[8]
        
        if (useChrPrefix == "True"):
            useChrPrefixBool = True
        else:
            useChrPrefixBool = False
        
        # yield all the information about the current coordinate
        yield (id, chr, startCoordinate, stopCoordinate, bamFile, refFastaFile, baseQuality, mappingQuality, useChrPrefixBool)

    inputFileHandler.close()
    return


def main():
    # python2.7 inspectPileups.py ../data/test/inspectPileups.txt -l DEBUG
    
    # InputFile: aTcgaId, aChrom, aStartCoordinate, aStopCoordinate, (n * (aBamFile, aFastaFile, aBaseQuality, aMappingQuality, aUseChrPrefix)), anIsDebug
    # OutputFile: aTcgaId, aChrom, aStartCoordinate, aStopCoordinate, (n * (OriginalPileup, OrigianlBQs, AfterBaqPileup, AfterBaqBQs, AfterCovertRawPileup, AfterConvertRawBQs, AfterFilterBQPileup, AferFilterBQBQs))
    
    usage = "usage: python2.7 %prog inputFile [options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    i_cmdLineParser.add_option("-m", "--minBaseQuality", type="int", default=int(10), dest="minBaseQuality", metavar="MIN_BASE_QUAL", help="the minimum final base quality, %default by default")
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, STDOUT by default")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(2,5,2)
    i_argLength = len(sys.argv)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)
    
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_inputFilename = str(i_cmdLineArgs[0])
    
    # get the optional params with default values
    i_logLevel = i_cmdLineOptions.logLevel 
    i_minBaseQual = i_cmdLineOptions.minBaseQuality
        
    i_dirList = []
    i_readFilenameList = [i_inputFilename]
    i_writeFilenameList = []
    i_outputFilename = None
    i_outputFileHandler = None
    i_logFilename = None
    if (i_cmdLineOptions.outputFilename != None):
        i_outputFilename = str(i_cmdLineOptions.outputFilename)
        i_outputFileHandler = open(i_outputFilename, "w")
        i_writeFilenameList += [i_outputFilename]
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
        i_writeFilenameList += [i_logFilename]
        
    # check for any errors
    if (not rnaEditingUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)
        
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
        logging.debug("inputFilename=%s" % i_inputFilename)
        logging.debug("outputFilename=%s" % i_outputFilename)
        logging.debug("logLevel=%s" % i_logLevel)
        logging.debug("logFilename=%s" % i_logFilename)
        
    # get a generator for the lines in the input file
    i_inputGenerator = get_input(i_inputFilename, i_debug)
    
    for (id, chr, startCoordinate, stopCoordinate, bamFile, refFastaFile, baseQuality, mappingQuality, useChrPrefix) in i_inputGenerator:
        if (i_debug):
            logging.debug("Input Data: %s %s %s %s %s %s %s %s %s", id, chr, startCoordinate, stopCoordinate, bamFile, refFastaFile, baseQuality, mappingQuality, useChrPrefix)
            
        '''
        if not os.path.exists(bamFile):
            logging.warning('bam file: %s  does not exist.' % bamFile)
            continue
        
        if not os.path.exists(refFastaFile):
            logging.warning('reference fasta file: %s  does not exist.' % refFastaFile)
            continue
        '''
            
        outputList = [id, chr, startCoordinate, stopCoordinate]
        
        # get the pileups without BAQ
        (ref, originalPileup, originalBQs) = get_pileup(bamFile, refFastaFile, baseQuality, mappingQuality, False, chr, startCoordinate, stopCoordinate, useChrPrefix, i_debug)
        outputList += [originalPileup, originalBQs]
        
        # get the pileups with BAW
        (ref, baqPileup, baqBQs) = get_pileup(bamFile, refFastaFile, baseQuality, mappingQuality, True, chr, startCoordinate, stopCoordinate, useChrPrefix, i_debug)
        outputList += [baqPileup, baqBQs]
        
        # convert the raw reads
        (convertedPileup, covertedBQs) = convert_raw_reads(baqPileup, baqBQs, ref, i_debug)
        outputList += [convertedPileup, covertedBQs]
        
        # filter the raw reads
        (filteredPileup, filteredBQs) = filter_by_base_quality(convertedPileup, covertedBQs, i_minBaseQual, i_debug)
        outputList += [filteredPileup, filteredBQs]
        
        # output the info
        if (i_outputFilename is None):
            print >> sys.stdout, "\t".join(outputList)
        else:
            i_outputFileHandler.write(outputList)
        
    if (i_outputFilename is not None):
        i_outputFileHandler.close()
            
    return

main()
sys.exit(0)
