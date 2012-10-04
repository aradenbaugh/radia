#!/usr/bin/env python2.7

import sys
import time
import re
import subprocess
import datetime
import logging
#from argparse import ArgumentParser
from optparse import OptionParser
from itertools import izip
import rnaEditingUtil
import collections
#import cProfile
#import pstats


'''
'    Amie Radenbaugh - 09/01/2011
'    UCSC - RNA and DNA Integrated Analysis (RADIA)
'    Program name: "radia.py"
'
'    This program identifies RNA and DNA variants in .bam files.  The program is designed
'    to take in 4 .bam files:  DNA Normal, RNA Normal, DNA Tumor, and RNA Tumor.  For the 
'    normal DNA, the program outputs any differences compared to the reference which could 
'    be potential Germline mutations.  For the normal RNA, the program outputs any differences 
'    compared to the reference and the normal DNA which could be potential normal RNA-Editing
'    events.  For the tumor DNA, the program outputs any difference compared to the reference, 
'    normal DNA and normal RNA which could be potential Somatic mutations.  For the tumor
'    RNA, the program outputs any difference compared to the reference, normal DNA, normal 
'    RNA and tumor DNA which could be potential RNA-editing events.
'
'    The program is designed for 4 .bam files, but the user can also specify just one or two.
'    The program will report RNA and DNA variants.
'
'    This could be used at the top of BED files that are uploaded to the browser:
'        "browser position chr12:7055740-7070479"
'        "track name=rnaEditing description=RNA-Editing visibility=2 itemRgb=On"
'''

# this regular expression is used to remove insertions and deletions from raw reads
# a read could look like:  "T$TT+3AGGGT+2AG+2AG.-2AGGG..-1A"
# insertions start with a "+", deletions with a "-"
# in theory, there could be multiple digits
i_numOfIndelsRegEx = re.compile("[+-]{1}(\\d+)")

# this regular expression will match any number of valid cDNA strings
i_cDNARegEx = re.compile("[ACGTNacgtn]+")

# this regular expression will match full TCGA sample Ids, e.g. TCGA-AG-A016-01A-01R or TCGA-37-4133-10A-01D
i_tcgaNameRegEx = re.compile("TCGA-(\\w){2}-(\\w){4}-(\\w){3}-(\\w){3}")


def get_chrom_size(aChrom, anInputStream, anIsDebug):
    '''
    ' This function reads from a .tab file that contains the sizes for each chromosome.
    ' The file has two columns:  the first column has the chromosome identifier, the
    ' second column has the size of the chromosome.  The columns are separated by tabs.
    ' Here is an example:
    ' chr1    249250621
    ' chr2    243199373
    ' ...       ...
    '
    ' aChrom: The chrom size to return
    ' anInputStream: The input stream for the chrom sizes file
    '''
     
    for line in anInputStream:
          
        # if it is an empty line, then just continue
        if (line.isspace() or line.startswith("#")):
            continue;

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        if (anIsDebug):
            logging.debug("get_chrom_size(): looking for size of chrom %s, line=%s", aChrom, line)	
        
        # split the line on the tab
        splitLine = re.split("\t", line)

        # the coordinate is the second element
        chr = splitLine[0]
        size = int(splitLine[1])
        
        # the chromosome from the file has the following format "chr1"
        # the aChrom variable only has the integer 1
        if (chr == "chr" + str(aChrom)):
            if (anIsDebug):
                logging.debug("get_chrom_size(): found size of chrom %s, size=%s", aChrom, size)
            return size
    return


def get_batch_end_coordinate(aStartCoordinate, anEndCoordinate, aBatchSize):
    '''
    ' This function takes a start coordinate, an end coordinate, and a batch size and
    ' returns the next appropriate batch end coordinate which is either the start coordinate
    ' plus the batch size if this is less than the final end coordinate otherwise the end
    ' coordinate.
    '
    ' aStartCoordinate:  A start coordinate
    ' anEndCoordinate:  A stop coordinate
    ' aBatchSize:  A batch size
    '''
    if ((aStartCoordinate + aBatchSize) < anEndCoordinate):
        return (aStartCoordinate + aBatchSize)
    else:
        return (anEndCoordinate)
    

def get_sam_data(aSamFile, anIsDebug):
    '''
    ' This function is used during testing to read data from a .sam input file.  It uses the python 
    ' generator to yield the information for one coordinate position at a time.  It ignores empty lines 
    ' and strips trailing \r\n characters.  This function yields the chromosome, coordinate, reference base, 
    ' number of reads, raw reads, and quality scores.
    '
    ' aSamFile:  A .sam file
    '''
    
    # open the sam file
    samFileHandler = open(aSamFile, "r")
     
    for line in samFileHandler:
          
        # if the samtools select statement returns no reads which can happen when the batch size is
        # small and the selection is done in an area with no reads, then a warning message will be
        # returned that starts with "[mpileup]".  We can ignore the message and move on to the next
        # select statement.
        if (line.isspace() or line.startswith("[mpileup]")):
            continue;

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        if (anIsDebug):
            logging.debug("Original SAM pileup: %s", line)	
        
        # split the .sam line on the tab
        splitLine = re.split("\t", line)

        # the coordinate is the second element
        chr = splitLine[0]
        coordinate = int(splitLine[1])
        reference = splitLine[2].upper()
        numOfReads = int(splitLine[3])
        reads = splitLine[4]
        qualScores = splitLine[5]
       
        # yield all the information about the current coordinate
        yield (chr, coordinate, reference, numOfReads, reads, qualScores)

    samFileHandler.close()
    return
    

def get_bam_data(aBamFile, aFastaFile, aBaseQual, aMappingQual, aChrom, aStartCoordinate, aStopCoordinate, aBatchSize, aUseChrPrefix, anIsRna, anIsDebug):
    '''
    ' This function uses the python generator to yield the information for one coordinate at a time.
    ' In order to reduce the time and memory overhead of loading the entire .bam file into memory at
    ' once, this function reads in chunks of data at a time.  The number of coordinates that should be
    ' read into memory at a given time is determined by the "aBatchSize" parameter.  This function uses the 
    ' samtools "mpileup" command to make a selection. 
    '
    ' The original start and end coordinates are specified by the "aStartCoordinate" and "anEndCoordinate" 
    ' parameters and are typically initialized to 0 and the size of the chromosome respectively. This function
    ' will loop over the .bam file, selecting "aBatchSize" number of coordinates into memory at once.  Each line that
    ' is selected will be processed and yielded using the python generator.  When all lines from the current batch 
    ' are processed, the start and end coordinates will be incremented, and the next selection will be made from the
    ' .bam file.  This process continues until the end of the chromosome has been reached.
    '
    ' This function yields the chromosome, coordinate, reference base, number of reads, raw reads, and the quality scores.
    '
    ' aBamFile:            A .bam file to be read from
    ' aFastaFile:          The FASTA file that should be used in the samtools command which is needed for the reference base.
    ' aBaseQual:           The base quality score that should be used in the samtools command
    ' aMappingQual:        The mapping quality score that should be used in the samtools command
    ' aChrom:              The chromosome that should be used in the samtools command
    ' aStartCoordinate:    The initial start coordinate (typically zero)
    ' aStopCoordinate:     The initial stop coordinate (typically the size of the chromosome)
    ' aBatchSize:          The number of coordinates to load into memory at one time
    ' aUseChrPrefix:       Whether the 'chr' should be used in the region parameter of the samtools command
    '''
    
    # initialize the first start and stop coordinates
    # the stop coordinate is calculated according to the "aBatchSize" param
    currentStartCoordinate = aStartCoordinate
    currentStopCoordinate = get_batch_end_coordinate(currentStartCoordinate, aStopCoordinate, aBatchSize)

    # while we still have coordinates to select from the .bam file
    while (currentStartCoordinate < aStopCoordinate):
        # execute the samtools command
        #pileupsBatch = execute_samtools_cmd_new(aBamFile, aFastaFile, aBaseQual, aMappingQual, aChrom, aUseChrPrefix, currentStartCoordinate, currentStopCoordinate, anIsDebug) 
        #pileups = re.split("\n", pileupsBatch)
        pileups = execute_samtools_cmd(aBamFile, aFastaFile, aBaseQual, aMappingQual, aChrom, aUseChrPrefix, currentStartCoordinate, currentStopCoordinate, anIsRna, anIsDebug)
        
        if (anIsDebug):        
            logging.debug("samtools number of lines selected from %s to %s = %s", currentStartCoordinate, currentStopCoordinate, len(pileups))

        # if there was data for these coordinates
        if (len(pileups) > 0):
            # for each line representing one coordinate
            for line in pileups:
                # if the samtools select statement returns no reads which can happen when the batch size is
                # small and the selection is done in an area with no reads, then a warning message will be
                # returned that starts with "[mpileup]".  We can ignore the message and move on to the next
                # select statement.
                if (line.isspace() or line.startswith("[mpileup]") or line.startswith("<mpileup>")):
                    continue;

                # strip the carriage return and newline characters
                line = line.rstrip("\r\n")

                # split the line on the tab
                splitLine = re.split("\t", line)
            
                if (anIsDebug):    
                    logging.debug("Original BAM pileup: %s", line)

                #if (len(splitLine) > 1):
                # the coordinate is the second element
                chr = splitLine[0]
                coordinate = int(splitLine[1])
                reference = splitLine[2].upper()
                numOfReads = int(splitLine[3])
                reads = splitLine[4]
                qualScores = splitLine[5]
                    
                # yield all the information about the current coordinate
                yield (chr, coordinate, reference, numOfReads, reads, qualScores)
                
        # calculate a new start and stop coordinate for the next select statement
        currentStartCoordinate = currentStartCoordinate + aBatchSize
        currentStopCoordinate = get_batch_end_coordinate(currentStartCoordinate, aStopCoordinate, aBatchSize)
    
    return

'''
def execute_samtools_cmd_new(aBamFile, aFastaFile, aBaseQuality, aMappingQuality, aChrom, aUseChrPrefix, aStartCoordinate, aStopCoordinate, anIsDebug):
    
    ' This function executes an external command.  The command is the "samtools mpileup" command which returns all 
    ' the information about the sequencing reads for specific coordinates.  There are two things to be careful about
    ' when using the samtools mpileup command.  Some .bam files use the 'chr' prefix when specifying the region to 
    ' select with the -r argument.  If the 'chr' prefix is required, then specify the --useChrPrefix argument and also
    ' make sure that the fasta file that is specified also has the 'chr' prefix.  Here are some examples of the commands
    ' that can be copy/pasted to the command line to view the output:
    '
    ' samtools mpileup -f /inside/home/dzerbino/data/hg19.fa -Q 20 -q 10 -r chr1:10000-20000 /inside/depot/wlee/bccagsc_070610/TCGA-AB-2995*/all_reads.bam
    ' samtools mpileup -f /inside/depot/fa/all_sequences.fasta -Q 20 -q 10 -r chr1:855155-1009900 /inside/depot/lusc/exome/TCGA-21-1070-01A-01W-0782-08.bam
    ' samtools mpileup -f /inside/depot/fa/all_sequences.fasta -Q 10 -q 10 -r Y:5721100-5721100 /inside/home/cwilks/bb_pipeline/runs/mutcomp_ex3_usrerun/bams/TCGA-AG-A008*10A*IlluminaGA*.bam
    ' samtools mpileup -f /inside/depot/fa/all_sequences.fasta -Q 10 -q 10 -r Y:5721100-5721100 /inside/home/cwilks/bb_pipeline/runs/mutcomp_ex3_usrerun/bams/TCGA-AG-A008*01A*IlluminaGA*.bam
    '
    ' aBamFile:            A .bam file to be read from
    ' aFastaFile:          The FASTA file which is needed for the reference base.
    ' aBaseQuality:        The base quality score for the samtools command
    ' aMappingQuality:     The mapping quality score for the samtools command
    ' aChrom:              The chromosome that we are selecting from
    ' aStartCoordinate:    The start coordinate of the selection
    ' aStopCoordinate:     The stop coordinate of the selection
    ' aUseChrPrefix:       Whether the 'chr' should be used in the samtools command
    
    
    if (aUseChrPrefix):
        # create the samtools command
        samtoolsSelectStatement = "samtools mpileup -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r chr" + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
        #print >> sys.stderr, samtoolsSelectStatement
    else:
        # create the samtools command
        samtoolsSelectStatement = "samtools mpileup -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r " + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
    
    # output the samtools command
    if (anIsDebug):
        logging.debug(samtoolsSelectStatement)
    
    # execute the samtools command
    timeSamtoolsStart = time.time()
    samtoolsCall = subprocess.Popen(samtoolsSelectStatement, shell=True, bufsize=-1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    
    # communicate() waits for the process to finish
    (pileups, samtoolsStdErr) = samtoolsCall.communicate()
    
    timeSamtoolsEnd = time.time()
    timeSpent = timeSamtoolsEnd-timeSamtoolsStart
    
    if (anIsDebug):
        logging.debug("Time spent executing samtools command: %s hrs, %s mins, %s secs", (timeSpent/3600), (timeSpent/60), (timeSpent))  
    
    # if the return code is None, then  the process is not yet finished
    # communicate() waits for the process to finish, poll() does not
    if (samtoolsCall.returncode == None):
        logging.warning("The samtools mpileup command is not done, and you are moving on...possibly without all the data?")
    # if samtools returned a return code != 0, then an error occurred
    # warning: previous versions of samtools did not return a return code!
    elif (samtoolsCall.returncode != 0):
        logging.warning("The return code of '%s' from the samtools mpileup command indicates an error.", samtoolsCall.returncode)
        logging.warning("Warning/error from %s:\n%s", samtoolsSelectStatement, samtoolsStdErr)
    
    return pileups
'''


def execute_samtools_cmd(aBamFile, aFastaFile, aBaseQuality, aMappingQuality, aChrom, aUseChrPrefix, aStartCoordinate, aStopCoordinate, anIsRna, anIsDebug):
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
        if (anIsRna):
            # add the -A option
            #samtoolsSelectStatement = "samtools mpileup -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -A -r chr" + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
            samtoolsSelectStatement = "samtools mpileup -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r chr" + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
        else:
            samtoolsSelectStatement = "samtools mpileup -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r chr" + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
            
    else:
        if (anIsRna):
            # add the -A option
            #samtoolsSelectStatement = "samtools mpileup -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -A -r " + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
            samtoolsSelectStatement = "samtools mpileup -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r " + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
        else:
            samtoolsSelectStatement = "samtools mpileup -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r " + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
    
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
    
    # initialize some counts
    indelCount = 0
    
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
            indelCount += 1
            # get the pattern that matched the reg ex, i.e. +3 or -2
            indel = match.group()
            # the length of the indel is the number following the + or - sign
            lengthOfIndel = indel[1:len(indel)]
            # this reg ex will specifically match indels with the specified length, i.e. +3AGG, -2AG
            indelRegEx = re.compile("\\" + indel + "[ACGTNacgtn=]{" + lengthOfIndel + "}")
            
            # we can simply remove the indels and replace them with empty strings for now
            # there are no base quality scores for indels that need to be removed
            aStringOfRawReads = indelRegEx.sub("", aStringOfRawReads) 
            
        if (indelCount > 0):
            logging.debug("%s indels found in %s", indelCount, aStringOfRawReads)
            
    # count starts and stops
    starts = aStringOfRawReads.count("^")
    stops = aStringOfRawReads.count("$")
        
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
     
    return (finalPileups, finalQuals, lenFinalPileups, starts, stops, indelCount)     


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
    numBasesDict = collections.defaultdict(int)
    sumBaseQualsDict = collections.defaultdict(int)
    numPlusStrandDict = collections.defaultdict(int)
                
    # loop through the reads and the corresponding quality scores
    for (base, rawScore) in izip(aStringOfReads, aStringOfQualScores):
        convertedScore = ord(rawScore)-33
        # the scores are in ascii, so convert them to integers
        if (convertedScore >= aMinBaseQualityScore):
            
            # count the ones on the plus strand
            if (base in "ACGTN"):
                numPlusStrandDict[base] += 1
            # convert all to plus strand after counting
            else:
                base = base.upper()
            
            # keep the base and quality
            pileups += base
            qualScores += rawScore
            
            # keep track of the number of each base and the corresponding qual score
            numBasesDict[base] += 1
            sumBaseQualsDict[base] += convertedScore
                
    return (pileups, qualScores, len(pileups), numBasesDict, sumBaseQualsDict, numPlusStrandDict)               


def format_bam_output(aChrom, aRefList, anAltList, anAltCountsDict, anAltPerDict, aStringReads, aStringQualScores, aNumBases, aStartsCount, aStopsCount, anIndelCount, aBaseCountsDict, aQualitySumsOfBasesDict, aPlusStrandCountsDict, anIsDebug):
    '''
    ' This function converts information from a .bam mpileup coordinate into a format that can be output to a VCF formatted file.
    ' This function calculates the average overall base quality score, strand bias, and fraction of reads supporting the alternative.
    ' It also calculates the allele specific depth, average base quality score, strand bias, and fraction of reads supporting the alternative.
    ' The format for the output in VCF is:  GT:DP:INDEL:START:STOP:AD:AF:BQ:SB.
    '
    ' aDnaSet:  A set of dna found at this position
    ' anAltList: A list of alternative alleles found thus far
    ' aStringReads:  A string of reads that have been converted from raw format and filtered
    ' aStringQualScores: A string of quality scores for the reads
    ' aStartsCount:  The number of bases that were at the start of the read
    ' aStopsCount:  The number of bases that were at the stop of the read
    ' anIndelCount:  The number of indels at this position
    ' aBaseCountsDict:  A dictionary with the number of bases of each type
    ' aQualitySumsOfBasesDict:  A dictionary with the sum of all quality scores for each type of base
    ' aPlusStrandCountsDict:  The number of bases that occurred on the plus strand
    '''
    
    # initialize the return variables
    returnString = "."
    sumAltReadSupport = 0
    
    # if we have reads at this position
    if (aNumBases > 0):
        
        #format = "GT:DP:INDEL:START:STOP:AD:AF:BQ:SB"
        
        #vcfHeader += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        #vcfHeader += "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n"
        #vcfHeader += "##FORMAT=<ID=INDEL,Number=1,Type=Integer,Description=\"Number of indels\">\n"
        #vcfHeader += "##FORMAT=<ID=START,Number=1,Type=Integer,Description=\"Number of reads starting at this position\">\n"
        #vcfHeader += "##FORMAT=<ID=STOP,Number=1,Type=Integer,Description=\"Number of reads stopping at this position\">\n"
        #vcfHeader += "##FORMAT=<ID=AD,Number=.,Type=Float,Description=\"Depth of reads supporting alleles 0/1/2/3\">\n"
        #vcfHeader += "##FORMAT=<ID=AF,Number=.,Type=Float,Description=\"Fraction of reads supporting alleles 0/1/2/3\">\n"
        #vcfHeader += "##FORMAT=<ID=BQ,Number=.,Type=Float,Description=\"Avg base quality for reads supporting alleles 0/1/2/3\">\n"
        #vcfHeader += "##FORMAT=<ID=SB,Number=.,Type=Float,Description=\"Strand Bias for reads supporting alleles 0/1/2/3\">\n"
            
        # initialize some lists
        depths = list()
        readSupports = list()
        baseQuals = list()
        strandBias = list()
        refBaseCountsDict = {}
        altBaseCountsDict = {}
        
        # for each base in the ref list and alt list
        # the order matters for the output
        for base in (aRefList + anAltList):
            # get the number of times the base occurs
            count = aBaseCountsDict[base]
            depths.append(count)

            # calculate the allele specific fraction of read support
            readSupport = round(count/float(aNumBases), 2)
            readSupports.append(readSupport)
            
            # if the base is an alt, then count it for the overall read support
            if (base in anAltList):
                sumAltReadSupport += count
                anAltCountsDict[base] += count
                # keep track of the read support fractions for all alts    
                anAltPerDict[int(readSupport*100)] += 1
                altBaseCountsDict[base] = count
            else:
                refBaseCountsDict[base] = count
            
            # calculate the allele specific avg quality and plus strand scores
            if (count > 0):
                avgBaseQuality = round(aQualitySumsOfBasesDict[base]/float(count),2)
                avgPlusStrandBias = round(aPlusStrandCountsDict[base]/float(count),2)
            else:
                avgBaseQuality = 0.0
                avgPlusStrandBias = 0.0
                
            baseQuals.append(avgBaseQuality)
            strandBias.append(avgPlusStrandBias)
            
        # get the genotype:
        #    if chrom X or Y
        #        then genotype = the ref or alt with the max read depth
        #    if there are only reads for ref
        #        then genotype = 0/0
        #    if there are only reads for alt
        #        then genotype = 1/1
        #    if there are reads for both ref and alt, then pick the ones with max counts
        #        then genotype = 0/1
        
        genotypes = list()
        refAltList = aRefList + anAltList
        singleGenotypeChroms = ["chrX", "chrY", "chrM", "chrMT", "X", "Y", "M", "MT"]
        if (aChrom in singleGenotypeChroms):
            if aBaseCountsDict:
                # get the overall max count and max base, and assign the genotype
                (maxOverallBase, maxOverallCount) = max(aBaseCountsDict.iteritems(), key=lambda x:x[1])
                genotypes.append(refAltList.index(maxOverallBase))
            else:
                genotypes.append(0)
        else:
            if refBaseCountsDict:
                # get the max ref count
                (maxRefBase, maxRefCount) = max(refBaseCountsDict.iteritems(), key=lambda x:x[1])
                maxRefIndex = refAltList.index(maxRefBase)
            else:
                maxRefCount = 0
                
            if altBaseCountsDict:
                # get the max of alt count
                (maxAltBase, maxAltCount) = max(altBaseCountsDict.iteritems(), key=lambda x:x[1])
                maxAltIndex = refAltList.index(maxAltBase)
            else:
                maxAltCount = 0
            
            # if we have some ref reads
            if (maxRefCount > 0):
                # assign the ref index
                genotypes.append(maxRefIndex)
                
                # if we have some alts
                if (maxAltCount > 0):
                    # assign the alt index
                    genotypes.append(maxAltIndex)
                else:
                    # no alts = 0/0 for genotype
                    genotypes.append(maxRefIndex)
                    depths.append(0)
                    readSupports.append(0)
                    baseQuals.append(0.0)
                    strandBias.append(0.0)
            # we have no refs
            else:
                # remove the max alt base
                del altBaseCountsDict[maxAltBase]
                
                if altBaseCountsDict:
                    # get the second max alt
                    (secondMaxAltBase, secondMaxAltCount) = max(altBaseCountsDict.iteritems(), key=lambda x:x[1])
                    secondMaxAltIndex = refAltList.index(secondMaxAltBase) 
                else:
                    secondMaxAltCount = 0
                    
                # if both counts are above zero, then add them both
                if (maxAltCount > 0 and secondMaxAltCount > 0):
                    genotypes.append(maxAltIndex)
                    genotypes.append(secondMaxAltIndex)
                # if only one is, add them twice
                elif (maxAltCount > 0):
                    genotypes.append(maxAltIndex)
                    genotypes.append(maxAltIndex)
                else: 
                    genotypes.append(secondMaxAltIndex)
                    genotypes.append(secondMaxAltIndex)
                            
        # create a list of each of the elements, then join them by colon
        outputList = ("/".join(map(str, genotypes)), str(aNumBases), str(anIndelCount), str(aStartsCount), str(aStopsCount), ",".join(map(str, depths)), ",".join(map(str, readSupports)), ",".join(map(str, baseQuals)), ",".join(map(str, strandBias)))
        returnString = ":".join(outputList)
        
    # return the string representation and overall calculations       
    return (returnString, anAltCountsDict, anAltPerDict, sumAltReadSupport)


def get_next_pileup(aGenerator):
    '''
    ' This function returns the next pileup from a generator that yields pileups.  If the user doesn't
    ' specify three .bam files, then the generator will be "None", so just return some default values.
    ' If we reach the end of a file, the generator will throw a StopIteration, so just catch it and 
    ' return some default values.  Otherwise, return the appropriate pileup information.
    '
    ' aGenerator:  A .bam mpileup generator that yields the next pileup
    '''
    if (aGenerator == None):
        return False, "", -1, "", 0, "", ""
    else:
        try:
            # get the next line
            (chr, coordinate, refBase, numReads, reads, qualScores) = aGenerator.next() 
            return True, chr, int(coordinate), refBase, int(numReads), reads, qualScores                     
        except StopIteration:
            return False, "", -1, "", 0, "", ""


def find_variants(aChr, aCoordinate, aRefBase, aNumBases, aReads, aBaseQuals, aPreviousUniqueBases, aPreviousBaseCounts, aReadDepthDict, anAltPerDict, aCoordinateWithData, aDnaSet, aRefList, anAltList, anAltCountsDict, aHasValidData, aShouldOutput, aGainModCount, aLossModCount, aGainModType, aLossModType, anInfoDict, aMinTotalNumBases, aMinAltNumBases, aBaseQual, aBaseQualsList, aSourcePrefix, anIsDebug):
    '''
    ' This function finds variants in the pileups that are provided.  The user specifies the reads and quality scores at a given coordinate
    ' and aDnaSet that the reads should be compared to.  This function first converts the samtools pileup of reads into human-readable reads
    ' and counts the number of bases on the plus and minus strands, the number of bases at the start and end of reads, and the number of indels.
    ' This function then ensures that the bases in the reads pass the minimum base quality score.  If the number of remaining bases is greater
    ' than or equal to the minimum total of bases specified by 'aMinTotalNumBases', then this function looks to see if there are any variants
    ' in the data.
    '
    ' The user specifies 'aDnaSet' that is empty when processing normal DNA.  This function automatically adds the reference base, so no
    ' pre-processing of aDnaSet is needed.  After the reference has been added, the function looks for variants in the reads.  If a base
    ' is not in aDnaSet and there are at least 'aMinAltNumBases' of them, then this function adds the variant to 'aModTypesSet'.  After 
    ' all the variants have been processed, the unique reads at this position are added to 'aDnaSet' which is used in the next steps: 
    ' looking for somatic variations and rna-editing events.
    '
    ' This function returns:
    ' bamOutputString - The concatenated form of the pileup data at this position
    ' aDnaSet - A set of 'parent' DNA {Ref for germline variants, ref + normal for somatic mutations, ref + normal + tumor for rna-editing}
    ' aHasValidData - If there was valid data at this position
    ' aShouldOutput - If there were any variants found to be output
    ' aModCount - The number of variants found
    ' aModTypesSet - The set of variants found
    '''
    
    # default outputs
    bamOutputString = "."
    sumOfBaseQuals = 0
    sumOfStrandBiases = 0
    sumOfAltReads = 0
    uniqueBases = ""
    oneAboveMinAltBasesFlag = False
    setBelowMinAltBasesFlag = False 
    baseCountsDict = collections.defaultdict(int)
    
    # convert the raw reads into human-readable reads
    (convertedReads, convertedBaseQuals, aNumBases, starts, stops, indels) = convert_raw_reads(aReads, aBaseQuals, aRefBase, anIsDebug)
    
    # keep track of the read coverage
    aReadDepthDict[aNumBases] += 1
    if (aNumBases > 0):
        aCoordinateWithData += 1

    if (anIsDebug):    
        logging.debug("After convert_raw_reads(): %s %s %s %s %s %s %s %s %s", aChr, aCoordinate, aRefBase, aNumBases, convertedReads, convertedBaseQuals, starts, stops, indels)
    
    # if we still have some bases
    if (aNumBases > 0):
               
        # filter out the bases that don't meet the minimum base quality scores
        (convertedReads, convertedBaseQuals, aNumBases, baseCountsDict, qualitySumsOfBasesDict, plusStrandCountsDict) = filter_by_base_quality(convertedReads, convertedBaseQuals, aBaseQual, anIsDebug)  
        
        if (anIsDebug):
            logging.debug("After filter_by_base_quality(): %s %s %s %s %s %s %s %s %s %s %s %s", aChr, aCoordinate, aRefBase, aNumBases, convertedReads, convertedBaseQuals, starts, stops, indels, baseCountsDict, qualitySumsOfBasesDict, plusStrandCountsDict)
        
        # if we still have some bases
        if (aNumBases > 0):
            # keep track of the base quals for later filtering
            #aBaseQualsList.append([aSourcePrefix, aChr, str(aCoordinate), aReads, aBaseQuals, convertedReads, convertedBaseQuals])
                
            # add the reference base
            aDnaSet.add(aRefBase)
                    
            # for each unique base
            for base in set(convertedReads):
                
                # keep track of every ALT base in the order that it's found
                if (base not in aDnaSet and base not in anAltList):
                    anAltList.append(base)
                
                # if we have enough total bases
                if (aNumBases >= aMinTotalNumBases):
                
                    aHasValidData = True
                    
                    # keep track of every unique base
                    uniqueBases += base
                    
                    # if there is a base that doesn't match the reference
                    if (base not in aDnaSet):            
                        
                        # if we have enough ALT bases
                        if (baseCountsDict[base] >= aMinAltNumBases):
                            oneAboveMinAltBasesFlag = True
                            aShouldOutput = True
                            aGainModCount += 1
                            
                            if (anIsDebug):
                                logging.debug("Modification found!  Base '%s' from the reads does not exist in the parent DNA set %s", base, aDnaSet)
                            
                            # if this is the first modification found, then record it in the "SS" field
                            if (len(anInfoDict["SS"]) == 0):
                                if (aGainModType == "GERM"):
                                    anInfoDict["SS"].append("1")
                                elif (aGainModType == "SOM"):
                                    anInfoDict["SS"].append("2")
                                    anInfoDict["SOMATIC"].append("True")
                                elif (aGainModType.find("EDIT") != -1):
                                    anInfoDict["SS"].append("5")
                                # other
                                else:
                                    anInfoDict["SS"].append("4")
                            
                            # it didn't match anything so far, so output all possible combinations
                            for dna in aDnaSet:
                                anInfoDict["MT"].append(aGainModType)
                                anInfoDict["MC"].append(dna + ">" + base)
                        else:
                            setBelowMinAltBasesFlag = True
                                
                    # check for LOH's
                    if ("GERM" in anInfoDict["MT"]):
                        for base in aPreviousUniqueBases:
                            # if we had enough previous germline bases, and they are now lost
                            if ((aPreviousBaseCounts[base] >= aMinAltNumBases) and (base not in convertedReads)):
                                aShouldOutput = True
                                aLossModCount += 1
                                
                                if (anIsDebug):
                                    logging.debug("Loss found!  Base '%s' from the parent DNA set no longer exists in the child reads %s", base, convertedReads)
                                
                                # remove the previous germline SS
                                if ("1" in anInfoDict["SS"]):
                                    germIndex = anInfoDict["SS"].index("1")
                                    anInfoDict["SS"].pop(germIndex)
                                
                                # set the SS to LOH
                                if (aLossModType == "LOH"):
                                    anInfoDict["SS"].append("3")
                                else:
                                    anInfoDict["SS"].append("4")
                            
                                # there could be more than one germline, so just remove them all
                                while ("GERM" in anInfoDict["MT"]):
                                    # remove the previous germline mutations       
                                    germIndex = anInfoDict["MT"].index("GERM")
                                    anInfoDict["MT"].pop(germIndex)
                                    anInfoDict["MC"].pop(germIndex)
                                    
                                # add the LOH
                                anInfoDict["MT"].append(aLossModType)
                                anInfoDict["MC"].append(aPreviousUniqueBases + ">" + uniqueBases)
                        
            # add the unique reads for the next step
            aDnaSet = aDnaSet.union(set(convertedReads))
            
            # get the summary output for the pileups at this position
            (bamOutputString, anAltCountsDict, anAltPerDict, sumOfAltReads) = format_bam_output(aChr, aRefList, anAltList, anAltCountsDict, anAltPerDict, convertedReads, convertedBaseQuals, aNumBases, starts, stops, indels, baseCountsDict, qualitySumsOfBasesDict, plusStrandCountsDict, anIsDebug)
            sumOfBaseQuals = sum(qualitySumsOfBasesDict.itervalues())
            sumOfStrandBiases = sum(plusStrandCountsDict.itervalues())
                 
    return (bamOutputString, uniqueBases, baseCountsDict, aReadDepthDict, anAltPerDict, aCoordinateWithData, aDnaSet, anAltList, anAltCountsDict, aHasValidData, aShouldOutput, (aNumBases < aMinTotalNumBases), (not oneAboveMinAltBasesFlag and setBelowMinAltBasesFlag), aGainModCount, aLossModCount, anInfoDict, aNumBases, indels, starts, stops, sumOfBaseQuals, sumOfStrandBiases, sumOfAltReads, aBaseQualsList)
        
        
def get_vcf_header(aVCFFormat, aRefId, aRefURL, aRefFilename, aRadiaVersion, aPatientId, aChrom, aParamDict, aFilenameList, aLabelList, aDescList, aPlatformList, aSourceList, anIsDebug):
    # initialize the column headers
    columnHeaders = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    
    # get the initial fields
    vcfHeader = ""
    vcfHeader += "##fileformat=" + aVCFFormat + "\n"
    vcfHeader += "##tcgaversion=1.0\n"
    vcfHeader += "##fileDate=" + datetime.date.today().strftime("%Y%m%d") + "\n"
    vcfHeader += "##center=UCSC\n"
    vcfHeader += "##source=\"RADIA pipeline " + aRadiaVersion + "\"\n"
    if (aRefFilename != None):
        vcfHeader += "##reference=file://" + aRefFilename + "\n"
    else:
        vcfHeader += "##reference=<ID=" + aRefId + ",source=\"" + aRefURL + "\">\n"
        
    vcfHeader += "##phasing=none\n"
        
    # add RADIA param info
    aParamDict["algorithm"] = "RADIA"
    aParamDict["version"] = "1.0"
    #vcfHeader += "##vcfProcessLog=<"
    vcfHeader += "##vcfGenerator=<"
    for (paramName) in sorted(aParamDict.iterkeys()):
        paramValue = aParamDict[paramName]
        if (paramValue != None):
            # don't output the defaults for files that aren't specified
            if (paramName.startswith("dnaNormal") and "DNA_NORMAL" not in aLabelList):
                continue;
            elif (paramName.startswith("rnaNormal") and "RNA_NORMAL" not in aLabelList):
                continue;
            elif (paramName.startswith("dnaTumor") and "DNA_TUMOR" not in aLabelList):
                continue;
            elif (paramName.startswith("rnaTumor") and "RNA_TUMOR" not in aLabelList):
                continue;
            else:
                vcfHeader += paramName + "=<" + str(paramValue) + ">,"
    vcfHeader = vcfHeader.rstrip(",")
    vcfHeader += ">\n"
    
    vcfHeader += "##INDIVIDUAL=" + aPatientId + "\n"
    
    # get the sample fields
    for (filename, label, description, platform, source) in izip(aFilenameList, aLabelList, aDescList, aPlatformList, aSourceList):
        # determine the SampleName from the filename
        # i_tcgaNameRegEx = re.compile("TCGA-(\\w){2}-(\\w){4}-(\\w){3}-(\\w){3}")
        matchObj = i_tcgaNameRegEx.search(filename)
        if (matchObj == None):
            logging.warning("Can't determine the 'SampleName' from the filename %s for the VCF header.", filename)
            vcfHeader += "##SAMPLE=<ID=" + label + ",SampleName=" + aPatientId + ",Individual=" + aPatientId + ",Description=\"" + description + "\",File=\"" + filename + "\",Platform=\"" + platform + "\",Source=\"" + source + "\">\n"
        else:
            vcfHeader += "##SAMPLE=<ID=" + label + ",SampleName=" + matchObj.group() + ",Individual=" + aPatientId + ",Description=\"" + description + "\",File=\"" + filename + "\",Platform=\"" + platform + "\",Source=\"" + source + "\">\n"
        columnHeaders.append(label)

    # get the info fields
    vcfHeader += "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n"
    vcfHeader += "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of unique alleles across all samples\">\n"
    vcfHeader += "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n" 
    vcfHeader += "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele frequency, for each ALT allele, in the same order as listed\">\n" 
    vcfHeader += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth for all samples\">\n"
    vcfHeader += "##INFO=<ID=INDEL,Number=1,Type=Integer,Description=\"Number of indels for all samples\">\n"
    vcfHeader += "##INFO=<ID=START,Number=1,Type=Integer,Description=\"Number of reads starting at this position across all samples\">\n"
    vcfHeader += "##INFO=<ID=STOP,Number=1,Type=Integer,Description=\"Number of reads stopping at this position across all samples\">\n"
    vcfHeader += "##INFO=<ID=BQ,Number=1,Type=Float,Description=\"Overall average base quality\">\n"
    #vcfHeader += "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Overall average mapping quality\">\n"
    vcfHeader += "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Overall strand bias\">\n"
    vcfHeader += "##INFO=<ID=FA,Number=1,Type=Float,Description=\"Overall fraction of reads supporting ALT\">\n"
    vcfHeader += "##INFO=<ID=MT,Number=.,Type=String,Description=\"Modification types at this position\">\n"
    vcfHeader += "##INFO=<ID=MC,Number=.,Type=String,Description=\"Modification base changes at this position\">\n"
    vcfHeader += "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Indicates if record is a somatic mutation\">\n"
    vcfHeader += "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type, can be SNP, INS or DEL\">\n" 
    #vcfHeader += "##INFO=<ID=DEL,Number=1,Type=Integer,Description=\"Number of deletions in all samples at this position\">\n"
    #vcfHeader += "##INFO=<ID=INS,Number=1,Type=Integer,Description=\"Number of insertions in all samples at this position\">\n"
    #vcfHeader += "##INFO=<ID=VC,Number=1,Type=String,Description=\"Somatic variant classification (Intergenic, DEL, INS)\">\n"
    #vcfHeader += "##INFO=<ID=ProtCh,Number=1,Type=String,Description=\"Protein change due to somatic variant\">\n"

    # get the filter fields
    vcfHeader += "##FILTER=<ID=noref,Description=\"Position skipped, reference=N\">\n"
    vcfHeader += "##FILTER=<ID=diffref,Description=\"Position skipped, different references in files\">\n"
    vcfHeader += "##FILTER=<ID=mbt,Description=\"Minimum total bases is less than user-specified cut-off\">\n"
    vcfHeader += "##FILTER=<ID=mba,Description=\"Minimum ALT bases is less than user-specified cut-off\">\n"
    
    # get the format fields
    # these fields are sample specific
    vcfHeader += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    vcfHeader += "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position in the sample\">\n"
    vcfHeader += "##FORMAT=<ID=INDEL,Number=1,Type=Integer,Description=\"Number of indels\">\n"
    vcfHeader += "##FORMAT=<ID=START,Number=1,Type=Integer,Description=\"Number of reads starting at this position\">\n"
    vcfHeader += "##FORMAT=<ID=STOP,Number=1,Type=Integer,Description=\"Number of reads stopping at this position\">\n"
    vcfHeader += "##FORMAT=<ID=AD,Number=.,Type=Float,Description=\"Depth of reads supporting alleles\">\n"
    vcfHeader += "##FORMAT=<ID=AF,Number=.,Type=Float,Description=\"Fraction of reads supporting alleles\">\n"
    vcfHeader += "##FORMAT=<ID=BQ,Number=.,Type=Float,Description=\"Avg base quality for reads supporting alleles\">\n"
    #vcfHeader += "##FORMAT=<ID=MQ,Number=.,Type=Float,Description=\"Avg mapping quality for reads supporting alleles\">\n"
    vcfHeader += "##FORMAT=<ID=SB,Number=.,Type=Float,Description=\"Strand Bias for reads supporting alleles\">\n"
    vcfHeader += "##FORMAT=<ID=SS,Number=1,Type=Integer,Description=\"Variant status relative to non-adjacent Normal, 0=wildtype,1=germline,2=somatic,3=LOH,4=unknown,5=rnaEditing\">\n"
    vcfHeader += "##FORMAT=<ID=SSC,Number=1,Type=Integer,Description=\"Somatic score between 0 and 255\">\n"
    
    '''
    ##INFO=<ID=fa20,Number=0,Type=Flag,Description="Fraction of ALT below 20% of reads">
    vcfHeader += "##INFO=<ID=FA,Number=0,Type=Flag,Description=\"Fraction of ALT below " + anAltFilterCutoff + "% of reads\">\n"
    '''
    
    # get the filter fields
    ##FILTER=<ID=q10,Description="Genotype Quality < 10">
    ##FILTER=<ID=ma,Description="Position in germline has 2+ support for 2+ alleles">
    ##FILTER=<ID=idl10,Description="Position is within 10 bases of an indel">
    ##FILTER=<ID=idls5,Description="Less than 5 reads supporting indel in appropriate tissue">
    ##FILTER=<ID=pbias,Description="Positional bias, all reads supporting ALT are in first or last third of read">
    ##FILTER=<ID=sbias,Description="Strand bias, majority of reads supporting ALT are on forward OR reverse strand">
    ##FILTER=<ID=mc3,Description="Greater than 3 reads of somatic allele in germline">
    
    '''
    vcfHeader += "##FILTER=<ID=dbsnp,Description=\"dbSNP membership, build " + aDbSnpVersion + "\">\n"
    vcfHeader += "##FILTER=<ID=blqual,Description=\"Position overlaps 1000 Genomes Project mapping quality blacklist\">\n"
    vcfHeader += "##FILTER=<ID=bldp,Description=\"Position overlaps 1000 Genomes Project depth blacklist\">\n"
    vcfHeader += "##FILTER=<ID=blinternal,Description=\"Position overlaps internal blacklist\">\n"
    vcfHeader += "##FILTER=<ID=repeats,Description=\"Position overlaps with repeat regions identified with RepeatMasker\">\n"
    vcfHeader += "##FILTER=<ID=minps,Description=\"Min patient support is not > " + aMinPSCutoff + "\">\n"
    vcfHeader += "##FILTER=<ID=minpspercent,Description=\"Min patient support is not > " + aMinPSPercent + "\">\n"
    vcfHeader += "##FILTER=<ID=mintotaldna,Description=\"Min  " + aMinPSCutoff + "\">\n"
    vcfHeader += "##FILTER=<ID=mintotalrna,Description=\"Min  " + aMinPSPercent + "\">\n"
    vcfHeader += "##FILTER=<ID=minaltdna,Description=\"Min  " + aMinPSCutoff + "\">\n"
    vcfHeader += "##FILTER=<ID=minaltrna,Description=\"Min  " + aMinPSPercent + "\">\n"
    vcfHeader += "##FILTER=<ID=minaltsupport,Description=\"Min  " + aMinPSCutoff + "\">\n"
    vcfHeader += "##FILTER=<ID=mincoverage,Description=\"Min  " + aMinPSPercent + "\">\n"
    '''
    
    vcfHeader += "#" + "\t".join(columnHeaders)
    return vcfHeader
            
               
def main():
    
    # command for running this on a small test case: 
    #python2.7 radia.py TCGA-AB-2995 12 --normalUseChr --tumorUseChr --rnaUseChr -n ../data/test/TCGA-AB-2995_normal.sam -t ../data/test/TCGA-AB-2995_tumor.sam -r ../data/test/TCGA-AB-2995_rna.sam
    
    # commands for running this on real data:
    #python2.7 radia.py SF4454 Y -n /inside/grotto/astrocytoma/hg18Files/A02404_81MJ7ABXX_2.bam -t /inside/grotto/astrocytoma/hg18Files/A02397_81MH7ABXX_5.bam -r /inside/grotto/astrocytoma/hg18Files/A02406_819T0ABXX_2.bam -f /inside/depot/fa/all_sequences.fasta -c /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes_sorted.tab -o /inside/grotto/astrocytoma/tripleBam/SF4454_chrY.rad
    #python2.7 radia.py TCGA-BH-A18M 10 -n /inside/depot3/cwilks/brca_from_kolossus/TCGA-BH-A18M-11A-33D-A12B-09_IlluminaGA-DNASeq_exome.bam -t /inside/depot3/cwilks/brca_from_kolossus/TCGA-BH-A18M-01A-11D-A12B-09_IlluminaGA-DNASeq_exome.bam -r /inside/depot3/cwilks/brca_from_kolossus/UNCID_368568.TCGA-BH-A18M-11A-33R-A12D-07.110714_UNC13-SN749_0082_AD0DGMABXX.3.trimmed.annotated.translated_to_genomic.bam -f /inside/depot/fa/GRCh37-lite.fa --rnaFasta=/inside/depot/fa/hg19.fasta --rnaUseChr -c /inside/home/aradenba/rnaEditing/data/hg19/hg19_chromSizes_sorted.tab -o /inside/grotto/users/aradenba/data/hg19/brca/TCGA-BH-A18M_chr10.rad
    
    #python2.7 radia.py TCGA-AG-A008 22 -n /inside/home/cwilks/bb_pipeline/runs/mutcomp_ex3_usrerun/bams/TCGA-AG-A008-10A-01W-A00K-09_IlluminaGA-DNASeq_exome.bam_HOLD_QC_PENDING.bam
    #                                   -t /inside/home/cwilks/bb_pipeline/runs/mutcomp_ex3_usrerun/bams/TCGA-AG-A008-01A-01W-A00K-09_IlluminaGA-DNASeq_exome.bam
    #                                   -f /inside/depot/fa/all_sequences.fasta 
    #                                   -c /inside/home/aradenba/rnaEditing/data/hg18/hg18_chromSizes_sorted.tab 
    #                                   -o /inside/grotto/users/aradenba/data/hg18/benchmark/radia/
    # -e hg18 -u /inside/depot/fa/all_sequences.fasta
     
    i_radiaVersion = "v1.0"
    
    # create the usage statement
    usage = "usage: python2.7 %prog id chrom [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    # add the optional parameters
    i_cmdLineParser.add_option("-b", "--batchSize", type="int", dest="batchSize", default=int(1000000), metavar="BATCH_SIZE", help="the size of the samtool selections that are loaded into memory at one time, %default by default")
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, STDOUT by default")
    i_cmdLineParser.add_option("-c", "--chromSizesFilename", dest="chromSizesFilename", metavar="CHROM_SIZES_FILE", help="the name of the file with the chromosome sizes")
    i_cmdLineParser.add_option("-f", "--fastaFilename", dest="fastaFilename", metavar="FASTA_FILE", help="the name of the fasta file that can be used on all .bams, see below for specifying individual fasta files for each .bam file")
    i_cmdLineParser.add_option("-p", "--useChrPrefix", action="store_true", default=False, dest="useChrPrefix", help="include this argument if the 'chr' prefix should be used in the samtools command for all .bams, see below for specifying the prefix for individual .bam files")
    i_cmdLineParser.add_option("-v", "--vcfFileFormat", dest="vcfFileFormat", default="VCFv4.1", help="the vcf file format that should be used, %default by default")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    i_cmdLineParser.add_option("-e", "--refId", default="ncbi-build-36.1", dest="refId", metavar="REF_ID", help="the reference Id - used in the reference VCF meta tag, %default by default")
    i_cmdLineParser.add_option("-u", "--refUrl", default="http://www.hgsc.bcm.tmc.edu/collaborations/human-reference/hsap36.1-hg18.fasta", dest="refUrl", metavar="REF_URL", help="the URL for the reference - used in the reference VCF meta tag, %default by default")
    i_cmdLineParser.add_option("-m", "--refFilename", default="/inside/depot/fa/all_sequences.fasta", dest="refFilename", metavar="REF_FILE", help="the location of the reference - used in the reference VCF meta tag, %default by default")
    i_cmdLineParser.add_option("-a", "--startCoordinate", type="int", default=int(1), dest="startCoordinate", metavar="START_COORDINATE", help="the start coordinate for testing small regions, %default by default")
    i_cmdLineParser.add_option("-z", "--stopCoordinate", type="int", default=int(0), dest="stopCoordinate", metavar="STOP_COORDINATE", help="the stop coordinate for testing small regions, %default by default")
    i_cmdLineParser.add_option("-d", "--dataSource", default="cgHub", dest="dataSource", metavar="DATA_SOURCE", help="the source of the data - used in the sample VCF meta tag, %default by default")
    i_cmdLineParser.add_option("-q", "--sequencingPlatform", default="Illumina", dest="sequencingPlatform", metavar="SEQ_PLATFORM", help="the sequencing platform - used in the sample VCF meta tag, %default by default")
    i_cmdLineParser.add_option("-s", "--statsDir", dest="statsDir", metavar="STATS_DIR", help="a stats directory where some basic stats can be output")
    
    # params for normal DNA
    i_cmdLineParser.add_option("-n", "--dnaNormalFilename", dest="dnaNormalFilename", metavar="DNA_NORMAL_FILE", help="the name of the normal DNA .bam file")
    i_cmdLineParser.add_option("", "--dnaNormalMinTotalBases", type="int", default=int(4), dest="dnaNormalMinTotalNumBases", metavar="DNA_NOR_MIN_TOTAL_BASES", help="the minimum number of overall normal DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinAltBases", type="int", default=int(2), dest="dnaNormalMinAltNumBases", metavar="DNA_NOR_MIN_ALT_BASES", help="the minimum number of alternative normal DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalBaseQual", type="int", default=int(10), dest="dnaNormalMinBaseQuality", metavar="DNA_NOR_BASE_QUAL", help="the minimum normal DNA base quality, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMapQual", type="int", default=int(10), dest="dnaNormalMinMappingQuality", metavar="DNA_NOR_MAP_QUAL", help="the minimum normal DNA mapping quality, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalUseChr", action="store_true", default=False, dest="dnaNormalUseChrPrefix", help="include this argument if the 'chr' prefix should be used in the samtools command for the normal DNA .bam file")
    i_cmdLineParser.add_option("", "--dnaNormalFasta", dest="dnaNormalFastaFilename", metavar="DNA_NOR_FASTA_FILE", help="the name of the fasta file for the normal DNA .bam file")
    i_cmdLineParser.add_option("", "--dnaNormalMitochon", default = "M", dest="dnaNormalMitochon", metavar="DNA_NOR_MITOCHON", help="the short name for the mitochondrial DNA ('M' or 'MT'), %default by default")
    
    # params for tumor DNA
    i_cmdLineParser.add_option("-t", "--dnaTumorFilename", dest="dnaTumorFilename", metavar="DNA_TUMOR_FILE", help="the name of the tumor DNA .bam file")
    i_cmdLineParser.add_option("", "--dnaTumorMinTotalBases", type="int", default=int(4), dest="dnaTumorMinTotalNumBases", metavar="DNA_TUM_MIN_TOTAL_BASES", help="the minimum number of overall tumor DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinAltBases", type="int", default=int(2), dest="dnaTumorMinAltNumBases", metavar="DNA_TUM_MIN_ALT_BASES", help="the minimum number of alternative tumor DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorBaseQual", type="int", default=int(10), dest="dnaTumorMinBaseQuality", metavar="DNA_TUM_BASE_QUAL", help="the minimum tumor DNA base quality, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMapQual", type="int", default=int(10), dest="dnaTumorMinMappingQuality", metavar="DNA_TUM_MAP_QUAL", help="the minimum tumor DNA mapping quality, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorUseChr", action="store_true", default=False, dest="dnaTumorUseChrPrefix", help="include this argument if the 'chr' prefix should be used in the samtools command for the tumor DNA .bam file")
    i_cmdLineParser.add_option("", "--dnaTumorFasta", dest="dnaTumorFastaFilename", metavar="DNA_TUM_FASTA_FILE", help="the name of the fasta file for the tumor DNA .bam file")
    i_cmdLineParser.add_option("", "--dnaTumorMitochon", default = "M", dest="dnaTumorMitochon", metavar="DNA_TUM_MITOCHON", help="the short name for the mitochondrial DNA ('M' or 'MT'), %default by default")
    
    # params for normal RNA
    i_cmdLineParser.add_option("-x", "--rnaNormalFilename", dest="rnaNormalFilename", metavar="RNA_NORMAL_FILE", help="the name of the normal RNA-Seq .bam file")
    i_cmdLineParser.add_option("", "--rnaNormalMinTotalBases", type="int", default=int(4), dest="rnaNormalMinTotalNumBases", metavar="RNA_NOR_MIN_TOTAL_BASES", help="the minimum number of overall normal RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinAltBases", type="int", default=int(2), dest="rnaNormalMinAltNumBases", metavar="RNA_NOR_MIN_ALT_BASES", help="the minimum number of alternative normal RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalBaseQual", type="int", default=int(10), dest="rnaNormalMinBaseQuality", metavar="RNA_NOR_BASE_QUAL", help="the minimum normal RNA-Seq base quality, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMapQual", type="int", default=int(10), dest="rnaNormalMinMappingQuality", metavar="RNA_NOR_MAP_QUAL", help="the minimum normal RNA-Seq mapping quality, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalUseChr", action="store_true", default=False, dest="rnaNormalUseChrPrefix", help="include this argument if the 'chr' prefix should be used in the samtools command for the normal RNA .bam file")
    i_cmdLineParser.add_option("", "--rnaNormalFasta", dest="rnaNormalFastaFilename", metavar="RNA_NOR_FASTA_FILE", help="the name of the fasta file for the normal RNA .bam file")    
    i_cmdLineParser.add_option("", "--rnaNormalMitochon", default = "M", dest="rnaNormalMitochon", metavar="RNA_NOR_MITOCHON", help="the short name for the mitochondrial RNA ('M' or 'MT'), %default by default")
    
    # params for tumor RNA
    i_cmdLineParser.add_option("-r", "--rnaTumorFilename", dest="rnaTumorFilename", metavar="RNA_TUMOR_FILE", help="the name of the tumor RNA-Seq .bam file")
    i_cmdLineParser.add_option("", "--rnaTumorMinTotalBases", type="int", default=int(4), dest="rnaTumorMinTotalNumBases", metavar="RNA_TUM_MIN_TOTAL_BASES", help="the minimum number of overall tumor RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinAltBases", type="int", default=int(2), dest="rnaTumorMinAltNumBases", metavar="RNA_TUM_MIN_ALT_BASES", help="the minimum number of alternative tumor RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorBaseQual", type="int", default=int(10), dest="rnaTumorMinBaseQuality", metavar="RNA_TUM_BASE_QUAL", help="the minimum tumor RNA-Seq base quality, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMapQual", type="int", default=int(10), dest="rnaTumorMinMappingQuality", metavar="RNA_TUM_MAP_QUAL", help="the minimum tumor RNA-Seq mapping quality, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorUseChr", action="store_true", default=False, dest="rnaTumorUseChrPrefix", help="include this argument if the 'chr' prefix should be used in the samtools command for the tumor RNA .bam file")
    i_cmdLineParser.add_option("", "--rnaTumorFasta", dest="rnaTumorFastaFilename", metavar="RNA_TUM_FASTA_FILE", help="the name of the fasta file for the tumor RNA .bam file")    
    i_cmdLineParser.add_option("", "--rnaTumorMitochon", default = "M", dest="rnaTumorMitochon", metavar="RNA_TUM_MITOCHON", help="the short name for the mitochondrial RNA ('M' or 'MT'), %default by default")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(3,35,1)
    i_argLength = len(sys.argv)
    
    #print sys.stderr, ("The number of arguments given by the user is not expected.  User-length %s is not in %s." % i_argLength, i_possibleArgLengths)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        #print "usage: python2.7 rnaDnaAnalyzer.py id chrom -b batchSize -m outputFormat -o outputFilename -c chromSizesFilename -f fastaFilename -p
        #              --dnaNormalFilename=dnaNormalFilename --dnaNormalMinTotalBases=10 --dnaNormalMinAltBases=2 --dnaNormalBaseQual=20 --dnaNormalMapQual=20 --dnaNormalUseChr --dnaNormalFasta=dnaNormalFastaFilename 
        #              --dnaTumorFilename=dnaTumorFilename --dnaTumorMinTotalBases=10 --dnaTumorMinAltBases=2 --dnaTumorBaseQual=20 --dnaTumorMapQual=20 --dnaTumorUseChr --dnaTumorFasta=dnaTumorFastaFilename 
        #              --rnaNormalFilename=rnaNormalFilename --rnaNormalMinTotalBases=10 --rnaNormalMinTotalBases=2 --rnaNormalBaseQual=20 --rnaNormalMapQual=20 --rnaNormalUseChr --rnaNormalFasta=rnaNormalFastaFilename"
        #              --rnaTumorFilename=rnaTumorFilename --rnaTumorMinTotalBases=10 --rnaTumorMinTotalBases=2 --rnaTumorBaseQual=20 --rnaTumorMapQual=20 --rnaTumorUseChr --rnaTumorFasta=rnaTumorFastaFilename"
        i_cmdLineParser.print_help()
        sys.exit(1)
    
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_cmdLineOptionsDict = vars(i_cmdLineOptions)
    i_id = str(i_cmdLineArgs[0])
    i_chrom = str(i_cmdLineArgs[1])
    
    # get the optional params with default values
    i_batchSize = i_cmdLineOptions.batchSize
    i_useChrPrefix = i_cmdLineOptions.useChrPrefix
    i_vcfFormat = i_cmdLineOptions.vcfFileFormat
    i_logLevel = i_cmdLineOptions.logLevel 
    i_startCoordinate = i_cmdLineOptions.startCoordinate
    i_stopCoordinate = i_cmdLineOptions.stopCoordinate
    i_refId = i_cmdLineOptions.refId
    i_refUrl = i_cmdLineOptions.refUrl
    i_refFilename = i_cmdLineOptions.refFilename
    i_dataSource = i_cmdLineOptions.dataSource
    i_sequencingPlatform = i_cmdLineOptions.sequencingPlatform
    
    i_dnaNormMinTotalNumBases = i_cmdLineOptions.dnaNormalMinTotalNumBases
    i_dnaNormMinAltNumBases = i_cmdLineOptions.dnaNormalMinAltNumBases
    i_dnaNormBaseQual = i_cmdLineOptions.dnaNormalMinBaseQuality
    i_dnaNormMapQual = i_cmdLineOptions.dnaNormalMinMappingQuality
    i_dnaNormUseChr = i_cmdLineOptions.dnaNormalUseChrPrefix
    i_dnaNormMitochon = i_cmdLineOptions.dnaNormalMitochon
    
    i_dnaTumMinTotalNumBases = i_cmdLineOptions.dnaTumorMinTotalNumBases
    i_dnaTumMinAltNumBases = i_cmdLineOptions.dnaTumorMinAltNumBases
    i_dnaTumBaseQual = i_cmdLineOptions.dnaTumorMinBaseQuality
    i_dnaTumMapQual = i_cmdLineOptions.dnaTumorMinMappingQuality
    i_dnaTumUseChr = i_cmdLineOptions.dnaTumorUseChrPrefix
    i_dnaTumMitochon = i_cmdLineOptions.dnaTumorMitochon
    
    i_rnaNormMinTotalNumBases = i_cmdLineOptions.rnaNormalMinTotalNumBases
    i_rnaNormMinAltNumBases = i_cmdLineOptions.rnaNormalMinAltNumBases
    i_rnaNormBaseQual = i_cmdLineOptions.rnaNormalMinBaseQuality
    i_rnaNormMapQual = i_cmdLineOptions.rnaNormalMinMappingQuality
    i_rnaNormUseChr = i_cmdLineOptions.rnaNormalUseChrPrefix
    i_rnaNormMitochon = i_cmdLineOptions.rnaNormalMitochon
    
    i_rnaTumMinTotalNumBases = i_cmdLineOptions.rnaTumorMinTotalNumBases
    i_rnaTumMinAltNumBases = i_cmdLineOptions.rnaTumorMinAltNumBases
    i_rnaTumBaseQual = i_cmdLineOptions.rnaTumorMinBaseQuality
    i_rnaTumMapQual = i_cmdLineOptions.rnaTumorMinMappingQuality
    i_rnaTumUseChr = i_cmdLineOptions.rnaTumorUseChrPrefix
    i_rnaTumMitochon = i_cmdLineOptions.rnaTumorMitochon
    
    # the user can specify that the prefix should be used on all bams with one param
    if (i_useChrPrefix):
        i_dnaNormUseChr = True
        i_dnaTumUseChr = True
        i_rnaNormUseChr = True
        i_rnaTumUseChr = True
    
    # try to get any optional parameters with no defaults    
    i_readFilenameList = []
    i_writeFilenameList = []
    i_dirList = []
    filenames = []
    labels = []
    descriptions = []
    
    i_outputFilename = None
    i_logFilename = None
    i_dnaNormalFilename = None
    i_dnaNormalGenerator = None
    i_dnaTumorFilename = None
    i_dnaTumorGenerator = None
    i_rnaNormalFilename = None
    i_rnaNormalGenerator = None
    i_rnaTumorFilename = None
    i_rnaTumorGenerator = None
    i_dnaNormalFastaFilename = None
    i_dnaTumorFastaFilename = None
    i_rnaNormalFastaFilename = None
    i_rnaTumorFastaFilename = None
    i_chromSizesFilename = None
    i_statsDir = None
    if (i_cmdLineOptions.dnaNormalFilename != None):
        i_dnaNormalFilename = str(i_cmdLineOptions.dnaNormalFilename)
        i_readFilenameList += [i_dnaNormalFilename]   
        filenames += [i_dnaNormalFilename]
        labels += ["DNA_NORMAL"]
        descriptions += ["Normal DNA Sample"]
    if (i_cmdLineOptions.rnaNormalFilename != None):
        i_rnaNormalFilename = str(i_cmdLineOptions.rnaNormalFilename)
        i_readFilenameList += [i_rnaNormalFilename] 
        filenames += [i_rnaNormalFilename] 
        labels += ["RNA_NORMAL"]
        descriptions += ["Normal RNA Sample"]
    if (i_cmdLineOptions.dnaTumorFilename != None):
        i_dnaTumorFilename = str(i_cmdLineOptions.dnaTumorFilename)
        i_readFilenameList += [i_dnaTumorFilename] 
        filenames += [i_dnaTumorFilename]  
        labels += ["DNA_TUMOR"]
        descriptions += ["Tumor DNA Sample"]
    if (i_cmdLineOptions.rnaTumorFilename != None):
        i_rnaTumorFilename = str(i_cmdLineOptions.rnaTumorFilename)
        i_readFilenameList += [i_rnaTumorFilename]  
        filenames += [i_rnaTumorFilename]
        labels += ["RNA_TUMOR"]
        descriptions += ["Tumor RNA Sample"]
    if (i_cmdLineOptions.outputFilename != None):
        i_outputFilename = str(i_cmdLineOptions.outputFilename)
        i_writeFilenameList += [i_outputFilename]
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
        i_writeFilenameList += [i_logFilename]
    if (i_cmdLineOptions.chromSizesFilename != None):
        i_chromSizesFilename = str(i_cmdLineOptions.chromSizesFilename)
        i_readFilenameList += [i_chromSizesFilename]
    if (i_cmdLineOptions.statsDir != None):
        i_statsDir = str(i_cmdLineOptions.statsDir)
        i_dirList += [i_statsDir]
        
    # if a universal fasta file is specified, then use it
    if (i_cmdLineOptions.fastaFilename != None):
        i_dnaNormalFastaFilename = str(i_cmdLineOptions.fastaFilename)
        i_dnaTumorFastaFilename = str(i_cmdLineOptions.fastaFilename)
        i_rnaNormalFastaFilename = str(i_cmdLineOptions.fastaFilename)
        i_rnaTumorFastaFilename = str(i_cmdLineOptions.fastaFilename)
        
    # if individual fasta files are specified, they over-ride the universal one
    if (i_cmdLineOptions.dnaNormalFastaFilename != None):
        i_dnaNormalFastaFilename = str(i_cmdLineOptions.dnaNormalFastaFilename)
    if (i_cmdLineOptions.dnaTumorFastaFilename != None):
        i_dnaTumorFastaFilename = str(i_cmdLineOptions.dnaTumorFastaFilename)
    if (i_cmdLineOptions.rnaNormalFastaFilename != None):
        i_rnaNormalFastaFilename = str(i_cmdLineOptions.rnaNormalFastaFilename)
    if (i_cmdLineOptions.rnaTumorFastaFilename != None):
        i_rnaTumorFastaFilename = str(i_cmdLineOptions.rnaTumorFastaFilename)
        
    if (i_dnaNormalFastaFilename != None):
        i_readFilenameList += [i_dnaNormalFastaFilename]  
    if (i_dnaTumorFastaFilename != None):
        i_readFilenameList += [i_dnaTumorFastaFilename]
    if (i_rnaNormalFastaFilename != None):
        i_readFilenameList += [i_rnaNormalFastaFilename]
    if (i_rnaTumorFastaFilename != None):
        i_readFilenameList += [i_rnaTumorFastaFilename]
    if (i_refFilename != None):
        i_readFilenameList += [i_refFilename]
        
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
        logging.debug("outputFilename=%s" % i_outputFilename)
        logging.debug("logLevel=%s" % i_logLevel)
        logging.debug("logFilename=%s" % i_logFilename)
        logging.debug("batchSize=%s" % i_batchSize)
        logging.debug("chromSizeFile=%s" % i_chromSizesFilename)
        logging.debug("vcfFormat=%s" % i_vcfFormat)
        logging.debug("i_startCoordinate=%s" % i_startCoordinate)
        logging.debug("i_stopCoordinate=%s" % i_stopCoordinate)
        logging.debug("i_refId=%s" % i_refId)
        logging.debug("i_refUrl=%s" % i_refUrl)
        logging.debug("i_refFilename=%s" % i_refFilename)
        logging.debug("i_statsDir=%s" % i_statsDir)
        
        logging.debug("dnaNormal=%s" % i_dnaNormalFilename)
        logging.debug("dna normal fasta File: %s" % i_dnaNormalFastaFilename)
        logging.debug("dna normal baseQual: %s" % i_dnaNormBaseQual)
        logging.debug("dna normal mappingQual: %s" % i_dnaNormMapQual)
        logging.debug("dna normal minTotalBases: %s" % i_dnaNormMinTotalNumBases)
        logging.debug("dna normal minAltBases: %s" % i_dnaNormMinAltNumBases)
        logging.debug("dna normal usePrefix? %s" % i_dnaNormUseChr)
        logging.debug("dna normal mitochon %s" % i_dnaNormMitochon)
    
        logging.debug("dnaTumor=%s" % i_dnaTumorFilename)
        logging.debug("dna tumor fasta File: %s" % i_dnaTumorFastaFilename)
        logging.debug("dna tumor baseQual: %s" % i_dnaTumBaseQual)
        logging.debug("dna tumor mappingQual: %s" % i_dnaTumMapQual)
        logging.debug("dna tumor minTotalBases: %s" % i_dnaTumMinTotalNumBases)
        logging.debug("dna tumor minAltBases: %s" % i_dnaTumMinAltNumBases)
        logging.debug("dna tumor usePrefix? %s" % i_dnaTumUseChr)
        logging.debug("dna tumor mitochon %s" % i_dnaTumMitochon)
        
        logging.debug("rnaNormal=%s" % i_rnaNormalFilename)
        logging.debug("rna normal fasta File: %s" % i_rnaNormalFastaFilename)
        logging.debug("rna normal baseQual: %s" % i_rnaNormBaseQual)
        logging.debug("rna normal mappingQual: %s" % i_rnaNormMapQual)
        logging.debug("rna normal minTotalBases: %s" % i_rnaNormMinTotalNumBases)
        logging.debug("rna normal minAltBases: %s" % i_rnaNormMinAltNumBases)
        logging.debug("rna normal usePrefix? %s" % i_rnaNormUseChr)
        logging.debug("rna normal mitochon %s" % i_rnaNormMitochon)
        
        logging.debug("rnaTumor=%s" % i_rnaTumorFilename)
        logging.debug("rna tumor fasta File: %s" % i_rnaTumorFastaFilename)
        logging.debug("rna tumor baseQual: %s" % i_rnaTumBaseQual)
        logging.debug("rna tumor mappingQual: %s" % i_rnaTumMapQual)
        logging.debug("rna tumor minTotalBases: %s" % i_rnaTumMinTotalNumBases)
        logging.debug("rna tumor minAltBases: %s" % i_rnaTumMinAltNumBases)
        logging.debug("rna tumor usePrefix? %s" % i_rnaTumUseChr)
        logging.debug("rna tumor mitochon %s" % i_rnaTumMitochon)
                    
    # check for any errors
    if (not rnaEditingUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)
        
    # the user must specify at least one .bam file
    if (i_dnaNormalFilename == None and i_dnaTumorFilename == None and i_rnaNormalFilename == None and i_rnaTumorFilename == None):
        #print >> sys.stderr, "You must specify at least one .bam file."
        logging.critical("You must specify at least one .bam file.")
        sys.exit(1)
            
    # get the generators that will yield the pileups
    # Note:  Use the "get_sam_data" when testing locally on a .sam file
    #        Use the "get_bam_data" when running on real .bam file data
    if ((i_dnaNormalFilename != None and i_dnaNormalFilename.endswith(".sam")) or 
        (i_dnaTumorFilename != None and i_dnaTumorFilename.endswith(".sam")) or 
        (i_rnaNormalFilename != None and i_rnaNormalFilename.endswith(".sam")) or 
        (i_rnaTumorFilename != None and i_rnaTumorFilename.endswith(".sam"))):
        
        if (i_dnaNormalFilename != None):
            i_dnaNormalGenerator = get_sam_data(i_dnaNormalFilename, i_debug)
            
        if (i_dnaTumorFilename != None):
            i_dnaTumorGenerator = get_sam_data(i_dnaTumorFilename, i_debug)
            
        if (i_rnaNormalFilename != None):
            i_rnaNormalGenerator = get_sam_data(i_rnaNormalFilename, i_debug)
        
        if (i_rnaTumorFilename != None):
            i_rnaTumorGenerator = get_sam_data(i_rnaTumorFilename, i_debug)
            
        i_startCoordinate = 150100
        i_stopCoordinate = 160000
        i_batchSize = 100
        
        if (i_debug):
            logging.debug("test i_startCoordinate: %s" % i_startCoordinate)
            logging.debug("test i_stopCoordinate: %s" % i_stopCoordinate)
            logging.debug("test i_batchSize: %s" % i_batchSize)

    else:    
        # make sure the user specified the necessary files
        if ((i_dnaNormalFilename != None and i_dnaNormalFastaFilename == None) or 
            (i_dnaTumorFilename != None and i_dnaTumorFastaFilename == None) or 
            (i_rnaNormalFilename != None and i_rnaNormalFastaFilename == None) or 
            (i_rnaTumorFilename != None and i_rnaTumorFastaFilename == None) or 
            i_chromSizesFilename == None):
            #print >> sys.stderr, "You must specify the fasta files and a file with the chrom sizes when calculating RNA and DNA variants on .bam files."
            logging.critical("You must specify the fasta files and a file with the chrom sizes when calculating DNA and RNA variants on .bam files.")
            sys.exit(1)
        
        # get the stop coordinate if it hasn't been specified
        if (i_stopCoordinate == 0):
            i_chromSizeFileHandler = open(i_chromSizesFilename, "r") 
            i_stopCoordinate = get_chrom_size(i_chrom, i_chromSizeFileHandler, i_debug)
            i_chromSizeFileHandler.close()
            
        #i_startCoordinate = 150100
        #i_stopCoordinate = 160000
        #i_batchSize = 100
        
        if (i_debug):
            logging.debug("original i_startCoordinate: %s" % i_startCoordinate)
            logging.debug("original i_stopCoordinate: %s" % i_stopCoordinate)
            logging.debug("original i_batchSize: %s" % i_batchSize)
                    
        if (i_dnaNormalFilename != None):
            # some bams/reference use "M", some use "MT"
            if (i_chrom == "M"):
                # get the DNA matched-normal generator
                i_dnaNormalGenerator = get_bam_data(i_dnaNormalFilename, i_dnaNormalFastaFilename, i_dnaNormBaseQual, i_dnaNormMapQual, i_dnaNormMitochon, i_startCoordinate, i_stopCoordinate, i_batchSize, i_dnaNormUseChr, False, i_debug)                      
            else:
                # get the DNA matched-normal generator
                i_dnaNormalGenerator = get_bam_data(i_dnaNormalFilename, i_dnaNormalFastaFilename, i_dnaNormBaseQual, i_dnaNormMapQual, i_chrom, i_startCoordinate, i_stopCoordinate, i_batchSize, i_dnaNormUseChr, False, i_debug)                      
    
        if (i_dnaTumorFilename != None):
            # some bams/reference use "M", some use "MT"
            if (i_chrom == "M"):
                # get the DNA Tumor generator
                i_dnaTumorGenerator = get_bam_data(i_dnaTumorFilename, i_dnaTumorFastaFilename, i_dnaTumBaseQual, i_dnaTumMapQual, i_dnaTumMitochon, i_startCoordinate, i_stopCoordinate, i_batchSize, i_dnaTumUseChr, False, i_debug)
            else:
                # get the DNA Tumor generator
                i_dnaTumorGenerator = get_bam_data(i_dnaTumorFilename, i_dnaTumorFastaFilename, i_dnaTumBaseQual, i_dnaTumMapQual, i_chrom, i_startCoordinate, i_stopCoordinate, i_batchSize, i_dnaTumUseChr, False, i_debug)   
        if (i_rnaNormalFilename != None):
            # some bams/reference use "M", some use "MT"
            if (i_chrom == "M"):
                # get the RNA-Seq normal generator
                i_rnaNormalGenerator = get_bam_data(i_rnaNormalFilename, i_rnaNormalFastaFilename, i_rnaNormBaseQual, i_rnaNormMapQual, i_rnaNormMitochon, i_startCoordinate, i_stopCoordinate, i_batchSize, i_rnaNormUseChr, True, i_debug)
            else:
                # get the RNA-Seq normal generator
                i_rnaNormalGenerator = get_bam_data(i_rnaNormalFilename, i_rnaNormalFastaFilename, i_rnaNormBaseQual, i_rnaNormMapQual, i_chrom, i_startCoordinate, i_stopCoordinate, i_batchSize, i_rnaNormUseChr, True, i_debug)
    
        if (i_rnaTumorFilename != None):
            # some bams/reference use "M", some use "MT"
            if (i_chrom == "M"):
                # get the RNA-Seq tumor generator
                i_rnaTumorGenerator = get_bam_data(i_rnaTumorFilename, i_rnaTumorFastaFilename, i_rnaTumBaseQual, i_rnaTumMapQual, i_rnaTumMitochon, i_startCoordinate, i_stopCoordinate, i_batchSize, i_rnaTumUseChr, True, i_debug)
            else:
                # get the RNA-Seq tumor generator
                i_rnaTumorGenerator = get_bam_data(i_rnaTumorFilename, i_rnaTumorFastaFilename, i_rnaTumBaseQual, i_rnaTumMapQual, i_chrom, i_startCoordinate, i_stopCoordinate, i_batchSize, i_rnaTumUseChr, True, i_debug)
    
    # open the output stream
    platforms = [i_sequencingPlatform] * len(filenames)
    sources = [i_dataSource] * len(filenames) 
    i_cmdLineOptionsDict["startCoordinate"] = i_startCoordinate
    i_cmdLineOptionsDict["stopCoordinate"] = i_stopCoordinate
    vcfHeader = get_vcf_header(i_vcfFormat, i_refId, i_refUrl, i_refFilename, i_radiaVersion, i_id, i_chrom, i_cmdLineOptionsDict, filenames, labels, descriptions, platforms, sources, i_debug)
    i_outputFileHandler = None
    if (i_outputFilename != None):
        i_outputFileHandler = open(i_outputFilename, "w")
        # output vcf meta information
        i_outputFileHandler.write(vcfHeader + "\n")
        #i_outputFileHandler.write("# RNA and DNA (RAD) Analyzer run on:" + str(timestamp) + "\n")
        #i_outputFileHandler.write("# User-specified Parameters:" + str(sys.argv) + "\n")
    else:
        # output vcf meta information
        print >> sys.stdout, vcfHeader
        #print >> sys.stdout, "# RNA and DNA (RAD) Analyzer run on:", str(timestamp)
        #print >> sys.stdout, "# User-specified Parameters:", sys.argv

    
    startTime = time.time()

    '''       
        # for each coordinate
            # if we have normal dna
                # compare to reference -> germline mutations
                
            # if we have normal rna-seq
                # characterize normal muts
                # identify normal rna-editing
                
            # if we have tumor dna
                # compare to reference and normal -> somatic mutations
            
            # if we have tumor rna-seq
                # characterize tumor muts
                # identify tumor rna-editing
    '''
    
    # get the first pileup from each file
    # if a file is not specified, then the "moreLines" flags will be set to false and initial values will be returned
    (moreDnaNormalLines, dnaNormalChr, dnaNormalCoordinate, dnaNormalRefBase, dnaNormalNumBases, dnaNormalReads, dnaNormalQualScores) = get_next_pileup(i_dnaNormalGenerator)
    (moreRnaNormalLines, rnaNormalChr, rnaNormalCoordinate, rnaNormalRefBase, rnaNormalNumBases, rnaNormalReads, rnaNormalQualScores) = get_next_pileup(i_rnaNormalGenerator)
    (moreDnaTumorLines, dnaTumorChr, dnaTumorCoordinate, dnaTumorRefBase, dnaTumorNumBases, dnaTumorReads, dnaTumorQualScores) = get_next_pileup(i_dnaTumorGenerator)
    (moreRnaTumorLines, rnaTumorChr, rnaTumorCoordinate, rnaTumorRefBase, rnaTumorNumBases, rnaTumorReads, rnaTumorQualScores) = get_next_pileup(i_rnaTumorGenerator)
    
    # initialize some variables
    format = "GT:DP:INDEL:START:STOP:AD:AF:BQ:SB"
    countRnaDnaCoordinateOverlap = 0
    totalGerms = 0
    totalSoms = 0
    totalNormEdits = 0
    totalTumEdits = 0
    totalNoRef = 0
    totalLohs = 0
    totalNormNotExp = 0
    totalTumNotExp = 0
    countRefMismatches = 0
    dnaSet = set()
    altList = list()
    refList = list()
    filterList = list()
    #baseQualsList = list()
    coordinateBaseQualsList = list()
    previousBaseCounts = collections.defaultdict(int)
    dnaNormalPreviousBaseCounts = collections.defaultdict(int)
    altCountsDict = collections.defaultdict(int)
    infoDict = collections.defaultdict(list)
    dnaNormalReadDPDict = collections.defaultdict(int)
    rnaNormalReadDPDict = collections.defaultdict(int)
    dnaTumorReadDPDict = collections.defaultdict(int)
    rnaTumorReadDPDict = collections.defaultdict(int)
    dnaNormalAltPercentDict = collections.defaultdict(int)
    rnaNormalAltPercentDict = collections.defaultdict(int)
    dnaTumorAltPercentDict = collections.defaultdict(int)
    rnaTumorAltPercentDict = collections.defaultdict(int)
    dnaNormalCoordinateWithData = 0
    dnaTumorCoordinateWithData = 0
    rnaNormalCoordinateWithData = 0
    rnaTumorCoordinateWithData = 0
    
    # for each coordinate that we'd like to investigate
    for currentCoordinate in range(i_startCoordinate, i_stopCoordinate):
        
        if (i_debug):
            logging.debug("currentCoordinate: %s", currentCoordinate)
            logging.debug("Initial NormalDNAData: %s %s %s %s %s %s", dnaNormalChr, dnaNormalCoordinate, dnaNormalRefBase, dnaNormalNumBases, dnaNormalReads, dnaNormalQualScores)
            logging.debug("Initial NormalRNAData: %s %s %s %s %s %s", rnaNormalChr, rnaNormalCoordinate, rnaNormalRefBase, rnaNormalNumBases, rnaNormalReads, rnaNormalQualScores)
            logging.debug("Initial TumorDNAData: %s %s %s %s %s %s", dnaTumorChr, dnaTumorCoordinate, dnaTumorRefBase, dnaTumorNumBases, dnaTumorReads, dnaTumorQualScores)
            logging.debug("Initial TumorRNAData: %s %s %s %s %s %s", rnaTumorChr, rnaTumorCoordinate, rnaTumorRefBase, rnaTumorNumBases, rnaTumorReads, rnaTumorQualScores)
            
        # empty the set of DNA for each new coordinate
        dnaSet.clear()
        altCountsDict.clear()
        infoDict.clear()
        del altList[:]
        del refList[:]
        del filterList[:]
        del coordinateBaseQualsList[:]
        
        setMinTotalBasesFlag = True
        setMinAltBasesFlag = True
        shouldOutput = False
        hasDNA = False
        hasRNA = False
        totalSamples = 0
        totalReadDepth = 0
        totalIndels = 0
        totalStarts = 0
        totalStops = 0
        totalSumBaseQual = 0
        totalSumStrandBias = 0
        totalAltReadDepth = 0
        
        # create some default output in case there are no reads for one dataset but there are for others
        #columnHeaders = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        vcfOutputList = [i_chrom, str(currentCoordinate), "."]
        dnaNormalOutputString = "."
        dnaTumorOutputString = "."
        rnaNormalOutputString = "."
        rnaTumorOutputString = "."
        previousUniqueBases = ""
        dnaNormalPreviousBases = ""
    
        # create the ref list for this coordinate
        if (dnaNormalCoordinate == currentCoordinate and dnaNormalRefBase not in refList):
            refList.append(dnaNormalRefBase)
        
        if (rnaNormalCoordinate == currentCoordinate and rnaNormalRefBase not in refList):
            refList.append(rnaNormalRefBase)
        
        if (dnaTumorCoordinate == currentCoordinate and dnaTumorRefBase not in refList):
            refList.append(dnaTumorRefBase)
        
        if (rnaTumorCoordinate == currentCoordinate and rnaTumorRefBase not in refList):
            refList.append(rnaTumorRefBase)
        
        # if we aren't debugging and we have an "N" in the ref or more than one ref, then just ignore this coordinate and move on to the next
        if (not i_debug and ("N" in refList or len(refList) > 1)):
            # if there are more lines, and the coordinate is <= the current coordinate, then get the next pileup
            if (moreDnaNormalLines and dnaNormalCoordinate <= currentCoordinate):
                (moreDnaNormalLines, dnaNormalChr, dnaNormalCoordinate, dnaNormalRefBase, dnaNormalNumBases, dnaNormalReads, dnaNormalQualScores) = get_next_pileup(i_dnaNormalGenerator)                      
            
            if (moreRnaNormalLines and rnaNormalCoordinate <= currentCoordinate):
                (moreRnaNormalLines, rnaNormalChr, rnaNormalCoordinate, rnaNormalRefBase, rnaNormalNumBases, rnaNormalReads, rnaNormalQualScores) = get_next_pileup(i_rnaNormalGenerator)                   
            
            if (moreDnaTumorLines and dnaTumorCoordinate <= currentCoordinate):
                (moreDnaTumorLines, dnaTumorChr, dnaTumorCoordinate, dnaTumorRefBase, dnaTumorNumBases, dnaTumorReads, dnaTumorQualScores) = get_next_pileup(i_dnaTumorGenerator)                      
                   
            if (moreRnaTumorLines and rnaTumorCoordinate <= currentCoordinate):
                (moreRnaTumorLines, rnaTumorChr, rnaTumorCoordinate, rnaTumorRefBase, rnaTumorNumBases, rnaTumorReads, rnaTumorQualScores) = get_next_pileup(i_rnaTumorGenerator)
            
            # continue to the next coordinate
            continue;
                        
        # if we have normal reads at the current position
        if (dnaNormalCoordinate == currentCoordinate):
            
            # specify the normal constants
            gainModType = "GERM"
            lossModType = "NOREF"
        
            # process the normal DNA
            (dnaNormalOutputString, dnaNormalPreviousBases, dnaNormalPreviousBaseCounts, dnaNormalReadDPDict, dnaNormalAltPercentDict, dnaNormalCoordinateWithData, dnaSet, altList, altCountsDict, hasDNA, shouldOutput, numTotalBasesFilter, numAltBasesFilter, totalGerms, totalNoRef, infoDict, numBases, indels, starts, stops, totalBaseQual, totalStrandBias, totalAltReadSupport, coordinateBaseQualsList) = find_variants(dnaNormalChr, dnaNormalCoordinate, dnaNormalRefBase, dnaNormalNumBases, dnaNormalReads, dnaNormalQualScores, previousUniqueBases, previousBaseCounts, dnaNormalReadDPDict, dnaNormalAltPercentDict, dnaNormalCoordinateWithData, dnaSet, refList, altList, altCountsDict, hasDNA, shouldOutput, totalGerms, totalNoRef, gainModType, lossModType, infoDict, i_dnaNormMinTotalNumBases, i_dnaNormMinAltNumBases, i_dnaNormBaseQual, coordinateBaseQualsList, "DNA_NORMAL", i_debug)
            
            #if (hasDNA):
            if (numBases > 0):
                totalSamples += 1
                totalReadDepth += numBases
                totalIndels += indels
                totalStarts += starts
                totalStops += stops
                totalSumBaseQual += totalBaseQual
                totalSumStrandBias += totalStrandBias
                totalAltReadDepth += totalAltReadSupport
                setMinTotalBasesFlag = (setMinTotalBasesFlag and numTotalBasesFilter)
                setMinAltBasesFlag = (setMinAltBasesFlag and numAltBasesFilter)
                
        # if we have normal rna-seq reads at the current position
        if (rnaNormalCoordinate == currentCoordinate):
            
            # if either a normal of tumor file is specified, we will label them as edits
            # if neither a normal file nor a tumor file is specified, we will label them as variants
            if (i_dnaNormalFilename == None and i_dnaTumorFilename == None):
                gainModType = "RNA_NOR_VAR"
            else:
                gainModType = "NOR_EDIT"
                
            lossModType = "NOTEXP"
            previousUniqueBases = ""
            
            (rnaNormalOutputString, previousUniqueBases, previousBaseCounts, rnaNormalReadDPDict, rnaNormalAltPercentDict, rnaNormalCoordinateWithData, dnaSet, altList, altCountsDict, hasRNA, shouldOutput, numTotalBasesFilter, numAltBasesFilter, totalNormEdits, totalNormNotExp, infoDict, numBases, indels, starts, stops, totalBaseQual, totalStrandBias, totalAltReadSupport, coordinateBaseQualsList) = find_variants(rnaNormalChr, rnaNormalCoordinate, rnaNormalRefBase, rnaNormalNumBases, rnaNormalReads, rnaNormalQualScores, previousUniqueBases, previousBaseCounts, rnaNormalReadDPDict, rnaNormalAltPercentDict, rnaNormalCoordinateWithData, dnaSet, refList, altList, altCountsDict, hasRNA, shouldOutput, totalNormEdits, totalNormNotExp, gainModType, lossModType, infoDict, i_rnaNormMinTotalNumBases, i_rnaNormMinAltNumBases, i_rnaNormBaseQual, coordinateBaseQualsList, "RNA_NORMAL", i_debug)    
            
            #if (hasRNA):
            if (numBases > 0):
                totalSamples += 1
                totalReadDepth += numBases
                totalIndels += indels
                totalStarts += starts
                totalStops += stops
                totalSumBaseQual += totalBaseQual
                totalSumStrandBias += totalStrandBias
                totalAltReadDepth += totalAltReadSupport
                setMinTotalBasesFlag = (setMinTotalBasesFlag and numTotalBasesFilter)
                setMinAltBasesFlag = (setMinAltBasesFlag and numAltBasesFilter)
                
        # if we have tumor reads at the current position
        if (dnaTumorCoordinate == currentCoordinate):
                
            # if a normal file is specified, we will label them as somatic mutations
            # otherwise, we will just call them variants
            if (i_dnaNormalFilename != None):
                gainModType = "SOM"
            else:
                gainModType = "DNA_TUM_VAR"
            lossModType = "LOH"
            
            # process the tumor DNA
            (dnaTumorOutputString, previousUniqueBases, previousBaseCounts, dnaTumorReadDPDict, dnaTumorAltPercentDict, dnaTumorCoordinateWithData, dnaSet, altList, altCountsDict, hasDNA, shouldOutput, numTotalBasesFilter, numAltBasesFilter, totalSoms, totalLohs, infoDict, numBases, indels, starts, stops, totalBaseQual, totalStrandBias, totalAltReadSupport, coordinateBaseQualsList) = find_variants(dnaTumorChr, dnaTumorCoordinate, dnaTumorRefBase, dnaTumorNumBases, dnaTumorReads, dnaTumorQualScores, dnaNormalPreviousBases, dnaNormalPreviousBaseCounts, dnaTumorReadDPDict, dnaTumorAltPercentDict, dnaTumorCoordinateWithData, dnaSet, refList, altList, altCountsDict, hasDNA, shouldOutput, totalSoms, totalLohs, gainModType, lossModType, infoDict, i_dnaTumMinTotalNumBases, i_dnaTumMinAltNumBases, i_dnaTumBaseQual, coordinateBaseQualsList, "DNA_TUMOR", i_debug)
            
            #if (hasDNA):
            if (numBases > 0):
                totalSamples += 1
                totalReadDepth += numBases
                totalIndels += indels
                totalStarts += starts
                totalStops += stops
                totalSumBaseQual += totalBaseQual
                totalSumStrandBias += totalStrandBias
                totalAltReadDepth += totalAltReadSupport
                setMinTotalBasesFlag = (setMinTotalBasesFlag and numTotalBasesFilter)
                setMinAltBasesFlag = (setMinAltBasesFlag and numAltBasesFilter)
            
        # if we have tumor rna-seq reads at the current position
        if (rnaTumorCoordinate == currentCoordinate):
            
            # if either a normal of tumor file is specified, we will label them as edits
            # if neither a normal file nor a tumor file is specified, we will label them as variants
            if (i_dnaNormalFilename == None and i_dnaTumorFilename == None):
                gainModType = "RNA_TUM_VAR"
            else:
                gainModType = "TUM_EDIT"
            lossModType = "NOTEXP"
            previousUniqueBases = ""
            
            (rnaTumorOutputString, previousUniqueBases, previousBaseCounts, rnaTumorReadDPDict, rnaTumorAltPercentDict, rnaTumorCoordinateWithData, dnaSet, altList, altCountsDict, hasRNA, shouldOutput, numTotalBasesFilter, numAltBasesFilter, totalTumEdits, totalTumNotExp, infoDict, numBases, indels, starts, stops, totalBaseQual, totalStrandBias, totalAltReadSupport, coordinateBaseQualsList) = find_variants(rnaTumorChr, rnaTumorCoordinate, rnaTumorRefBase, rnaTumorNumBases, rnaTumorReads, rnaTumorQualScores, previousUniqueBases, previousBaseCounts, rnaTumorReadDPDict, rnaTumorAltPercentDict, rnaTumorCoordinateWithData, dnaSet, refList, altList, altCountsDict, hasRNA, shouldOutput, totalTumEdits, totalTumNotExp, gainModType, lossModType, infoDict, i_rnaTumMinTotalNumBases, i_rnaTumMinAltNumBases, i_rnaTumBaseQual, coordinateBaseQualsList, "RNA_TUMOR", i_debug)    
            
            # if (hasRNA)
            if (numBases > 0):
                totalSamples += 1
                totalReadDepth += numBases
                totalIndels += indels
                totalStarts += starts
                totalStops += stops
                totalSumBaseQual += totalBaseQual
                totalSumStrandBias += totalStrandBias
                totalAltReadDepth += totalAltReadSupport
                setMinTotalBasesFlag = (setMinTotalBasesFlag and numTotalBasesFilter)
                setMinAltBasesFlag = (setMinAltBasesFlag and numAltBasesFilter)
            
        # count the number of ref mismatches
        if (len(refList) > 1):
            countRefMismatches += 1
            
        # if we should output
        if (shouldOutput or i_debug):
            # now that we know we should output these coordinate base quals, add them to the final list to be output
            #baseQualsList += coordinateBaseQualsList
            
            # the chrom, position, and Id columns have been filled
            #columnHeaders = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            
            # add the ref, alt, and score
            vcfOutputList.append(",".join(refList))
            vcfOutputList.append(",".join(altList))
            vcfOutputList.append("0")
            
            # add filters
            # if one of the references is "N", then set the filter and log it
            if ("N" in refList):
                filterList.append("noref")
            # if there is more than one reference, then set the filter and log it
            if (len(refList) > 1):   
                filterList.append("diffref")  
            # if there aren't enough total bases, then set the filter and log it    
            if (setMinTotalBasesFlag):
                filterList.append("mbt")    
            # if there aren't enough ALT bases, then set the filter and log it    
            if (setMinAltBasesFlag):
                filterList.append("mba")    
            # if there are no filters thus far, then pass it    
            if (len(filterList) == 0):
                filterList.append("PASS")
            
            # if we pass the basic filters, or if we are debugging    
            if (("PASS" in filterList and shouldOutput) or (i_debug)):
                vcfOutputList.append(";".join(filterList))
                
                #vcfHeader += "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n"
                #vcfHeader += "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of unique alleles across all samples\">\n"
                #vcfHeader += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth for all samples\">\n"
                #vcfHeader += "##INFO=<ID=INDEL,Number=1,Type=Integer,Description=\"Number of indels for all samples\">\n"
                #vcfHeader += "##INFO=<ID=START,Number=1,Type=Integer,Description=\"Number of reads that started at this position across all samples\">\n"
                #vcfHeader += "##INFO=<ID=STOP,Number=1,Type=Integer,Description=\"Number of reads that stopped at this position across all samples\">\n"
                #vcfHeader += "##INFO=<ID=INDEL,Number=1,Type=Integer,Description=\"Number of indels for all samples\">\n"
                #vcfHeader += "##INFO=<ID=BQ,Number=1,Type=Float,Description=\"Overall average base quality\">\n"
                #vcfHeader += "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Overall average reads on plus strand\">\n"
                #vcfHeader += "##INFO=<ID=FA,Number=1,Type=Float,Description=\"Overall fraction of reads supporting ALT\">\n"
                #vcfHeader += "##INFO=<ID=MT,Number=.,Type=String,Description=\"Modification types at this position\">\n"
                #vcfHeader += "##INFO=<ID=MC,Number=.,Type=String,Description=\"Modification base changes at this position\">\n"
                
                # add the alt counts and frequencies in the same order as the alt list 
                for base in altList:
                    infoDict["AC"].append(str(altCountsDict[base]))
                    infoDict["AF"].append(str(round(altCountsDict[base]/float(totalReadDepth),2)))
                
                # add modTypes to info
                infoDict["NS"].append(str(totalSamples))
                infoDict["AN"].append(str(len(dnaSet)))
                infoDict["DP"].append(str(totalReadDepth))
                infoDict["INDEL"].append(str(totalIndels))
                infoDict["START"].append(str(totalStarts))
                infoDict["STOP"].append(str(totalStops))
                infoDict["VT"].append("SNP")
                if (totalReadDepth > 0):
                    infoDict["BQ"].append(str(round(totalSumBaseQual/float(totalReadDepth),2)))
                    infoDict["SB"].append(str(round(totalSumStrandBias/float(totalReadDepth),2)))
                    infoDict["FA"].append(str(round(totalAltReadDepth/float(totalReadDepth),2)))
                
                # add info
                infoField = ""
                for key in sorted(infoDict.iterkeys()):
                    if ("True" in infoDict[key]):
                        infoField += key + ";"
                    else:    
                        infoField += key + "=" + ",".join(infoDict[key]) + ";"
                
                vcfOutputList.append(infoField.rstrip(";"))
                
                # add format
                vcfOutputList.append(format)
                
                # add the sample specific data
                if (i_dnaNormalFilename != None):
                    vcfOutputList.append(dnaNormalOutputString)
                if (i_rnaNormalFilename != None):
                    vcfOutputList.append(rnaNormalOutputString)
                if (i_dnaTumorFilename != None):
                    vcfOutputList.append(dnaTumorOutputString)
                if (i_rnaTumorFilename != None):
                    vcfOutputList.append(rnaTumorOutputString)
                
                # output
                if ("PASS" not in filterList):
                    logging.info("\t".join(vcfOutputList))
                elif (shouldOutput):
                    if (i_outputFileHandler != None):
                        i_outputFileHandler.write("\t".join(vcfOutputList) + "\n")
                    else:
                        print >> sys.stdout, "\t".join(vcfOutputList)
                        
                    if (i_debug):
                        logging.debug("finalOutput: %s", "\t".join(vcfOutputList))    
             
        # count coordinates when we have both DNA and RNA
        if (hasDNA and hasRNA):
            countRnaDnaCoordinateOverlap += 1
        
        # if there are more lines, and the coordinate is <= the current coordinate, then get the next pileup
        if (moreDnaNormalLines and dnaNormalCoordinate <= currentCoordinate):
            (moreDnaNormalLines, dnaNormalChr, dnaNormalCoordinate, dnaNormalRefBase, dnaNormalNumBases, dnaNormalReads, dnaNormalQualScores) = get_next_pileup(i_dnaNormalGenerator)                      
        if (moreRnaNormalLines and rnaNormalCoordinate <= currentCoordinate):
            (moreRnaNormalLines, rnaNormalChr, rnaNormalCoordinate, rnaNormalRefBase, rnaNormalNumBases, rnaNormalReads, rnaNormalQualScores) = get_next_pileup(i_rnaNormalGenerator)                   
        if (moreDnaTumorLines and dnaTumorCoordinate <= currentCoordinate):
            (moreDnaTumorLines, dnaTumorChr, dnaTumorCoordinate, dnaTumorRefBase, dnaTumorNumBases, dnaTumorReads, dnaTumorQualScores) = get_next_pileup(i_dnaTumorGenerator)                      
        if (moreRnaTumorLines and rnaTumorCoordinate <= currentCoordinate):
            (moreRnaTumorLines, rnaTumorChr, rnaTumorCoordinate, rnaTumorRefBase, rnaTumorNumBases, rnaTumorReads, rnaTumorQualScores) = get_next_pileup(i_rnaTumorGenerator)   
    
    if (i_statsDir != None):
        # output the read depth info
        i_readDPFileHandler = open(i_statsDir + i_id + "_" + i_chrom + "_readDepth.tab", "w")
        i_readDPFileHandler.write("ReadDepth\tDNANormDP\tRNANormDP\tDNATumDP\tRNATumDP\n")
        
        # get the overall max
        dnaNormalMaxDP = 0
        rnaNormalMaxDP = 0
        dnaTumorMaxDP = 0
        rnaTumorMaxDP = 0
        if dnaNormalReadDPDict:
            dnaNormalMaxDP = max(dnaNormalReadDPDict.iterkeys())
        if rnaNormalReadDPDict:
            rnaNormalMaxDP = max(rnaNormalReadDPDict.iterkeys())
        if dnaTumorReadDPDict:
            dnaTumorMaxDP = max(dnaTumorReadDPDict.iterkeys())
        if rnaTumorReadDPDict:
            rnaTumorMaxDP = max(rnaTumorReadDPDict.iterkeys())
        
        maxDP = max(dnaNormalMaxDP, dnaTumorMaxDP, rnaNormalMaxDP, rnaTumorMaxDP)
        for index in range(1,(maxDP+1)):
            i_readDPFileHandler.write(str(index) + "\t" + str(dnaNormalReadDPDict[index]) + "\t" + str(rnaNormalReadDPDict[index]) + "\t" + str(dnaTumorReadDPDict[index]) + "\t" + str(rnaTumorReadDPDict[index]) +"\n")
        i_readDPFileHandler.close()
        
        # output the alt percent info
        i_altPerFileHandler = open(i_statsDir + i_id + "_" + i_chrom + "_altPercent.tab", "w")
        i_altPerFileHandler.write("AltPercent\tDNANormAP\tRNANormAP\tDNATumAP\RNATumAP\n")
        
        for index in range(1,101):
            i_altPerFileHandler.write(str(index) + "\t" + str(dnaNormalAltPercentDict[index]) + "\t" + str(rnaNormalAltPercentDict[index]) + "\t" + str(dnaTumorAltPercentDict[index]) + "\t" + str(rnaTumorAltPercentDict[index])+"\n")
        i_altPerFileHandler.close()
        
        # output all of the base quality info for later filtering
        #i_baseQualsFileHandler = open(i_statsDir + i_id + "_" + i_chrom + "_baseQuals.tab", "w")
        #i_baseQualsFileHandler.write("Source\tChrom\tCoordinate\tOriginalReads\tOriginalBaseQuals\tReads\tBaseQuals\n")
        #for baseQualOutput in baseQualsList:
        #    i_baseQualsFileHandler.write("\t".join(baseQualOutput) + "\n")
        #i_baseQualsFileHandler.close()
        
        #i_variantCountsFileHandler = open(i_statsDir + "variantCounts.tab", "w")
        #i_variantCountsFileHandler.write("PatientId\tChrom\tGERM\tSOM\tNormEdits\tTumEdits\tLOHs\n")
        #i_variantCountsFileHandler.close()
        
        # output the variant counts
        i_variantCountsFileHandler = open(i_statsDir + "variantCounts.tab", "a")
        i_variantCountsFileHandler.write(i_id + "\t" + i_chrom + "\t" + str(totalGerms) + "\t" + str(totalSoms) + "\t" + str(totalNormEdits) + "\t" + str(totalTumEdits) + "\t" + str(totalLohs) + "\n")
        i_variantCountsFileHandler.close()
        
        #i_genStatsFileHandler = open(i_statsDir + "genStats.tab", "w")
        #i_genStatsFileHandler.write("PatientId\tChrom\tDNANormCoordinatesWithData\tDNATumCoordinatesWithData\n")
        #i_genStatsFileHandler.close()
        
        # output the coordinates with data
        i_genStatsFileHandler = open(i_statsDir + "genStats.tab", "a")
        i_genStatsFileHandler.write(i_id + "\t" + i_chrom + "\t" + str(dnaNormalCoordinateWithData) + "\t" + str(rnaNormalCoordinateWithData) + "\t" + str(dnaTumorCoordinateWithData) + "\t" + str(rnaTumorCoordinateWithData) + "\n")
        i_genStatsFileHandler.close()
    
    stopTime = time.time()                
    logging.warning("Chrom %s and Id %s: Total RNA and DNA coordinates that overlapped=%s, Total GERMs=%s, Total SOMs=%s, Total Normal EDITs=%s, Total Tumor EDITs=%s, Total LOHs=%s, Total Ref Mismatches=%s", i_chrom, i_id, countRnaDnaCoordinateOverlap, totalGerms, totalSoms, totalNormEdits, totalTumEdits, totalLohs, countRefMismatches)
    logging.warning("Chrom %s and Id %s: Total time=%s hrs, %s mins, %s secs, binSize=%s", i_chrom, i_id, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime), i_batchSize) 
           
    # close the files 
    if (i_outputFilename != None):
        i_outputFileHandler.close()
    return
 
#cProfile.run('main()', 'perfStats')
#p = pstats.Stats('perfStats')
#p.sort_stats('name')
#p.print_stats()
#p.sort_stats('cumulative').print_stats(50)

 
main()    
sys.exit(0)
