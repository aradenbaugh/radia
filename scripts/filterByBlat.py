#!/usr/bin/env python

import sys
import time
from optparse import OptionParser
from itertools import izip
import radiaUtil
import collections
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


def get_vcf_data(aVcfFile, aPassOnlyFlag, anIsDebug):
    '''
    ' This function reads from a VCF input file and uses the python generator
    ' to yield the information one line at a time.  It ignores empty lines and
    ' strips trailing \r\n characters.  This function yields all the
    ' information from the VCF file.
    '
    ' aVcfFile:         A VCF file
    ' aPassOnlyFlag:    If all calls should be processed or only those calls
    '                   that passed the filters thus far
    ' anIsDebug:         A flag for outputting debug messages to STDERR
    '''

    # open the VCF file
    fileHandler = radiaUtil.get_read_fileHandler(aVcfFile)

    for line in fileHandler:

        # if it is an empty line, then just continue
        # if is is a header line, then just continue
        if (line.isspace() or line.startswith("#")):
            continue

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        # if (anIsDebug):
        #    logging.debug("VCF: %s", line)

        # if we are only suppose to process the passed calls
        # and this call has not passed, then skip it
        if (aPassOnlyFlag and "PASS" not in line):
            continue

        yield line

    fileHandler.close()
    return


def parse_blat_output(aBlatFile, anOutputFormat, anIsDebug):
    '''
    ' This function parses the output from BLAT.  Two formats are supported:
    ' BLAST NCBI-8 and PSL.  It groups all of the information from one query
    ' sequence and uses the python generator to yield the information.  It
    ' ignores empty lines and strips trailing \r\n characters.
    '
    ' aBlatFile:         A output file from BLAT
    ' anOutputFormat:    BLAST or PSL
    ' anIsDebug:         A flag for outputting debug messages to STDERR
    '''

    # open the file
    fileHandler = radiaUtil.get_read_fileHandler(aBlatFile)
    blatHitsDict = collections.defaultdict(list)
    previousPrefix = ""

    for line in fileHandler:

        # if it is an empty line, then just continue
        # if is is a header line, then just continue
        if (line.isspace() or line.startswith("#")):
            continue

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        # if (anIsDebug):
        #    logging.debug("BLAT: %s", line)

        # split the line on the tab
        splitLine = line.split("\t")

        # get the coordinate data =
        # rnaTumor_7_55196749_HS2144:2:1108:17342:164248
        if (anOutputFormat == "PSL"):
            # the PSL output has a bunch of header lines that we want to skip
            # if the first column can't be converted into an int, then skip
            try:
                int(splitLine[0])
            except ValueError:
                continue
            blatId = splitLine[9]
        elif (anOutputFormat == "BLAST"):
            blatId = splitLine[0]

        prefix = "_".join(blatId.split("_")[0:3])
        readId = "_".join(blatId.split("_")[0:4])

        # this catches all of the matches except the first one
        if (prefix == previousPrefix):
            blatHitsDict[readId].append(line)
            '''
            if (anIsDebug):
                logging.debug("prefixes match, current=%s, prev=%s",
                              prefix, previousPrefix)
            '''
        # if the prefixes don't match and the blatHitsDict is not empty:
        # we've reached a new set of blat hits, so yield the previous ones
        elif blatHitsDict:
            '''
            if (anIsDebug):
                logging.debug("new prefix=%s, prev=%s", prefix, previousPrefix)
                logging.debug("yielding len blatHits=%s", len(blatHitsDict))
            '''
            # yield the blat hits for this prefix
            yield blatHitsDict
            # clear the blat hits dict for the next matches
            blatHitsDict.clear()
            # set the previous prefix and start filling
            # the dict with the first prefix
            previousPrefix = prefix
            blatHitsDict[readId].append(line)
            '''
            if (anIsDebug):
                logging.debug("after yield current=%s, prev=%s",
                              prefix, previousPrefix)
            '''
        # if the prefixes don't match, and the blatHitsDict is empty:
        # this is the first line of the VCF, set the previous prefix
        # and add it to the blatHitsDict
        else:
            blatHitsDict[readId].append(line)
            previousPrefix = prefix

    # this one is needed to yield the very last blatHitsDict when all
    # lines of the VCF have been processed
    yield blatHitsDict
    return


def is_blat_hit_overlap(aBlatHitChrom, aBlatHitStart, aBlatHitStop,
                        aChromList, aPosList, anIsDebug):
    '''
    if (anIsDebug):
        logging.debug("chromList=%s, posList=%s, blatChr=%s, blatStart=%s, " +
                      "blatStop=%s", aChromList, aPosList, aBlatHitChrom,
                      aBlatHitStart, aBlatHitStop)
    '''
    for (chrom, pos) in izip(aChromList, list(map(int, aPosList))):
        # vcfs are 1-based and pysam requires a 0-based coordinate
        pos = pos - 1

        # if the current chrom and position overlaps with the
        # blat hit chrom, start and stop positions
        # the blat results report the coordinates based on the strand
        # so check both directions here
        if (chrom == aBlatHitChrom or ("chr" + chrom) == aBlatHitChrom):
            if ((pos >= aBlatHitStart and pos <= aBlatHitStop) or
                (pos >= aBlatHitStop and pos <= aBlatHitStart)):
                return True
    return False


def is_valid_read_blast_format(aBlatHitsList, aChromList, aPosList,
                               anRnaIncludeSecondaryAlignmentsFlag,
                               anOrderMagnitude, anIsDebug):
    '''
    ' This method determines if the read is valid.  It compares all of the
    ' BLAT hits to determine if the read mapped to another location in the
    ' genome with a better score.
    '
    ' aBlatHitsList:
    '    The list of blat hits from the output
    ' aChromList:
    '    The DNA chrom or RNA transcript names
    ' aPosList:
    '    The DNA coordinate or RNA coordinates
    ' anRnaIncludeSecondaryAlignmentsFlag:
    '    A flag to include the RNA secondary alignments
    ' anOrderMagnitude:
    '    The order of magnitude
    ' anIsDebug:
    '    The debug flag
    '''
    # set defaults
    validRead = ""
    # the minOverlapEvalue is the minimum evalue from a blat hit
    # that overlaps with the target location
    minOverlapEValue = sys.float_info.max
    # the minOtherEvalue is the minimum evalue from any blat hit
    # that may overlap or not overlap with the target location
    minOtherEValue = sys.float_info.max
    # the maxOverlapIdentity is the maximum identity from a blat hit
    # that overlaps with the target location
    maxOverlapIdentity = sys.float_info.min
    # the maxOverlapIdentity is the maximum identity from any blat hit
    # that may overlap or not overlap with the target location
    maxOtherIdentity = sys.float_info.min
    readCount = 0

    # for each read, investigate the blat hits to see if this read is valid
    for blatHit in aBlatHitsList:
        '''
        if (anIsDebug):
            logging.debug("blatHit: %s", blatHit)
        '''
        # If the RNA reads are mapped to the transcriptome,
        # a blat hit could look something like this:
        # rnaTumor_3_79482763_readName_A_41_11_355_133
        # NM_181742       95.49   133     6       0       1       133
        # 4688    4556    1.2e-63 240.0

        # split the line on the tab
        splitLine = blatHit.split("\t")

        blatReadIdFull = splitLine[0]
        blatReadIdSplit = blatReadIdFull.split("_")
        prefix = blatReadIdSplit[0]
        # chrom = blatReadIdSplit[1]
        # coordinate = blatReadIdSplit[2]
        # readName = blatReadIdSplit[3]
        # strandedBase = blatReadIdSplit[4]
        # baseQualConverted = blatReadIdSplit[5]
        # readMapQual = blatReadIdSplit[6]
        readFlag = int(blatReadIdSplit[7])
        readLength = int(blatReadIdSplit[8])
        readLengthHalf = readLength/2

        isRNA = False
        isDNA = False
        if prefix.startswith("rna"):
            isRNA = True
        else:
            isDNA = True

        # MapSplice doesn't set the proper pair flag for RNA-Seq reads,
        # so don't require proper pairs for RNA-Seq data
        # Also, we want to allow secondary alignments for RNA that has been
        # aligned to the transcriptome to allow for multiple isoforms
        paired = int(readFlag) & 0x1
        properPaired = int(readFlag) & 0x2
        secondaryAlignment = int(readFlag) & 0x100
        # if this is a DNA read, then it needs to be properly paired
        # and the primary alignment
        # if this is an RNA read and it is the primary alignment, or
        # if this is an RNA read and secondary alignments are allowed
        if ((isDNA and paired and properPaired and not secondaryAlignment) or
            (isRNA and not secondaryAlignment) or
            (isRNA and anRnaIncludeSecondaryAlignmentsFlag)):

            readCount += 1

            blatChrom = splitLine[1]
            blatIdentity = float(splitLine[2])
            blatAlignmentLength = int(splitLine[3])
            # blatNumMismatches = int(splitLine[4])
            # blatGapOpenings = int(splitLine[5])
            # blatQueryStart = int(splitLine[6])
            # blatQueryStop = int(splitLine[7])
            blatRefStart = int(splitLine[8])
            blatRefStop = int(splitLine[9])
            blatEValue = float(splitLine[10])
            # blatBitScore = float(splitLine[11])

            # if the blat hit overlaps with the target location
            if is_blat_hit_overlap(blatChrom, blatRefStart, blatRefStop,
                                   aChromList, aPosList, anIsDebug):
                '''
                if (anIsDebug):
                    logging.debug("blat hit overlaps!!")
                    logging.debug("EVAL_BEF: blatEval=%s, " +
                                  "minOverlapEval=%s, " +
                                  "minOtherEValue=%s",
                                  blatEValue,
                                  minOverlapEValue,
                                  minOtherEValue)
                '''
                # if the blat alignment length is greater than half the read
                # if the current evalue is less than the min overlap evalue,
                if (blatAlignmentLength > readLengthHalf and
                    blatEValue < minOverlapEValue):
                    # if the minOverlapEvalue is less than the minOtherEvalue
                    # then set the minOtherEvalue to the overlapping evalue
                    if (minOverlapEValue < minOtherEValue):
                        minOtherEValue = minOverlapEValue
                    # set the minOverlapEvalue to the lower blat hit evalue
                    minOverlapEValue = blatEValue
                    validRead = blatHit
                    '''
                    if (anIsDebug and blatAlignmentLength < readLengthHalf):
                        logging.warning("Setting the minOverlapEvalue even " +
                                        "though the blatAlignLen is < half " +
                                        "the read length!! blatAlignLen=%s, " +
                                        "readLenHalf=%s", blatAlignmentLength,
                                        readLengthHalf)
                    '''
                # the current evalue is not less than the overlapping eval,
                # so check if it is less than the min other evalue
                elif (blatAlignmentLength > readLengthHalf and
                      blatEValue < minOtherEValue):
                    minOtherEValue = minOverlapEValue

                '''
                if (anIsDebug):
                    logging.debug("EVAL_AFT: blatEval=%s, " +
                                  "minOverlapEval=%s, minOtherEValue=%s",
                                  blatEValue, minOverlapEValue,
                                  minOtherEValue)

                    logging.debug("IDENT_BEF: blatAlignLen=%s, " +
                                  "readLenHalf=%s, blatIdentity=%s, " +
                                  "maxOverlapIdentity=%s, maxOtherIdentity=%s",
                                  blatAlignmentLength, readLengthHalf,
                                  blatIdentity, maxOverlapIdentity,
                                  maxOtherIdentity)
                '''

                # if the blat alignment length is greater than half the read
                # and the current blat percent identity is greater than
                # the current overlapping max identity
                if (blatAlignmentLength > readLengthHalf and
                    blatIdentity > maxOverlapIdentity):

                    # if the maxOverlapIdentity is greater than the
                    # maxOtherIdentity, then set the maxOtherIdentity to the
                    # overlapping evalue
                    if (maxOverlapIdentity > maxOtherIdentity):
                        maxOtherIdentity = maxOverlapIdentity
                    # set the maxOverlapIdentity to the greater blat identity
                    maxOverlapIdentity = blatIdentity
                    validRead = blatHit
                # the current blat maxIdentity is not greater than the
                # max overlapping identity, so check if it is greater than the
                # max other identity
                elif (blatAlignmentLength > readLengthHalf and
                      blatIdentity > maxOtherIdentity):
                    maxOtherIdentity = blatIdentity
                '''
                if (anIsDebug):
                    logging.debug("IDENT_AFT: blatAlignLen=%s, " +
                                  "readLenHalf=%s, blatIdentity=%s, " +
                                  "maxOverlapIdentity=%s, " +
                                  "maxOtherIdentity=%s",
                                  blatAlignmentLength, readLengthHalf,
                                  blatIdentity, maxOverlapIdentity,
                                  maxOtherIdentity)
                '''
            # else the blat hit does not overlap with the target location
            else:
                '''
                if (anIsDebug):
                    logging.debug("blat hit doesn't overlap position")
                    logging.debug("EVAL_BEF: blatEval=%s, " +
                                  "minOverlapEval=%s, " +
                                  "minOtherEValue=%s",
                                  blatEValue, minOverlapEValue,
                                  minOtherEValue)

                    logging.debug("IDENT_BEF: blatAlignLen=%s, " +
                                  "readLenHalf=%s, blatIdentity=%s, " +
                                  "maxOverlapIdentity=%s, " +
                                  "maxOtherIdentity=%s",
                                  blatAlignmentLength, readLengthHalf,
                                  blatIdentity, maxOverlapIdentity,
                                  maxOtherIdentity)
                '''
                # if the blat alignment length is greater than half the read
                # if the current blat evalue is less than the min other eval
                if (blatAlignmentLength > readLengthHalf and
                    blatEValue < minOtherEValue):
                    minOtherEValue = blatEValue
                # if the blat alignment length is greater than half the read,
                # and the current blat percent identity is greater than the
                # max other identity
                if (blatAlignmentLength > readLengthHalf and
                    blatIdentity > maxOtherIdentity):
                    maxOtherIdentity = blatIdentity
                '''
                if (anIsDebug):
                    logging.debug("IDENT_AFT: blatAlignLen=%s, " +
                                  "readLenHalf=%s, blatIdentity=%s, " +
                                  "maxOverlapIdentity=%s, " +
                                  "maxOtherIdentity=%s",
                                  blatAlignmentLength, readLengthHalf,
                                  blatIdentity, maxOverlapIdentity,
                                  maxOtherIdentity)

                    logging.debug("EVAL_AFT: blatEval=%s, " +
                                  "minOverlapEval=%s, " +
                                  "minOtherEValue=%s",
                                  blatEValue, minOverlapEValue,
                                  minOtherEValue)
                '''
            '''
            if (anIsDebug):
                logging.debug("BLAT_HIT_FINAL: " +
                              "minEval=%s, minOtherEValue=%s, " +
                              "maxIdentity=%s, maxOtherIdentity=%s",
                              minOverlapEValue, minOtherEValue,
                              maxOverlapIdentity, maxOtherIdentity)
            '''
    '''
    if (anIsDebug):
        logging.debug("READ_FINAL: minOverlapEval=%s, minOtherEValue=%s",
                      minOverlapEValue, minOtherEValue)
        logging.debug("READ_FINAL: maxOverlapIdentity=%s, maxOtherIdentity=%s",
                      maxOverlapIdentity, maxOtherIdentity)
    '''
    # if none of the hits were processed, then just return false
    if (readCount == 0):
        return (False, validRead)
    # if no read covered the desired location
    if (minOverlapEValue == sys.float_info.max):
        return (False, validRead)

    # if the overlapping read has the lowest e-value (within a certain order
    # of magnitude) or the highest identity, then return true, otherwise our
    # read is not mapped to the best spot in the genome, so return false
    if (anOrderMagnitude == 0):
        if (minOverlapEValue < minOtherEValue or
            maxOverlapIdentity > maxOtherIdentity):
            '''
            if (anIsDebug):
                logging.debug("found valid read!!")
            '''
            return (True, validRead)
        else:
            '''
            if (anIsDebug):
                logging.debug("no valid read!!")
            '''
            return (False, validRead)
    else:
        if ((minOverlapEValue/(1/(10.0**anOrderMagnitude))) < minOtherEValue or
            (maxOverlapIdentity > maxOtherIdentity)):
            return (True, validRead)
        else:
            return (False, validRead)

    return


def is_valid_read_psl_format(aBlatHitsList, aChromList, aPosList,
                             anRnaIncludeSecondaryAlignmentsFlag, anIsDebug):
    '''
    ' This method determines if the read is valid.  It compares all of the
    ' BLAT hits to determine if the read mapped to another location in the
    ' genome with a better score.
    '
    ' aBlatHitsList:
    '    The list of blat hits from the output
    ' aChromList:
    '    The DNA chrom or RNA transcript names
    ' aPosList:
    '    The DNA coordinate or RNA coordinates
    ' anRnaIncludeSecondaryAlignmentsFlag:
    '    A flag to include the RNA secondary alignments
    ' anIsDebug:
    '    The debug flag
    '''
    # set defaults
    validRead = ""
    # the maxOverlappingScore is the maximum score from a blat hit
    # that overlaps with the target location
    maxOverlappingScore = -sys.maxint - 1
    # the maxOtherScore is the maximum identity from any blat hit
    # that may overlap or not overlap with the target location
    maxOtherScore = -sys.maxint - 1
    readCount = 0

    # for each read, investigate the blat hits to see if this read is valid
    for blatHit in aBlatHitsList:
        if (anIsDebug):
            logging.debug("blatHit: %s", blatHit)

        # If the RNA reads are mapped to the transcriptome,
        # a blat hit could look something like this:
        # 109 9 0 0 0 0 0 0 + rnaTumor_3_79482763_readName_A_41_1_355_118
        # 118 0 118 NR_126435 4671 2294 2412 1 118, 0, 2294,

        # split the line on the tab
        splitLine = blatHit.split("\t")

        blatReadIdFull = splitLine[9]
        blatReadIdSplit = blatReadIdFull.split("_")
        prefix = blatReadIdSplit[0]
        # chrom = blatReadIdSplit[1]
        # coordinate = blatReadIdSplit[2]
        # readName = blatReadIdSplit[3]
        # strandedBase = blatReadIdSplit[4]
        # baseQualConverted = blatReadIdSplit[5]
        # readMapQual = blatReadIdSplit[6]
        readFlag = int(blatReadIdSplit[7])
        # readLength = int(blatReadIdSplit[8])

        isRNA = False
        isDNA = False
        if prefix.startswith("rna"):
            isRNA = True
        else:
            isDNA = True

        # MapSplice doesn't set the proper pair flag for RNA-Seq reads,
        # so don't require proper pairs for RNA-Seq data
        # Also, we want to allow secondary alignments for RNA that has been
        # aligned to the transcriptome to allow for multiple isoforms
        paired = int(readFlag) & 0x1
        properPaired = int(readFlag) & 0x2
        secondaryAlignment = int(readFlag) & 0x100
        # if this is a DNA read, then it needs to be properly paired
        # and the primary alignment
        # if this is an RNA read and it is the primary alignment, or
        # if this is an RNA read and secondary alignments are allowed
        if ((isDNA and paired and properPaired and not secondaryAlignment) or
            (isRNA and not secondaryAlignment) or
            (isRNA and anRnaIncludeSecondaryAlignmentsFlag)):

            readCount += 1

            matches = int(splitLine[0])
            mismatches = int(splitLine[1])
            repeatMatches = int(splitLine[2])
            # numNs = int(splitLine[3])
            qNumInserts = int(splitLine[4])
            # qNumInsertBases = int(splitLine[5])
            tNumInserts = int(splitLine[6])
            # tNumInsertBases = int(splitLine[7])
            # strand = splitLine[8]
            # blatReadIdFill = splitLine[9]
            # qSize = int(splitLine[10])
            # qStart = int(splitLine[11])
            # qEnd = int(splitLine[12])
            tName = splitLine[13]
            # tSize = int(splitLine[14])
            tStart = int(splitLine[15])
            tEnd = int(splitLine[16])
            # blockCount = int(splitLine[17])
            # blockSizes = splitLine[18].split(",")
            # qStarts = splitLine[19].split(",")
            # tStarts = splitLine[20].split(",")

            # this is the blat score on the web interface
            matchScore = (1 * (matches + repeatMatches))
            mismatchScore = (1 * mismatches)
            blatScore = matchScore - mismatchScore - qNumInserts - tNumInserts

            '''
            if (anIsDebug):
                logging.debug("blatScore=%s, matches=%s, repeatMatches=%s, " +
                              "mismatches=%s, qNumInserts=%s, tNumInserts=%s",
                              blatScore, matches, repeatMatches,
                              mismatches, qNumInserts, tNumInserts)
            '''

            # if the blat hit overlaps with the target location
            if is_blat_hit_overlap(tName, tStart, tEnd, aChromList,
                                   aPosList, anIsDebug):
                '''
                if (anIsDebug):
                    logging.debug("blat hit overlaps!!")
                    logging.debug("SCORE_BEF: blatScore=%s, " +
                                  "maxOverlappingScore=%s, " +
                                  "maxOtherScore=%s",
                                  blatScore,
                                  maxOverlappingScore,
                                  maxOtherScore)
                '''

                # if the current score is greater than the max overlap score
                if (blatScore > maxOverlappingScore):
                    # if the maxOverlappingScore is greater than the
                    # maxOtherScore, then set the maxOtherScore to the
                    # overlapping score
                    if (maxOverlappingScore > maxOtherScore):
                        maxOtherScore = maxOverlappingScore
                    # set the maxOverlappingScore to the greater blat score
                    maxOverlappingScore = blatScore
                    validRead = blatHit

                # the current blat score is not greater than the max
                # overlapping score, so check if it is greater than the
                # max other score
                elif (blatScore > maxOtherScore):
                    maxOtherScore = blatScore

                '''
                if (anIsDebug):
                    logging.debug("SCORE_AFT: blatScore=%s, " +
                                  "maxOverlappingScore=%s, " +
                                  "maxOtherScore=%s",
                                  blatScore,
                                  maxOverlappingScore,
                                  maxOtherScore)
                '''

            # else the blat hit does not overlap with the target location
            else:
                '''
                if (anIsDebug):
                    logging.debug("blat hit doesn't overlap position")
                    logging.debug("SCORE_BEF: blatScore=%s, " +
                                  "maxOverlappingScore=%s, " +
                                  "maxOtherScore=%s",
                                  blatScore,
                                  maxOverlappingScore,
                                  maxOtherScore)
                '''
                # if the current blat score is greater than
                # the max other identity
                if (blatScore > maxOtherScore):
                    maxOtherScore = blatScore

                '''
                if (anIsDebug):
                    logging.debug("SCORE_AFT: blatScore=%s, " +
                                  "maxOverlappingScore=%s, " +
                                  "maxOtherScore=%s",
                                  blatScore,
                                  maxOverlappingScore,
                                  maxOtherScore)
                '''
            '''
            if (anIsDebug):
                logging.debug("BLAT_HIT_FINAL: " +
                              "maxOverlappingScore=%s, maxOtherScore=%s",
                              maxOverlappingScore, maxOtherScore)
            '''

    # if none of the hits were processed, then just return false
    if (readCount == 0):
        return (False, validRead)
    # if no read covered the desired location
    elif (maxOverlappingScore == -sys.maxint - 1):
        return (False, validRead)

    '''
    if (anIsDebug):
        logging.debug("READ_FINAL: maxOverlappingScore=%s, maxOtherScore=%s",
                      maxOverlappingScore, maxOtherScore)
    '''

    # if another read had a better score, then return false
    # otherwise, our read is mapped to the best spot in the genome
    if (maxOtherScore > maxOverlappingScore):
        '''
        if (anIsDebug):
            logging.debug("no valid read!!")
        '''
        return (False, validRead)
    else:
        '''
        if (anIsDebug):
            logging.debug("found valid read!!")
        '''
        return (True, validRead)


def main():

    # command for running this on a small test case:
    # python filterByBlat.py TCGA-00-4454 7 ../data/test/TCGA-00-4454_EGFR.vcf
    # ../data/test/TCGA-00-4454_EGFR.fa ../data/test/TCGA-00-4454_EGFR.psl
    # --dnaNormalFilename=../data/test/TCGA-00-4454_EGFR.reads

    startTime = time.time()

    # create the usage statement
    usage = "usage: python %prog id chrom vcfFile blatOutputFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)

    # add the optional parameters
    i_cmdLineParser.add_option(
        "-c", "--allVCFCalls", action="store_false", default=True,
        dest="passedVCFCallsOnly",
        help="by default only the VCF calls that have passed all filters " +
             "thus far are processed, include this argument if all of the " +
             "VCF calls should be processed")
    i_cmdLineParser.add_option(
        "-k", "--keepPreviousFilters", action="store_true", default=False,
        dest="keepPreviousFilters",
        help="by default the previous filters are overwritten with the blat " +
             "filter, include this argument if the previous filters should " +
             "be kept")

    i_cmdLineParser.add_option(
        "-o", "--outputFilename",
        dest="outputFilename", metavar="OUTPUT_FILE", default=sys.stdout,
        help="the name of the output file, STDOUT by default")
    i_cmdLineParser.add_option(
        "-b", "--blatOutputFormat",
        dest="blatOutputFormat", metavar="OUTPUT_FORMAT", default="BLAST",
        help="the BLAT output format, BLAST by default")
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
        "", "--transcriptNameTag",
        dest="transcriptNameTag",
        help="the INFO key where the original transcript name can be found")
    i_cmdLineParser.add_option(
        "", "--transcriptCoordinateTag",
        dest="transcriptCoordinateTag",
        help="the INFO key where the original transcript" +
             "coordinate can be found")
    i_cmdLineParser.add_option(
        "", "--transcriptStrandTag",
        dest="transcriptStrandTag",
        help="the INFO key where the original transcript strand can be found")
    i_cmdLineParser.add_option(
        "", "--rnaIncludeSecondaryAlignments",
        action="store_true", default=False,
        dest="rnaIncludeSecondaryAlignments",
        help="if you align the RNA to transcript isoforms, then you may " +
             "want to include RNA secondary alignments in the pileup")

    i_cmdLineParser.add_option(
        "-n", "--blatDnaNormalReads", action="store_true", default=False,
        dest="blatDnaNormalReads",
        help="include this argument if the normal DNA reads " +
             "should be processed")
    i_cmdLineParser.add_option(
        "-x", "--blatRnaNormalReads", action="store_true", default=False,
        dest="blatRnaNormalReads",
        help="include this argument if the normal RNA reads " +
             "should be processed")
    i_cmdLineParser.add_option(
        "-t", "--blatDnaTumorReads", action="store_true", default=False,
        dest="blatDnaTumorReads",
        help="include this argument if the tumor DNA reads " +
             "should be processed")
    i_cmdLineParser.add_option(
        "-r", "--blatRnaTumorReads", action="store_true", default=False,
        dest="blatRnaTumorReads",
        help="include this argument if the tumor RNA reads " +
             "should be processed")

    i_cmdLineParser.add_option(
        "-d", "--minReadDepth", type="int", default=int(4),
        dest="minReadDepth", metavar="MIN_READ_DP",
        help="the minimum number of valid reads that are necessary, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "-p", "--minReadPercent", type="float", default=float(0.10),
        dest="minReadPercent", metavar="MIN_READ_PCT",
        help="the minimum percentage of valid reads that are necessary, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "-m", "--minOrderMagnitude", type="int", default=float(0),
        dest="minOrderMagnitude", metavar="MIN_ORDER_MAGNITUDE",
        help="the minimum order of magnitude difference between the blat " +
             "hit at the query position vs. the next best blat hit in order " +
             "for the read to be valid, %default by default")

    '''
    i_cmdLineParser.add_option(
        "-e", "--minEValue", type="float", default=float(10e-6),
        dest="minEValue", metavar="MIN_EVALUE",
        help="the minimum e-value needed for a blat hit to be significant, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "-u", "--maxIdentity", type="float", default=float(0.95),
        dest="maxIdentity", metavar="MAX_IDENTITY",
        help="the maximum match length adjusted identity for a blat hit to " +
             "be significant, %default by default")
    i_cmdLineParser.add_option(
        "-l", "--minIdentity", type="float", default=float(0.5),
        dest="minIdentity", metavar="MIN_IDENTITY",
        help="the minimum match length adjusted identity for a blat hit to " +
             "be significant, %default by default")
    '''

    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(5, 27, 1)
    i_argLength = len(sys.argv)

    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)

    # get the required parameters
    (cmdLineOpts, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_id = str(i_cmdLineArgs[0])
    i_vcfFilename = str(i_cmdLineArgs[1])
    i_blatOutputFilename = str(i_cmdLineArgs[2])

    # get the optional params with default values
    i_passedVCFCallsOnlyFlag = cmdLineOpts.passedVCFCallsOnly
    i_keepPreviousFiltersFlag = cmdLineOpts.keepPreviousFilters
    i_blatOutputFormat = cmdLineOpts.blatOutputFormat
    i_logLevel = cmdLineOpts.logLevel
    i_rnaIncludeSecondaryAlignments = cmdLineOpts.rnaIncludeSecondaryAlignments
    i_minReadDepth = cmdLineOpts.minReadDepth
    i_minReadPercent = cmdLineOpts.minReadPercent
    i_minOrderMagnitude = cmdLineOpts.minOrderMagnitude
    # i_minEValue = cmdLineOpts.minEValue
    # i_maxIdentity = cmdLineOpts.maxIdentity
    # i_minIdentity = cmdLineOpts.minIdentity

    i_blatDnaNormalReads = cmdLineOpts.blatDnaNormalReads
    i_blatDnaTumorReads = cmdLineOpts.blatDnaTumorReads
    i_blatRnaNormalReads = cmdLineOpts.blatRnaNormalReads
    i_blatRnaTumorReads = cmdLineOpts.blatRnaTumorReads

    # try to get any optional parameters with no defaults
    i_outputFilename = None
    i_logFilename = None
    i_transcriptNameTag = None
    i_transcriptCoordinateTag = None
    i_transcriptStrandTag = None
    if (cmdLineOpts.outputFilename is not None):
        i_outputFilename = cmdLineOpts.outputFilename
    if (cmdLineOpts.logFilename is not None):
        i_logFilename = cmdLineOpts.logFilename
    if (cmdLineOpts.transcriptNameTag is not None):
        i_transcriptNameTag = cmdLineOpts.transcriptNameTag
    if (cmdLineOpts.transcriptCoordinateTag is not None):
        i_transcriptCoordinateTag = cmdLineOpts.transcriptCoordinateTag
    if (cmdLineOpts.transcriptStrandTag is not None):
        i_transcriptStrandTag = cmdLineOpts.transcriptStrandTag

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
        logging.debug("vcfFilename=%s", i_vcfFilename)
        logging.debug("blatOutputFilename=%s", i_blatOutputFilename)
        logging.debug("passedCallsOnly? %s", i_passedVCFCallsOnlyFlag)
        logging.debug("keepPreviousFiltersFlag? %s", i_keepPreviousFiltersFlag)
        logging.debug("blatOutputFormat=%s", i_blatOutputFormat)

        logging.debug("transcriptNameTag %s", i_transcriptNameTag)
        logging.debug("transcriptCoordinateTag %s", i_transcriptCoordinateTag)
        logging.debug("transcriptStrandTag %s", i_transcriptStrandTag)
        logging.debug("rnaInclSecAlign=%s" % i_rnaIncludeSecondaryAlignments)

        logging.debug("blatDnaNormal? %s", i_blatDnaNormalReads)
        logging.debug("blatDnaTumor? %s", i_blatDnaTumorReads)
        logging.debug("blatRnaNormal? %s", i_blatRnaNormalReads)
        logging.debug("blatRnaTumor? %s", i_blatRnaTumorReads)

        logging.debug("minReadDepth=%s", i_minReadDepth)
        logging.debug("minReadPercent=%s", i_minReadPercent)
        logging.debug("minOrderMagnitude=%s", i_minOrderMagnitude)

    # check for any errors
    i_writeFilenameList = []
    if (cmdLineOpts.outputFilename is not sys.stdout):
        i_writeFilenameList = [i_outputFilename]
    if (i_logFilename is not None):
        i_writeFilenameList = [i_logFilename]

    i_readFilenameList = [i_vcfFilename, i_blatOutputFilename]

    if (not radiaUtil.check_for_argv_errors(None,
                                            i_readFilenameList,
                                            i_writeFilenameList)):
        sys.exit(1)

    # open the output stream
    if i_outputFilename is not sys.stdout:
        i_outputFileHandler = radiaUtil.get_write_fileHandler(i_outputFilename)
    else:
        i_outputFileHandler = i_outputFilename

    # get the VCF generator
    i_vcfGenerator = get_vcf_data(i_vcfFilename,
                                  i_passedVCFCallsOnlyFlag,
                                  i_debug)

    # get the blat hits generator
    i_blatGenerator = parse_blat_output(i_blatOutputFilename,
                                        i_blatOutputFormat,
                                        i_debug)

    for (vcfLine, blatHitsDict) in izip(i_vcfGenerator, i_blatGenerator):

        if (i_debug):
            logging.debug("VCF Line=%s", vcfLine)
            logging.debug("Len Blat Hits=%s", len(blatHitsDict))

        # parse the VCF line
        splitLine = vcfLine.split("\t")

        # the coordinate is the second element
        vcfChr = splitLine[0]
        vcfStopCoordinate = int(splitLine[1])
        vcfIds = splitLine[2]
        vcfRef = splitLine[3]
        vcfAlts = splitLine[4]
        vcfScore = splitLine[5]
        vcfFilterSet = set(splitLine[6].split(";"))
        vcfInfoList = splitLine[7].split(";")
        vcfInfoDict = collections.defaultdict(list)
        for info in vcfInfoList:
            keyValueList = info.split("=")
            # some keys are just singular without a value (e.g. DB, etc.)
            if (len(keyValueList) == 1):
                vcfInfoDict[keyValueList[0]] = ["True"]
            else:
                # the value can be a comma separated list
                vcfInfoDict[keyValueList[0]] = keyValueList[1].split(",")
        vcfRestOfLine = "\t".join(splitLine[8:])

        modTypes = vcfInfoDict["MT"]
        modTypeFilters = dict()
        atLeastOnePass = False
        for modType in modTypes:

            blatOverallReadDepth = 0
            numValidReads = 0

            prefix = ""
            if (modType == "GERM" and i_blatDnaNormalReads):
                prefix = "dnaNormal"
            elif (modType == "NOR_EDIT" and i_blatRnaNormalReads):
                prefix = "rnaNormal"
            elif (modType == "SOM" and i_blatDnaTumorReads):
                prefix = "dnaTumor"
            elif ((modType == "SOM" or modType == "TUM_EDIT") and
                  i_blatRnaTumorReads):
                prefix = "rnaTumor"

            # get the expected prefix
            vcfKey = "_".join([prefix, vcfChr, str(vcfStopCoordinate)])

            # for each read, investigate the blat
            # hits to see if this read is valid
            for (readId, blatHitList) in blatHitsDict.iteritems():
                if (i_debug):
                    logging.debug("num of blat hits for read %s=%s",
                                  readId, len(blatHitList))

                # if the readId does not start with the vcfKey,
                # then something is wrong. the VCF and blat hits
                # need to be in sync...
                if (not readId.startswith(vcfKey)):
                    logging.error("The blat query seems to be out of sync " +
                                    "with the blat hits.")
                    logging.error("VCF Line=%s", vcfLine)
                    logging.error("readId=%s, blatHitsDict=%s",
                                    readId, blatHitsDict[readId][1])
                    sys.exit(1)

                blatOverallReadDepth += 1

                # find out if the read is valid or if it
                # maps to other places in the genome
                if (i_blatOutputFormat == "PSL"):
                    # if we should process the transcripts
                    if ((i_transcriptNameTag is not None) and
                        (i_transcriptNameTag in vcfInfoDict)):
                        (isValidRead, validRead) = is_valid_read_psl_format(
                                        blatHitList,
                                        vcfInfoDict[i_transcriptNameTag],
                                        vcfInfoDict[i_transcriptCoordinateTag],
                                        i_rnaIncludeSecondaryAlignments,
                                        i_debug)
                    else:
                        (isValidRead, validRead) = is_valid_read_psl_format(
                                        blatHitList,
                                        [vcfChr],
                                        [vcfStopCoordinate],
                                        i_rnaIncludeSecondaryAlignments,
                                        i_debug)

                elif (i_blatOutputFormat == "BLAST"):
                    # if we should process the transcripts
                    if ((i_transcriptNameTag is not None) and
                        (i_transcriptNameTag in vcfInfoDict)):
                        (isValidRead, validRead) = is_valid_read_blast_format(
                                        blatHitList,
                                        vcfInfoDict[i_transcriptNameTag],
                                        vcfInfoDict[i_transcriptCoordinateTag],
                                        i_rnaIncludeSecondaryAlignments,
                                        i_minOrderMagnitude,
                                        i_debug)
                    else:
                        (isValidRead, validRead) = is_valid_read_blast_format(
                                        blatHitList,
                                        [vcfChr],
                                        [vcfStopCoordinate],
                                        i_rnaIncludeSecondaryAlignments,
                                        i_minOrderMagnitude,
                                        i_debug)

                # if we have only one valid blat hit, then the read doesn't
                # map to other places in the genome very well, so let's use it
                if (isValidRead):
                    numValidReads += 1

                    if (i_debug):
                        logging.debug("ValidRead: %s", validRead)
                elif (i_debug):
                    logging.debug("not a valid read")

            if (blatOverallReadDepth > 0):
                tmpAltPct = numValidReads/float(blatOverallReadDepth)
                altPercent = round(tmpAltPct, 2)
            else:
                altPercent = 0.0

            if (numValidReads < i_minReadDepth or
                altPercent < i_minReadPercent):
                modTypeFilters[modType] = "blat"
            else:
                modTypeFilters[modType] = "PASS"
                atLeastOnePass = True

            if (i_debug):
                logging.debug("blatOverallReadDepth=%s, numValidReads=%s, " +
                              "altPercent=%s", str(blatOverallReadDepth),
                              str(numValidReads), str(altPercent))
                logging.debug("modType=%s, passed? %s", modType,
                              modTypeFilters[modType])
                logging.debug("blatFilter originalDepth=%s, validBlatDepth=%s",
                              str(blatOverallReadDepth), str(numValidReads))

        # make a copy of the list to manipulate
        modTypesTmpList = list(modTypes)
        modChanges = vcfInfoDict["MC"]
        # if at least one passed, then remove the ones that didn't
        for (modType, modChange) in izip(modTypes, modChanges):
            # if at least one passed, then remove the ones that didn't
            if (atLeastOnePass):
                if (modTypeFilters[modType] == "blat"):
                    modTypesTmpList.remove(modType)
                    modChanges.remove(modChange)

        # set the modTypes and modChanges
        vcfInfoDict["MT"] = modTypesTmpList
        vcfInfoDict["MC"] = modChanges

        # if at least one passed, then set pass
        if (atLeastOnePass):
            vcfFilterSet = ["PASS"]
        else:
            # if the user wants to keep the previous filters
            if (i_keepPreviousFiltersFlag):
                # if the call previous passed, then just set blat
                if (len(vcfFilterSet) == 1 and "PASS" in vcfFilterSet):
                    vcfFilterSet = ["blat"]
                # otherwise, add it to the previous filters
                else:
                    vcfFilterSet.add("blat")
            # otherwise, just set the blat filter
            else:
                vcfFilterSet = ["blat"]

            # update the mod filters
            modTypes = vcfInfoDict["MT"]
            modChanges = vcfInfoDict["MC"]
            origins = vcfInfoDict["ORIGIN"]
            if vcfInfoDict["MF"] is None:
                modFilters = []
            else:
                modFilters = vcfInfoDict["MF"]
            if vcfInfoDict["MFT"] is None:
                modFilterTypes = []
            else:
                modFilterTypes = vcfInfoDict["MFT"]

            for origin in origins:
                for (modType, modChange) in izip(modTypes, modChanges):
                    modFilterTypes.append("_".join([origin,
                                                    modType,
                                                    modChange]))
                    modFilters.append("_".join(vcfFilterSet))

            vcfInfoDict["MF"] = modFilters
            vcfInfoDict["MFT"] = modFilterTypes

        output = [vcfChr, str(vcfStopCoordinate), vcfIds, vcfRef,
                  vcfAlts, vcfScore, ";".join(vcfFilterSet)]

        # add the modified info dict
        infoField = ""
        for key in sorted(vcfInfoDict.iterkeys()):
            if (len(vcfInfoDict[key]) == 0):
                continue
            elif ("True" in vcfInfoDict[key]):
                infoField += key + ";"
            else:
                infoField += key + "=" + ",".join(vcfInfoDict[key]) + ";"

        output.append(infoField.rstrip(";"))
        output.append(vcfRestOfLine)

        i_outputFileHandler.write("\t".join(output) + "\n")

    stopTime = time.time()
    logging.info("filterByBlat.py for Id %s: Total time=%s hrs, %s mins, " +
                 "%s secs", i_id, ((stopTime-startTime)/(3600)),
                 ((stopTime-startTime)/60), (stopTime-startTime))

    # close the files
    if (i_outputFilename is not sys.stdout):
        i_outputFileHandler.close()

    return


main()
sys.exit(0)
