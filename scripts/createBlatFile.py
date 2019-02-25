#!/usr/bin/env python

__requires__ = ['pysam>=0.8.1']

import pysam
import sys
import time
from optparse import OptionParser
import radiaUtil
import collections
import logging
from itertools import izip
import os


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

i_reverseCompDict = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def get_vcf_data(aVcfFile, aHeaderFile, aPassOnlyFlag, anIsDebug):
    '''
    ' This function reads from a VCF input file and uses the python generator
    ' to yield the information one line at a time.  It ignores empty lines and
    ' strips trailing \r\n characters.  This function yields all the
    ' information from the VCF file.
    '
    ' aVcfFile:       A VCF file
    ' aHeaderFile:    A VCF file to get the header lines
    ' aPassOnlyFlag:  If all calls should be processed or only those calls
    '                 that passed the filters thus far
    ' anIsDebug:      A flag for outputting debug messages to STDERR
    '''

    # open the header file
    fileHandler = radiaUtil.get_read_fileHandler(aHeaderFile)

    for line in fileHandler:

        # if it is an empty line, then just continue
        if (line.isspace()):
            continue

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        # if (anIsDebug):
        #    logging.debug("VCF Header: %s", line)

        # if we find the column headers
        if ("#CHROM" in line):
            columnsLine = line.lstrip("#")
            columnsList = columnsLine.split("\t")
            columnsList = columnsList[9:len(columnsList)]
            continue

        # if we find the vcfGenerator line, then create the dict of params
        elif ("vcfGenerator" in line):
            # generatorLine = line.rstrip(">")
            # generatorLine = generatorLine.lstrip("##vcfGenerator=<")
            generatorLine = line[0:(len(line)-1)]
            # print "generatorLine: %s", generatorLine
            generatorLine = generatorLine[16:len(generatorLine)]
            # print "generatorLine: %s", generatorLine
            generatorParamsList = generatorLine.split(",")
            generatorParamsDict = {}

            # create a dictionary of existing params
            for param in generatorParamsList:
                (key, value) = param.split("=")
                value = value.rstrip(">")
                value = value.lstrip("<")
                generatorParamsDict[key] = value
            continue

        # if we are done with the header, then stop
        elif (not line.startswith("#")):
            break

    fileHandler.close()

    # open the VCF file
    fileHandler = radiaUtil.get_read_fileHandler(aVcfFile)

    for line in fileHandler:

        # if it is an empty line, then just continue
        if (line.isspace()):
            continue

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        # if (anIsDebug):
        #    logging.debug("VCF: %s", line)

        # if we find the column headers
        if ("#CHROM" in line):
            columnsLine = line.lstrip("#")
            columnsList = columnsLine.split("\t")
            columnsList = columnsList[9:len(columnsList)]
            continue

        # if we find the vcfGenerator line, then create the dict of params
        elif ("vcfGenerator" in line):
            # generatorLine = line.rstrip(">")
            # generatorLine = generatorLine.lstrip("##vcfGenerator=<")
            generatorLine = line[0:(len(line)-1)]
            # print "generatorLine: %s", generatorLine
            generatorLine = generatorLine[16:len(generatorLine)]
            # print "generatorLine: %s", generatorLine
            generatorParamsList = generatorLine.split(",")
            generatorParamsDict = {}

            # create a dictionary of existing params
            for param in generatorParamsList:
                (key, value) = param.split("=")
                value = value.rstrip(">")
                value = value.lstrip("<")
                generatorParamsDict[key] = value
            continue

        # these are header lines, so just continue
        elif (line.startswith("#")):
            continue

        # if we are only suppose to process the passed calls
        # and this call has not passed, then skip it
        elif (aPassOnlyFlag and "PASS" not in line):
            continue

        # split the line on the tab
        splitLine = line.split("\t")

        # the coordinate is the second element
        chrom = splitLine[0]
        stopCoordinate = int(splitLine[1])
        idList = splitLine[2].split(";")
        refList = splitLine[3].split(",")
        altList = splitLine[4].split(",")
        score = float(splitLine[5])
        filterSet = set(splitLine[6].split(";"))
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

        # yield all the information about the current coordinate
        yield (chrom, stopCoordinate, idList, refList, altList,
               score, filterSet, infoDict, "\t".join(splitLine[8:]),
               generatorParamsDict)
    fileHandler.close()
    return


def rev_comp_nucleotide(aNucleotide):
    '''
    ' This function returns the reverse complement of the parameter aNucleotide
    '
    ' aNucleotide:    A nucleotide to reverse complement
    '''
    if aNucleotide in i_reverseCompDict:
        return i_reverseCompDict[aNucleotide]
    else:
        logging.error("Trying to reverse complement an unknown nucleotide: %s",
                      aNucleotide)
        sys.exit(1)

    return None


def group_reads_by_name(aChromList, aPosList, aTranscriptStrandList,
                        anAlleleSet, aBamFile, aFastaFile, aBamOrigin,
                        anRnaIncludeSecondaryAlignmentsFlag, anIsDebug):

        readsDict = collections.defaultdict(list)

        # loop through all of the transcripts
        for (chrom, pos, strand) in izip(aChromList,
                                         list(map(int, aPosList)),
                                         aTranscriptStrandList):

            # vcfs are 1-based and pysam requires a 0-based coordinate
            pos = pos - 1

            # get the reference base from the fasta
            refBase = aFastaFile.fetch(chrom, pos, pos+1).upper()

            if (anIsDebug):
                logging.debug("getting pileups for chrom=%s, pos=%s, " +
                              "strand=%s", chrom, pos, strand)

            # get the pileups
            for pileupColumn in aBamFile.pileup(chrom, pos, pos+1,
                                                stepper="nofilter"):

                # move through pileup until at the correct position
                if pileupColumn.pos < pos:
                    '''
                    if (anIsDebug):
                        logging.debug("continue pileupColumn.pos=%s, pos=%s",
                                      pileupColumn.pos, pos)
                    '''
                    continue
                if pileupColumn.pos > pos:
                    '''
                    if (anIsDebug):
                        logging.debug("break out pileupColumn.pos=%s, pos=%s",
                                      pileupColumn.pos, pos)
                    '''
                    break

                totalReads = 0
                keptReads = 0

                # loop through the reads and create a dictionary
                for pileupRead in pileupColumn.pileups:

                    totalReads += 1
                    alignedRead = pileupRead.alignment

                    # skip over reads with problems
                    if (alignedRead.is_qcfail or
                        alignedRead.is_unmapped or
                        alignedRead.is_duplicate or
                        pileupRead.is_del or
                        pileupRead.is_refskip):
                        '''
                        if (anIsDebug):
                            logging.debug("read is unmapped, duplicate, " +
                                          "qcfail, del, or refskip at %s:%s",
                                          chrom, pos)
                        '''
                        continue
                    # no longer needed?? indel length for the position
                    # following the current pileup site
                    if pileupRead.indel != 0:
                        '''
                        if (anIsDebug):
                            logging.debug("base is an indel at %s:%s",
                                          chrom, pos)
                        '''
                        continue
                    # skip over secondary mappings for RNA if the
                    # param is not set to include them
                    if (aBamOrigin == "RNA" and
                        alignedRead.is_secondary and
                        not anRnaIncludeSecondaryAlignmentsFlag):
                        '''
                        if (anIsDebug):
                            logging.debug("read is secondary alignment but " +
                                          "flag to include secondary " +
                                          "alignments for RNA is false %s:%s",
                                          chrom, pos)
                        '''
                        continue

                    # keep a dict of all reads, using the readName as the key
                    # due to the inclusion of secondary alignments for RNA-Seq,
                    # there could be more than 2 reads that are paired
                    oneReadDict = {}
                    keptReads += 1
                    qPos = pileupRead.query_position
                    oneReadDict["alignedRead"] = alignedRead
                    oneReadDict["pileupRead"] = pileupRead
                    # old names vs. new names
                    # qname            queryName
                    # rname            reference_name
                    # pos              reference_start
                    # mapq             mapping_quality
                    # qlen             query_length
                    # pnext or mpos    next_reference_start
                    # isize or tlen    template_length

                    # oneReadDict["qname"] = alignedRead.query_name
                    oneReadDict["flag"] = alignedRead.flag
                    oneReadDict["name"] = alignedRead.query_name
                    oneReadDict["start"] = alignedRead.reference_start
                    oneReadDict["mapQual"] = alignedRead.mapping_quality
                    # oneReadDict["cigar"] = alignedRead.cigar
                    # oneReadDict["mateName"] = alignedRead.next_reference_name
                    oneReadDict["mateStart"] = alignedRead.next_reference_start
                    oneReadDict["insertSize"] = alignedRead.template_length
                    oneReadDict["sequence"] = alignedRead.seq
                    # oneReadDict["qualities"] = alignedRead.qual
                    oneReadDict["qlen"] = len(alignedRead.query_sequence)
                    oneReadDict["base"] = alignedRead.seq[qPos]
                    oneReadDict["baseQual"] = alignedRead.qual[qPos]
                    oneReadDict["sequenceIndex"] = pileupRead.query_position
                    oneReadDict["chrom"] = chrom
                    oneReadDict["pos"] = pos
                    oneReadDict["strand"] = strand
                    oneReadDict["refBase"] = refBase

                    # the RSEM fasta for transcripts has reads in the 5' to 3'
                    # direction. the RNA bam files have reads aligned to the
                    # RSEM fasta in the 5' to 3' direction. if the transcript
                    # is on the genomic "-" strand, then we need to reverse
                    # complement the base. otherwise, we just use the forward
                    # strand base
                    strandedBase = alignedRead.seq[qPos]
                    if (strand is not None and strand == "-"):
                        strandedBase = rev_comp_nucleotide(strandedBase)
                    oneReadDict["strandedBase"] = strandedBase

                    '''
                    if (anIsDebug):
                        logging.debug("name=%s, oneReadDict=%s",
                                      alignedRead.query_name, oneReadDict)
                    '''

                    if (strandedBase in anAlleleSet):
                        '''
                        if (anIsDebug):
                            logging.debug("strandedBase=%s, is in " +
                                          "alleleSet=%s: name=%s " +
                                          "oneReadDict=%s",
                                          strandedBase, anAlleleSet,
                                          alignedRead.query_name, oneReadDict)
                        '''

                        # add it to the dictionary of reads
                        readsDict[alignedRead.query_name].append(oneReadDict)

                    elif anIsDebug:
                        logging.debug("throwing out read where " +
                                      "strandedBase=%s is not in the " +
                                      "alleleSet=%s",
                                      strandedBase, anAlleleSet)

                if (anIsDebug):
                    logging.debug("group_reads_by_name(): %s:%s added %s " +
                                  "out of %s reads to initial reads dict, " +
                                  "readsDictLen=%s", chrom, pos, keptReads,
                                  totalReads, len(readsDict.keys()))

        return readsDict


def find_non_overlapping_reads(aReadsDict, aMinBaseQual, anIsDebug):
    '''
    ' This function loops through the reads with the same name.  Due to the
    ' inclusion of secondary alignments for RNA-Seq, there could be more
    ' than 2 reads.  If the reads overlap and the bases agree, then we keep
    ' the read with the highest base quality. If the reads overlap and the
    ' bases disagree, then we only keep the read with the highest base quality
    ' if the other reads have a low quality.  If the reads don't overlap,
    ' then keep all the reads.
    '
    ' aReadsDict:    The dictionary of reads grouped by their read name
    ' aMinBaseQual:  A minimum base quality score for the base
    ' anIsDebug:     A flag for outputting debug messages to STDERR
    '''
    # loop through reads
    #     if reads overlap:
    #         if the bases agree:
    #             keep the one with the highest base quality
    #         else bases don't agree:
    #            keep the read with the high quality only if
    #            the other reads have a low quality
    #     else:
    #        keep all non-overlaps
    nonOverlapReadsList = []
    for (readName, readList) in aReadsDict.iteritems():
        '''
        if (anIsDebug):
            logging.debug("readName=%s, len=%s, readList=%s", readName,
                          len(readList), readList)
        '''
        # if there is only one read, then no overlapping reads exist
        if (len(readList) == 1):
            # if (anIsDebug):
            #    logging.debug("only one read %s", readName)
            nonOverlapReadsList.append(readList[0])
        # check for overlapping reads
        else:
            maxBaseQual = -1
            maxBaseQualIndex = -1
            nextMaxBaseQual = -1
            basesSet = set()

            # for each read
            for index in range(0, len(readList)):
                readStart = readList[index]["start"]
                readSeqIndex = readList[index]["sequenceIndex"]
                readMateStart = readList[index]["mateStart"]
                readLength = readList[index]["qlen"]
                # abs(core.pos + p->qpos - core.mpos) < (core.l_qseq)
                # samtools view pos = reference_start and
                # samtools view pos = 1-based leftmost coordinate
                # pysam pos = reference_start = 0-based leftmost coordinate
                # samtools view pnext (1-based) = pysam mpos (0-based) =
                #    next_reference_start = the pos of the mate/next read
                # sequenceIndex pysam qpos = query_position = position of
                #    the read base at the pileup site, 0-based.
                #    None if is_del or is_refskip is set.
                # len(sequence) or pysam qlen = query_length =
                #    length of the query sequence

                '''
                if (anIsDebug):
                    logging.debug("readStart(pos)=%s, " +
                                  "readSeqIndex(qpos)=%s, " +
                                  "readMateStart(mpos)=%s, " +
                                  "readLength(l_qseq)=%s", readStart,
                                  readSeqIndex, readMateStart, readLength)
                    logging.debug("%s (readStart + readSeqIndex): " +
                                  "abs(readStart + readSeqIndex - " +
                                  "readMateStart) %s <? " +
                                  "readAlignedLength %s",
                                  (readStart + readSeqIndex),
                                  abs(readStart + readSeqIndex -
                                      readMateStart),
                                  readLength)
                '''

                # if abs(core.pos + p->qpos - core.mpos) < core.l_qseq)
                readSpan = readStart + readSeqIndex - readMateStart
                if (abs(readSpan) < readLength):
                    '''
                    if (anIsDebug):
                        logging.debug("reads overlap readStart=%s, " +
                                      "sequenceIndex=%s, " +
                                      "readMateStart=%s, " +
                                      "abs=%s, alignedLen=%s",
                                      str(readStart), str(readSeqIndex),
                                      str(readMateStart),
                                      str(abs(readSpan)),
                                      str(readLength))
                    '''

                    base = readList[index]["base"]
                    strandedBase = readList[index]["strandedBase"]
                    qual = ord(readList[index]["baseQual"])-33
                    basesSet.add(base)
                    if (qual > maxBaseQual):
                        # if the current max is higher than the next max,
                        # then propagate the current max down
                        if (maxBaseQual > nextMaxBaseQual):
                            nextMaxBaseQual = maxBaseQual
                        maxBaseQual = qual
                        maxBaseQualIndex = index
                    elif (qual > nextMaxBaseQual):
                        nextMaxBaseQual = qual

                    '''
                    if (anIsDebug):
                        logging.debug("base=%s, strandedBase=%s, qual=%s, " +
                                      "maxBaseQual=%s, " +
                                      "maxBaseQualIndex=%s, " +
                                      "basesSet=%s", base, strandedBase,
                                      str(qual), str(maxBaseQual),
                                      str(maxBaseQualIndex), basesSet)
                    '''
                # if they don't overlap, then just add the reads
                else:
                    '''
                    if (anIsDebug):
                        logging.debug("reads don't overlap " +
                                      "readStart=%s, sequenceIndex=%s, " +
                                      "readMateStart=%s, abs=%s, " +
                                      "alignedLen=%s", str(readStart),
                                      str(readSeqIndex),
                                      str(readMateStart),
                                      str(abs(readSpan)),
                                      str(readLength))
                    '''
                    nonOverlapReadsList.append(readList[index])

            # if all the bases agree, then keep the one with the
            # highest base quality and avoid double counting
            if (len(basesSet) == 1):
                '''
                if (anIsDebug):
                    logging.debug("all bases in overlapping reads " +
                                  "agree, keeping=%s, %s %s",
                                  str(maxBaseQualIndex), readName,
                                  readList)
                '''
                nonOverlapReadsList.append(readList[maxBaseQualIndex])
            # if the bases disagree, then it's likely a sequencing error.
            # only keep the read with a high qual if the other reads
            # have a low qual
            elif(maxBaseQual >= aMinBaseQual and
                 nextMaxBaseQual < aMinBaseQual):
                '''
                if (anIsDebug):
                    logging.debug("not all bases in overlapping reads " +
                                  "agree, but keeping a high qual one, " +
                                  "maxBaseQual=%s, nextMaxBaseQual=%s, " +
                                  "minBaseQual=%s, keeping=%s, %s %s",
                                  str(maxBaseQual), str(nextMaxBaseQual),
                                  str(aMinBaseQual), str(maxBaseQualIndex),
                                  readName, readList)
                '''
                nonOverlapReadsList.append(readList[maxBaseQualIndex])
            '''
            elif anIsDebug:
                logging.debug("bases in overlapping reads don't agree " +
                              "and no obvious read to select based on " +
                              "qual scores, maxBaseQual=%s, " +
                              "nextMaxBaseQual=%s, minBaseQual=%s, %s %s",
                              str(maxBaseQual), str(nextMaxBaseQual),
                              str(aMinBaseQual), readName, readList)
            '''
    if (anIsDebug):
        logging.debug("find_non_overlapping_reads(): " +
                      "nonOverlapReadsLen=%s", len(nonOverlapReadsList))

    return nonOverlapReadsList


def write_to_blat_file(aBlatFileHandler, aGenomicChr, aGenomicCoordinate,
                       aChromList, aPosList, aTxStrandList,
                       aParamsDict, anInfoDict, aPrefix, anAltOnlyFlag,
                       anRnaIncludeSecondaryAlignmentsFlag, anIsDebug):
    '''
    ' This function gets all of the reads at a specific
    ' coordinate and creates a BLAT input file.
    '
    ' aBlatFileHandler:
    '    A file handler where all the BLAT query data is written
    ' aGenomicChr:
    '    The genomic chromosome
    ' aGenomicCoordinate:
    '    The genomic stop coordinate
    ' aChromList:
    '    A list of chromosomes or transcript names to process
    ' aPosList:
    '    A list of coordinates to process
    ' aTxStrandList:
    '    A list of transcript strands
    ' aParamsDict:
    '    The parameters from the header
    ' anInfoDict:
    '    A dictionary of the INFO column
    ' aPrefix:
    '    The prefix for the parameters dictionary
    ' anAltOnlyFlag:
    '    If all reads should be processed or only those
    '    that have the alternate allele
    ' anRnaIncludeSecondaryAlignmentsFlag:
    '    A flag to include the RNA secondary alignments
    ' anIsDebug:
    '    A flag for outputting debug messages to STDERR
    '''

    # get the info for executing the samtools command
    bamFilename = aParamsDict[aPrefix + "Filename"]
    fastaFilename = aParamsDict[aPrefix + "FastaFilename"]
    minBaseQual = aParamsDict[aPrefix + "MinBaseQuality"]
    bamOrigin = anInfoDict["ORIGIN"]

    if (not os.path.isfile(bamFilename)):
        logging.critical("The BAM file specified in the VCF header " +
                         "does not exist: %s", bamFilename)
        sys.exit(1)

    if (not os.path.isfile(fastaFilename)):
        logging.critical("The FASTA file specified in the VCF header " +
                         "does not exist: %s", fastaFilename)
        sys.exit(1)

    bamFile = pysam.Samfile(bamFilename, 'rb')
    fastaFile = pysam.Fastafile(fastaFilename)

    # get the set of alleles for this call
    modChanges = anInfoDict["MC"]
    modTypes = anInfoDict["MT"]
    altSet = set()
    refSet = set()
    for (modChange, modType) in izip(modChanges, modTypes):
        (ref, alt) = modChange.split(">")
        refSet.add(ref)
        altSet.add(alt)

    alleleSet = refSet.union(altSet)

    # group all of the reads by name
    readsDict = group_reads_by_name(aChromList, aPosList,
                                    aTxStrandList, alleleSet,
                                    bamFile, fastaFile, bamOrigin,
                                    anRnaIncludeSecondaryAlignmentsFlag,
                                    anIsDebug)

    # get all of the non-overlapping reads
    nonOverlapReadsList = find_non_overlapping_reads(readsDict,
                                                     minBaseQual,
                                                     anIsDebug)

    numRefsFound = 0
    numAltsFound = 0
    numNonRefAltBasesFound = 0
    pysamBases = ""
    pysamQuals = ""

    for readDict in nonOverlapReadsList:
        readBase = readDict["base"]
        strandedBase = readDict["strandedBase"]
        readBaseQual = readDict["baseQual"]
        readSequence = readDict["sequence"]
        readName = readDict["name"]
        readMapQual = readDict["mapQual"]
        readLength = readDict["qlen"]
        readFlag = readDict["flag"]

        pysamBases += strandedBase
        pysamQuals += readBaseQual

        # count the number of times a base is not the ref nor the alt
        if (strandedBase not in refSet and strandedBase not in altSet):
            if (anIsDebug):
                logging.debug("base=%s (orgBase=%s) is not a ref=%s, nor an " +
                              "alt=%s", strandedBase, readBase, refSet, altSet)
            numNonRefAltBasesFound += 1

        # count the refs
        if (strandedBase in refSet):
            if (anIsDebug):
                logging.debug("base=%s (orgBase=%s) is a ref=%s",
                              strandedBase, readBase, refSet)

            baseQualityConverted = ord(readBaseQual)-33
            if (baseQualityConverted >= int(minBaseQual)):
                numRefsFound += 1

        # if only the reads with the alts should be processed and
        # the base is an alt, or if all reads should be processed
        elif ((anAltOnlyFlag and strandedBase in altSet) or
              (not anAltOnlyFlag)):

            baseQualityConverted = ord(readBaseQual)-33

            if (anIsDebug):
                logging.debug("baseQualOrg=%s, baseQualityOrgOrd=%s, " +
                              "baseQualConverted=%s, minBaseQual=%s, " +
                              "greaterThan? %s", readBaseQual,
                              str(ord(readBaseQual)),
                              str(baseQualityConverted),
                              str(minBaseQual),
                              str(baseQualityConverted >= int(minBaseQual)))

            if (baseQualityConverted >= int(minBaseQual)):
                numAltsFound += 1
                if (anIsDebug):
                    logging.debug("found an alt for %s:%s, " +
                                  "base=%s, orgBase=%s, baseQual=%s",
                                  aChromList, aPosList, strandedBase,
                                  readBase, str(baseQualityConverted))

                # we want to use the genomic chr and stop coordinate in
                # the output instead of the transcript coordinate

                # rnaTumor_3_79482763_readName_A_41_11_355_133
                # NM_181742       95.49   133     6       0       1       133
                # 4688    4556    1.2e-63 240.0

                outputList = [aPrefix, aGenomicChr, str(aGenomicCoordinate),
                              readName.replace("_", ""), strandedBase,
                              str(baseQualityConverted), str(readMapQual),
                              str(readFlag), str(readLength)]
                aBlatFileHandler.write("> " + "_".join(outputList) + "\n")
                aBlatFileHandler.write(readSequence + "\n")

    if (anIsDebug):
        numTotal = numRefsFound + numAltsFound + numNonRefAltBasesFound
        logging.debug("For %s:%s, numNonRefAltBasesFound=%s, " +
                      "numAltsFound=%s, numRefsFound=%s, numTotal=%s",
                      aChromList, aPosList, str(numNonRefAltBasesFound),
                      str(numAltsFound), str(numRefsFound), str(numTotal))
        logging.debug("pysam bases=%s", pysamBases)
        logging.debug("pysam quals=%s", pysamQuals)

    return


def main():

    # command for running this on a small test case:
    # python createBlatFile.py TCGA-00-4454 7
    # ../data/test/TCGA-00-4454_EGFR.vcf
    # ../data/test/tmp/
    # --dnaNormalFilename=../data/test/TCGA-00-4454_EGFR.reads

    startTime = time.time()

    # create the usage statement
    usage = "usage: python %prog id vcfFile headerFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)

    # add the optional parameters
    i_cmdLineParser.add_option(
        "-c", "--allVCFCalls", action="store_false", default=True,
        dest="passedVCFCallsOnly",
        help="by default only the VCF calls that have passed all filters " +
             "thus far are processed, include this argument if all of the " +
             "VCF calls should be processed")
    i_cmdLineParser.add_option(
        "-b", "--allReadBases", action="store_false", default=True,
        dest="altBasesOnly",
        help="by default only the reads with the alternate base are " +
             "processed, include this argument if all of the reads " +
             "should be processed")

    i_cmdLineParser.add_option(
        "-o", "--outputFilename",
        dest="outputFilename", metavar="OUTPUT_FILE", default=sys.stdout,
        help="the name of the output file, STDOUT by default")
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

    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(3, 22, 1)
    i_argLength = len(sys.argv)

    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)

    # get the required parameters
    (cmdLineOpts, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_id = i_cmdLineArgs[0]
    i_vcfFilename = i_cmdLineArgs[1]
    i_headerFilename = i_cmdLineArgs[2]

    # get the optional params with default values
    i_passedVCFCallsOnlyFlag = cmdLineOpts.passedVCFCallsOnly
    i_altBasesOnlyFlag = cmdLineOpts.altBasesOnly
    i_logLevel = cmdLineOpts.logLevel
    i_rnaIncludeSecondaryAlignments = cmdLineOpts.rnaIncludeSecondaryAlignments

    i_blatDnaNormalReads = cmdLineOpts.blatDnaNormalReads
    i_blatDnaTumorReads = cmdLineOpts.blatDnaTumorReads
    i_blatRnaNormalReads = cmdLineOpts.blatRnaNormalReads
    i_blatRnaTumorReads = cmdLineOpts.blatRnaTumorReads

    # try to get any optional parameters with no defaults
    i_readFilenameList = [i_vcfFilename, i_headerFilename]
    i_writeFilenameList = []

    i_logFilename = None
    i_outputFilename = None
    i_transcriptNameTag = None
    i_transcriptCoordinateTag = None
    i_transcriptStrandTag = None
    if (cmdLineOpts.logFilename is not None):
        i_logFilename = cmdLineOpts.logFilename
        i_writeFilenameList += [i_logFilename]
    if (cmdLineOpts.outputFilename is not None):
        i_outputFilename = cmdLineOpts.outputFilename
        if (cmdLineOpts.outputFilename is not sys.stdout):
            i_writeFilenameList += [i_outputFilename]
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

    # set the debug flag
    i_debug = (i_numericLogLevel == logging.DEBUG)

    # output some debug info
    if (i_debug):
        logging.debug("id=%s", i_id)
        logging.debug("vcfFilename=%s", i_vcfFilename)
        logging.debug("headerFilename=%s", i_headerFilename)
        logging.debug("outputFilename=%s", i_outputFilename)

        logging.debug("passedCallsOnly? %s", i_passedVCFCallsOnlyFlag)
        logging.debug("altBasesOnlyFlag? %s", i_altBasesOnlyFlag)

        logging.debug("transcriptNameTag %s", i_transcriptNameTag)
        logging.debug("transcriptCoordinateTag %s", i_transcriptCoordinateTag)
        logging.debug("transcriptStrandTag %s", i_transcriptStrandTag)
        logging.debug("rnaInclSecAlign=%s" % i_rnaIncludeSecondaryAlignments)

        logging.debug("blatDnaNormal? %s", i_blatDnaNormalReads)
        logging.debug("blatDnaTumor? %s", i_blatDnaTumorReads)
        logging.debug("blatRnaNormal? %s", i_blatRnaNormalReads)
        logging.debug("blatRnaTumor? %s", i_blatRnaTumorReads)

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
    vcfGenerator = get_vcf_data(i_vcfFilename,
                                i_headerFilename,
                                i_passedVCFCallsOnlyFlag,
                                i_debug)

    # for each VCF call that should be investigated
    for (vcfChr, vcfStopCoordinate, vcfId, vcfRefList, vcfAltList, vcfScore,
         vcfFilterSet, vcfInfoDict, restOfLine, vcfParamsDict) in vcfGenerator:
        if (i_debug):
            logging.debug("VCF Data: %s %s %s %s %s %s %s %s %s", vcfChr,
                          vcfStopCoordinate, vcfId, vcfRefList,
                          vcfAltList, vcfScore, vcfFilterSet,
                          vcfInfoDict, restOfLine)

        modTypes = vcfInfoDict["MT"]
        for modType in modTypes:

            # get the reads contributing to a call and
            # put them in a blat query file
            if (modType == "GERM" and i_blatDnaNormalReads):
                write_to_blat_file(i_outputFileHandler,
                                   vcfChr,
                                   vcfStopCoordinate,
                                   [vcfChr],
                                   [vcfStopCoordinate],
                                   [None],
                                   vcfParamsDict,
                                   vcfInfoDict,
                                   "dnaNormal",
                                   i_altBasesOnlyFlag,
                                   False,
                                   i_debug)

            elif (modType == "NOR_EDIT" and i_blatRnaNormalReads):
                # if we should process the transcripts
                if ((i_transcriptNameTag is not None) and
                    (i_transcriptNameTag in vcfInfoDict)):
                    write_to_blat_file(i_outputFileHandler,
                                       vcfChr,
                                       vcfStopCoordinate,
                                       vcfInfoDict[i_transcriptNameTag],
                                       vcfInfoDict[i_transcriptCoordinateTag],
                                       vcfInfoDict[i_transcriptStrandTag],
                                       vcfParamsDict,
                                       vcfInfoDict,
                                       "rnaNormal",
                                       i_altBasesOnlyFlag,
                                       i_rnaIncludeSecondaryAlignments,
                                       i_debug)
                else:
                    write_to_blat_file(i_outputFileHandler,
                                       vcfChr,
                                       vcfStopCoordinate,
                                       [vcfChr],
                                       [vcfStopCoordinate],
                                       [None],
                                       vcfParamsDict,
                                       vcfInfoDict,
                                       "rnaNormal",
                                       i_altBasesOnlyFlag,
                                       i_rnaIncludeSecondaryAlignments,
                                       i_debug)

            elif (modType == "SOM" and i_blatDnaTumorReads):
                write_to_blat_file(i_outputFileHandler,
                                   vcfChr,
                                   vcfStopCoordinate,
                                   [vcfChr],
                                   [vcfStopCoordinate],
                                   [None],
                                   vcfParamsDict,
                                   vcfInfoDict,
                                   "dnaTumor",
                                   i_altBasesOnlyFlag,
                                   False,
                                   i_debug)

            elif ((modType == "SOM" or modType == "TUM_EDIT") and
                (i_blatRnaTumorReads)):
                # if we should process the transcripts
                if ((i_transcriptNameTag is not None) and
                    (i_transcriptNameTag in vcfInfoDict)):
                    write_to_blat_file(i_outputFileHandler,
                                       vcfChr,
                                       vcfStopCoordinate,
                                       list(vcfInfoDict[i_transcriptNameTag]),
                                       vcfInfoDict[i_transcriptCoordinateTag],
                                       vcfInfoDict[i_transcriptStrandTag],
                                       vcfParamsDict,
                                       vcfInfoDict,
                                       "rnaTumor",
                                       i_altBasesOnlyFlag,
                                       i_rnaIncludeSecondaryAlignments,
                                       i_debug)
                else:
                    write_to_blat_file(i_outputFileHandler,
                                       vcfChr,
                                       vcfStopCoordinate,
                                       [vcfChr],
                                       [vcfStopCoordinate],
                                       [None],
                                       vcfParamsDict,
                                       vcfInfoDict,
                                       "rnaTumor",
                                       i_altBasesOnlyFlag,
                                       i_rnaIncludeSecondaryAlignments,
                                       i_debug)

    stopTime = time.time()
    logging.info("createBlatFile.py Id %s: Total time=%s hrs, %s mins, " +
                 "%s secs", i_id, ((stopTime-startTime)/(3600)),
                 ((stopTime-startTime)/60), (stopTime-startTime))

    # close the files
    if (i_outputFilename is not sys.stdout):
        i_outputFileHandler.close()

    return


main()
sys.exit(0)
