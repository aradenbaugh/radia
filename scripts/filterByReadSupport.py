#!/usr/bin/env python

__requires__ = ['pysam>=0.8.1']

import pysam
import sys
import os
from optparse import OptionParser
import myVCF
import radiaUtil
import logging
from itertools import izip
import collections
import math

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

i_cigarDict = {}
i_cigarDict[0] = "match"
i_cigarDict[1] = "insertion"
i_cigarDict[2] = "deletion"
i_cigarDict[3] = "refskipped"
i_cigarDict[4] = "softclipped"
i_cigarDict[5] = "hardclipped"
i_cigarDict[6] = "padding"
i_cigarDict[7] = "seqmatch"
i_cigarDict[8] = "seqmismatch"


def get_passing_germline_alts(aCurrData):
    alts = []
    # for passing germline calls,
    # there should only be one,
    # but double-check anyway
    for modChange in aCurrData.infoDict["MC"]:
        (ref, alt) = modChange.split(">")

        # we can use the genomic forward strand alts because
        # mismatch_counts() will convert the reads to the forward strand
        alts.append(alt)

    return alts


def parse_vcf(aVCFFilename, aTxNameTag, aTxCoordinateTag, anIsDebug):

    vcf = radiaUtil.get_read_fileHandler(aVCFFilename)
    currVCF = myVCF.VCF()

    dnaNormBam = None
    dnaTumBam = None
    dnaTumFasta = None
    dnaTumChrPrefix = False
    rnaNormBam = None
    rnaNormFasta = None
    rnaNormChrPrefix = False
    rnaTumBam = None
    rnaTumFasta = None
    rnaTumChrPrefix = False

    germDict = {}
    txGermDict = {}
    mutationsDict = {}
    filterDict = {}

    # loop through and categorize all calls
    for line in vcf:

        # if we find the vcfGenerator line, then parse out the file info
        if ("vcfGenerator" in line):
            generatorLine = line[0:(len(line)-1)]
            generatorLine = generatorLine[16:len(generatorLine)]
            generatorParamsList = generatorLine.split(",")
            generatorParamsDict = {}

            # create a dictionary of existing params
            for param in generatorParamsList:
                (key, value) = param.split("=")
                value = value.rstrip(">")
                value = value.lstrip("<")
                generatorParamsDict[key] = value

            if ("dnaNormalFilename" in generatorParamsDict):
                dnaNormBam = generatorParamsDict["dnaNormalFilename"]
            if ("rnaNormalFilename" in generatorParamsDict):
                rnaNormBam = generatorParamsDict["rnaNormalFilename"]
                rnaNormFasta = generatorParamsDict["rnaNormalFastaFilename"]
                rnaNormChrPrefix = generatorParamsDict["rnaNormalUseChrPrefix"]
            if ("dnaTumorFilename" in generatorParamsDict):
                dnaTumBam = generatorParamsDict["dnaTumorFilename"]
                dnaTumFasta = generatorParamsDict["dnaTumorFastaFilename"]
                dnaTumChrPrefix = generatorParamsDict["dnaTumorUseChrPrefix"]
            if ("rnaTumorFilename" in generatorParamsDict):
                rnaTumBam = generatorParamsDict["rnaTumorFilename"]
                rnaTumFasta = generatorParamsDict["rnaTumorFastaFilename"]
                rnaTumChrPrefix = generatorParamsDict["rnaTumorUseChrPrefix"]
            continue
        # we found the chrom header line
        elif line.startswith("#CHROM"):
            currVCF.set_headers(line[1:].strip().split("\t"))
            continue
        # elif skip over the other header lines
        elif (line.startswith("#")):
            continue

        # now we are to the data
        dataAsList = line.strip().split("\t")
        # if the number of VCF columns doesn't
        # equal the number of VCF header columns
        if len(dataAsList) == len(currVCF.headers):

            # parse each VCF line
            currData = currVCF.make_data(dataAsList)

            # keep track of the passing germline and LOH calls
            if (("SNP" in currData.infoDict["VT"] and
                 "PASS" in currData.filterSet) and
                ("1" in currData.infoDict["SS"] or
                 "3" in currData.infoDict["SS"])):
                # initialize the dict
                if (currData.chrom not in germDict):
                    germDict[currData.chrom] = {}

                # add the germline alleles with the
                # genomic coordinates to the germDict
                altsList = get_passing_germline_alts(currData)
                germDict[currData.chrom][str(currData.pos-1)] = altsList

                # if we have transcript names and coordinates for these calls
                if (aTxNameTag is not None and aTxCoordinateTag is not None):
                    # for each transcript name and coordinate
                    txNames = currData.infoDict[aTxNameTag]
                    txCoordinates = currData.infoDict[aTxCoordinateTag]
                    for (txName, txCoordinate) in izip(txNames,
                                                       txCoordinates):
                        # initialize the dict
                        if txName not in txGermDict:
                            txGermDict[txName] = {}

                        # add the germline alleles with the
                        # transcript coordinates to the txGermDict
                        altsList = get_passing_germline_alts(currData)
                        txGermDict[txName][str(int(txCoordinate)-1)] = altsList
            '''
            # keep track of the passing somatic calls
            elif ("SNP" in currData.infoDict["VT"] and
                  "2" in currData.infoDict["SS"] and
                  "PASS" in currData.filterSet):
                if (currData.chrom not in mutationsDict):
                    mutationsDict[currData.chrom] = {}
                mutationsDict[currData.chrom][str(currData.pos-1)] = currData
            # keep track of the filtered calls
            elif ("SNP" in currData.infoDict["VT"] and
                  "2" in currData.infoDict["SS"] and
                  "PASS" not in currData.filterSet):
                if (currData.chrom not in filterDict):
                    filterDict[currData.chrom] = {}
                filterDict[currData.chrom][str(currData.pos-1)] = currData
            '''
        else:
            logging.warning("The number of VCF columns (%s) doesn't equal " +
                            "the number of VCF header columns (%s).",
                            len(dataAsList), len(currVCF.headers))
            logging.warning("Here are the header columns: %s", currVCF.headers)
            logging.warning("Here is the VCF line: %s", line.strip())

    # if (anIsDebug):
    #    logging.debug("germDict=%s", germDict)
    #    logging.debug("txGermDict=%s", txGermDict)
    #    logging.debug("mutationsDict=%s", mutationsDict)
    #    logging.debug("filterDict=%s", filterDict)

    return (dnaNormBam, rnaNormBam, rnaNormFasta, rnaNormChrPrefix,
            dnaTumBam, dnaTumFasta, dnaTumChrPrefix,
            rnaTumBam, rnaTumFasta, rnaTumChrPrefix,
            germDict, txGermDict, mutationsDict, filterDict)


def low_base_or_map_quals(aPileupread, aParamsDict, anIsDebug):

    # mapping quality scores are already converted to ints
    if aPileupread.alignment.mapq < aParamsDict["minMapQual"]:
        if (anIsDebug):
            logging.debug("low_base_or_map_quals MQ=%s < minMQ=%s",
                          aPileupread.alignment.mapq,
                          aParamsDict["minMapQual"])
        return True

    # the pysam documentation says: "base quality scores are unsigned chars but
    # they are *not* the ASCII encoded values, so no offset of 33 needs to be
    # subtracted", but this is only true for the 'query_qualities' and
    # 'query_alignment_qualities'. we need to subtract the offset for the
    # 'qual' field
    alignedRead = aPileupread.alignment
    qPos = aPileupread.query_position
    baseQualConverted = ord(alignedRead.qual[qPos])-33
    # queryQual = alignedRead.query_qualities[qPos]
    # baseQualConverted = ord(alignedRead.query_alignment_qualities[qPos])-33

    '''
    if (anIsDebug):
        logging.debug("alignedRead alignment.qual=%s", alignedRead.qual)
        logging.debug("alignedRead alignment.query_qualities=%s",
                      alignedRead.query_qualities)
        logging.debug("alignedRead alignment.query_alignment_qualities=%s",
                      alignedRead.query_alignment_qualities)
        logging.debug("low_base_or_map_quals baseQuality qpos=%s, " +
                      "orgQual=%s, ordQual=%s, convertedQual=%s, minBQ=%s",
                      qPos, alignedRead.qual[qPos],
                      ord(alignedRead.qual[qPos]),
                      baseQualConverted, aParamsDict["minBaseQual"])
    '''

    if baseQualConverted < aParamsDict["minBaseQual"]:
        if (anIsDebug):
            logging.debug("low_base_or_map_quals base BQ=%s < minBQ=%s",
                          baseQualConverted, aParamsDict["minBaseQual"])
        return True

    if (anIsDebug):
        logging.debug("low_base_or_map_quals found nothing")

    return False


def low_neighbor_base_quals(aPileupread, aParamsDict, anIsDebug):

    alignedRead = aPileupread.alignment
    qPos = aPileupread.query_position

    # calculate indices to check (e.g. 5 before and 5 after)
    start = qPos - aParamsDict["numNeighborBases"]
    if (start < 0):
        start = 0

    stop = qPos + aParamsDict["numNeighborBases"]
    if stop > alignedRead.rlen:
        stop = alignedRead.rlen
        '''
        if (anIsDebug):
            logging.debug("neighbors past the length query_position=%s, " +
                          "rlen=%s, start=%s, stop=%s",
                          qPos, alignedRead.rlen, start, stop)
        '''

    '''
    if (anIsDebug and (stop - start < 2 * aParamsDict["numNeighborBases"])):
        logging.debug("neighbor bases not 2x numNeighborBases: " +
                      "query_position=%s, rlen=%s, start=%s, stop=%s",
                      qPos, alignedRead.rlen, start, stop)
    '''

    # if the region examined is somehow negative
    if start > stop:
        '''
        if (anIsDebug):
            logging.debug("checking negative neighbor bases " +
                          "query_position=%s, rlen=%s, start=%s, stop=%s",
                          qPos, alignedRead.rlen, start, stop)
        '''
        return False

    # if anything in neighborhood has too low of base quality
    index = start
    for qualScore in alignedRead.qual[start:(stop + 1)]:
        # the pysam documentation says: "base quality scores are unsigned chars
        # but they are *not* the ASCII encoded values, so no offset of 33 needs
        # to be subtracted", but this is only true for the 'query_qualities'
        # and 'query_alignment_qualities'. we need to subtract the offset for
        # the 'qual' field
        baseQualConverted = ord(qualScore)-33
        '''
        if (anIsDebug):
            logging.debug("checking neighbor base query_pos=%s, start=%s, " +
                          "stop=%s, index=%s, qual=%s, ordQual=%s, " +
                          "convertedQual=%s", aPileupread.query_position,
                          start, stop, index, qualScore, ord(qualScore),
                          baseQualConverted)
        '''

        if baseQualConverted < aParamsDict["minNeighborBaseQual"]:
            if (anIsDebug):
                logging.debug("low_neighbor_base_quals neighboring base " +
                              "BQ=%s < minNBQ=%s", baseQualConverted,
                              aParamsDict["minNeighborBaseQual"])
            return True
        index += 1

    if (anIsDebug):
        logging.debug("low_neighbor_base_quals found nothing")

    return False


def reverse_complement_nucleotide(aNucleotide):
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


def is_mutation(aReadDict, aRef, anAlt, aBamOrigin, anMMPDict, anIsDebug):

    orgReadBase = aReadDict["base"]
    readBase = aReadDict["base"]
    orgRefBase = aReadDict["refBase"]
    refBase = aReadDict["refBase"]
    alignedRead = aReadDict["alignedRead"]

    # the RSEM fasta for transcripts has reads in the 5' to 3' direction
    # the RNA bam files have reads aligned to the RSEM fasta in the
    # 5' to 3' direction.  if the transcript is on the genomic "-" strand,
    # then we need to reverse complement the reads and the fasta ref
    # otherwise, we can just compare the forward strand base and fasta ref
    if (aBamOrigin == "RNA" and
        aReadDict["strand"] is not None and
        aReadDict["strand"] == "-"):
        readBase = reverse_complement_nucleotide(readBase)
        refBase = reverse_complement_nucleotide(refBase)

    if readBase not in anMMPDict:
        anMMPDict[readBase] = collections.defaultdict(int)

    if alignedRead.is_secondary:
        anMMPDict[readBase]["secondary"] += 1
        anMMPDict[readBase]["total"] += 1
    else:
        anMMPDict[readBase]["total"] += 1

    if (readBase == aRef):
        if (anIsDebug):
            logging.debug("is_mutation() base matches reference, " +
                          "queryPos=%s, chrom=%s, pos=%s, orgBase=%s, " +
                          "orgRef=%s, base=%s, ref=%s, vcfAlt=%s, vcfRef=%s",
                          aReadDict["sequenceIndex"], aReadDict["chrom"],
                          aReadDict["pos"], readBase, refBase, orgReadBase,
                          orgRefBase, anAlt, aRef)
        return False
    elif readBase == anAlt:
        if (anIsDebug):
            logging.debug("is_mutation() base matches alt, queryPos=%s, " +
                          "chrom=%s, pos=%s, orgBase=%s, orgRef=%s, " +
                          "base=%s, ref=%s, vcfAlt=%s, vcfRef=%s",
                          aReadDict["sequenceIndex"], aReadDict["chrom"],
                          aReadDict["pos"], readBase, refBase, orgReadBase,
                          orgRefBase, anAlt, aRef)
        return True

    if (anIsDebug):
        logging.debug("is_mutation() found nothing, queryPos=%s, chrom=%s, " +
                      "pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s, " +
                      "vcfAlt=%s, vcfRef=%s", aReadDict["sequenceIndex"],
                      aReadDict["chrom"], aReadDict["pos"], readBase, refBase,
                      orgReadBase, orgRefBase, anAlt, aRef)
    return False


def mismatch_counts(aCigarNum, aChrom, aRefIndex, aQueryIndex, aMutSS,
                    aTxStrand, anAlignedRead, aGermDict,
                    aTxGermDict, aFastafile, aBamOrigin,
                    aParamsDict, anIsDebug):
    refs = 0
    germs = 0
    muts = 0
    edits = 0

    for offset in range(aCigarNum):
        # get the ref base
        refPos = aRefIndex + offset
        orgRefBase = aFastafile.fetch(aChrom, refPos, refPos+1).upper()
        refBase = aFastafile.fetch(aChrom, refPos, refPos+1).upper()

        # get the query base
        orgReadBase = anAlignedRead.seq[aQueryIndex + offset]
        readBase = anAlignedRead.seq[aQueryIndex + offset]

        # the RSEM fasta for transcripts has reads in the 5' to 3' direction
        # the RNA bam files have reads aligned to the RSEM fasta in the
        # 5' to 3' direction.  if the transcript is on the genomic "-" strand,
        # then we need to reverse complement the reads and the fasta ref
        # otherwise, we just compare the forward strand base and fasta ref
        if (aBamOrigin == "RNA" and
            aTxStrand is not None and
            aTxStrand == "-"):
            readBase = reverse_complement_nucleotide(readBase)
            refBase = reverse_complement_nucleotide(refBase)

        if readBase == refBase:
            refs += 1
            if (anIsDebug):
                logging.debug("mismatch_counts() base matches ref:  " +
                              "aRefIndex=%s, aQueryIndex=%s, offset=%s, " +
                              "chrom=%s, pos=%s, orgBase=%s, orgRef=%s, " +
                              "base=%s, ref=%s, txStrand=%s, r1=%s, r2=%s, " +
                              "rev=%s, mateRev=%s",
                              aRefIndex, aQueryIndex, offset, aChrom, refPos,
                              orgReadBase, orgRefBase, readBase, refBase,
                              aTxStrand, anAlignedRead.is_read1,
                              anAlignedRead.is_read2, anAlignedRead.is_reverse,
                              anAlignedRead.mate_is_reverse)
        elif (aChrom in aGermDict and
              str(refPos) in aGermDict[aChrom] and
              readBase in aGermDict[aChrom][str(refPos)]):
            germs += 1
            if (anIsDebug):
                logging.debug("mismatch_counts() germline found:  " +
                              "aRefIndex=%s, aQueryIndex=%s, offset=%s, " +
                              "chrom=%s, pos=%s, orgBase=%s, orgRef=%s, " +
                              "base=%s, ref=%s, txStrand=%s, r1=%s, r2=%s, " +
                              "rev=%s, mateRev=%s",
                              aRefIndex, aQueryIndex, offset, aChrom, refPos,
                              orgReadBase, orgRefBase, readBase, refBase,
                              aTxStrand, anAlignedRead.is_read1,
                              anAlignedRead.is_read2, anAlignedRead.is_reverse,
                              anAlignedRead.mate_is_reverse)
        elif (aChrom in aTxGermDict and
              str(refPos) in aTxGermDict[aChrom] and
              readBase in aTxGermDict[aChrom][str(refPos)]):
            germs += 1
            if (anIsDebug):
                logging.debug("mismatch_counts() transcript germline found: " +
                              "aRefIndex=%s, aQueryIndex=%s, offset=%s, " +
                              "chrom=%s, pos=%s, orgBase=%s, orgRef=%s, " +
                              "base=%s, ref=%s, germlineAlts=%s, " +
                              "txStrand=%s, r1=%s, r2=%s, " +
                              "rev=%s, mateRev=%s", aRefIndex,
                              aQueryIndex, offset, aChrom, refPos, orgReadBase,
                              orgRefBase, readBase, refBase,
                              aTxGermDict[aChrom][str(refPos)],
                              aTxStrand, anAlignedRead.is_read1,
                              anAlignedRead.is_read2, anAlignedRead.is_reverse,
                              anAlignedRead.mate_is_reverse)
        # if this is an A>G RNA editing event which occur in clusters
        # use the original stranded bases instead of the genomic bases
        elif (aMutSS == "4" and orgRefBase == "A" and orgReadBase == "G"):
            edits += 1
            if (anIsDebug):
                logging.debug("mismatch_counts() A>G editing found: " +
                              "aRefIndex=%s, aQueryIndex=%s, offset=%s, " +
                              "chrom=%s, pos=%s, orgBase=%s, orgRef=%s, " +
                              "base=%s, ref=%s, txStrand=%s, r1=%s, r2=%s, " +
                              "rev=%s, mateRev=%s", aRefIndex, aQueryIndex,
                              offset, aChrom, refPos, orgReadBase, orgRefBase,
                              readBase, refBase, aTxStrand,
                              anAlignedRead.is_read1, anAlignedRead.is_read2,
                              anAlignedRead.is_reverse,
                              anAlignedRead.mate_is_reverse)
        else:
            muts += 1
            if (anIsDebug):
                logging.debug("mismatch_counts() mutation found:  " +
                              "aRefIndex=%s, aQueryIndex=%s, offset=%s, " +
                              "chrom=%s, pos=%s, orgBase=%s, orgRef=%s, " +
                              "base=%s, ref=%s, txStrand=%s, r1=%s, r2=%s, " +
                              "rev=%s, mateRev=%s", aRefIndex, aQueryIndex,
                              offset, aChrom, refPos, orgReadBase, orgRefBase,
                              readBase, refBase, aTxStrand,
                              anAlignedRead.is_read1, anAlignedRead.is_read2,
                              anAlignedRead.is_reverse,
                              anAlignedRead.mate_is_reverse)

        # as soon as the muts count is >= the maxMutsPerRead, return
        if (muts >= aParamsDict["maxMutsPerRead"]):
            return refs, germs, muts, edits

    return refs, germs, muts, edits


class Club():
    def __init__(self, aVCFFilename, aTxNameTag, aTxCoordinateTag, anIsDebug):

        (self.dnaNormBamFilename,
         self.rnaNormBamFilename, self.rnaNormFastaFilename,
         self.rnaNormChrPrefix, self.dnaTumBamFilename,
         self.dnaTumFastaFilename, self.dnaTumChrPrefix,
         self.rnaTumBamFilename, self.rnaTumFastaFilename,
         self.rnaTumChrPrefix, self.germDict, self.txGermDict,
         self.mutationsDict, self.filterDict) = parse_vcf(aVCFFilename,
                                                          aTxNameTag,
                                                          aTxCoordinateTag,
                                                          anIsDebug)

        if (self.dnaTumBamFilename is not None):
            if (not os.path.isfile(self.dnaTumBamFilename)):
                logging.error("The BAM file specified in the header does " +
                              "not exist: %s", self.dnaTumBamFilename)
                sys.exit(1)

            if (self.dnaTumFastaFilename is None or
                not os.path.isfile(self.dnaTumFastaFilename)):
                logging.error("The FASTA file specified in the header " +
                              "does not exist: %s", self.dnaTumFastaFilename)
                sys.exit(1)

            self.dnaTumBamFile = pysam.Samfile(self.dnaTumBamFilename, 'rb')
            self.dnaTumFastaFile = pysam.Fastafile(self.dnaTumFastaFilename)

        if (self.rnaNormBamFilename is not None):
            if (not os.path.isfile(self.rnaNormBamFilename)):
                logging.error("The BAM file specified in the header " +
                              "does not exist: %s", self.rnaNormBamFilename)
                sys.exit(1)

            if (self.rnaNormFastaFilename is None or
                not os.path.isfile(self.rnaNormFastaFilename)):
                logging.error("The FASTA file specified in the header " +
                              "does not exist: %s", self.rnaNormFastaFilename)
                sys.exit(1)

            self.rnaNormBamFile = pysam.Samfile(self.rnaNormBamFilename, 'rb')
            self.rnaNormFastaFile = pysam.Fastafile(self.rnaNormFastaFilename)

        if (self.rnaTumBamFilename is not None):
            if (not os.path.isfile(self.rnaTumBamFilename)):
                logging.error("The BAM file specified in the header " +
                              "does not exist: %s", self.rnaTumBamFilename)
                sys.exit(1)

            if (self.rnaTumFastaFilename is None or
                not os.path.isfile(self.rnaTumFastaFilename)):
                logging.error("The FASTA file specified in the header " +
                              "does not exist: %s", self.rnaTumFastaFilename)
                sys.exit(1)

            self.rnaTumBamFile = pysam.Samfile(self.rnaTumBamFilename, 'rb')
            self.rnaTumFastaFile = pysam.Fastafile(self.rnaTumFastaFilename)

        # initialize the factorial list
        self.factLogList = collections.defaultdict(float)
        self.maxSize = 0
        self.init_factorial(1000)

        return

    def rev_comp_nucleotide(self, aNucleotide):
        '''
        ' This function returns the reverse complement
        ' of the parameter aNucleotide
        '
        ' aNucleotide:    A nucleotide to reverse complement
        '''
        if aNucleotide in i_reverseCompDict:
            return i_reverseCompDict[aNucleotide]
        else:
            logging.error("Trying to reverse complement an "
                          "unknown nucleotide: %s", aNucleotide)
            sys.exit(1)
            return None

    def group_reads_by_name(self, aChromList, aPosList, aTxStrandList,
                            anAlleleSet, aBamOrigin, aParamsDict,
                            aMutSS, aMutType, anIsDebug):

        readsDict = collections.defaultdict(list)

        # loop through all of the transcripts
        for (chrom, pos, strand) in izip(aChromList,
                                         list(map(int, aPosList)),
                                         aTxStrandList):

            if (aMutSS == "Somatic" or aMutSS == "2"):
                if (aBamOrigin == "RNA"):
                    bamFile = self.rnaTumBamFile
                    fastaFile = self.rnaTumFastaFile
                    if (self.rnaTumChrPrefix == "True"):
                        chrom = "chr" + chrom
                else:
                    bamFile = self.dnaTumBamFile
                    fastaFile = self.dnaTumFastaFile
                    if (self.dnaTumChrPrefix == "True"):
                        chrom = "chr" + chrom
            elif (aMutSS == "4"):
                if (aMutType == "NOR_EDIT"):
                    bamFile = self.rnaNormBamFile
                    fastaFile = self.rnaNormFastaFile
                    if (self.rnaNormChrPrefix == "True"):
                        chrom = "chr" + chrom
                elif (aMutType == "TUM_EDIT"):
                    bamFile = self.rnaTumBamFile
                    fastaFile = self.rnaTumFastaFile
                    if (self.rnaTumChrPrefix == "True"):
                        chrom = "chr" + chrom

            # vcfs are 1-based and pysam requires a 0-based coordinate
            pos = pos - 1

            # get the reference base from the fasta
            refBase = fastaFile.fetch(chrom, pos, pos+1).upper()

            if (anIsDebug):
                logging.debug("getting pileups for chrom=%s, pos=%s, " +
                              "strand=%s", chrom, pos, strand)

            # get the pileups
            for pileupColumn in bamFile.pileup(chrom, pos, pos+1,
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
                        not aParamsDict["rnaIncludeSecondaryAlignments"]):
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
                    # there could be more than 2 reads
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
                    # oneReadDict["flag"] = alignedRead.flag
                    # oneReadDict["rname"] = alignedRead.reference_name
                    oneReadDict["start"] = alignedRead.reference_start
                    # oneReadDict["mapQual"] = alignedRead.mapping_quality
                    # oneReadDict["cigar"] = alignedRead.cigar
                    # oneReadDict["mateName"] = alignedRead.next_reference_name
                    oneReadDict["mateStart"] = alignedRead.next_reference_start
                    # oneReadDict["insertSize"] = alignedRead.template_length
                    # oneReadDict["sequence"] = alignedRead.seq
                    # oneReadDict["qualities"] = alignedRead.qual
                    oneReadDict["qlen"] = len(alignedRead.query_sequence)
                    oneReadDict["base"] = alignedRead.seq[qPos]
                    oneReadDict["baseQual"] = alignedRead.qual[qPos]
                    oneReadDict["sequenceIndex"] = pileupRead.query_position
                    oneReadDict["read1"] = alignedRead.is_read1
                    oneReadDict["read2"] = alignedRead.is_read2
                    oneReadDict["reverse"] = alignedRead.is_reverse
                    oneReadDict["mateReverse"] = alignedRead.mate_is_reverse
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
                        strandedBase = self.rev_comp_nucleotide(strandedBase)
                    oneReadDict["strandedBase"] = strandedBase

                    '''
                    if (anIsDebug):
                        logging.debug("name=%s, oneReadDict=%s",
                                      alignedRead.query_name, oneReadDict)
                    '''

                    if (strandedBase in anAlleleSet):
                        if (anIsDebug):
                            logging.debug("strandedBase=%s, is in " +
                                          "alleleSet=%s: name=%s " +
                                          "oneReadDict=%s",
                                          strandedBase, anAlleleSet,
                                          alignedRead.query_name, oneReadDict)

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

    def find_non_overlapping_reads(self, aReadsDict, aMinBaseQual, anIsDebug):
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

                    # if the reads overlap
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
                            logging.debug("base=%s, qual=%s, " +
                                          "maxBaseQual=%s, " +
                                          "maxBaseQualIndex=%s, " +
                                          "basesSet=%s", base, str(qual),
                                          str(maxBaseQual),
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

    def is_perfect(self, aPileupRead, aBamOrigin, aFastaFile, aMutSS,
                   aReadDict, aParamsDict, anIsDebug):

        # check to see if this read is 'perfect'
        alignedRead = aPileupRead.alignment

        # MapSplice doesn't set the proper pair flag for RNA-Seq reads,
        # so only do this for DNA reads
        if (aBamOrigin == "DNA" and not alignedRead.is_proper_pair):
            return False, "improperpair"
        # if it has an insertion, it's not perfect
        if ("I" in alignedRead.cigarstring):
            return False, "insertion"
        # if it has a deletion, it's not perfect
        if ("D" in alignedRead.cigarstring):
            return False, "deletion"
        # if the neighbor bases have low quals, it's not perfect
        if low_neighbor_base_quals(aPileupRead, aParamsDict, anIsDebug):
            return False, "neighborbasequals"

        # alignedRead.pos or alignedRead.reference_start is
        # the 0-based leftmost coordinate
        refIndex = alignedRead.pos
        # alignedRead.qstart or alignedRead.query_alignment_start is
        # the start index of the aligned query portion of the
        # sequence (0-based, inclusive) i.e. the index in the query
        # sequence where the first ref is consumed
        queryIndex = 0

        refCount = 0
        mutCount = 0
        germCount = 0
        editCount = 0
        softClippedCount = 0

        # Op  Description                                Consumes  Consumes
        #                                                Query     Ref
        # M   alignment match (match or mismatch)        yes       yes
        # I   insertion to the ref                       yes       no
        # D   deletion from the ref                      no        yes
        # N   skipped region from the ref                no        yes
        # S   soft-clipping (clipped seqs not in query)  yes       no
        # H   hard-clipping (clipped seqs not in query)  no        no
        # P   padding (silent deletion from padded ref)  no        no
        # =   sequence match                             yes       yes
        # X   sequence mismatch                          yes       yes

        for cigarTuple in alignedRead.cigar:

            (cigarOp, cigarNum) = cigarTuple

            '''
            if (anIsDebug):
                logging.debug("is_perfect() alignedread.cigar=%s, " +
                              "cigarstring=%s, cigarOpCode=%s, cigarOp=%s, " +
                              "cigarNum=%s", alignedRead.cigar,
                              alignedRead.cigarstring, cigarOp,
                              i_cigarDict[cigarOp], cigarNum)
                logging.debug("refIndex=%s, queryIndex=%s",
                              refIndex, queryIndex)
            '''
            '''
            # this i_cigarDict deciphers the cigarOp integer code
            i_cigarDict[0] = "match"
            i_cigarDict[1] = "insertion"
            i_cigarDict[2] = "deletion"
            i_cigarDict[3] = "refskipped"
            i_cigarDict[4] = "softclipped"
            i_cigarDict[5] = "hardclipped"
            i_cigarDict[6] = "padding"
            i_cigarDict[7] = "seqmatch"
            i_cigarDict[8] = "seqmismatch"
            '''

            # if it aligns
            if i_cigarDict[cigarOp] == "match":
                # count the number of mismatches across the
                # aligned portion of the read
                refs, germs, muts, edits = mismatch_counts(cigarNum,
                                                           aReadDict["chrom"],
                                                           refIndex,
                                                           queryIndex,
                                                           aMutSS,
                                                           aReadDict["strand"],
                                                           alignedRead,
                                                           self.germDict,
                                                           self.txGermDict,
                                                           aFastaFile,
                                                           aBamOrigin,
                                                           aParamsDict,
                                                           anIsDebug)
                refCount += refs
                germCount += germs
                mutCount += muts
                editCount += edits

                refIndex += cigarNum
                queryIndex += cigarNum

            elif i_cigarDict[cigarOp] == "insertion":
                queryIndex += cigarNum
            elif i_cigarDict[cigarOp] == "deletion":
                refIndex += cigarNum
            elif i_cigarDict[cigarOp] == "refskipped":
                refIndex += cigarNum
            elif i_cigarDict[cigarOp] == "softclipped":
                queryIndex += cigarNum
                softClippedCount += cigarNum
            elif i_cigarDict[cigarOp] == "hardclipped":
                continue
            elif i_cigarDict[cigarOp] == "padded":
                continue
            # the base in the read matches the reference
            elif i_cigarDict[cigarOp] == "seqmatch":
                refIndex += cigarNum
                queryIndex += cigarNum
                refCount += cigarNum
            # the base in the read does not match the reference
            elif i_cigarDict[cigarOp] == "seqmismatch":
                # count the number of mismatches across the
                # aligned portion of the read
                refs, germs, muts, edits = mismatch_counts(cigarNum,
                                                           aReadDict["chrom"],
                                                           refIndex,
                                                           queryIndex,
                                                           aMutSS,
                                                           aReadDict["strand"],
                                                           alignedRead,
                                                           self.germDict,
                                                           self.txGermDict,
                                                           aFastaFile,
                                                           aBamOrigin,
                                                           aParamsDict,
                                                           anIsDebug)
                refCount += refs
                germCount += germs
                mutCount += muts
                editCount += edits

                refIndex += cigarNum
                queryIndex += cigarNum
            else:
                logging.error("Unexpected value in the read cigar string " +
                              "%s at position %s", alignedRead.cigar,
                              alignedRead.tid)

            if (anIsDebug):
                logging.debug("counts for this read: refs=%s, muts=%s, " +
                              "germs=%s, edits=%s, softclipped=%s", refCount,
                              mutCount, germCount, editCount, softClippedCount)

            # if the number of muts in this read exceeds the limit
            # then apply the maxMutsPerRead filter
            if (mutCount >= aParamsDict["maxMutsPerRead"]):
                return False, "maxMuts"

            # softClippedCount is number of soft clipped bases across the read
            # if the percent of soft clipped bases is greater than the param,
            # this read is not perfect
            if (len(alignedRead.query_sequence) > 0):
                readLen = len(alignedRead.query_sequence)
                softClipPct = round(softClippedCount/float(readLen), 2)
                if (softClipPct >= aParamsDict["maxReadSoftClipPct"]):
                    if (anIsDebug):
                        '''
                        logging.debug("softClipped=%s, qlen=%s, " +
                                      "infer_qlen=%s, len(seq)=%s, ",
                                      softClippedCount,
                                      alignedRead.query_length,
                                      alignedRead.infer_query_length(),
                                      len(alignedRead.query_sequence))
                        '''
                        logging.debug("read soft clipped too much, " +
                                      "softClipped=%s, len(seq)=%s, " +
                                      "softClipPct=%s, maxSFP=%s",
                                      softClippedCount,
                                      len(alignedRead.query_sequence),
                                      softClipPct,
                                      aParamsDict["maxReadSoftClipPct"])
                    return False, "maxSoftClips"

        if (anIsDebug):
            logging.debug("this perfect read has refs=%s, muts=%s, " +
                          "germs=%s, edits=%s, softClip=%s", refCount,
                          mutCount, germCount, editCount, softClippedCount)

        return True, "perfect"

    def set_score(self, aCurrData, aScorePassingOnlyFlag, anIsDebug):

        # calculating the score is computationally expensive,
        # so only do it for all passing calls by default
        pvalue = 0.98
        phred = 0
        if ((aScorePassingOnlyFlag and "PASS" in aCurrData.filterSet) or
            (not aScorePassingOnlyFlag)):

            if ("GERM" in aCurrData.infoDict["MT"]):
                for (modType, modChange) in izip(aCurrData.infoDict["MT"],
                                                 aCurrData.infoDict["MC"]):
                    if (anIsDebug):
                        logging.debug("modType=%s, modChange=%s",
                                      modType, modChange)

                    if (modType == "GERM"):
                        ref, alt = modChange.split(">")
                        refIndex = aCurrData.allelesList.index(ref)
                        altIndex = aCurrData.allelesList.index(alt)

                        totalRefReads = 0
                        totalAltReads = 0
                        if (aCurrData.dnaNormalDict is not None):
                            refDp = aCurrData.dnaNormalDict["AD"][refIndex]
                            altDp = aCurrData.dnaNormalDict["AD"][altIndex]
                            if (refDp != "."):
                                totalRefReads += int(refDp)
                            if (altDp != "."):
                                totalAltReads += int(altDp)
                        if (aCurrData.rnaNormalDict is not None):
                            refDp = aCurrData.rnaNormalDict["AD"][refIndex]
                            altDp = aCurrData.rnaNormalDict["AD"][altIndex]
                            if (refDp != "."):
                                totalRefReads += int(refDp)
                            if (altDp != "."):
                                totalAltReads += int(altDp)
                        if (aCurrData.dnaTumorDict is not None):
                            refDp = aCurrData.dnaTumorDict["AD"][refIndex]
                            altDp = aCurrData.dnaTumorDict["AD"][altIndex]
                            if (refDp != "."):
                                totalRefReads += int(refDp)
                            if (altDp != "."):
                                totalAltReads += int(altDp)
                        if (aCurrData.rnaTumorDict is not None):
                            refDp = aCurrData.rnaTumorDict["AD"][refIndex]
                            altDp = aCurrData.rnaTumorDict["AD"][altIndex]
                            if (refDp != "."):
                                totalRefReads += int(refDp)
                            if (altDp != "."):
                                totalAltReads += int(altDp)

                        totalCoverage = totalRefReads + totalAltReads

                        newMaxSize = totalCoverage * 2
                        if (newMaxSize > self.maxSize):
                            self.init_factorial(newMaxSize)

                        pvalue, phred = self.get_score(totalCoverage,
                                                       0,
                                                       totalRefReads,
                                                       totalAltReads,
                                                       self.maxSize,
                                                       self.factLogList,
                                                       i_debug)
                        if (anIsDebug):
                            logging.debug("pval=%s, phred=%s", pvalue, phred)

                        aCurrData.qual = str(phred)
                        aCurrData.infoDict["SSC"] = ["0"]
                        aCurrData.infoDict["PVAL"] = [str(pvalue)]

            elif ("SOM" in aCurrData.infoDict["MT"]):
                for (modType, modChange) in izip(aCurrData.infoDict["MT"],
                                                 aCurrData.infoDict["MC"]):

                    if (anIsDebug):
                        logging.debug("modType=%s, modChange=%s",
                                      modType, modChange)
                    if (modType == "SOM"):
                        ref, alt = modChange.split(">")
                        refIndex = aCurrData.allelesList.index(ref)
                        altIndex = aCurrData.allelesList.index(alt)

                        normalRefReads = 0
                        normalAltReads = 0
                        tumorRefReads = 0
                        tumorAltReads = 0

                        if (aCurrData.dnaNormalDict is not None):
                            normRefDp = aCurrData.dnaNormalDict["AD"][refIndex]
                            normAltDp = aCurrData.dnaNormalDict["AD"][altIndex]
                            if (normRefDp != "."):
                                normalRefReads += int(normRefDp)
                            if (normAltDp != "."):
                                normalAltReads += int(normAltDp)
                        if (aCurrData.rnaNormalDict is not None):
                            normRefDp = aCurrData.rnaNormalDict["AD"][refIndex]
                            normAltDp = aCurrData.rnaNormalDict["AD"][altIndex]
                            if (normRefDp != "."):
                                normalRefReads += int(normRefDp)
                            if (normAltDp != "."):
                                normalAltReads += int(normAltDp)
                        if (aCurrData.dnaTumorDict is not None):
                            tumRefDp = aCurrData.dnaTumorDict["AD"][refIndex]
                            tumAltDp = aCurrData.dnaTumorDict["AD"][altIndex]
                            if (tumRefDp != "."):
                                tumorRefReads += int(tumRefDp)
                            if (tumAltDp != "."):
                                tumorAltReads += int(tumAltDp)
                        if (aCurrData.rnaTumorDict is not None):
                            tumRefDp = aCurrData.rnaTumorDict["AD"][refIndex]
                            tumAltDp = aCurrData.rnaTumorDict["AD"][altIndex]
                            if (tumRefDp != "."):
                                tumorRefReads += int(tumRefDp)
                            if (tumAltDp != "."):
                                tumorAltReads += int(tumAltDp)

                        totalCoverage = (normalRefReads + normalAltReads +
                                         tumorRefReads + tumorAltReads)

                        newMaxSize = totalCoverage
                        if (newMaxSize > self.maxSize):
                            self.init_factorial(newMaxSize)

                        pvalue, phred = self.get_score(normalRefReads,
                                                       normalAltReads,
                                                       tumorRefReads,
                                                       tumorAltReads,
                                                       self.maxSize,
                                                       self.factLogList,
                                                       i_debug)
                        if (anIsDebug):
                            logging.debug("pval=%s, phred=%s", pvalue, phred)

                        aCurrData.qual = str(phred)
                        aCurrData.infoDict["SSC"] = [str(phred)]
                        aCurrData.infoDict["PVAL"] = [str(pvalue)]

            elif ("TUM_EDIT" in aCurrData.infoDict["MT"]):
                for (modType, modChange) in izip(aCurrData.infoDict["MT"],
                                                 aCurrData.infoDict["MC"]):

                    if (anIsDebug):
                        logging.debug("modType=%s, modChange=%s",
                                      modType, modChange)

                    if (modType == "TUM_EDIT"):
                        ref, alt = modChange.split(">")
                        refIndex = aCurrData.allelesList.index(ref)
                        altIndex = aCurrData.allelesList.index(alt)

                        totalRefReads = 0
                        totalAltReads = 0
                        editingRefReads = 0
                        editingAltReads = 0
                        if (aCurrData.dnaNormalDict is not None):
                            refDp = aCurrData.dnaNormalDict["AD"][refIndex]
                            altDp = aCurrData.dnaNormalDict["AD"][altIndex]
                            if (refDp != "."):
                                totalRefReads += int(refDp)
                            if (altDp != "."):
                                totalAltReads += int(altDp)
                        if (aCurrData.rnaNormalDict is not None):
                            refDp = aCurrData.rnaNormalDict["AD"][refIndex]
                            altDp = aCurrData.rnaNormalDict["AD"][altIndex]
                            editRefDp = aCurrData.rnaNormalDict["AD"][refIndex]
                            editAltDp = aCurrData.rnaNormalDict["AD"][altIndex]
                            if (refDp != "."):
                                totalRefReads += int(refDp)
                            if (altDp != "."):
                                totalAltReads += int(altDp)
                            if (editRefDp != "."):
                                editingRefReads += int(editRefDp)
                            if (editAltDp != "."):
                                editingAltReads += int(editAltDp)
                        if (aCurrData.dnaTumorDict is not None):
                            refDp = aCurrData.dnaTumorDict["AD"][refIndex]
                            altDp = aCurrData.dnaTumorDict["AD"][altIndex]
                            if (refDp != "."):
                                totalRefReads += int(refDp)
                            if (altDp != "."):
                                totalAltReads += int(altDp)
                        if (aCurrData.rnaTumorDict is not None):
                            refDp = aCurrData.rnaTumorDict["AD"][refIndex]
                            altDp = aCurrData.rnaTumorDict["AD"][altIndex]
                            editRefDp = aCurrData.rnaTumorDict["AD"][refIndex]
                            editAltDp = aCurrData.rnaTumorDict["AD"][altIndex]
                            if (refDp != "."):
                                totalRefReads += int(refDp)
                            if (altDp != "."):
                                totalAltReads += int(altDp)
                            if (editRefDp != "."):
                                editingRefReads += int(editRefDp)
                            if (editAltDp != "."):
                                editingAltReads += int(editAltDp)

                        totalCoverage = totalRefReads + totalAltReads

                        newMaxSize = (totalCoverage +
                                      editingRefReads +
                                      editingAltReads)
                        if (newMaxSize > self.maxSize):
                            self.init_factorial(newMaxSize)

                        pvalue, phred = self.get_score(totalCoverage,
                                                       0,
                                                       editingRefReads,
                                                       editingAltReads,
                                                       self.maxSize,
                                                       self.factLogList,
                                                       i_debug)
                        if (anIsDebug):
                            logging.debug("pval=%s, phred=%s", pvalue, phred)

                        aCurrData.qual = str(phred)
                        aCurrData.infoDict["SSC"] = ["0"]
                        aCurrData.infoDict["PVAL"] = [str(pvalue)]
            else:
                # these are RNA_TUM_VAR lines that have little or no DNA
                aCurrData.qual = str(phred)
                aCurrData.infoDict["SSC"] = [str(phred)]
                aCurrData.infoDict["PVAL"] = [str(pvalue)]
        else:
            # calculating the score is computationally expensive,
            # so only do it for all passing calls by default
            # these are non-passing lines
            aCurrData.qual = str(phred)
            aCurrData.infoDict["SSC"] = [str(phred)]
            currData.infoDict["PVAL"] = [str(pvalue)]

        if (anIsDebug):
            logging.debug("pos=%s, qual=%s, pval=%s, phred=%s",
                          currData.pos, currData.qual, pvalue, phred)

        return

    def get_score(self, anExpRefCount, anExpMutCount, anObsRefCount,
                  anObsMutCount, aMaxSize, aFactLogList, anIsDebug):

        if (anIsDebug):
            logging.debug("anExpRefCount=%s, anExpMutCount=%s, " +
                          "anObsRefCount=%s, anObsMutCount=%s",
                          anExpRefCount, anExpMutCount,
                          anObsRefCount, anObsMutCount)
            logging.debug("1st R-tail maxSize=%s", aMaxSize)

        pValue = self.get_right_tail(anExpRefCount, anExpMutCount,
                                     anObsRefCount, anObsMutCount,
                                     aMaxSize, aFactLogList, anIsDebug)

        if (anIsDebug):
            logging.debug("after R-tail p-value=%s", pValue)

        if math.isnan(pValue):
            logging.warning("Warning: unable to calculate R-tail p-value: " +
                            "expRef=%s, expAlt=%s, obsRef=%s, obsAlt=%s",
                            anExpRefCount, anExpMutCount,
                            anObsRefCount, anObsMutCount)

        # If p-value is 1, do left-sided test
        if(pValue >= 0.999):

            if (anIsDebug):
                logging.debug("1st L-tail maxSize=%s", aMaxSize)

            pValue = self.get_left_tail(anExpRefCount, anExpMutCount,
                                        anObsRefCount, anObsMutCount,
                                        aMaxSize, aFactLogList, anIsDebug)

            if math.isnan(pValue):
                logging.warning("Warning: unable to calculate L-tail " +
                                "p-value: expRef=%s, expAlt=%s, obsRef=%s, " +
                                "obsAlt=%s", anExpRefCount, anExpMutCount,
                                anObsRefCount, anObsMutCount)

        # if pValue is nan, then an error occurred
        # warning messages were written
        if math.isnan(pValue):
            pValue = 1

        # calculate the phred score
        if (pValue > 0):
            phred = 0 - (10 * math.log(pValue, 10))
        elif (pValue < 0):
            # pValue = 1 - pValue
            # phred = 0 - (10 * math.log(pValue, 10))
            phred = 0
        # the pvalue is 0
        else:
            phred = 255

        if (phred > 255):
            phred = 255

        return '{0:1.2e}'.format(pValue), int(round(phred))

    def get_pvalue(self, a, b, c, d, aFactLogList, anIsDebug):
        n = a + b + c + d

        if (anIsDebug):
            logging.debug("a=%s, b=%s, c=%s, d=%s, n=%s, a+b=%s, c+d=%s, " +
                          "a+c=%s, b+d=%s", a, b, c, d, n, a+b, c+d, a+c, b+d)
            logging.debug("f[a+b]=%s, f[c+d]=%s, f[a+c]=%s, f[b+d]=%s, add=%s",
                          aFactLogList[a+b], aFactLogList[c+d],
                          aFactLogList[a+c], aFactLogList[b+d],
                          (aFactLogList[a + b] + aFactLogList[c + d] +
                           aFactLogList[a + c] + aFactLogList[b + d]))
            logging.debug("f[a]=%s, f[b]=%s, f[c]=%s, f[d]=%s, f[n]=%s, " +
                          "minuses=%s", aFactLogList[a], aFactLogList[b],
                          aFactLogList[c], aFactLogList[d], aFactLogList[n],
                          (aFactLogList[a] + aFactLogList[b] +
                           aFactLogList[c] + aFactLogList[d] +
                           aFactLogList[n]))

        p = ((aFactLogList[a + b] + aFactLogList[c + d] +
              aFactLogList[a + c] + aFactLogList[b + d]) -
             (aFactLogList[a] + aFactLogList[b] + aFactLogList[c] +
              aFactLogList[d] + aFactLogList[n]))

        if (anIsDebug):
            logging.debug("factorial=%s, pvalue=%s", p, math.exp(p))

        return math.exp(p)

    def get_right_tail(self, a, b, c, d, aMaxSize, aFactLogList, anIsDebug):
        n = a + b + c + d
        if (n > aMaxSize):
            logging.warning("Problem calculating R-tail: n=%s > maxSize=%s",
                            n, aMaxSize)
            return float("nan")

        # get the first p-value
        p = self.get_pvalue(a, b, c, d, aFactLogList, anIsDebug)

        if (anIsDebug):
            logging.debug("doing R-tail: p=%s, a=%s, b=%s, c=%s, d=%s",
                          p, a, b, c, d)

        if (c < b):
            minRight = c
        else:
            minRight = b

        for i in xrange(0, minRight):
            # if we've reached the pvalue giving us the
            # max phred score of 255, just return
            if (p < 0.00000000000000000000000001):
                return p

            a += 1
            b -= 1
            c -= 1
            d += 1
            pTemp = self.get_pvalue(a, b, c, d, aFactLogList, anIsDebug)
            p += pTemp

            if (anIsDebug):
                logging.debug("doing round %s:", i)
                logging.debug("\tpTemp = %s", pTemp)
                logging.debug("\ta=%s, b=%s, c=%s, d=%s",  a, b, c, d)
        return p

    def get_left_tail(self, a, b, c, d, aMaxSize, aFactLogList, anIsDebug):
        n = a + b + c + d
        if (n > aMaxSize):
            logging.warning("Problem calculating L-tail: n=%s > maxSize=%s",
                            n, aMaxSize)
            return float("nan")

        # get the first p-value
        p = self.get_pvalue(a, b, c, d, aFactLogList, anIsDebug)

        if (anIsDebug):
            logging.debug("doing L-tail: p=%s, a=%s, b=%s, c=%s, d=%s",
                          p, a, b, c, d)

        if (a < d):
            minLeft = a
        else:
            minLeft = d

        for i in xrange(0, minLeft):
            # if we've reached the pvalue giving us the
            # max phred score of 255, just return
            if (p < 0.00000000000000000000000001):
                logging.warning("L-tail breaking out of round %s of %s",
                                i, minLeft)
                return p

            a -= 1
            b += 1
            c += 1
            d -= 1
            pTemp = self.get_pvalue(a, b, c, d, aFactLogList, anIsDebug)
            p += pTemp

            if (anIsDebug):
                logging.debug("doing round %s:", i)
                logging.debug("\tpTemp=%s", pTemp)
                logging.debug("\ta=%s, b=%s, c=%s, d=%s",  a, b, c, d)

        return p

    def init_factorial(self, aNewMaxSize):

        if (self.maxSize == 0):
            self.factLogList[0] = 0.0
            self.maxSize = 1

        for index in xrange(self.maxSize, aNewMaxSize+1):
            self.factLogList[index] = (self.factLogList[index - 1] +
                                       math.log(index))

        self.maxSize = aNewMaxSize

        return

    def parse_info_field(self, anInfoField):

        # parse the info field and create a dict
        infoDict = collections.defaultdict(list)
        infoFieldList = anInfoField.split(";")
        for info in infoFieldList:
            keyValueList = info.split("=")
            # some keys are just singular without a value (e.g. DB, etc.)
            if (len(keyValueList) == 1):
                infoDict[keyValueList[0]] = ["True"]
            else:
                # the value can be a comma separated list
                infoDict[keyValueList[0]] = keyValueList[1].split(",")

        return infoDict

    def filter_by_read_support(self, aCurrData, aTxNameTag, aTxCoordinateTag,
                               aTxStrandTag, aMutType, aModChange, aBamOrigin,
                               aParamsDict, anIsDebug):

        # expects a 0-based pos
        refCount = 0
        mutCountReads = 0
        mutCountQualReads = 0
        lowQualCount = 0
        starts = 0
        middles = 0
        ends = 0
        numPerfect = 0
        perfectStarts = 0
        perfectMiddles = 0
        perfectEnds = 0
        perfectForStrand = 0
        perfectRevStrand = 0
        readsWithMaxMuts = 0
        readsWithMaxSoftClips = 0
        readsWithDels = 0
        readsWithIns = 0
        readsWithNeighborBaseQuals = 0
        readsWithImproperPairs = 0
        mmpDict = dict()

        bases = ""
        baseQuals = ""

        # if we should use the transcript information
        if (aBamOrigin == "RNA" and
            aTxNameTag is not None and
            aTxCoordinateTag is not None and
            aTxStrandTag is not None and
            aCurrData.infoDict[aTxNameTag] is not None and
            currData.infoDict[aTxCoordinateTag] is not None and
            currData.infoDict[aTxStrandTag] is not None):
            chromList = currData.infoDict[aTxNameTag]
            posList = currData.infoDict[aTxCoordinateTag]
            strandList = currData.infoDict[aTxStrandTag]
        # if we should use the genomic information
        else:
            chromList = [aCurrData.chrom]
            posList = [aCurrData.pos]
            strandList = [None]

        mutSS = currData.infoDict["SS"][0]

        if (anIsDebug):
            logging.debug("begin filter for %s and %s and %s, mutSS=%s, " +
                          "mutType=%s, origin=%s", chromList, posList,
                          strandList, mutSS, aMutType, aBamOrigin)
            logging.debug("parmsDict: %s", aParamsDict)

        if (mutSS == "Somatic" or mutSS == "2"):
            if (aBamOrigin == "RNA"):
                fastaFile = self.rnaTumFastaFile
                if "MMP" not in currData.rnaTumorDict:
                    currData.rnaTumorDict["MMP"] = list()
                mmpList = currData.rnaTumorDict["MMP"]
            else:
                fastaFile = self.dnaTumFastaFile
                if "MMP" not in currData.dnaTumorDict:
                    currData.dnaTumorDict["MMP"] = list()
                mmpList = currData.dnaTumorDict["MMP"]
        elif (mutSS == "4"):
            if (aMutType == "TUM_EDIT"):
                fastaFile = self.rnaTumFastaFile
                if "MMP" not in currData.rnaTumorDict:
                    currData.rnaTumorDict["MMP"] = list()
                mmpList = currData.rnaTumorDict["MMP"]
            elif (aMutType == "NOR_EDIT"):
                fastaFile = self.rnaNormFastaFile
                if "MMP" not in currData.rnaNormalDict:
                    currData.rnaNormalDict["MMP"] = list()
                mmpList = currData.rnaNormalDict["MMP"]

        # add the previous and next reference base to the INFO
        chrom = chromList[0]
        # vcfs are 1-based and pysam requires a 0-based coordinate
        pos = int(posList[0]) - 1
        prevRefBase = fastaFile.fetch(chrom, pos-1, pos).upper()
        # refBase = fastaFile.fetch(chrom, pos, pos+1).upper()
        nextRefBase = fastaFile.fetch(chrom, pos+1, pos+2).upper()
        currData.infoDict["PN"] = [prevRefBase]
        currData.infoDict["NN"] = [nextRefBase]

        # get the set of alleles for this call
        alleleSet = set()
        (ref, alt) = modChange.split(">")
        alleleSet.add(ref)
        alleleSet.add(alt)

        # group all of the reads by name
        readsDict = self.group_reads_by_name(chromList, posList, strandList,
                                             alleleSet, aBamOrigin,
                                             aParamsDict, mutSS,
                                             aMutType, anIsDebug)

        # get all of the non-overlapping reads
        minBaseQual = aParamsDict["minBaseQual"]
        nonOverlapReadsList = self.find_non_overlapping_reads(readsDict,
                                                              minBaseQual,
                                                              anIsDebug)

        # loop over the non-overlapping reads
        for readDict in nonOverlapReadsList:
            alignedRead = readDict["alignedRead"]
            pileupRead = readDict["pileupRead"]

            readBase = readDict["base"]
            readBaseQual = readDict["baseQual"]
            # readSequence = readDict["sequence"]
            # readSeqIndex = readDict["sequenceIndex"]
            # readName = readDict["queryName"]
            # readMapQual = readDict["mapQual"]
            # readInsertSize = readDict["insertSize"]
            # readFlag = readDict["flag"]

            '''
            if (anIsDebug):
                logging.debug("found aligned read at: %s:%s = %s",
                              readDict["chrom"], readDict["pos"], alignedRead)
            '''

            # if this read supports the alternative allele
            if is_mutation(readDict, ref, alt, aBamOrigin, mmpDict, anIsDebug):

                '''
                if (anIsDebug):
                    logging.debug("found aligned read supporting variant " +
                                  "at: %s:%s = %s", readDict["chrom"],
                                  readDict["pos"], alignedRead)
                '''

                mutCountReads += 1

                if (anIsDebug):
                    logging.debug("found read with mutation, number of " +
                                  "reads with mutations=%s", mutCountReads)

                # if the base and map quals are good enough
                if (not low_base_or_map_quals(pileupRead,
                                              aParamsDict,
                                              anIsDebug)):

                    mutCountQualReads += 1

                    bases += readBase
                    baseQuals += readBaseQual

                    if (anIsDebug):
                        logging.debug("mutation read passed base and map " +
                                      "quals, qualReads=%s", mutCountQualReads)

                    # count the positions in the reads
                    readLength = len(alignedRead.query_sequence)
                    if (pileupRead.query_position/float(readLength) <= 0.33):
                        starts += 1
                    elif (pileupRead.query_position/float(readLength) <= 0.66):
                        middles += 1
                    else:
                        ends += 1

                    if (anIsDebug):
                        logging.debug("pbias query_pos=%s, readLength=%s, " +
                                      "starts=%s, middles=%s, ends=%s",
                                      pileupRead.query_position, readLength,
                                      starts, middles, ends)

                    # see if this read is perfect
                    isPerfectFlag, reasonNotPerf = self.is_perfect(pileupRead,
                                                                   aBamOrigin,
                                                                   fastaFile,
                                                                   mutSS,
                                                                   readDict,
                                                                   aParamsDict,
                                                                   anIsDebug)

                    if (anIsDebug):
                        logging.debug("isPerfect?=%s, reason=%s",
                                      isPerfectFlag, reasonNotPerf)

                    # keep track of the counts
                    if reasonNotPerf == "improperpair":
                        readsWithImproperPairs += 1
                    elif reasonNotPerf == "insertion":
                        readsWithIns += 1
                    elif reasonNotPerf == "deletion":
                        readsWithDels += 1
                    elif reasonNotPerf == "neighborbasequals":
                        readsWithNeighborBaseQuals += 1
                    elif reasonNotPerf == "maxMuts":
                        readsWithMaxMuts += 1
                    elif reasonNotPerf == "maxSoftClips":
                        readsWithMaxSoftClips += 1

                    # if we found a perfect read
                    if isPerfectFlag:

                        numPerfect += 1

                        # determine the strand
                        if alignedRead.is_reverse:
                            perfectRevStrand += 1
                        else:
                            perfectForStrand += 1

                        # count the perfect positions in the reads
                        readLen = len(alignedRead.query_sequence)
                        qPos = pileupRead.query_position
                        if (qPos/float(readLen) <= 0.33):
                            perfectStarts += 1
                        elif (qPos/float(readLen) <= 0.66):
                            perfectMiddles += 1
                        else:
                            perfectEnds += 1

                        if (anIsDebug):
                            logging.debug("perfectpbias: query_pos=%s, " +
                                          "readLength=%s, perfectStarts=%s, " +
                                          "perfectMiddles=%s, perfectEnds=%s",
                                          pileupRead.query_position,
                                          readLength,
                                          perfectStarts,
                                          perfectMiddles,
                                          perfectEnds)

                    if (anIsDebug):
                        logging.debug("perfect read counts for %s:%s, " +
                                      "mutSS=%s, mutType=%s, numPerfect=%s",
                                      readDict["chrom"], readDict["pos"],
                                      mutSS, aMutType, numPerfect)
                else:
                    lowQualCount += 1
            else:
                refCount += 1
                bases += readBase
                baseQuals += readBaseQual

        if (anIsDebug):
            logging.debug("pysam bases=%s", bases)
            logging.debug("pysam quals=%s", baseQuals)
            logging.debug("final at pos %s:%s, mutSS=%s, mutType=%s, " +
                          "refCt=%s, mutCountReads=%s, mutCountQualReads=%s",
                          chromList, posList, mutSS, aMutType, refCount,
                          mutCountReads, mutCountQualReads)
            logging.debug("final at pos %s:%s, mutSS=%s, mutType=%s, " +
                          "lowQualCount=%s, improperPairs=%s, ins=%s, " +
                          "dels=%s, neighborBQ=%s, maxMuts=%s, " +
                          "maxSoftClips=%s, numPerfect=%s", chromList,
                          posList, mutSS, aMutType, lowQualCount,
                          readsWithImproperPairs, readsWithIns, readsWithDels,
                          readsWithNeighborBaseQuals, readsWithMaxMuts,
                          readsWithMaxSoftClips, numPerfect)

        # keep track of all filters
        filters = set()

        # all of the counts have been done for the genomic forward strand
        # so we can use aModChange here to access the MMP
        multiMappingPct = 0.0
        # multiple mapping of the mutCountReads
        (source, target) = aModChange.split(">")
        if ((target in mmpDict) and
            (mmpDict[target]["total"] > 0) and
            (mmpDict[target]["total"] >= aParamsDict["minMultiMapDepth"])):
            # get the multi mapping percent
            secondary = mmpDict[target]["secondary"]
            total = mmpDict[target]["total"]
            multiMappingPct = round(secondary/float(total), 2)
            if (multiMappingPct >= aParamsDict["maxMultiMapPct"]):
                filters.add("mxmmp")
                if (anIsDebug):
                    logging.debug("checkfilter multiMap applied for %s:%s, " +
                                  "mutSS=%s, mutType=%s, " +
                                  "secondaryTargetReads=%s, " +
                                  "totalTargetReads=%s, mmp=%s",
                                  chromList, posList, mutSS, aMutType,
                                  mmpDict[target]["secondary"],
                                  mmpDict[target]["total"], multiMappingPct)
        elif (anIsDebug and target in mmpDict):
            logging.debug("checkfilter multiMap no minDepth for %s:%s, " +
                          "mutSS=%s, mutType=%s, secondaryTargetReads=%s, " +
                          "totalTargetReads=%s", chromList, posList, mutSS,
                          aMutType, mmpDict[target]["secondary"],
                          mmpDict[target]["total"])
        elif (anIsDebug):
            logging.debug("checkfilter multiMap no minDepth for %s:%s, " +
                          "mutSS=%s, mutType=%s", chromList, posList,
                          mutSS, aMutType)

        # only apply strand bias to the perfect reads with
        # mutations if we have enough reads
        if ((numPerfect > 0) and
            (numPerfect >= aParamsDict["minStrandBiasDepth"])):
            sbias = round(perfectForStrand/float(numPerfect), 2)
            # if sbias < 0.1 or sbias > 0.9:
            if ((sbias > (aParamsDict["maxStrandBias"])) or
                (sbias < (1.0 - aParamsDict["maxStrandBias"]))):
                filters.add("perfsbias")
                if (anIsDebug):
                    logging.debug("checkfilter perfsbias for %s:%s, " +
                                  "numPerfect=%s, filters=%s, forstrand=%s, " +
                                  "revstrand=%s, sbias=%s", chromList, posList,
                                  numPerfect, filters, perfectForStrand,
                                  perfectRevStrand, sbias)
        elif (anIsDebug):
                logging.debug("checkfilter sbias no minDepth for %s:%s, " +
                              "numPerfect=%s", chromList, posList, numPerfect)

        # only apply positional bias to the reads with
        # mutations if we have enough reads
        if ((mutCountQualReads > 0) and
            (mutCountQualReads >= aParamsDict["minPositionBiasDepth"])):
            pbiasStarts = round(starts/float(mutCountQualReads), 2)
            pbiasEnds = round(ends/float(mutCountQualReads), 2)
            '''
            if (anIsDebug):
                logging.debug("mutCountQualReads=%s, pbiasStarts=%s, " +
                              "pbiasEnds=%s", mutCountQualReads, pbiasStarts,
                              pbiasEnds)
            '''
            if (pbiasStarts >= aParamsDict["maxPositionBias"]):
                filters.add("pbias")
                '''
                if (anIsDebug):
                    logging.debug("pbias from starts starts=%s, middles=%s, " +
                                  "ends=%s, mutCountQualReads=%s", starts,
                                  middles, ends, mutCountQualReads)
                '''
            elif (pbiasEnds >= aParamsDict["maxPositionBias"]):
                filters.add("pbias")
                '''
                if (anIsDebug):
                    logging.debug("pbias from ends starts=%s, middles=%s, " +
                                  "ends=%s, mutCountQualReads=%s", starts,
                                  middles, ends, mutCountQualReads)
                '''
            if (anIsDebug and "pbias" in filters):
                logging.debug("checkfilter pbias for %s:%s, " +
                              "mutCountQualReads=%s, filters=%s, " +
                              "starts=%s (%s), middles=%s (%s), ends=%s (%s)",
                              chromList, posList, mutCountQualReads, filters,
                              starts, starts/float(mutCountQualReads),
                              middles, middles/float(mutCountQualReads),
                              ends, ends/float(mutCountQualReads))
        elif (anIsDebug):
            logging.debug("checkfilter pbias no minDepth for %s:%s, " +
                          "mutCountQualReads=%s", chromList, posList,
                          mutCountQualReads)

        # only apply positional bias to the perfect reads with
        # mutations if we have enough reads
        if ((numPerfect > 0) and
            (numPerfect >= aParamsDict["minPositionBiasDepth"])):
            perfPbiasStarts = round(perfectStarts/float(numPerfect), 2)
            perfPbiasEnds = round(perfectEnds/float(numPerfect), 2)
            '''
            if (anIsDebug):
                logging.debug("numPerfect=%s, perfectpbiasStarts=%s, " +
                              "perfectpbiasEnds=%s", numPerfect,
                              perfPbiasStarts, perfPbiasEnds)
            '''
            if (perfPbiasStarts >= aParamsDict["maxPositionBias"]):
                filters.add("perfpbias")
                '''
                if (anIsDebug):
                    logging.debug("perfectpbias from starts " +
                                  "perfectStarts=%s, perfectMiddles=%s, " +
                                  "perfectEnds=%s, numPerfect=%s",
                                  perfectStarts, perfectMiddles,
                                  perfectEnds, numPerfect)
                '''
            elif (perfPbiasEnds >= aParamsDict["maxPositionBias"]):
                filters.add("perfpbias")
                '''
                if (anIsDebug):
                    logging.debug("perfectpbias from ends perfectStarts=%s, " +
                                  "perfectMiddles=%s, perfectEnds=%s, " +
                                  "numPerfect=%s", perfectStarts,
                                  perfectMiddles, perfectEnds, numPerfect)
                '''
            if (anIsDebug and "perfpbias" in filters):
                logging.debug("checkfilter perfpbias for %s:%s, " +
                              "numPerfect=%s, filters=%s, " +
                              "perfectStarts=%s (%s), " +
                              "perfectMiddles=%s (%s), " +
                              "perfectEnds=%s (%s)", chromList, posList,
                              numPerfect, filters,
                              perfectStarts, perfectStarts/float(numPerfect),
                              perfectMiddles, perfectMiddles/float(numPerfect),
                              perfectEnds, perfectEnds/float(numPerfect))
        elif (anIsDebug):
            logging.debug("checkfilter perfectpbias no minDepth for %s:%s, " +
                          "numPerfect=%s", chromList, posList, numPerfect)

        # if we don't have enough perfect reads
        if numPerfect < aParamsDict["minPerfectReads"]:
            filters.add("perfmnad")
        # check if the percent of perfect reads is high enough
        if (mutCountReads + refCount > 0):
            perfectPct = round(numPerfect/float(mutCountReads + refCount), 2)
            if perfectPct < aParamsDict["minPerfectReadsPct"]:
                filters.add("perfmnaf")
            elif (anIsDebug):
                logging.debug("perfectPct=%s >= minPerfectPct=%s",
                              perfectPct, aParamsDict["minPerfectReadsPct"])

        # all of the counts have been done for the genomic forward strand
        # so we can use the vcf ref and altList here to set the MMP
        del mmpList[:]
        for base in ([aCurrData.ref] + aCurrData.altList):
            if ((base in mmpDict) and (mmpDict[base]["total"] > 0)):
                secondary = mmpDict[base]["secondary"]
                total = mmpDict[base]["total"]
                mmpList.append(str(round(secondary/float(total), 2)))
            else:
                mmpList.append("0.0")

        # return the superset of filters applied
        if (len(filters) > 0):
            return filters
        # if no filters have been returned thus far, this call passes
        else:
            return set(["PASS"])


if __name__ == '__main__':

    usage = "usage: python %prog vcfFile [Options]"
    cmdLineParser = OptionParser(usage=usage)
    cmdLineParser.add_option(
        "-c", "--allVCFCalls", action="store_false", default=True,
        dest="passedVCFCallsOnly",
        help="by default only the VCF calls that have passed all filters " +
             "thus far are processed, include this argument if all of the " +
             "VCF calls should be processed")
    cmdLineParser.add_option(
        "-s", "--scoreAllVCFCalls", action="store_false", default=True,
        dest="scorePassingVCFCallsOnly",
        help="by default the score will only be calculated for all passing " +
             "VCF calls, include this argument if the score should be " +
             "calculated for all VCF calls")
    cmdLineParser.add_option(
        "-o", "--outputFilename", default=sys.stdout,
        dest="outputFilename", metavar="OUTPUT_FILE",
        help="the name of the output file, STDOUT by default")
    cmdLineParser.add_option(
        "-l", "--log", default="WARNING",
        dest="logLevel", metavar="LOG",
        help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), " +
             "%default by default")
    cmdLineParser.add_option(
        "-g", "--logFilename",
        dest="logFilename", metavar="LOG_FILE",
        help="the name of the log file, STDERR by default")

    cmdLineParser.add_option(
        "", "--transcriptNameTag",
        dest="transcriptNameTag",
        help="the INFO key for the original transcript name")
    cmdLineParser.add_option(
        "", "--transcriptCoordinateTag",
        dest="transcriptCoordinateTag",
        help="the INFO key for the original transcript coordinate")
    cmdLineParser.add_option(
        "", "--transcriptStrandTag",
        dest="transcriptStrandTag",
        help="the INFO key for the original transcript strand")
    cmdLineParser.add_option(
        "", "--rnaIncludeSecondaryAlignments",
        action="store_true", default=False,
        dest="rnaIncludeSecondaryAlignments",
        help="if you align the RNA to transcript isoforms, then you may " +
             "want to include RNA secondary alignments in the " +
             "samtools mpileups")

    cmdLineParser.add_option(
        "", "--minPerfectReads",
        type="int", default=int(4),
        dest="minPerfectReads", metavar="MIN_PERFECT_READS",
        help="the minimum number of perfect reads needed for a call to " +
             "pass, %default by default")
    cmdLineParser.add_option(
        "", "--minPerfectReadsPct",
        type="float", default=float(0.10),
        dest="minPerfectReadsPct", metavar="MIN_PERFECT_READS_PCT",
        help="the minimum percentage of perfect reads supporting a variant " +
             "at a position, %default by default")
    cmdLineParser.add_option(
        "", "--maxMutsPerRead",
        type="int", default=int(4),
        dest="maxMutsPerRead", metavar="MAX_MUTS_PER_READ",
        help="the maximum number of mutations allowed in a perfect read, " +
             "%default by default")
    cmdLineParser.add_option(
        "", "--maxReadSoftClipPct",
        type="float", default=float(0.10),
        dest="maxReadSoftClipPct", metavar="MAX_READ_SOFT_CLIP_PCT",
        help="the maximum percentage of bases that can be soft-clipped " +
             "for a perfect read, %default by default")
    cmdLineParser.add_option(
        "", "--minBaseQual",
        type="int", default=int(10),
        dest="minBaseQual", metavar="MIN_BASE_QUAL",
        help="the minimum base quality for bases supporting the ALT, " +
             "%default by default")
    cmdLineParser.add_option(
        "", "--minMapQual",
        type="int", default=int(10),
        dest="minMapQual", metavar="MIN_MAP_QUAL",
        help="the minimum mapping quality for reads supporting the ALT, " +
             "%default by default")
    cmdLineParser.add_option(
        "", "--numNeighborBases",
        type="int", default=int(5),
        dest="numNeighborBases", metavar="NUM_NEIGHBOR_BASES",
        help="the number of neighboring bases to check for base quality, " +
             "%default by default")
    cmdLineParser.add_option(
        "", "--minNeighborBaseQual",
        type="int", default=int(10),
        dest="minNeighborBaseQual", metavar="MIN_NEIGHBOR_BASE_QUAL",
        help="the minimum neighboring bases quality, %default by default")
    cmdLineParser.add_option(
        "", "--maxStrandBias",
        type="float", default=float(0.90),
        dest="maxStrandBias", metavar="MAX_STRAND_BIAS",
        help="the maximum percentage of strand bias on reads that support " +
             "the ALT, %default by default")
    cmdLineParser.add_option(
        "", "--minStrandBiasDepth",
        type="int", default=int(4),
        dest="minStrandBiasDepth", metavar="MIN_STRAND_BIAS_DP",
        help="the minimum total depth needed for the strand bias filter " +
             "to be applied, %default by default")
    cmdLineParser.add_option(
        "", "--maxPositionBias",
        type="float", default=float(0.95),
        dest="maxPositionBias", metavar="MAX_POSITION_BIAS",
        help="the maximum percentage of positional bias for reads that " +
             "support the ALT, %default by default")
    cmdLineParser.add_option(
        "", "--minPositionBiasDepth",
        type="int", default=int(4),
        dest="minPositionBiasDepth", metavar="MIN_POSITION_BIAS_DP",
        help="the minimum total depth needed for the positional bias filter " +
             "to be applied, %default by default")
    cmdLineParser.add_option(
        "", "--maxMultiMapPct",
        type="float", default=float(0.90),
        dest="maxMultiMapPct", metavar="MAX_MULTI_MAP_PCT",
        help="the maximum percentage of secondary mapping reads that " +
             "support the ALT, %default by default")
    cmdLineParser.add_option(
        "", "--minMultiMapDepth",
        type="int", default=int(4),
        dest="minMultiMapDepth", metavar="MIN_MULTI_MAP_DP",
        help="the minimum total depth needed for the multi map filter " +
             "to be applied, %default by default")

    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(1, 46, 1)
    i_argLength = len(sys.argv)

    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        cmdLineParser.print_help()
        sys.exit(1)

    (cmdLineOpts, cmdLineArgs) = cmdLineParser.parse_args()

    i_readFilenameList = []
    i_writeFilenameList = []
    i_dirList = []

    # get the required params
    vcfFilename = cmdLineArgs[0]
    i_readFilenameList += [vcfFilename]

    # get the optional parameters with defaults
    i_passedVCFCallsOnlyFlag = cmdLineOpts.passedVCFCallsOnly
    i_scorePassingVCFCallsOnly = cmdLineOpts.scorePassingVCFCallsOnly
    i_logLevel = cmdLineOpts.logLevel
    i_rnaIncludeSecondaryAlignments = cmdLineOpts.rnaIncludeSecondaryAlignments
    i_maxReadSoftClipPct = cmdLineOpts.maxReadSoftClipPct
    i_minPerfectReads = cmdLineOpts.minPerfectReads
    i_minPerfectReadsPct = cmdLineOpts.minPerfectReadsPct
    i_minBaseQual = cmdLineOpts.minBaseQual
    i_minMapQual = cmdLineOpts.minMapQual
    i_numNeighborBases = cmdLineOpts.numNeighborBases
    i_minNeighborBaseQual = cmdLineOpts.minNeighborBaseQual
    i_maxStrandBias = cmdLineOpts.maxStrandBias
    i_minStrandBiasDepth = cmdLineOpts.minStrandBiasDepth
    i_maxPositionBias = cmdLineOpts.maxPositionBias
    i_minPositionBiasDepth = cmdLineOpts.minPositionBiasDepth
    i_maxMutsPerRead = cmdLineOpts.maxMutsPerRead
    i_maxMultiMapPct = cmdLineOpts.maxMultiMapPct
    i_minMultiMapDepth = cmdLineOpts.minMultiMapDepth

    # get the optional parameters with out defaults
    i_outputFilename = None
    i_logFilename = None
    i_txNameTag = None
    i_txCoordinateTag = None
    i_txStrandTag = None
    if (cmdLineOpts.outputFilename is not None):
        i_outputFilename = cmdLineOpts.outputFilename
        if (cmdLineOpts.outputFilename is not sys.stdout):
            i_writeFilenameList += [i_outputFilename]
    if (cmdLineOpts.logFilename is not None):
        i_logFilename = str(cmdLineOpts.logFilename)
        i_writeFilenameList += [i_logFilename]
    if (cmdLineOpts.transcriptNameTag is not None):
        i_txNameTag = cmdLineOpts.transcriptNameTag
    if (cmdLineOpts.transcriptCoordinateTag is not None):
        i_txCoordinateTag = cmdLineOpts.transcriptCoordinateTag
    if (cmdLineOpts.transcriptStrandTag is not None):
        i_txStrandTag = cmdLineOpts.transcriptStrandTag

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

    # do some debugging
    if (i_debug):
        logging.debug("vcfFilename=%s", vcfFilename)
        logging.debug("outputFilename=%s", i_outputFilename)
        logging.debug("logLevel=%s", i_logLevel)
        logging.debug("logFile=%s", i_logFilename)
        logging.debug("passingCallsOnlyFlag=%s", i_passedVCFCallsOnlyFlag)
        logging.debug("scorePassingCallsOnly=%s", i_scorePassingVCFCallsOnly)
        logging.debug("txNameTag=%s", i_txNameTag)
        logging.debug("txCoordinateTag=%s", i_txCoordinateTag)
        logging.debug("txStrandTag=%s", i_txStrandTag)
        logging.debug("rnaInclSecAlign=%s" % i_rnaIncludeSecondaryAlignments)
        logging.debug("maxReadSoftClipPct=%s" % i_maxReadSoftClipPct)
        logging.debug("minPerfectReads=%s" % i_minPerfectReads)
        logging.debug("minPerfectReadsPct=%s" % i_minPerfectReadsPct)
        logging.debug("minBaseQual=%s" % i_minBaseQual)
        logging.debug("minMapQual=%s" % i_minMapQual)
        logging.debug("numNeighborBases=%s" % i_numNeighborBases)
        logging.debug("minNeighborBaseQual=%s" % i_minNeighborBaseQual)
        logging.debug("maxStrandBias=%s" % i_maxStrandBias)
        logging.debug("minStrandBiasDepth=%s" % i_minStrandBiasDepth)
        logging.debug("maxPositionBias=%s" % i_maxPositionBias)
        logging.debug("minPositionBiasDepth=%s" % i_minPositionBiasDepth)
        logging.debug("maxMutsPerRead=%s" % i_maxMutsPerRead)
        logging.debug("maxMultiMapPct=%s" % i_maxMultiMapPct)
        logging.debug("minMultiMapDepth=%s" % i_minMultiMapDepth)

    params = {}
    params["rnaIncludeSecondaryAlignments"] = i_rnaIncludeSecondaryAlignments
    params["maxReadSoftClipPct"] = i_maxReadSoftClipPct
    params["minPerfectReads"] = i_minPerfectReads
    params["minPerfectReadsPct"] = i_minPerfectReadsPct
    params["maxMutsPerRead"] = i_maxMutsPerRead
    params["minBaseQual"] = i_minBaseQual
    params["minMapQual"] = i_minMapQual
    params["numNeighborBases"] = i_numNeighborBases
    params["minNeighborBaseQual"] = i_minNeighborBaseQual
    params["maxStrandBias"] = i_maxStrandBias
    params["minStrandBiasDepth"] = i_minStrandBiasDepth
    params["maxPositionBias"] = i_maxPositionBias
    params["minPositionBiasDepth"] = i_minPositionBiasDepth
    params["maxMultiMapPct"] = i_maxMultiMapPct
    params["minMultiMapDepth"] = i_minMultiMapDepth

    # check for any errors
    if (not radiaUtil.check_for_argv_errors(i_dirList,
                                            i_readFilenameList,
                                            i_writeFilenameList)):
        sys.exit(1)

    vcf = radiaUtil.get_read_fileHandler(vcfFilename)
    currVCF = myVCF.VCF()

    club = Club(vcfFilename, i_txNameTag, i_txCoordinateTag, i_debug)

    if i_outputFilename is not sys.stdout:
        i_outputFileHandler = radiaUtil.get_write_fileHandler(i_outputFilename)
    else:
        i_outputFileHandler = i_outputFilename

    hasAddedFilterHeader = False

    # loop through and categorize all calls
    for line in vcf:

        # add the filters to the header
        if ((not hasAddedFilterHeader) and line.startswith("##FILTER")):
            hasAddedFilterHeader = True
            i_outputFileHandler.write(
                "##FILTER=<ID=perfmnad,Description=\"The number of perfect " +
                "ALT reads is less than the minimum\">\n")
            i_outputFileHandler.write(
                "##FILTER=<ID=perfmnaf,Description=\"The percentage of " +
                "perfect ALT reads is less than the minimum\">\n")
            i_outputFileHandler.write(
                "##FILTER=<ID=perfsbias,Description=\"A strand bias " +
                "exists on the perfect ALT reads\">\n")
            i_outputFileHandler.write(
                "##FILTER=<ID=perfpbias,Description=\"A positional " +
                "bias exists on the perfect ALT reads\">\n")
            i_outputFileHandler.write(
                "##FILTER=<ID=mxmmp,Description=\"The percentage of " +
                "ALT multi-mapping reads is greater than the maximum\">\n")
            i_outputFileHandler.write(
                "##FILTER=<ID=pbias,Description=\"A positional bias " +
                "exists\">\n")

            i_outputFileHandler.write(line)
            continue
        # set the header columns
        elif line.startswith("#CHROM"):
            i_outputFileHandler.write(line)
            currVCF.set_headers(line[1:].strip().split("\t"))
            continue
        # output the header lines
        elif (line.startswith("#")):
            i_outputFileHandler.write(line)
            continue

        # now we are to a data line
        dataAsList = line.strip().split("\t")
        if len(dataAsList) == len(currVCF.headers):
            # get the current line
            currData = currVCF.make_data(dataAsList)

            # if we should only process passing calls and this call passes
            # or we should process all calls
            if ((i_passedVCFCallsOnlyFlag and "PASS" in currData.filterSet) or
                (not i_passedVCFCallsOnlyFlag)):

                # if this is a somatic mutation or an RNA editing event
                if (("SNP" in currData.infoDict["VT"]) and
                    ("2" in currData.infoDict["SS"] or
                     "4" in currData.infoDict["SS"])):

                    dnaFilters = set()
                    rnaFilters = set()

                    # Somatic calls can be made by the DNA or RNA
                    # RNA editing calls can be made by the normal or tumor RNA
                    for origin in currData.infoDict["ORIGIN"]:
                        if (origin == "DNA"):
                            # a somatic call made by the DNA
                            modType = currData.infoDict["MT"][0]
                            modChange = currData.infoDict["MC"][0]
                            dnaReadFilters = club.filter_by_read_support(
                                currData,
                                i_txNameTag,
                                i_txCoordinateTag,
                                i_txStrandTag,
                                modType,
                                modChange,
                                "DNA",
                                params,
                                i_debug)
                            dnaFilters = dnaFilters.union(dnaReadFilters)

                        # if we already passed using the DNA,
                        # then don't bother checking the RNA
                        elif ("PASS" not in dnaFilters):
                            # for RNA editing events, a call can have both
                            # normal and tumor editing, loop through them both
                            modTypes = currData.infoDict["MT"]
                            modChanges = currData.infoDict["MC"]
                            for (modType, modChange) in izip(modTypes,
                                                             modChanges):
                                rnaReadFilters = club.filter_by_read_support(
                                    currData,
                                    i_txNameTag,
                                    i_txCoordinateTag,
                                    i_txStrandTag,
                                    modType,
                                    modChange,
                                    "RNA",
                                    params,
                                    i_debug)
                                rnaFilters = rnaFilters.union(rnaReadFilters)

                    # if it passed by DNA or RNA, then it passed
                    if ("PASS" in dnaFilters or "PASS" in rnaFilters):
                        currData.filterSet = set(["PASS"])
                    # if this call was passing until now,
                    # then just add the new filters from here
                    elif ("PASS" in currData.filterSet):
                        currData.filterSet = dnaFilters.union(rnaFilters)
                    # if this call was not passing until now,
                    # then add the filters to the previous filters
                    else:
                        filters = dnaFilters.union(rnaFilters)
                        currData.filterSet = currData.filterSet.union(filters)

                    if (i_debug):
                        logging.debug("dnaFilter=%s, rnaFilter=%s, " +
                                      "currData.filter=%s", dnaFilters,
                                      rnaFilters, currData.filterSet)

            # set the score
            club.set_score(currData, i_scorePassingVCFCallsOnly, i_debug)

            # output the final line
            i_outputFileHandler.write(str(currData) + "\n")
        else:
            logging.warning("The number of VCF columns (%s) doesn't equal " +
                            "the number of VCF header columns (%s).",
                            len(dataAsList), len(currVCF.headers))
            logging.warning("Here are the header columns: %s", currVCF.headers)
            logging.warning("Here is the VCF line: %s", line.strip())
