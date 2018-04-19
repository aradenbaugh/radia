#!/usr/bin/env python
__requires__=['pysam>=0.8.1']
import pkg_resources
import pysam
import sys
import os
import time
from optparse import OptionParser
import myvcf
import gzip
import radiaUtil
import logging
from itertools import izip
import collections

'''
'    RNA and DNA Integrated Analysis (RADIA) identifies RNA and DNA variants in NGS data.
'    Copyright (C) 2010-2018  Amie Radenbaugh
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

cigardict = {}
cigardict[0] = "match"
cigardict[1] = "insertion"
cigardict[2] = "deletion"
cigardict[3] = "skipped"
cigardict[4] = "softclipped"
cigardict[5] = "hardclipped"
cigardict[6] = "padding"
cigardict[7] = "seqmatch"
cigardict[8] = "seqmismatch"


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


def getPassingGermlineAlts(aCurrData):
    alts = []
    # for passing germline calls, there should only be one, but double-check anyway
    modChanges = aCurrData.info["MC"].split(",")
    # keep track of germline alts
    for modChange in modChanges:
        (ref, alt) = modChange.split(">")
        alts.append(alt)
        
    return alts

    
def parsevcf(filename, aTranscriptNameTag, aTranscriptCoordinateTag, anIsDebug):
    vcf = get_read_fileHandler(filename)
    currVCF = myvcf.VCF()

    dnaNormalBam = None
    dnaTumorBam = None
    dnaTumorFasta = None
    dnaTumorChrPrefix = False
    rnaNormalBam = None
    rnaNormalFasta = None
    rnaNormalChrPrefix = False
    rnaTumorBam = None
    rnaTumorFasta = None
    rnaTumorChrPrefix = False
    
    germlineDict = {}
    transcriptGermlineDict = {}
    mutationsDict = {}
    filterDict = {}
    
    # loop through and categorize all calls
    # first process the header lines
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
                dnaNormalBam = generatorParamsDict["dnaNormalFilename"]
            if ("rnaNormalFilename" in generatorParamsDict):
                rnaNormalBam = generatorParamsDict["rnaNormalFilename"]
                rnaNormalFasta = generatorParamsDict["rnaNormalFastaFilename"]
                rnaNormalChrPrefix = generatorParamsDict["rnaNormalUseChrPrefix"]
            if ("dnaTumorFilename" in generatorParamsDict):
                dnaTumorBam = generatorParamsDict["dnaTumorFilename"]
                dnaTumorFasta = generatorParamsDict["dnaTumorFastaFilename"]
                dnaTumorChrPrefix = generatorParamsDict["dnaTumorUseChrPrefix"]
            if ("rnaTumorFilename" in generatorParamsDict):
                rnaTumorBam = generatorParamsDict["rnaTumorFilename"]
                rnaTumorFasta = generatorParamsDict["rnaTumorFastaFilename"]
                rnaTumorChrPrefix = generatorParamsDict["rnaTumorUseChrPrefix"]
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
        # if the number of VCF columns doesn't equal the number of VCF header columns
        if len(dataAsList) == len(currVCF.headers):
            # parse each VCF line
            curr_data = currVCF.make_data(dataAsList)
             
            # keep track of the passing germline and loh calls
            #if curr_data.info["VT"] == "SNP" and curr_data.filter == ["PASS"] and (curr_data.info["SS"] == "Germline" or curr_data.info["SS"] == "LOH"):
            if curr_data.info["VT"] == "SNP" and curr_data.filter == ["PASS"] and (curr_data.info["SS"] == "1" or curr_data.info["SS"] == "3"):
                # initialize the dict
                if (curr_data.chrom not in germlineDict):
                    germlineDict[curr_data.chrom] = {}
            
                # add the germline calls with the genomic coordinates to the germlineDict
                # RADIA calls will have the MC (Modification Change) in the info
                if ("MC" in curr_data.info):
                    germlineDict[curr_data.chrom][str(curr_data.pos-1)] = getPassingGermlineAlts(curr_data)
                else:
                    # if the call doesn't have the MC in the info, then just use the alt field
                    germlineDict[curr_data.chrom][str(curr_data.pos-1)] = curr_data.alt
                
                # if we have transcript names and coordinates for these calls
                if (aTranscriptNameTag != None and aTranscriptCoordinateTag != None):
                    # for each transcript name and coordinate 
                    for (transcriptName, transcriptCoordinate) in izip(curr_data.info[aTranscriptNameTag].split(","), curr_data.info[aTranscriptCoordinateTag].split(",")):
                        # initialize the dict
                        if transcriptName not in transcriptGermlineDict:
                            transcriptGermlineDict[transcriptName] = {}
                    
                        # RADIA calls will have the MC (Modification Change) in the info
                        if ("MC" in curr_data.info):
                            transcriptGermlineDict[transcriptName][str(int(transcriptCoordinate)-1)] = getPassingGermlineAlts(curr_data)
                        else:
                            # if the call doesn't have the MC in the info, then just use the alt field
                            transcriptGermlineDict[transcriptName][str(int(transcriptCoordinate)-1)] = curr_data.alt
                        
            # keep track of the passing somatic calls
            #elif curr_data.info["VT"] == "SNP" and curr_data.info["SS"] == "Somatic" and curr_data.filter == ["PASS"]:
            elif curr_data.info["VT"] == "SNP" and curr_data.info["SS"] == "2" and curr_data.filter == ["PASS"]:
                if (curr_data.chrom not in mutationsDict):
                    mutationsDict[curr_data.chrom] = {}
                mutationsDict[curr_data.chrom][str(curr_data.pos-1)] = curr_data
            # keep track of the filtered calls
            #elif curr_data.info["VT"] == "SNP" and curr_data.info["SS"] == "Somatic" and curr_data.filter != ["PASS"]:
            elif curr_data.info["VT"] == "SNP" and curr_data.info["SS"] == "2" and curr_data.filter != ["PASS"]:                
                if (curr_data.chrom not in filterDict):
                    filterDict[curr_data.chrom] = {}
                filterDict[curr_data.chrom][str(curr_data.pos-1)] = curr_data
        else:
            logging.error("The number of VCF columns (%s) doesn't equal the number of VCF header columns (%s).", len(dataAsList), len(currVCF.headers))
            logging.error("Here are the VCF header columns: %s", currVCF.headers)
            logging.error("Here is the VCF line: %s", line.strip())
            sys.exit(1)
    
    if (anIsDebug):
        logging.debug("germlineDict=%s", germlineDict)
        logging.debug("transcriptGermlineDict=%s", transcriptGermlineDict)
        logging.debug("mutationsDict=%s", mutationsDict)
        logging.debug("filterDict=%s", filterDict)
    
    return (dnaNormalBam, rnaNormalBam, rnaNormalFasta, rnaNormalChrPrefix, 
            dnaTumorBam, dnaTumorFasta, dnaTumorChrPrefix, 
            rnaTumorBam, rnaTumorFasta, rnaTumorChrPrefix, 
            germlineDict, transcriptGermlineDict, mutationsDict, filterDict)
    

def check_base_and_map_quals(aPileupread, aParamsDict, anIsDebug):
    
    # mapping quality scores are already converted to ints
    if aPileupread.alignment.mapq < aParamsDict["minMapQual"]: #MINMQ
        if (anIsDebug):
            logging.debug("checkread read MQ=%s < minMQ=%s", aPileupread.alignment.mapq, aParamsDict["minMapQual"])
        return False
    
    # the pysam documentation says: "base quality scores are unsigned chars but
    # they are *not* the ASCII encoded values, so no offset of 33 needs to be subtracted"
    # but this is only true for the 'query_qualities' and 'query_alignment_qualities'
    # we need to subtract the offset for the 'qual' field
    baseQualConverted = ord(aPileupread.alignment.qual[aPileupread.query_position])-33
    #queryQual = aPileupread.alignment.query_qualities[aPileupread.query_position]
    #baseQualConverted = ord(aPileupread.alignment.query_alignment_qualities[aPileupread.query_position])-33

    if (anIsDebug):
        #logging.debug("alignedRead alignment.qual=%s", aPileupread.alignment.qual)
        #logging.debug("alignedRead alignment.query_qualities=%s", aPileupread.alignment.query_qualities)
        #logging.debug("alignedRead alignment.query_alignment_qualities=%s", aPileupread.alignment.query_alignment_qualities)
        logging.debug("check_base_and_map_quals baseQuality qpos=%s, orgQual=%s, ordQual=%s, convertedQual=%s, minBQ=%s", aPileupread.query_position, aPileupread.alignment.qual[aPileupread.query_position], ord(aPileupread.alignment.qual[aPileupread.query_position]), baseQualConverted, aParamsDict["minBaseQual"])

    if baseQualConverted < aParamsDict["minBaseQual"]: #MINBQ
        if (anIsDebug):
            logging.debug("checkread base BQ=%s < minBQ=%s", baseQualConverted, aParamsDict["minBaseQual"])
        return False 
    
    # calculate indices to check (e.g. 5 before and 5 after)
    start = aPileupread.query_position - aParamsDict["numNeighborBases"]
    #if start < 2:
    #    start = 2
    if (start < 0):
        start = 0
    
    stop = aPileupread.query_position + aParamsDict["numNeighborBases"]
    #if stop > pileupread.alignment.rlen - 2:
    #    stop = pileupread.alignment.rlen - 2
    if stop > aPileupread.alignment.rlen:
        if (anIsDebug):
            logging.debug("neighbors past the length query_position=%s, rlen=%s, start=%s, stop=%s", aPileupread.query_position, aPileupread.alignment.rlen, start, stop)
        stop = aPileupread.alignment.rlen
    
    if (anIsDebug and (stop - start < 2 * aParamsDict["numNeighborBases"])):
        logging.debug("neighbor bases not 2x numNeighborBases: query_position=%s, rlen=%s, start=%s, stop=%s", aPileupread.query_position, aPileupread.alignment.rlen, start, stop)
        
    # if the region examined is somehow negative
    if start > stop:
        if (anIsDebug):
            logging.debug("checking negative neighbor bases query_position=%s, rlen=%s, start=%s, stop=%s", aPileupread.query_position, aPileupread.alignment.rlen, start, stop)
        return False

    # if anything in neighborhood has too low of base quality
    index = start
    for qualScore in aPileupread.alignment.qual[start:(stop + 1)]:
        # the pysam documentation says: "base quality scores are unsigned chars but they are *not* the ASCII encoded values, so no offset of 33 needs to be subtracted"
        # but the documentation seems to be wrong, we need to subtract the offset for samtools mpileup, samtools view, and pysam
        baseQualConverted = ord(qualScore)-33
        if (anIsDebug):
            logging.debug("checking neighbor base query_pos=%s, start=%s, stop=%s, index=%s, qual=%s, ordQual=%s, convertedQual=%s", aPileupread.query_position, start, stop, index, qualScore, ord(qualScore), baseQualConverted)
        index += 1
        if baseQualConverted < aParamsDict["minNeighborBaseQual"]: #MINNQS
            if (anIsDebug):
                logging.debug("checkread neighboring base BQ=%s < minNQS=%s", baseQualConverted, aParamsDict["minNeighborBaseQual"])
            return False

    if (anIsDebug):
        logging.debug("checkread found nothing")
     
    return True


def reverse_complement_nucleotide(aNucleotide):
    '''
    ' This function returns the reverse complement of the parameter aNucleotide
    '
    ' aNucleotide:    A nucleotide to reverse complement
    '''
    if aNucleotide in i_reverseCompDict:
        return i_reverseCompDict[aNucleotide]
    else:
        logging.error("Trying to reverse complement an unknown nucleotide: %s", aNucleotide)
        sys.exit(1)
        
    return None


def is_mutation(aReadDict, aRefList, anAltList, aBamOrigin, anIsDebug):
    
    orgReadBase = aReadDict["base"]
    readBase = aReadDict["base"]
    orgRefBase = aReadDict["refBase"]
    refBase = aReadDict["refBase"]
    
    # the RSEM fasta for transcripts has reads in 5' to 3' direction
    # the RNA bam files have reads aligned to the RSEM fasta in 5' to 3' direction
    # if the transcript is on the genomic "-" strand,
    # then we need to reverse complement the reads and the fasta
    # if this is an RNA call and the alignment is on the reverse strand, 
    # then reverse comp both the base and the fasta ref otherwise, 
    # just compare the forward strand base and fasta ref
    if (aBamOrigin == "RNA" and aReadDict["strand"] != None and aReadDict["strand"] == "-"):
        readBase = reverse_complement_nucleotide(readBase) 
        refBase = reverse_complement_nucleotide(refBase)
                
    if (readBase == refBase):
        if (anIsDebug):
            logging.debug("ismut() base matches reference, queryPos=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s, vcfRef=%s, vcfAlt=%s", aReadDict["sequenceIndex"], aReadDict["chrom"], aReadDict["pos"], readBase, refBase, orgReadBase, orgRefBase, aRefList, anAltList)
        return False
    elif readBase in anAltList:
        if (anIsDebug):
            logging.debug("ismut() base matches alt, queryPos=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s, vcfRef=%s, vcfAlt=%s", aReadDict["sequenceIndex"], aReadDict["chrom"], aReadDict["pos"], readBase, refBase, orgReadBase, orgRefBase, aRefList, anAltList)
        return True
    
    if (anIsDebug):
        logging.debug("ismut() found nothing, queryPos=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s, vcfRef=%s, vcfAlt=%s", aReadDict["sequenceIndex"], aReadDict["chrom"], aReadDict["pos"], readBase, refBase, orgReadBase, orgRefBase, aRefList, anAltList)
    return False


def mismatch_counts(aCigarOffset, aChrom, aCurrentPos, aCurrentIndex, aTranscriptStrand, anAlignedread, aGermlineDict, aTranscriptGermlineDict, aFastafile, aBamOrigin, anIsDebug):
    refs = 0
    germs = 0
    muts = 0
    
    for offset in range(aCigarOffset):
        pos = aCurrentPos + offset
        orgReadBase = anAlignedread.seq[aCurrentIndex + offset]
        orgRefBase = aFastafile.fetch(aChrom, pos, pos+1).upper()
        readBase = anAlignedread.seq[aCurrentIndex + offset]
        refBase = aFastafile.fetch(aChrom, pos, pos+1).upper()
        
        # the RSEM fasta for transcripts has reads in 5' to 3' direction
        # the RNA bam files have reads aligned to the RSEM fasta in 5' to 3' direction
        # if the transcript is on the genomic "-" strand,
        # then we need to reverse complement the reads and the fasta
        if (aBamOrigin == "RNA" and aTranscriptStrand != None and aTranscriptStrand == "-"):
            readBase = reverse_complement_nucleotide(readBase)
            refBase = reverse_complement_nucleotide(refBase)
        
        if readBase == refBase:
            refs += 1
            if (anIsDebug):
                logging.debug("nummut() base matches ref:  currentpos=%s, currentind=%s, offset=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s", aCurrentPos, aCurrentIndex, offset, aChrom, pos, orgReadBase, orgRefBase, readBase, refBase)
        elif aChrom in aGermlineDict and str(pos) in aGermlineDict[aChrom] and readBase in aGermlineDict[aChrom][str(pos)]:
            germs += 1
            if (anIsDebug):
                logging.debug("nummut() germline found:  currentpos=%s, currentind=%s, offset=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s", aCurrentPos, aCurrentIndex, offset, aChrom, pos, orgReadBase, orgRefBase, readBase, refBase)
        elif (aChrom in aTranscriptGermlineDict and str(pos) in aTranscriptGermlineDict[aChrom] and readBase in aTranscriptGermlineDict[aChrom][str(pos)]):
            germs += 1
            if (anIsDebug):
                logging.debug("nummut() transcript germline found:  currentpos=%s, currentind=%s, offset=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s, germlineAlts=%s", aCurrentPos, aCurrentIndex, offset, aChrom, pos, orgReadBase, orgRefBase, readBase, refBase, aTranscriptGermlineDict[aChrom][str(pos)])
        else:
            muts += 1
            if (anIsDebug):
                logging.debug("nummut() mutation found:  currentpos=%s, currentind=%s, offset=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s", aCurrentPos, aCurrentIndex, offset, aChrom, pos, orgReadBase, orgRefBase, readBase, refBase)

    return refs, germs, muts


def mutation_counts(anAlignedread, aGermlineDict, aTranscriptGermlineDict, aTranscriptStrand, aChrom, aFastafile, aBamOrigin, anIsDebug):
    
    # alignedread.pos or alignedread.reference_start is the 0-based leftmost coordinate
    currentpos = anAlignedread.pos
    # alignedread.qstart or alignedread.query_alignment_start is the start index of the aligned query portion of the sequence (0-based, inclusive)
    currentind = anAlignedread.qstart
    
    insCount = 0
    delCount = 0
    refCount = 0
    mutCount = 0
    germlineCount = 0
    softClippedCount = 0
    
    if (anIsDebug):
        logging.debug("mutation_counts() anAlignedread=%s", anAlignedread)
    
    for cig in anAlignedread.cigar:
        if (anIsDebug):
            logging.debug("alignedread.cigar=%s, cig=%s", anAlignedread.cigar, cig)
            logging.debug("currentpos=%s, currentind=%s", currentpos, currentind)
        
        # if it aligns
        if cigardict[cig[0]] == "match":
            
            # count the number of mismatches across the aligned portion of the read
            refs, germs, muts = mismatch_counts(cig[1], aChrom, currentpos, currentind, aTranscriptStrand, anAlignedread, aGermlineDict, aTranscriptGermlineDict, aFastafile, aBamOrigin, anIsDebug)
            refCount += refs
            germlineCount += germs
            mutCount += muts
            
            currentpos += cig[1]
            currentind += cig[1]
            
        elif cigardict[cig[0]] == "insertion":
            currentind += cig[1]
            insCount += cig[1]
        elif cigardict[cig[0]] == "deletion":
            currentpos += cig[1]
            delCount += cig[1]
        elif cigardict[cig[0]] == "skipped":
            continue
        elif cigardict[cig[0]] == "softclipped":
            softClippedCount += cig[1]
            continue
        elif cigardict[cig[0]] == "hardclipped":
            continue
        elif cigardict[cig[0]] == "padded":
            continue
        # the base in the read matches the reference
        elif cigardict[cig[0]] == "seqmatch":
            currentpos += cig[1]
            currentind += cig[1]
            refCount += cig[1]
        # the base in the read does not match the reference
        elif cigardict[cig[0]] == "seqmismatch":
            # count the number of mismatches across the aligned portion of the read
            refs, germs, muts = mismatch_counts(cig[1], aChrom, currentpos, currentind, aTranscriptStrand, anAlignedread, aGermlineDict, aTranscriptGermlineDict, aFastafile, aBamOrigin, anIsDebug)
            refCount += refs
            germlineCount += germs
            mutCount += muts
            
            currentpos += cig[1]
            currentind += cig[1]
        else:
            logging.error("Unexpected value in the read cigar string %s at position %s", anAlignedread.cigar, anAlignedread.tid)
        
        if (anIsDebug):    
            logging.debug("counts for this read: ins=%s, del=%s, refs=%s, muts=%s, germline=%s, softclipped=%s", insCount, delCount, refCount, mutCount, germlineCount, softClippedCount)
    
    return (insCount, delCount, mutCount, germlineCount, softClippedCount)


class Club():
    def __init__(self, vcffilename, aTranscriptNameTag, aTranscriptCoordinateTag, anIsDebug):
        (self.normalbamfilename, 
         self.rnanormalbamfilename, self.rnanormalfastafilename, self.rnanormalchrprefix, 
         self.tumorbamfilename, self.tumorfastafilename, self.tumorchrprefix, 
         self.rnatumorbamfilename, self.rnatumorfastafilename, self.rnatumorchrprefix, 
         self.germlineDict, self.transcriptGermlineDict, self.mutationsDict, self.filterDict) = parsevcf(
             vcffilename, aTranscriptNameTag, aTranscriptCoordinateTag, anIsDebug)
        
        if (self.tumorbamfilename != None):
            if (not os.path.isfile(self.tumorbamfilename)):
                print >> sys.stderr, "The BAM file specified in the header does not exist: ", self.tumorbamfilename
                sys.exit(1)
                
            self.tumorbamfile = pysam.Samfile(self.tumorbamfilename, 'rb')
    
        if (self.tumorfastafilename != None):
            if (not os.path.isfile(self.tumorfastafilename)):
                print >> sys.stderr, "The FASTA file specified in the header does not exist: ", self.tumorfastafilename
                sys.exit(1)

            self.fastafile = pysam.Fastafile(self.tumorfastafilename)
        
        if (self.rnanormalbamfilename != None):
            if (not os.path.isfile(self.rnanormalbamfilename)):
                print >> sys.stderr, "The BAM file specified in the header does not exist: ", self.rnanormalbamfilename
                sys.exit(1)
        
            if (self.rnanormalfastafilename == None or not os.path.isfile(self.rnanormalfastafilename)):
                print >> sys.stderr, "The FASTA file specified in the header does not exist: ", self.rnanormalfastafilename
                sys.exit(1)
            
            self.rnanormalbamfile = pysam.Samfile(self.rnanormalbamfilename, 'rb')
            self.rnanormalfastafile = pysam.Fastafile(self.rnanormalfastafilename)
        
        if (self.rnatumorbamfilename != None):
            if (not os.path.isfile(self.rnatumorbamfilename)):
                print >> sys.stderr, "The BAM file specified in the header does not exist: ", self.rnatumorbamfilename
                sys.exit(1)
        
            if (self.rnatumorfastafilename == None or not os.path.isfile(self.rnatumorfastafilename)):
                print >> sys.stderr, "The FASTA file specified in the header does not exist: ", self.rnatumorfastafilename
                sys.exit(1)
            
            self.rnatumorbamfile = pysam.Samfile(self.rnatumorbamfilename, 'rb')
            self.rnatumorfastafile = pysam.Fastafile(self.rnatumorfastafilename)
            
        
    def find_non_overlapping_reads(self, aReadsDict, aMinBaseQual, anIsDebug):
        '''
        ' This function loops through the reads with the same name.  Due to the inclusion of secondary
        ' alignments for RNA-Seq, there could be more than 2 reads that are paired.  If the read pairs overlap and
        ' the bases agree, then we keep the read with the highest base quality.  If the read pairs overlap and
        ' the bases disagree, then we only keep the read with the highest base quality if the other reads have a 
        ' low quality.  If the reads don't overlap, then keep all the reads. 
        '
        ' aReadsDict:    The dictionary of reads grouped by their read name
        ' aMinBaseQual:  A minimum base quality score for the base
        ' aCoordinate:   The genomic coordinate
        ' anIsDebug:     A flag for outputting debug messages to STDERR
        '''
        # loop through read pairs
        #     if read pairs overlap:
        #         if the bases agree:
        #             keep the one with the highest base quality
        #         else bases don't agree:
        #            keep the read with the high quality only if the other reads have a low quality
        #     else:
        #        keep all non-overlaps 
        nonOverlappingReadsList = []
        for (readName, readPairsList) in aReadsDict.iteritems():
            if (anIsDebug):
                logging.debug("readName=%s, len=%s, readList=%s", readName, len(readPairsList), readPairsList)
            
            # if there is only one read, then no overlapping reads exist
            if (len(readPairsList) == 1):
                if (anIsDebug):
                    logging.debug("only one read %s", readName)
                nonOverlappingReadsList.append(readPairsList[0])
            # check for overlapping reads
            else:
                maxBaseQual = -1
                maxBaseQualIndex = -1
                nextMaxBaseQual = -1
                basesSet = set()
                
                # for each read
                for index in range(0, len(readPairsList)):
                    readStart = readPairsList[index]["start"]
                    readSequenceIndex = readPairsList[index]["sequenceIndex"]
                    readMateStart = readPairsList[index]["mateStart"]
                    readLength = readPairsList[index]["qlen"]
                    # abs(p->b->core.pos + p->qpos - p->b->core.mpos) < p->b->core.l_qseq)
                    # samtools view pos = reference_start = 1-based leftmost coordinate
                    # pysam pos = reference_start = 0-based leftmost coordinate
                    # samtools view pnext (1-based) = pysam mpos (0-based) = next_reference_start = the position of the mate/next read.
                    # sequenceIndex pysam qpos = query_position = position of the read base at the pileup site, 0-based. None if is_del or is_refskip is set. 
                    # len(sequence) or pysam qlen = query_length length of the query sequence
                    
                    if (anIsDebug):
                        logging.debug("readStart(pos)=%s, readSequenceIndex(qpos)=%s, readMateStart(mpos)=%s, readLength(l_qseq)=%s", readStart, readSequenceIndex, readMateStart, readLength)
                        logging.debug("%s (readStart + readSequenceIndex): abs(readStart + readSequenceIndex - readMateStart) %s <? (readAlignedLength) %s", (readStart + readSequenceIndex), abs(readStart + readSequenceIndex - readMateStart), readLength)
                        
                    # if the reads overlap
                    # if abs(p->b->core.pos + p->qpos - p->b->core.mpos) < p->b->core.l_qseq)
                    if (abs(readStart + readSequenceIndex - readMateStart) < readLength):            
                        
                        if (anIsDebug):
                            logging.debug("reads overlap readStart=%s, sequenceIndex=%s, readMateStart=%s, abs=%s, alignedLen=%s", str(readStart), str(readSequenceIndex), str(readMateStart), str(abs(readStart + readSequenceIndex - readMateStart)), str(readLength))
                        base = readPairsList[index]["base"]
                        qual = ord(readPairsList[index]["baseQual"])-33
                        basesSet.add(base)
                        if (qual > maxBaseQual):
                            # if the current max is higher than the next max,
                            # then propagate the current max down
                            if (maxBaseQual > nextMaxBaseQual):
                                nextMaxBaseQual= maxBaseQual
                            maxBaseQual = qual
                            maxBaseQualIndex = index
                        elif (qual > nextMaxBaseQual):
                            nextMaxBaseQual = qual
                        if (anIsDebug):
                            logging.debug("base=%s, qual=%s, maxBaseQual=%s, maxBaseQualIndex=%s, basesSet=%s", base, str(qual), str(maxBaseQual), str(maxBaseQualIndex), basesSet)    
                    # if they don't overlap, then just add the reads
                    else:
                        if (anIsDebug):
                            logging.debug("reads don't overlap readStart=%s, sequenceIndex=%s, readMateStart=%s, abs=%s, alignedLen=%s", str(readStart), str(readSequenceIndex), str(readMateStart), str(abs(readStart + readSequenceIndex - readMateStart)), str(readLength))
                        nonOverlappingReadsList.append(readPairsList[index])
                        
                # if all the bases agree, then keep the one with the highest base quality and avoid double counting
                if (len(basesSet) == 1):
                    if (anIsDebug):
                        logging.debug("all bases in overlapping reads agree, keeping=%s, %s %s", str(maxBaseQualIndex), readName, readPairsList)
                    nonOverlappingReadsList.append(readPairsList[maxBaseQualIndex])
                # if the bases disagree, then it's likely a sequencing error. only keep the read with a high qual if the other reads have a low qual
                elif(maxBaseQual >= aMinBaseQual and nextMaxBaseQual < aMinBaseQual):
                    if (anIsDebug):
                        logging.debug("not all bases in overlapping reads agree, but keeping a high qual one, maxBaseQual=%s, nextMaxBaseQual=%s, minBaseQual=%s, keeping=%s, %s %s", str(maxBaseQual), str(nextMaxBaseQual), str(aMinBaseQual), str(maxBaseQualIndex), readName, readPairsList)
                    nonOverlappingReadsList.append(readPairsList[maxBaseQualIndex])
                elif anIsDebug:
                    logging.debug("bases in overlapping reads don't agree and no obvious read to select based on qual scores, maxBaseQual=%s, nextMaxBaseQual=%s, minBaseQual=%s, %s %s", str(maxBaseQual), str(nextMaxBaseQual), str(aMinBaseQual), readName, readPairsList)
        
        return nonOverlappingReadsList


    def checkfilter(self, aChromList, aPosList, aTranscriptStrandList, aRefList, anAltList, mutSS, mutType, aBamOrigin, aParamsDict, anIsDebug):
        
        if (anIsDebug):
            logging.debug("begin checkfilter for %s and %s, mutSS=%s, mutType=%s, origin=%s", aChromList, aPosList, mutSS, mutType, aBamOrigin)
            logging.debug("parmsDict: %s", aParamsDict)
        
        #expects a 0-based pos
        refCount = 0
        mutCountCheckedReads = 0
        mutCountReads = 0
        perfectStarts = 0
        perfectMiddles = 0
        perfectEnds = 0
        starts = 0
        middles = 0
        ends = 0
        readsWithMaxMuts = 0
        maxSoftClips = 0
        readsWithDels = 0
        readsWithIns = 0
        numPerfect = 0
        improperPairs = 0
        perfectForStrand = 0
        perfectRevStrand = 0
        readsDict = collections.defaultdict(list)
        bases = ""
        baseQuals = ""
        totalReads = 0
        keptReads = 0
        
        # loop through all of the transcripts
        #startTime = time.time()
        for (chrom, pos, strand) in izip(aChromList, list(map(int, aPosList)), aTranscriptStrandList):
        
            if (mutSS == "Somatic" or mutSS == "2"):
                if (aBamOrigin == "RNA"):
                    bamfile = self.rnatumorbamfile
                    fastafile = self.rnatumorfastafile
                    if (self.rnatumorchrprefix == "True"):
                        chrom = "chr" + chrom
                else:
                    bamfile = self.tumorbamfile
                    fastafile = self.fastafile
                    if (self.tumorchrprefix == "True"):
                        chrom = "chr" + chrom
            elif (mutSS == "4"):
                if (mutType == "NOR_EDIT"):
                    bamfile = self.rnanormalbamfile
                    fastafile = self.rnanormalfastafile
                    if (self.rnanormalchrprefix == "True"):
                        chrom = "chr" + chrom
                elif (mutType == "TUM_EDIT"):
                    bamfile = self.rnatumorbamfile
                    fastafile = self.rnatumorfastafile
                    if (self.rnatumorchrprefix == "True"):
                        chrom = "chr" + chrom
            
            # vcfs are 1-based and pysam requires a 0-based coordinate
            pos = pos - 1
            
            # get the reference base from the fasta
            refBase = fastafile.fetch(chrom, pos, pos+1).upper()
            
            if (anIsDebug):
                logging.debug("getting pileups for chrom=%s, pos=%s", chrom, pos)
            #logging.info("getting pileups for chrom=%s, pos=%s", chrom, pos)
            
            # future versions of pysam have some nice options
            #for pileupcolumn in bamfile.pileup(chrom, pos, pos+1, stepper="nofilter", 
            #                                   fastafile=fastafile, ignore_overlaps=True, 
            #                                   flag_filter=1540, ignore_orphans=True, 
            #                                   compute_baq=True, redo_baq=True):
            
            # get the pileups
            for pileupcolumn in bamfile.pileup(chrom, pos, pos+1, stepper="nofilter"):
                
                # move through pileup until at the correct position
                if pileupcolumn.pos < pos:
                    #if (anIsDebug):
                    #    logging.debug("continue pileupcolumn.pos=%s, pos=%s", pileupcolumn.pos, pos)
                    continue
                if pileupcolumn.pos > pos:
                    #if (anIsDebug):
                    #    logging.debug("break out pileupcolumn.pos=%s, pos=%s", pileupcolumn.pos, pos)
                    break
                
                # take care of over-lapping reads
                # loop through the reads and create a dictionary
                txTotalReads = 0
                txKeptReads = 0
                for pileupread in pileupcolumn.pileups:
                    txTotalReads += 1
                    alignedread = pileupread.alignment
                    
                    # skip over reads with problems
                    if (alignedread.is_qcfail or alignedread.is_unmapped or alignedread.is_duplicate or pileupread.is_del or pileupread.is_refskip):
                        #if (anIsDebug):
                        #    logging.debug("read is unmapped, duplicate, qcfail, del, or refskip at %s:%s", chrom, pos)
                        continue;
                    # no longer needed?? indel length for the position following the current pileup site
                    #if pileupread.indel != 0:
                    #    #if (anIsDebug):
                    #    #    logging.debug("base is an indel at %s:%s", chrom, pos)
                    #    return True
                    # skip over secondary mappings for RNA if the param is not set to include them
                    if (aBamOrigin == "RNA" and not aParamsDict["rnaIncludeSecondaryAlignments"] and pileupread.alignment.is_secondary):
                        #if (anIsDebug):
                        #    logging.debug("read is secondary alignment but flag to include secondary alignments for RNA is not set %s:%s", chrom, pos)
                        continue;
                    
                    # keep a dictionary of all reads, using the readName as the key
                    # due to the inclusion of secondary alignments for RNA-Seq, there could be more than 2 reads that are paired
                    oneReadDict = {}
                    txKeptReads += 1
                    oneReadDict["alignedRead"] = alignedread
                    oneReadDict["pileupRead"] = pileupread
                    #oneReadDict["qname"] = alignedread.query_name                          # qname
                    #oneReadDict["flag"] = alignedread.flag                                 # flag
                    #oneReadDict["rname"] = alignedread.reference_name                      # rname
                    oneReadDict["start"] = alignedread.reference_start                      # pos
                    #oneReadDict["mapQual"] = alignedread.mapping_quality                   # mapq
                    #oneReadDict["cigar"] = alignedread.cigar                               # cigar
                    #oneReadDict["mateName"] = alignedread.next_reference_name              # rnext
                    oneReadDict["mateStart"] = alignedread.next_reference_start             # pnext or mpos
                    #oneReadDict["insertSize"] = alignedread.template_length                # isize or tlen
                    #oneReadDict["sequence"] = alignedread.seq                              # seq
                    #oneReadDict["qualities"] = alignedread.qual                            # qual
                    oneReadDict["qlen"] = alignedread.query_length                          # qlen
                    oneReadDict["base"] = alignedread.seq[pileupread.query_position]
                    oneReadDict["baseQual"] = alignedread.qual[pileupread.query_position]
                    oneReadDict["sequenceIndex"] = pileupread.query_position
                    oneReadDict["chrom"] = chrom
                    oneReadDict["pos"] = pos
                    oneReadDict["strand"] = strand
                    oneReadDict["refBase"] = refBase
                    
                    #if (anIsDebug):
                    #    logging.debug("name=%s, oneReadDict=%s", alignedread.query_name, oneReadDict)
                    
                    # add it to the dictionary of reads
                    readsDict[alignedread.query_name].append(oneReadDict)
                    
                if (anIsDebug):
                    logging.debug("readsDictLen=%s", len(readsDict.keys()))
                    
            totalReads += txTotalReads
            keptReads += txKeptReads
            #logging.info("time_readSupport: %s:%s added %s out of %s reads to initial reads dict", chrom, pos, txKeptReads, txTotalReads)
        
        #stopTime = time.time()
        #logging.info("time_readSupport: %s:%s added %s out of %s reads to initial reads dict: Total time=%s hrs, %s mins, %s secs", aChromList, aPosList, keptReads, totalReads, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime))
        
        # get all of the non-overlapping reads
        #startTime = time.time()
        nonOverlappingReadsList = self.find_non_overlapping_reads(readsDict, aParamsDict["minBaseQual"], anIsDebug)
        #stopTime = time.time()
        #logging.info("time_readSupport: %s:%s nonoverlapping=%s out of %s: Total time=%s hrs, %s mins, %s secs", aChromList, aPosList, len(nonOverlappingReadsList), totalReads, ((stopTime-startTime)/(3600)), ((stopTime-startTime)/60), (stopTime-startTime))
        
        if (anIsDebug):
            logging.debug("nonOverlappingReadsListLen=%s", len(nonOverlappingReadsList))

        #startNonOverlappingTime = time.time()
        # loop over the non-overlapping reads
        for readDict in nonOverlappingReadsList:
            alignedread = readDict["alignedRead"]
            pileupread = readDict["pileupRead"]
            
            readBase = readDict["base"]
            readBaseQual = readDict["baseQual"]
            #readSequence = readDict["sequence"]
            #readSequenceIndex = readDict["sequenceIndex"]
            #readName = readDict["name"]
            #readMapQual = readDict["mapQual"]
            #readInsertSize = readDict["insertSize"]
            #readFlag = readDict["flag"]
            
            if (anIsDebug):
                logging.debug("found aligned read at: %s:%s = %s", readDict["chrom"], readDict["pos"], alignedread)
            
            # if this read supports the alternative allele
            if is_mutation(readDict, aRefList, anAltList, aBamOrigin, anIsDebug):
                
                mutCountReads += 1
                
                if (anIsDebug):
                    logging.debug("found read with mutation, number of reads with mutations=%s", mutCountReads)
                
                # if the base and map quals are good enough
                if check_base_and_map_quals(pileupread, aParamsDict, anIsDebug):
                    
                    mutCountCheckedReads +=1
                    
                    bases += readBase
                    baseQuals += readBaseQual
                    
                    if (anIsDebug):
                        logging.debug("mutation read passed base and map quals")
                    
                    # check to see how many are 'perfect' reads
                    perfect = True
                    alignedread = pileupread.alignment
                    
                    # count the number and types of mutations on this one read
                    insCount, delCount, mutCount, germCount, softClippedCount = mutation_counts(alignedread, self.germlineDict, self.transcriptGermlineDict, readDict["strand"], readDict["chrom"], fastafile, aBamOrigin, anIsDebug)
                    
                    if (anIsDebug):
                        logging.debug("this one read has ins=%s, del=%s, muts=%s, germline=%s, softClip=%s", insCount, delCount, mutCount, germCount, softClippedCount)
                    
                    if insCount:
                        readsWithIns += 1
                        perfect = False
                    if delCount:
                        readsWithDels += 1
                        perfect = False
                    # MapSplice doesn't set the proper pair flag for RNA-Seq reads, so only do this for DNA reads
                    if (aBamOrigin == "DNA" and not alignedread.is_proper_pair):
                        improperPairs += 1
                        perfect = False
                    # RNA editing events with mutSS=4 occur in clusters, so only apply this filter when mutSS!=4
                    if mutSS != "4" and mutCount >= aParamsDict["maxMutsPerRead"]:
                        readsWithMaxMuts += 1
                        perfect = False
                    
                    # count the positions in the reads
                    readLength = len(alignedread.query_sequence)
                    if (pileupread.query_position/float(readLength) <= 0.33):
                        starts += 1
                    elif (pileupread.query_position/float(readLength) <= 0.66):
                        middles += 1
                    else:
                        ends += 1
                    
                    if (anIsDebug):
                        logging.debug("query_pos=%s, readLength=%s, starts=%s, middles=%s, ends=%s", pileupread.query_position, readLength, starts, middles, ends)
                    
                    # softClippedCount is the number of soft clipped bases across this read
                    # if the percent of soft clipped bases is greater than the param, this read is not perfect
                    if softClippedCount:
                        if alignedread.query_length > 0:
                            if ((round(softClippedCount/float(alignedread.query_length)), 2) > aParamsDict["maxReadSoftClipPct"]):
                                maxSoftClips += 1
                                perfect = False
                                #if (anIsDebug):
                                #    logging.debug("read soft clipped too much, softClipped=%s, query_len=%s", softClippedCount, alignedread.query_length)
                        elif (alignedread.infer_query_length() > 0):
                            if ((round(softClippedCount/float(alignedread.infer_query_length())), 2) > aParamsDict["maxReadSoftClipPct"]):
                                maxSoftClips += 1
                                perfect = False
                                #if (anIsDebug):
                                #    logging.debug("read soft clipped too much, softClipped=%s, infer_query_len=%s", softClippedCount, alignedread.infer_query_length())
                        elif (len(alignedread.query_sequence) > 0):
                            if ((round(softClippedCount/float(len(alignedread.query_sequence))), 2) > aParamsDict["maxReadSoftClipPct"]):
                                maxSoftClips += 1
                                perfect = False
                                #if (anIsDebug):
                                #    logging.debug("read soft clipped too much, softClipped=%s, len seq=%s", softClippedCount, len(alignedread.query_sequence))
                        else:
                            logging.warning("Couldn't determine the read length and therefore apply the soft-clip filter for read=%s", alignedread)
                    
                    # if we found a perfect read
                    if perfect:
                        if (anIsDebug):
                            logging.debug("found perfect read at %s:%s with query_pos=%s, qlen=%s, readPositionPct=%s", readDict["chrom"], readDict["pos"], pileupread.query_position, alignedread.qlen, pileupread.query_position/float(alignedread.qlen))
                        
                        # determine the strand
                        if alignedread.is_reverse:
                            perfectRevStrand += 1
                        else:
                            perfectForStrand += 1
                        numPerfect += 1
                        
                        # count the perfect positions in the reads
                        readLength = len(alignedread.query_sequence)
                        if (pileupread.query_position/float(readLength) <= 0.33):
                            perfectStarts += 1
                        elif (pileupread.query_position/float(readLength) <= 0.66):
                            perfectMiddles += 1
                        else:
                            perfectEnds += 1
                        
                        if (anIsDebug):
                            logging.debug("query_pos=%s, readLength=%s, perfectStarts=%s, perfectMiddles=%s, perfectEnds=%s", pileupread.query_position, readLength, perfectStarts, perfectMiddles, perfectEnds)
                    
                    if (anIsDebug):
                        logging.debug("counts on mutRead for %s:%s, mutSS=%s, mutType=%s, properPair?=%s, perfect?=%s, numPerfect=%s", readDict["chrom"], readDict["pos"], mutSS, mutType, alignedread.is_proper_pair, perfect, numPerfect)
            else:
                refCount +=1
        
        #stopNonOverlappingTime = time.time()
        #logging.info("time_readSupport: %s:%s processing nonoverlapping reads: Total time=%s hrs, %s mins, %s secs", aChromList, aPosList, ((stopNonOverlappingTime-startNonOverlappingTime)/(3600)), ((stopNonOverlappingTime-startNonOverlappingTime)/60), (stopNonOverlappingTime-startNonOverlappingTime))
        
        if (anIsDebug):
            logging.debug("pysam bases=%s", bases)
            logging.debug("pysam quals=%s", baseQuals)
            logging.debug("final at pos %s:%s, mutSS=%s, mutType=%s, number of reads with ins=%s, dels=%s, maxMuts=%s, maxSoftClips=%s, numImproperPair=%s, numPerfect=%s", readDict["chrom"], readDict["pos"], mutSS, mutType, readsWithIns, readsWithDels, readsWithMaxMuts, maxSoftClips, improperPairs, numPerfect)
        
        # keep track of all filters
        filters = []
        
        # only apply strand bias to the perfect reads with mutations if we have enough reads
        if (numPerfect > 0) and (numPerfect >= aParamsDict["minStrandBiasDepth"]):
            sbias = float(perfectForStrand) / numPerfect
            #if sbias < 0.1 or sbias > 0.9:
            if (sbias > (aParamsDict["maxStrandBias"]) or sbias < (1.0 - aParamsDict["maxStrandBias"])):
                filters.append("perfectsbias")
                if (anIsDebug):
                    logging.debug("sbias flag for %s:%s, mutSS=%s, mutType=%s, forstrand=%s, revstrand=%s, sbias=%s", aChromList, aPosList, mutSS, mutType, perfectForStrand, perfectRevStrand, sbias)
        
        if (anIsDebug):
            logging.debug("checkfilter sbias for %s:%s, filters=%s, numPerfect=%s, mutCountCheckedReads=%s, mutCountReads=%s, refCount=%s", aChromList, aPosList, filters, numPerfect, mutCountCheckedReads, mutCountReads, refCount)

        # only apply positional bias to the reads with mutations if we have enough reads
        if ((mutCountCheckedReads > 0) and (mutCountCheckedReads >= aParamsDict["minPositionBiasDepth"])):
            if (anIsDebug):
                logging.debug("mutCountCheckedReads=%s, maxPBiasStarts=%s, maxPBiasEnds=%s", mutCountCheckedReads, round(starts/float(mutCountCheckedReads),2), round(ends/float(mutCountCheckedReads),2))
            if (round(starts/float(mutCountCheckedReads),2) >= aParamsDict["maxPositionBias"]):
                if (anIsDebug):
                    logging.debug("pbias from starts starts=%s, middles=%s, ends=%s, mutCountCheckedReads=%s", starts, middles, ends, mutCountCheckedReads)
                filters.append("pnewbias")
            elif (round(ends/float(mutCountCheckedReads),2) >= aParamsDict["maxPositionBias"]):
                if (anIsDebug):
                    logging.debug("pbias from ends starts=%s, middles=%s, ends=%s, mutCountCheckedReads=%s", starts, middles, ends, mutCountCheckedReads)
                filters.append("pnewbias")
        
        if (anIsDebug and mutCountCheckedReads > 0):
            logging.debug("checkfilter pbias for %s:%s, filters=%s, starts=%s (%s), middles=%s (%s), ends=%s (%s)", aChromList, aPosList, filters, starts, starts/float(mutCountCheckedReads), middles, middles/float(mutCountCheckedReads), ends, ends/float(mutCountCheckedReads))
        
        # only apply positional bias to the perfect reads with mutations if we have enough reads
        if ((numPerfect > 0) and (numPerfect >= aParamsDict["minPositionBiasDepth"])):   
            #if (anIsDebug):
            #    logging.debug("numPerfect=%s, maxPBiasStarts=%s, maxPBiasEnds=%s", numPerfect, round(perfectStarts/float(numPerfect),2), round(perfectEnds/float(numPerfect),2))
            if (round(perfectStarts/float(numPerfect),2) >= aParamsDict["maxPositionBias"]):
                #if (anIsDebug):
                #    logging.debug("perfectbias from starts perfectStarts=%s, perfectMiddles=%s, perfectEnds=%s", perfectStarts, perfectMiddles, perfectEnds)
                filters.append("perfectpbias")
            elif (round(perfectEnds/float(numPerfect),2) >= aParamsDict["maxPositionBias"]):
                #if (anIsDebug):
                #    logging.debug("perfectbias from ends perfectStarts=%s, perfectMiddles=%s, perfectEnds=%s", perfectStarts, perfectMiddles, perfectEnds)
                filters.append("perfectpbias")
                 
        if (anIsDebug and numPerfect > 0):
            logging.debug("checkfilter perfectpbias filters=%s, perfectStarts=%s (%s), perfectMiddles=%s (%s), perfectEnds=%s (%s)", filters, perfectStarts, perfectStarts/float(numPerfect), perfectMiddles, perfectMiddles/float(numPerfect), perfectEnds, perfectEnds/float(numPerfect))
        
        # if we don't have enough perfect reads
        if numPerfect < aParamsDict["minPerfectReads"]:
            filters.append("perfectcount")
        # check if the percent of perfect reads is high enough
        if (mutCountReads + refCount > 0):
            perfectPct = round(numPerfect/float(mutCountReads + refCount), 2)
            if perfectPct < aParamsDict["minPerfectReadsPct"]:
                filters.append("perfectpct")
            elif (anIsDebug):
                logging.debug("perfectPct=%s >= minPerfectPct=%s",  perfectPct, aParamsDict["minPerfectReadsPct"])
        
        # return the superset of filters applied
        if (len(filters) > 0):
            return filters
        # if no filters have been returned thus far, this call passes
        else:    
            return ["PASS"]
        

if __name__ == '__main__':
    
    usage = "usage: python %prog vcfFile [Options]"
    cmdLineParser = OptionParser(usage=usage)
    cmdLineParser.add_option("-c", "--allVCFCalls", action="store_false", default=True, dest="passedVCFCallsOnly", help="by default only the VCF calls that have passed all filters thus far are processed, include this argument if all of the VCF calls should be processed")
    cmdLineParser.add_option("-o", "--outputFilename", default=sys.stdout, dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, STDOUT by default")
    cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    
    cmdLineParser.add_option("", "--transcriptNameTag", dest="transcriptNameTag", help="the INFO key where the original transcript name can be found")
    cmdLineParser.add_option("", "--transcriptCoordinateTag", dest="transcriptCoordinateTag", help="the INFO key where the original transcript coordinate can be found")
    cmdLineParser.add_option("", "--transcriptStrandTag", dest="transcriptStrandTag", help="the INFO key where the original transcript strand can be found")
    cmdLineParser.add_option("", "--rnaIncludeSecondaryAlignments", action="store_true", default=False, dest="rnaIncludeSecondaryAlignments", help="if you align the RNA to transcript isoforms, then you may want to include RNA secondary alignments in the samtools mpileups")
    
    cmdLineParser.add_option("", "--minPerfectReads", type="int", default=int(4), dest="minPerfectReads", metavar="MIN_PERFECT_READS", help="the minimum number of perfect reads needed for a call to pass, %default by default")
    cmdLineParser.add_option("", "--minPerfectReadsPct", type="float", default=float(0.10), dest="minPerfectReadsPct", metavar="MIN_PERFECT_READS_PCT", help="the minimum percentage of perfect reads supporting a variant at a position, %default by default")
    cmdLineParser.add_option("", "--maxMutsPerRead", type="int", default=int(4), dest="maxMutsPerRead", metavar="MAX_MUTS_PER_READ", help="the maximum number of mutations allowed in a perfect read, %default by default")
    cmdLineParser.add_option("", "--maxReadSoftClipPct", type="float", default=float(0.10), dest="maxReadSoftClipPct", metavar="MAX_READ_SOFT_CLIP_PCT", help="the maximum percentage of bases that can be soft-clipped for a perfect read, %default by default")
    cmdLineParser.add_option("", "--minBaseQual", type="int", default=int(10), dest="minBaseQual", metavar="MIN_BASE_QUAL", help="the minimum base quality for bases supporting the ALT, %default by default")
    cmdLineParser.add_option("", "--minMapQual", type="int", default=int(10), dest="minMapQual", metavar="MIN_MAP_QUAL", help="the minimum mapping quality for reads supporting the ALT, %default by default")
    cmdLineParser.add_option("", "--numNeighborBases", type="int", default=int(5), dest="numNeighborBases", metavar="NUM_NEIGHBOR_BASES", help="the number of neighboring bases to check for base quality, %default by default")
    cmdLineParser.add_option("", "--minNeighborBaseQual", type="int", default=int(10), dest="minNeighborBaseQual", metavar="MIN_NEIGHBOR_BASE_QUAL", help="the minimum neighboring bases quality, %default by default")
    cmdLineParser.add_option("", "--maxStrandBias", type="float", default=float(0.90), dest="maxStrandBias", metavar="MAX_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    cmdLineParser.add_option("", "--minStrandBiasDepth", type="int", default=int(4), dest="minStrandBiasDepth", metavar="MIN_STRAND_BIAS_DP", help="the minimum total depth needed for the strand bias filter to be applied, %default by default")
    cmdLineParser.add_option("", "--maxPositionBias", type="float", default=float(0.95), dest="maxPositionBias", metavar="MAX_POSITION_BIAS", help="the maximum percentage of positional bias for reads that support the ALT, %default by default")
    cmdLineParser.add_option("", "--minPositionBiasDepth", type="int", default=int(4), dest="minPositionBiasDepth", metavar="MIN_POSITION_BIAS_DP", help="the minimum total depth needed for the positional bias filter to be applied, %default by default")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(1,38,1)
    i_argLength = len(sys.argv)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        cmdLineParser.print_help()
        sys.exit(1)
            
    (cmdLineOptions, cmdLineArgs) = cmdLineParser.parse_args()
    
    i_readFilenameList = []
    i_writeFilenameList = []
    i_dirList = []
    
    # get the required params
    vcfFilename = cmdLineArgs[0]
    i_readFilenameList += [vcfFilename]
    
    # get the optional parameters with defaults    
    i_passedVCFCallsOnlyFlag = cmdLineOptions.passedVCFCallsOnly
    i_logLevel = cmdLineOptions.logLevel
    i_rnaIncludeSecondaryAlignments = cmdLineOptions.rnaIncludeSecondaryAlignments
    i_maxReadSoftClipPct = cmdLineOptions.maxReadSoftClipPct
    i_minPerfectReads = cmdLineOptions.minPerfectReads
    i_minPerfectReadsPct = cmdLineOptions.minPerfectReadsPct
    i_minBaseQual = cmdLineOptions.minBaseQual
    i_minMapQual = cmdLineOptions.minMapQual
    i_numNeighborBases = cmdLineOptions.numNeighborBases
    i_minNeighborBaseQual = cmdLineOptions.minNeighborBaseQual
    i_maxStrandBias = cmdLineOptions.maxStrandBias
    i_minStrandBiasDepth = cmdLineOptions.minStrandBiasDepth
    i_maxPositionBias = cmdLineOptions.maxPositionBias
    i_minPositionBiasDepth = cmdLineOptions.minPositionBiasDepth
    i_maxMutsPerRead = cmdLineOptions.maxMutsPerRead
    
    # get the optional parameters with out defaults
    i_outputFilename = None
    i_logFilename = None
    i_transcriptNameTag = None
    i_transcriptCoordinateTag = None
    i_transcriptStrandTag = None
    if (cmdLineOptions.outputFilename != None):
        i_outputFilename = cmdLineOptions.outputFilename
        if (cmdLineOptions.outputFilename != sys.stdout):
            i_writeFilenameList += [i_outputFilename]
    if (cmdLineOptions.logFilename != None):
        i_logFilename = str(cmdLineOptions.logFilename)
        i_writeFilenameList += [i_logFilename]
    if (cmdLineOptions.transcriptNameTag != None):
        i_transcriptNameTag = cmdLineOptions.transcriptNameTag
    if (cmdLineOptions.transcriptCoordinateTag != None):
        i_transcriptCoordinateTag = cmdLineOptions.transcriptCoordinateTag
    if (cmdLineOptions.transcriptStrandTag != None):
        i_transcriptStrandTag = cmdLineOptions.transcriptStrandTag
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
    
    # do some debugging
    if (i_debug):
        logging.debug("vcfFilename=%s", vcfFilename)
        logging.debug("outputFilename=%s", i_outputFilename)
        logging.debug("logLevel=%s", i_logLevel)
        logging.debug("logFile=%s", i_logFilename)
        logging.debug("passingVCFCallsOnlyFlag=%s", i_passedVCFCallsOnlyFlag)
        logging.debug("transcriptNameTag=%s", i_transcriptNameTag)
        logging.debug("transcriptCoordinateTag=%s", i_transcriptCoordinateTag)
        logging.debug("transcriptStrandTag=%s", i_transcriptStrandTag)
        logging.debug("rnaIncludeSecondaryAlignments=%s" % i_rnaIncludeSecondaryAlignments)
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
        
    # check for any errors
    if (not radiaUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)
    
    vcf = get_read_fileHandler(vcfFilename)
    currVCF = myvcf.VCF()
    
    club = Club(vcfFilename, i_transcriptNameTag, i_transcriptCoordinateTag, i_debug)
    
    if i_outputFilename is not sys.stdout:
        i_outputFileHandler = get_write_fileHandler(i_outputFilename)
    else:
        i_outputFileHandler = i_outputFilename
    
    hasAddedFilterHeader = False
    headerLines = "##FILTER=<ID=perfectpct,Description=\"The percentage of perfect variant supporting reads was not above the minimum\">\n"
    headerLines += "##FILTER=<ID=perfectcount,Description=\"The number of perfect reads was not above the minimum\">\n"
    headerLines += "##FILTER=<ID=perfectsbias,Description=\"A strand bias exists on the perfect reads\">\n"
    headerLines += "##FILTER=<ID=perfectpbias,Description=\"A positional bias exists on the perfect reads\">\n"
    
    #loop through and categorize all calls
    for line in vcf:
        
        # add the filters to the header
        if ((not hasAddedFilterHeader) and line.startswith("##FILTER")):
            hasAddedFilterHeader = True       
            i_outputFileHandler.write(headerLines)
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
            curr_data = currVCF.make_data(dataAsList)
            
            # if this is a somatic mutation or an RNA editing event
            if curr_data.info["VT"] == "SNP" and curr_data.filter == ["PASS"] and (curr_data.info["SS"] == "Somatic" or curr_data.info["SS"] == "2" or curr_data.info["SS"] == "4"):
    
                # Pebbles doesn't have the ORIGIN flag
                if (curr_data.info["ORIGIN"] == None):
                    curr_data.filter = club.checkfilter([curr_data.chrom], [curr_data.pos], [None], curr_data.ref, curr_data.alt, curr_data.info["SS"], curr_data.info["MT"], "DNA", params, i_debug)
                # These are RADIA calls
                else:
                    dnaFilter = list()
                    rnaFilter = list()
                    originFlags = curr_data.info["ORIGIN"].split(",")
                    
                    # A call can be made by the DNA or RNA
                    for origin in originFlags:
                        if (origin == "DNA"):
                            if ((i_passedVCFCallsOnlyFlag and "PASS" in curr_data.filter) or (not i_passedVCFCallsOnlyFlag)):
                                # for RNA editing events, a call can have both normal and tumor editing, so loop through them both
                                modTypes = curr_data.info["MT"].split(",")
                                for modType in modTypes:
                                    dnaFilter += club.checkfilter([curr_data.chrom], [curr_data.pos], [None], curr_data.ref, curr_data.alt, curr_data.info["SS"], modType, "DNA", params, i_debug)
                        # if we already passed using the DNA, then don't bother checking the RNA
                        if ("PASS" not in dnaFilter):
                            if ((i_passedVCFCallsOnlyFlag and "PASS" in curr_data.filter) or (not i_passedVCFCallsOnlyFlag)):
                                # for RNA editing events, a call can have both normal and tumor editing, so loop through them both
                                modTypes = curr_data.info["MT"].split(",")
                                for modType in modTypes:
                                    try:
                                        # if the transcript name, coordinate, and strand should be used instead of the genomic chrom and coordinate
                                        if ((i_transcriptNameTag != None and i_transcriptCoordinateTag != None and i_transcriptStrandTag != None) and 
                                            (curr_data.info[i_transcriptNameTag] != None and curr_data.info[i_transcriptCoordinateTag] != None and curr_data.info[i_transcriptStrandTag] != None)):
                                            rnaFilter += club.checkfilter(curr_data.info[i_transcriptNameTag].split(","), curr_data.info[i_transcriptCoordinateTag].split(","), curr_data.info[i_transcriptStrandTag].split(","), curr_data.ref, curr_data.alt, curr_data.info["SS"], modType, "RNA", params, i_debug)
                                        else:
                                            rnaFilter += club.checkfilter([curr_data.chrom], [curr_data.pos], [None], curr_data.ref, curr_data.alt, curr_data.info["SS"], modType, "RNA", params, i_debug)
                                    except:
                                        logging.error("Problem with the following line: %s", line)
                                        raise
                    
                    # if it passed by DNA or RNA, then it passed          
                    if ("PASS" in dnaFilter or "PASS" in rnaFilter):
                        curr_data.filter = ["PASS"]
                    # if this call was passing until now, then just add the new filters from here
                    elif ("PASS" in curr_data.filter):
                        curr_data.filter = dnaFilter + rnaFilter
                    # if this call was not passing until now, then add the filters to the previous filters
                    else:
                        curr_data.filter += dnaFilter + rnaFilter
                        
                    logging.debug("dnaFilter=%s, rnaFilter=%s, curr_data.filter=%s", dnaFilter, rnaFilter, curr_data.filter)
            # output the final line
            i_outputFileHandler.write(str(curr_data) + "\n")
        else:
            logging.error("The number of VCF columns (%s) doesn't equal the number of VCF header columns (%s).", len(dataAsList), len(currVCF.headers))
            logging.error("Here are the VCF header columns: %s", currVCF.headers)
            logging.error("Here is the VCF line: %s", line.strip())
            sys.exit(1)
    

