#!/usr/bin/env python
__requires__=['pysam>=0.8.1']
import pkg_resources
import pysam
import sys
import os
from optparse import OptionParser
import myvcf
import gzip
import radiaUtil
import logging
from itertools import izip
import collections
import math

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

cigarDict = {}
cigarDict[0] = "match"
cigarDict[1] = "insertion"
cigarDict[2] = "deletion"
cigarDict[3] = "refskipped"
cigarDict[4] = "softclipped"
cigarDict[5] = "hardclipped"
cigarDict[6] = "padding"
cigarDict[7] = "seqmatch"
cigarDict[8] = "seqmismatch"


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


def get_passing_germline_alts(aCurrData):
    alts = []
    # for passing germline calls, there should only be one, but double-check anyway
    modChanges = aCurrData.info["MC"].split(",")
    # keep track of germline alts
    for modChange in modChanges:
        (ref, alt) = modChange.split(">")
        alts.append(alt)
        
    return alts

    
def parse_vcf(aVCFFilename, aTranscriptNameTag, aTranscriptCoordinateTag, anIsDebug):
    vcf = get_read_fileHandler(aVCFFilename)
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
                    germlineDict[curr_data.chrom][str(curr_data.pos-1)] = get_passing_germline_alts(curr_data)
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
                            transcriptGermlineDict[transcriptName][str(int(transcriptCoordinate)-1)] = get_passing_germline_alts(curr_data)
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
            logging.warning("The number of VCF columns (%s) doesn't equal the number of VCF header columns (%s).", len(dataAsList), len(currVCF.headers))
            logging.warning("Here are the VCF header columns: %s", currVCF.headers)
            logging.warning("Here is the VCF line: %s", line.strip())
    
    #if (anIsDebug):
    #    logging.debug("germlineDict=%s", germlineDict)
    #    logging.debug("transcriptGermlineDict=%s", transcriptGermlineDict)
    #    logging.debug("mutationsDict=%s", mutationsDict)
    #    logging.debug("filterDict=%s", filterDict)
    
    return (dnaNormalBam, rnaNormalBam, rnaNormalFasta, rnaNormalChrPrefix, 
            dnaTumorBam, dnaTumorFasta, dnaTumorChrPrefix, 
            rnaTumorBam, rnaTumorFasta, rnaTumorChrPrefix, 
            germlineDict, transcriptGermlineDict, mutationsDict, filterDict)
    

def low_base_or_map_quals(aPileupread, aParamsDict, anIsDebug):
    
    # mapping quality scores are already converted to ints
    if aPileupread.alignment.mapq < aParamsDict["minMapQual"]: #MINMQ
        if (anIsDebug):
            logging.debug("low_base_or_map_quals MQ=%s < minMQ=%s", aPileupread.alignment.mapq, aParamsDict["minMapQual"])
        return True
    
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
        logging.debug("low_base_or_map_quals baseQuality qpos=%s, orgQual=%s, ordQual=%s, convertedQual=%s, minBQ=%s", aPileupread.query_position, aPileupread.alignment.qual[aPileupread.query_position], ord(aPileupread.alignment.qual[aPileupread.query_position]), baseQualConverted, aParamsDict["minBaseQual"])
    
    if baseQualConverted < aParamsDict["minBaseQual"]: #MINBQ
        if (anIsDebug):
            logging.debug("low_base_or_map_quals base BQ=%s < minBQ=%s", baseQualConverted, aParamsDict["minBaseQual"])
        return True
    
    if (anIsDebug):
        logging.debug("low_base_or_map_quals found nothing")
    
    return False


def low_neighbor_base_quals(aPileupread, aParamsDict, anIsDebug):
    
    # calculate indices to check (e.g. 5 before and 5 after)
    start = aPileupread.query_position - aParamsDict["numNeighborBases"]
    if (start < 0):
        start = 0
    
    stop = aPileupread.query_position + aParamsDict["numNeighborBases"]
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
        # the pysam documentation says: "base quality scores are unsigned chars but
        # they are *not* the ASCII encoded values, so no offset of 33 needs to be subtracted"
        # but this is only true for the 'query_qualities' and 'query_alignment_qualities'
        # we need to subtract the offset for the 'qual' field
        baseQualConverted = ord(qualScore)-33
        #if (anIsDebug):
        #    logging.debug("checking neighbor base query_pos=%s, start=%s, stop=%s, index=%s, qual=%s, ordQual=%s, convertedQual=%s", aPileupread.query_position, start, stop, index, qualScore, ord(qualScore), baseQualConverted)
        
        if baseQualConverted < aParamsDict["minNeighborBaseQual"]: #MINNQS
            if (anIsDebug):
                logging.debug("low_neighbor_base_quals neighboring base BQ=%s < minNBQ=%s", baseQualConverted, aParamsDict["minNeighborBaseQual"])
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
            logging.debug("is_mutation() base matches reference, queryPos=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s, vcfRef=%s, vcfAlt=%s", aReadDict["sequenceIndex"], aReadDict["chrom"], aReadDict["pos"], readBase, refBase, orgReadBase, orgRefBase, aRefList, anAltList)
        return False
    elif readBase in anAltList:
        if (anIsDebug):
            logging.debug("is_mutation() base matches alt, queryPos=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s, vcfRef=%s, vcfAlt=%s", aReadDict["sequenceIndex"], aReadDict["chrom"], aReadDict["pos"], readBase, refBase, orgReadBase, orgRefBase, aRefList, anAltList)
        return True
    
    if (anIsDebug):
        logging.debug("is_mutation() found nothing, queryPos=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s, vcfRef=%s, vcfAlt=%s", aReadDict["sequenceIndex"], aReadDict["chrom"], aReadDict["pos"], readBase, refBase, orgReadBase, orgRefBase, aRefList, anAltList)
    return False


def mismatch_counts(aCigarNum, aChrom, aRefIndex, aQueryIndex, aTranscriptStrand, anAlignedread, aGermlineDict, aTranscriptGermlineDict, aFastafile, aBamOrigin, aParamsDict, anIsDebug):
    refs = 0
    germs = 0
    muts = 0
    
    for offset in range(aCigarNum):
        # get the ref base
        refPos = aRefIndex + offset
        orgRefBase = aFastafile.fetch(aChrom, refPos, refPos+1).upper()
        refBase = aFastafile.fetch(aChrom, refPos, refPos+1).upper()
        
        # get the query base
        orgReadBase = anAlignedread.seq[aQueryIndex + offset]
        readBase = anAlignedread.seq[aQueryIndex + offset]
        
        # the RSEM fasta for transcripts has reads in 5' to 3' direction
        # the RNA bam files have reads aligned to the RSEM fasta in 5' to 3' direction
        # if the transcript is on the genomic "-" strand,
        # then we need to reverse complement the reads and the fasta
        if (aBamOrigin == "RNA" and aTranscriptStrand != None and aTranscriptStrand == "-"):
            readBase = reverse_complement_nucleotide(readBase)
            refBase = reverse_complement_nucleotide(refBase)
        
        if readBase == refBase:
            refs += 1
            #if (anIsDebug):
            #    logging.debug("mismatch_counts() base matches ref:  aRefIndex=%s, aQueryIndex=%s, offset=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s", aRefIndex, aQueryIndex, offset, aChrom, refPos, orgReadBase, orgRefBase, readBase, refBase)
        elif aChrom in aGermlineDict and str(refPos) in aGermlineDict[aChrom] and readBase in aGermlineDict[aChrom][str(refPos)]:
            germs += 1
            #if (anIsDebug):
            #    logging.debug("mismatch_counts() germline found:  aRefIndex=%s, aQueryIndex=%s, offset=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s", aRefIndex, aQueryIndex, offset, aChrom, refPos, orgReadBase, orgRefBase, readBase, refBase)
        elif (aChrom in aTranscriptGermlineDict and str(refPos) in aTranscriptGermlineDict[aChrom] and readBase in aTranscriptGermlineDict[aChrom][str(refPos)]):
            germs += 1
            #if (anIsDebug):
            #    logging.debug("mismatch_counts() transcript germline found:  aRefIndex=%s, aQueryIndex=%s, offset=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s, germlineAlts=%s", aRefIndex, aQueryIndex, offset, aChrom, refPos, orgReadBase, orgRefBase, readBase, refBase, aTranscriptGermlineDict[aChrom][str(refPos)])
        else:
            muts += 1
            #if (anIsDebug):
            #    logging.debug("mismatch_counts() mutation found:  aRefIndex=%s, aQueryIndex=%s, offset=%s, chrom=%s, pos=%s, orgBase=%s, orgRef=%s, base=%s, ref=%s", aRefIndex, aQueryIndex, offset, aChrom, refPos, orgReadBase, orgRefBase, readBase, refBase)

        # as soon as the muts count is >= the maxMutsPerRead, return
        if (muts >= aParamsDict["maxMutsPerRead"]):
            return refs, germs, muts
        
    return refs, germs, muts


class Club():
    def __init__(self, aVCFFilename, aTranscriptNameTag, aTranscriptCoordinateTag, anIsDebug):
        (self.dnaNormalBamFilename,
         self.rnaNormalBamFilename, self.rnaNormalFastaFilename, self.rnaNormalChrPrefix,
         self.dnaTumorBamFilename, self.dnaTumorFastaFilename, self.dnaTumorChrPrefix,
         self.rnaTumorBamFilename, self.rnaTumorFastaFilename, self.rnaTumorChrPrefix,
         self.germlineDict, self.transcriptGermlineDict, self.mutationsDict, self.filterDict) = parse_vcf(
             aVCFFilename, aTranscriptNameTag, aTranscriptCoordinateTag, anIsDebug)
        
        if (self.dnaTumorBamFilename != None):
            if (not os.path.isfile(self.dnaTumorBamFilename)):
                print >> sys.stderr, "The BAM file specified in the header does not exist: ", self.dnaTumorBamFilename
                sys.exit(1)
                
            if (self.dnaTumorFastaFilename == None or not os.path.isfile(self.dnaTumorFastaFilename)):
                print >> sys.stderr, "The FASTA file specified in the header does not exist: ", self.dnaTumorFastaFilename
                sys.exit(1)
            
            self.dnaTumorBamFile = pysam.Samfile(self.dnaTumorBamFilename, 'rb')
            self.dnaTumorFastaFile = pysam.Fastafile(self.dnaTumorFastaFilename)
        
        if (self.rnaNormalBamFilename != None):
            if (not os.path.isfile(self.rnaNormalBamFilename)):
                print >> sys.stderr, "The BAM file specified in the header does not exist: ", self.rnaNormalBamFilename
                sys.exit(1)
        
            if (self.rnaNormalFastaFilename == None or not os.path.isfile(self.rnaNormalFastaFilename)):
                print >> sys.stderr, "The FASTA file specified in the header does not exist: ", self.rnaNormalFastaFilename
                sys.exit(1)
            
            self.rnaNormalBamFile = pysam.Samfile(self.rnaNormalBamFilename, 'rb')
            self.rnaNormalFastaFile = pysam.Fastafile(self.rnaNormalFastaFilename)
        
        if (self.rnaTumorBamFilename != None):
            if (not os.path.isfile(self.rnaTumorBamFilename)):
                print >> sys.stderr, "The BAM file specified in the header does not exist: ", self.rnaTumorBamFilename
                sys.exit(1)
        
            if (self.rnaTumorFastaFilename == None or not os.path.isfile(self.rnaTumorFastaFilename)):
                print >> sys.stderr, "The FASTA file specified in the header does not exist: ", self.rnaTumorFastaFilename
                sys.exit(1)
            
            self.rnaTumorBamFile = pysam.Samfile(self.rnaTumorBamFilename, 'rb')
            self.rnaTumorFastaFile = pysam.Fastafile(self.rnaTumorFastaFilename)
        
        return


    def group_reads_by_name(self, aChromList, aPosList, aTranscriptStrandList, aBamOrigin, aParamsDict, aMutSS, aMutType, anIsDebug):
        
        readsDict = collections.defaultdict(list)
        # loop through all of the transcripts
        for (chrom, pos, strand) in izip(aChromList, list(map(int, aPosList)), aTranscriptStrandList):
            
            if (aMutSS == "Somatic" or aMutSS == "2"):
                if (aBamOrigin == "RNA"):
                    bamFile = self.rnaTumorBamFile
                    fastaFile = self.rnaTumorFastaFile
                    if (self.rnaTumorChrPrefix == "True"):
                        chrom = "chr" + chrom
                else:
                    bamFile = self.dnaTumorBamFile
                    fastaFile = self.dnaTumorFastaFile
                    if (self.dnaTumorChrPrefix == "True"):
                        chrom = "chr" + chrom
            elif (aMutSS == "4"):
                if (aMutType == "NOR_EDIT"):
                    bamFile = self.rnaNormalBamFile
                    fastaFile = self.rnaNormalFastaFile
                    if (self.rnaNormalChrPrefix == "True"):
                        chrom = "chr" + chrom
                elif (aMutType == "TUM_EDIT"):
                    bamFile = self.rnaTumorBamFile
                    fastaFile = self.rnaTumorFastaFile
                    if (self.rnaTumorChrPrefix == "True"):
                        chrom = "chr" + chrom
            
            # vcfs are 1-based and pysam requires a 0-based coordinate
            pos = pos - 1
            
            # get the reference base from the fasta
            refBase = fastaFile.fetch(chrom, pos, pos+1).upper()
            
            if (anIsDebug):
                logging.debug("getting pileups for chrom=%s, pos=%s", chrom, pos)
            
            # get the pileups
            for pileupColumn in bamFile.pileup(chrom, pos, pos+1, stepper="nofilter"):
                
                # move through pileup until at the correct position
                if pileupColumn.pos < pos:
                    #if (anIsDebug):
                    #    logging.debug("continue pileupColumn.pos=%s, pos=%s", pileupColumn.pos, pos)
                    continue
                if pileupColumn.pos > pos:
                    #if (anIsDebug):
                    #    logging.debug("break out pileupColumn.pos=%s, pos=%s", pileupColumn.pos, pos)
                    break
                
                totalReads = 0
                keptReads = 0
                
                # loop through the reads and create a dictionary
                for pileupRead in pileupColumn.pileups:
                    
                    totalReads += 1
                    alignedRead = pileupRead.alignment
                    
                    # skip over reads with problems
                    if (alignedRead.is_qcfail or alignedRead.is_unmapped or alignedRead.is_duplicate or pileupRead.is_del or pileupRead.is_refskip):
                        #if (anIsDebug):
                        #    logging.debug("read is unmapped, duplicate, qcfail, del, or refskip at %s:%s", chrom, pos)
                        continue;
                    # no longer needed?? indel length for the position following the current pileup site
                    if pileupRead.indel != 0:
                    #    #if (anIsDebug):
                    #    #    logging.debug("base is an indel at %s:%s", chrom, pos)
                        continue;
                    # skip over secondary mappings for RNA if the param is not set to include them
                    if (aBamOrigin == "RNA" and not aParamsDict["rnaIncludeSecondaryAlignments"] and alignedRead.is_secondary):
                        #if (anIsDebug):
                        #    logging.debug("read is secondary alignment but flag to include secondary alignments for RNA is not set %s:%s", chrom, pos)
                        continue;
                    
                    # keep a dictionary of all reads, using the readName as the key
                    # due to the inclusion of secondary alignments for RNA-Seq, there could be more than 2 reads that are paired
                    oneReadDict = {}
                    keptReads += 1
                    oneReadDict["alignedRead"] = alignedRead
                    oneReadDict["pileupRead"] = pileupRead
                    #oneReadDict["queryName"] = alignedRead.query_name                      # qname
                    #oneReadDict["flag"] = alignedRead.flag                                 # flag
                    #oneReadDict["rname"] = alignedRead.reference_name                      # rname
                    oneReadDict["start"] = alignedRead.reference_start                      # pos
                    #oneReadDict["mapQual"] = alignedRead.mapping_quality                   # mapq
                    #oneReadDict["cigar"] = alignedRead.cigar                               # cigar
                    #oneReadDict["mateName"] = alignedRead.next_reference_name              # rnext
                    oneReadDict["mateStart"] = alignedRead.next_reference_start             # pnext or mpos
                    #oneReadDict["insertSize"] = alignedRead.template_length                # isize or tlen
                    #oneReadDict["sequence"] = alignedRead.seq                              # seq
                    #oneReadDict["qualities"] = alignedRead.qual                            # qual
                    oneReadDict["qlen"] = alignedRead.query_length                          # qlen
                    oneReadDict["base"] = alignedRead.seq[pileupRead.query_position]
                    oneReadDict["baseQual"] = alignedRead.qual[pileupRead.query_position]
                    oneReadDict["sequenceIndex"] = pileupRead.query_position
                    oneReadDict["chrom"] = chrom
                    oneReadDict["pos"] = pos
                    oneReadDict["strand"] = strand
                    oneReadDict["refBase"] = refBase
                    
                    #if (anIsDebug):
                    #    logging.debug("name=%s, oneReadDict=%s", alignedRead.query_name, oneReadDict)
                    
                    # add it to the dictionary of reads
                    readsDict[alignedRead.query_name].append(oneReadDict)
                    
                if (anIsDebug):
                    logging.debug("group_reads_by_name(): %s:%s added %s out of %s reads to initial reads dict, readsDictLen=%s", chrom, pos, keptReads, totalReads, len(readsDict.keys()))
        
        return readsDict
    
        
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
            #if (anIsDebug):
            #    logging.debug("readName=%s, len=%s, readList=%s", readName, len(readPairsList), readPairsList)
            
            # if there is only one read, then no overlapping reads exist
            if (len(readPairsList) == 1):
                #if (anIsDebug):
                #    logging.debug("only one read %s", readName)
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
                    
                    #if (anIsDebug):
                    #    logging.debug("readStart(pos)=%s, readSequenceIndex(qpos)=%s, readMateStart(mpos)=%s, readLength(l_qseq)=%s", readStart, readSequenceIndex, readMateStart, readLength)
                    #    logging.debug("%s (readStart + readSequenceIndex): abs(readStart + readSequenceIndex - readMateStart) %s <? (readAlignedLength) %s", (readStart + readSequenceIndex), abs(readStart + readSequenceIndex - readMateStart), readLength)
                        
                    # if the reads overlap
                    # if abs(p->b->core.pos + p->qpos - p->b->core.mpos) < p->b->core.l_qseq)
                    if (abs(readStart + readSequenceIndex - readMateStart) < readLength):            
                        
                        #if (anIsDebug):
                        #    logging.debug("reads overlap readStart=%s, sequenceIndex=%s, readMateStart=%s, abs=%s, alignedLen=%s", str(readStart), str(readSequenceIndex), str(readMateStart), str(abs(readStart + readSequenceIndex - readMateStart)), str(readLength))
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
                        #if (anIsDebug):
                        #    logging.debug("base=%s, qual=%s, maxBaseQual=%s, maxBaseQualIndex=%s, basesSet=%s", base, str(qual), str(maxBaseQual), str(maxBaseQualIndex), basesSet)
                    # if they don't overlap, then just add the reads
                    else:
                        #if (anIsDebug):
                        #    logging.debug("reads don't overlap readStart=%s, sequenceIndex=%s, readMateStart=%s, abs=%s, alignedLen=%s", str(readStart), str(readSequenceIndex), str(readMateStart), str(abs(readStart + readSequenceIndex - readMateStart)), str(readLength))
                        nonOverlappingReadsList.append(readPairsList[index])
                        
                # if all the bases agree, then keep the one with the highest base quality and avoid double counting
                if (len(basesSet) == 1):
                    #if (anIsDebug):
                    #    logging.debug("all bases in overlapping reads agree, keeping=%s, %s %s", str(maxBaseQualIndex), readName, readPairsList)
                    nonOverlappingReadsList.append(readPairsList[maxBaseQualIndex])
                # if the bases disagree, then it's likely a sequencing error. only keep the read with a high qual if the other reads have a low qual
                elif(maxBaseQual >= aMinBaseQual and nextMaxBaseQual < aMinBaseQual):
                    #if (anIsDebug):
                    #    logging.debug("not all bases in overlapping reads agree, but keeping a high qual one, maxBaseQual=%s, nextMaxBaseQual=%s, minBaseQual=%s, keeping=%s, %s %s", str(maxBaseQual), str(nextMaxBaseQual), str(aMinBaseQual), str(maxBaseQualIndex), readName, readPairsList)
                    nonOverlappingReadsList.append(readPairsList[maxBaseQualIndex])
                #elif anIsDebug:
                #    logging.debug("bases in overlapping reads don't agree and no obvious read to select based on qual scores, maxBaseQual=%s, nextMaxBaseQual=%s, minBaseQual=%s, %s %s", str(maxBaseQual), str(nextMaxBaseQual), str(aMinBaseQual), readName, readPairsList)
        
        if (anIsDebug):
            logging.debug("find_non_overlapping_reads(): nonOverlappingReadsListLen=%s", len(nonOverlappingReadsList))
        
        return nonOverlappingReadsList


    def is_perfect(self, aPileupRead, aBamOrigin, aFastaFile, aMutSS, aReadDict, aParamsDict, anIsDebug):
        
        # check to see if this read is 'perfect'
        alignedRead = aPileupRead.alignment
        
        # MapSplice doesn't set the proper pair flag for RNA-Seq reads, so only do this for DNA reads
        if (aBamOrigin == "DNA" and not alignedRead.is_proper_pair):
            #improperPairs += 1
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
        
        # alignedread.pos or alignedread.reference_start is the 0-based leftmost coordinate
        refIndex = alignedRead.pos
        # alignedread.qstart or alignedread.query_alignment_start is the start index of the aligned query portion of the sequence (0-based, inclusive)
        # the index in the query sequence where the first ref is consumed
        queryIndex = 0
        
        refCount = 0
        mutCount = 0
        germlineCount = 0
        softClippedCount = 0
        
        # Op    Description                                        Consumes Query    Consumes Ref
        # M     alignment match (can be match or mismatch)                yes            yes
        # I     insertion to the ref                                      yes            no
        # D     deletion from the ref                                     no             yes
        # N     skipped region from the ref                               no             yes
        # S     soft-clipping (clipped sequences not present in query)    yes            no
        # H     hard-clipping (clipped sequences not present in query)    no             no
        # P     padding (silent deletion from padded ref)                 no             no
        # =     sequence match                                            yes            yes
        # X     sequence mismatch                                         yes            yes
        
        for cigarTuple in alignedRead.cigar:
            
            (cigarOp, cigarNum) = cigarTuple
            
            #if (anIsDebug):
            #    logging.debug("is_perfect() alignedread.cigar=%s, cigarstring=%s, cigarOpCode=%s, cigarOp=%s, cigarNum=%s", alignedRead.cigar, alignedRead.cigarstring, cigarOp, cigarDict[cigarOp], cigarNum)
            #    logging.debug("refIndex=%s, queryIndex=%s", refIndex, queryIndex)
            
            '''
            # this cigardict deciphers the cigarOp integer code
            cigarDict[0] = "match"
            cigarDict[1] = "insertion"
            cigarDict[2] = "deletion"
            cigarDict[3] = "refskipped"
            cigarDict[4] = "softclipped"
            cigarDict[5] = "hardclipped"
            cigarDict[6] = "padding"
            cigarDict[7] = "seqmatch"
            cigarDict[8] = "seqmismatch"
            '''
            
            # if it aligns
            if cigarDict[cigarOp] == "match":
                # RNA editing events with mutSS=4 occur in clusters, so don't bother counting mismatches
                if (aMutSS != "4"):
                    # count the number of mismatches across the aligned portion of the read
                    refs, germs, muts = mismatch_counts(cigarNum, aReadDict["chrom"], refIndex, queryIndex, aReadDict["strand"], alignedRead, self.germlineDict, self.transcriptGermlineDict, aFastaFile, aBamOrigin, aParamsDict, anIsDebug)
                    refCount += refs
                    germlineCount += germs
                    mutCount += muts
                
                refIndex += cigarNum
                queryIndex += cigarNum
                
            elif cigarDict[cigarOp] == "insertion":
                queryIndex += cigarNum
            elif cigarDict[cigarOp] == "deletion":
                refIndex += cigarNum
            elif cigarDict[cigarOp] == "refskipped":
                refIndex += cigarNum
            elif cigarDict[cigarOp] == "softclipped":
                queryIndex += cigarNum
                softClippedCount += cigarNum
            elif cigarDict[cigarOp] == "hardclipped":
                continue
            elif cigarDict[cigarOp] == "padded":
                continue
            # the base in the read matches the reference
            elif cigarDict[cigarOp] == "seqmatch":
                refIndex += cigarNum
                queryIndex += cigarNum
                refCount += cigarNum
            # the base in the read does not match the reference
            elif cigarDict[cigarOp] == "seqmismatch":
                # RNA editing events with mutSS=4 occur in clusters, so don't bother counting mismatches
                if (aMutSS != "4"):
                    # count the number of mismatches across the aligned portion of the read
                    refs, germs, muts = mismatch_counts(cigarNum, aReadDict["chrom"], refIndex, queryIndex, aReadDict["strand"], alignedRead, self.germlineDict, self.transcriptGermlineDict, aFastaFile, aBamOrigin, aParamsDict, anIsDebug)
                    refCount += refs
                    germlineCount += germs
                    mutCount += muts
                
                refIndex += cigarNum
                queryIndex += cigarNum
            else:
                logging.error("Unexpected value in the read cigar string %s at position %s", alignedRead.cigar, alignedRead.tid)
            
            if (anIsDebug):
                logging.debug("counts for this read: refs=%s, muts=%s, germline=%s, softclipped=%s", refCount, mutCount, germlineCount, softClippedCount)
            
            # RNA editing events with mutSS=4 occur in clusters, so only apply this filter when mutSS!=4
            if aMutSS != "4" and mutCount >= aParamsDict["maxMutsPerRead"]:
                return False, "maxMuts"
            
            # softClippedCount is the number of soft clipped bases across this read
            # if the percent of soft clipped bases is greater than the param, this read is not perfect
            if (len(alignedRead.query_sequence) > 0):
                if (round(softClippedCount/float(len(alignedRead.query_sequence)), 2) >= aParamsDict["maxReadSoftClipPct"]):
                    if (anIsDebug):
                        #logging.debug("softClipped=%s, qlen=%s, infer_qlen=%s, len(seq)=%s, ", softClippedCount, alignedRead.query_length, alignedRead.infer_query_length(), len(alignedRead.query_sequence))
                        logging.debug("read soft clipped too much, softClipped=%s, len(seq)=%s, softClipPct=%s, maxSFP=%s", softClippedCount, len(alignedRead.query_sequence), round(softClippedCount/float(len(alignedRead.query_sequence)), 2), aParamsDict["maxReadSoftClipPct"])
                    return False, "maxSoftClips"
        
        if (anIsDebug):
            logging.debug("this perfect read has refs=%s, muts=%s, germline=%s, softClip=%s", refCount, mutCount, germlineCount, softClippedCount)
        
        return True, "perfect"


    def get_score(self, anExpRefCount, anExpMutCount, anObsRefCount, anObsMutCount, aMaxSize, aFactLogList, anIsDebug):
        
        if (anIsDebug):
            logging.debug("anExpRefCount=%s, anExpMutCount=%s, anObsRefCount=%s, anObsMutCount=%s", anExpRefCount, anExpMutCount, anObsRefCount, anObsMutCount)
            logging.debug("1st R-tail maxSize=%s", aMaxSize)
        
        pValue = self.get_right_tail(anExpRefCount, anExpMutCount, anObsRefCount, anObsMutCount, aMaxSize, aFactLogList, anIsDebug)
        
        if (anIsDebug):
            logging.debug("after R-tail p-value=%s", pValue)
        
        if math.isnan(pValue):
            logging.warning("Warning: unable to calculate R-tail p-value failure: expRef=%s, expAlt=%s, obsRef=%s, obsAlt=%s", anExpRefCount, anExpMutCount, anObsRefCount, anObsMutCount)
        
        # If p-value is 1, do left-sided test
        if(pValue >= 0.999):
            
            if (anIsDebug):
                logging.debug("1st L-tail maxSize=%s", aMaxSize)
            
            pValue = self.get_left_tail(anExpRefCount, anExpMutCount, anObsRefCount, anObsMutCount, aMaxSize, aFactLogList, anIsDebug)
            
            if math.isnan(pValue):
                logging.warning("Warning: unable to calculate L-tail p-value failure: expRef=%s, expAlt=%s, obsRef=%s, obsAlt=%s", anExpRefCount, anExpMutCount, anObsRefCount, anObsMutCount)
            
        # if pValue is nan, then an error occurred
        # warning messages were written
        if math.isnan(pValue):
            pValue = 1
        
        # calculate the phred score
        if (pValue > 0):
            phred = 0 - (10 * math.log(pValue, 10))
        elif (pValue < 0):
            #pValue = 1 - pValue
            #phred = 0 - (10 * math.log(pValue, 10))
            phred = 0
        # the pvalue is 0
        else:
            phred = 255
        
        if (phred > 255):
            phred = 255
        
        return '{0:1.2e}'.format(pValue), int(round(phred))


    def get_pvalue(self, a, b, c, d, aFactLogList, anIsDebug):
        try:
            n = a + b + c + d
            
            if (anIsDebug):
                logging.debug("a=%s, b=%s, c=%s, d=%s, n=%s, a+b=%s, c+d=%s, a+c=%s, b+d=%s", a, b, c, d, n, a+b, c+d, a+c, b+d)
                logging.debug("f[a+b]=%s, f[c+d]=%s, f[a+c]=%s, f[b+d]=%s, adds=%s", aFactLogList[a+b], aFactLogList[c+d], aFactLogList[a+c], aFactLogList[b+d], (aFactLogList[a + b] + aFactLogList[c + d] + aFactLogList[a + c] + aFactLogList[b + d]))
                logging.debug("f[a]=%s, f[b]=%s, f[c]=%s, f[d]=%s, f[n]=%s, minuses=%s", aFactLogList[a], aFactLogList[b], aFactLogList[c], aFactLogList[d], aFactLogList[n], (aFactLogList[a] + aFactLogList[b] + aFactLogList[c] + aFactLogList[d] + aFactLogList[n]))
            
            p = (aFactLogList[a + b] + aFactLogList[c + d] + aFactLogList[a + c] + aFactLogList[b + d]) - (aFactLogList[a] + aFactLogList[b] + aFactLogList[c] + aFactLogList[d] + aFactLogList[n])
            
            if (anIsDebug):
                logging.debug("factorial=%s, pvalue=%s", p, math.exp(p))
            
            return math.exp(p)
        except:
            return float("nan")


    def get_right_tail(self, a, b, c, d, aMaxSize, aFactLogList, anIsDebug):
        n = a + b + c + d
        if (n > aMaxSize):
            logging.warning("Problem calculating R-tail: n=%s > maxSize=%s", n, aMaxSize)
            return float("nan")
        
        # get the first p-value
        p = self.get_pvalue(a, b, c, d, aFactLogList, anIsDebug)
        
        if (anIsDebug):
            logging.debug("doing R-tail: p=%s, a=%s, b=%s, c=%s, d=%s", p, a, b, c, d)
        
        if (c < b):
            minRight = c
        else:
            minRight = b
        
        for i in xrange(0, minRight):
            # if we've reached the pvalue giving us the max phred score of 255, just return
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
            logging.warning("Problem calculating L-tail: n=%s > maxSize=%s", n, aMaxSize)
            return float("nan")
        
        # get the first p-value
        p = self.get_pvalue(a, b, c, d, aFactLogList, anIsDebug)
        
        if (anIsDebug):
            logging.debug("doing L-tail: p=%s, a=%s, b=%s, c=%s, d=%s", p, a, b, c, d)
        
        if (a < d):
            minLeft = a
        else:
            minLeft = d
        
        for i in xrange(0, minLeft):
            # if we've reached the pvalue giving us the max phred score of 255, just return
            if (p < 0.00000000000000000000000001):
                logging.warning("L-tail breaking out of round %s of %s", i, minLeft)
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
    
    
    def init_factorial(self, aFactLogList, aPreviousMaxSize, aNewMaxSize, anIsDebug):
        
        if (aPreviousMaxSize == 0):
            aFactLogList[0] = 0.0
            aPreviousMaxSize = 1
        
        for index in xrange(aPreviousMaxSize, aNewMaxSize+1):
            aFactLogList[index] = aFactLogList[index - 1] + math.log(index)
        
        return aFactLogList


    def parse_info_field(self, anInfoField, anIsDebug):
    
        # parse the info field and create a dict
        infoDict = collections.defaultdict(list)
        infoFieldList = anInfoField.split(";")
        for info in infoFieldList:
            keyValueList = info.split("=")
            # some keys are just singular without a value (e.g. DB, SOMATIC, etc.)
            if (len(keyValueList) == 1):
                infoDict[keyValueList[0]] = ["True"]
            else:
                # the value can be a comma separated list
                infoDict[keyValueList[0]] = keyValueList[1].split(",")
        
        return infoDict


    def parse_sample_data(self, aFormatField, aDataField, anAlleleList, anIsDebug):
        
        # GT:DP:AD:AF:INS:DEL:START:STOP:MQ0:MMQ:MQA:BQ:SB
        # 0/1:7:2,5:0.29,0.71:0:0:0:0:1,3:1,1:1,0:82,71:1.0,1.0
        
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth at this position in the sample">
        ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Depth of reads supporting allele (in order specified by GT)">
        ##FORMAT=<ID=AF,Number=.,Type=Float,Description="Fraction of reads supporting allele (in order specified by GT)">
        ##FORMAT=<ID=INS,Number=1,Type=Integer,Description="Number of small insertions at this location">
        ##FORMAT=<ID=DEL,Number=1,Type=Integer,Description="Number of small deletions at this location">
        ##FORMAT=<ID=START,Number=1,Type=Integer,Description="Number of reads starting at this position">
        ##FORMAT=<ID=STOP,Number=1,Type=Integer,Description="Number of reads stopping at this position">
        ##FORMAT=<ID=MQ0,Number=.,Type=Integer,Description="Number of mapping quality zero reads harboring allele (in order specified by GT)">
        ##FORMAT=<ID=MMQ,Number=.,Type=Integer,Description="Maximum mapping quality of read harboring allele (in order specified by GT)">
        ##FORMAT=<ID=MQA,Number=.,Type=Integer,Description="Avg mapping quality for reads supporting allele (in order specified by GT)">
        ##FORMAT=<ID=BQ,Number=.,Type=Integer,Description="Avg base quality for reads supporting allele (in order specified by GT)">
        ##FORMAT=<ID=SB,Number=.,Type=Float,Description="Strand Bias for reads supporting allele (in order specified by GT)">
        
        alleleSpecificFormats = ["AD", "AF", "BQ", "SB", "MQ0", "MMQ", "MQA"]
        dataDict = dict()
        formatFieldList = aFormatField.split(":")
        dataFieldList = aDataField.split(":")
        
        for (formatItem, dataItem) in izip(formatFieldList, dataFieldList):
            if (formatItem == "GT"):
                sep = "/"
            else:
                sep = ","
            
            # if this is an allele specific item
            if formatItem in alleleSpecificFormats:
                
                if (dataItem != "."):
                    tmpList = dataItem.split(sep)
                else:
                    tmpList = [0] * len(anAlleleList)
                
                alleleDict = dict()
                alleleIndex = 0
                for allele in anAlleleList:
                    alleleDict[allele] = tmpList[alleleIndex]
                    alleleIndex += 1
                
                dataDict[formatItem] = alleleDict
                
            # split on the separator
            else:
                if (dataItem != "."):
                    dataDict[formatItem] = dataItem.split(sep)
                else:
                    dataDict[formatItem] = "0"
        
        return dataDict


    def filter_by_read_support(self, aChromList, aPosList, aTranscriptStrandList, aRefList, anAltList, aMutSS, aMutType, aBamOrigin, aParamsDict, anIsDebug):
        
        if (anIsDebug):
            logging.debug("begin filter for %s and %s, mutSS=%s, mutType=%s, origin=%s", aChromList, aPosList, aMutSS, aMutType, aBamOrigin)
            logging.debug("parmsDict: %s", aParamsDict)
        
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
        mutCountReadsWithMultiMap = 0
        
        bases = ""
        baseQuals = ""
        
        if (aMutSS == "Somatic" or aMutSS == "2"):
            if (aBamOrigin == "RNA"):
                fastaFile = self.rnaTumorFastaFile
            else:
                fastaFile = self.dnaTumorFastaFile
        elif (aMutSS == "4"):
            if (aMutType == "TUM_EDIT"):
                fastaFile = self.rnaTumorFastaFile
            elif (aMutType == "NOR_EDIT"):
                fastaFile = self.rnaNormalFastaFile
        
        # group all of the reads by name
        readsDict = self.group_reads_by_name(aChromList, aPosList, aTranscriptStrandList, aBamOrigin, aParamsDict, aMutSS, aMutType, anIsDebug)

        # get all of the non-overlapping reads
        nonOverlappingReadsList = self.find_non_overlapping_reads(readsDict, aParamsDict["minBaseQual"], anIsDebug)
        
        # loop over the non-overlapping reads
        for readDict in nonOverlappingReadsList:
            alignedRead = readDict["alignedRead"]
            pileupRead = readDict["pileupRead"]
            
            readBase = readDict["base"]
            readBaseQual = readDict["baseQual"]
            #readSequence = readDict["sequence"]
            #readSequenceIndex = readDict["sequenceIndex"]
            #readName = readDict["queryName"]
            #readMapQual = readDict["mapQual"]
            #readInsertSize = readDict["insertSize"]
            #readFlag = readDict["flag"]
            
            #if (anIsDebug):
            #    logging.debug("found aligned read at: %s:%s = %s", readDict["chrom"], readDict["pos"], alignedRead)
            
            # if this read supports the alternative allele
            if is_mutation(readDict, aRefList, anAltList, aBamOrigin, anIsDebug):
                
                #if (anIsDebug):
                #    logging.debug("found aligned read supporting variant at: %s:%s = %s", readDict["chrom"], readDict["pos"], alignedRead)
                
                mutCountReads += 1
                if alignedRead.is_secondary:
                    mutCountReadsWithMultiMap += 1
                
                if (anIsDebug):
                    logging.debug("found read with mutation, number of reads with mutations=%s", mutCountReads)
                
                # if the base and map quals are good enough
                if (not low_base_or_map_quals(pileupRead, aParamsDict, anIsDebug)):
                    
                    mutCountQualReads +=1
                    bases += readBase
                    baseQuals += readBaseQual
                    
                    if (anIsDebug):
                        logging.debug("mutation read passed base and map quals, qualReads=%s", mutCountQualReads)
                    
                    # count the positions in the reads
                    readLength = len(alignedRead.query_sequence)
                    if (pileupRead.query_position/float(readLength) <= 0.33):
                        starts += 1
                    elif (pileupRead.query_position/float(readLength) <= 0.66):
                        middles += 1
                    else:
                        ends += 1
                    
                    if (anIsDebug):
                        logging.debug("pbias query_pos=%s, readLength=%s, starts=%s, middles=%s, ends=%s", pileupRead.query_position, readLength, starts, middles, ends)
                    
                    # see if this read is perfect
                    isPerfectFlag, reasonNotPerfect = self.is_perfect(pileupRead, aBamOrigin, fastaFile, aMutSS, readDict, aParamsDict, anIsDebug)
                    
                    if (anIsDebug):
                        logging.debug("isPerfect?=%s, reason=%s", isPerfectFlag, reasonNotPerfect)
                    
                    # keep track of the counts
                    if reasonNotPerfect == "improperpair":
                        readsWithImproperPairs += 1
                    elif reasonNotPerfect == "insertion":
                        readsWithIns += 1
                    elif reasonNotPerfect == "deletion":
                        readsWithDels += 1
                    elif reasonNotPerfect == "neighborbasequals":
                        readsWithNeighborBaseQuals += 1
                    elif reasonNotPerfect == "maxMuts":
                        readsWithMaxMuts += 1
                    elif reasonNotPerfect == "maxSoftClips":
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
                        readLength = len(alignedRead.query_sequence)
                        if (pileupRead.query_position/float(readLength) <= 0.33):
                            perfectStarts += 1
                        elif (pileupRead.query_position/float(readLength) <= 0.66):
                            perfectMiddles += 1
                        else:
                            perfectEnds += 1
                        
                        if (anIsDebug):
                            logging.debug("perfectpbias: query_pos=%s, readLength=%s, perfectStarts=%s, perfectMiddles=%s, perfectEnds=%s", pileupRead.query_position, readLength, perfectStarts, perfectMiddles, perfectEnds)
                    
                    if (anIsDebug):
                        logging.debug("perfect read counts for %s:%s, mutSS=%s, mutType=%s, numPerfect=%s", readDict["chrom"], readDict["pos"], aMutSS, aMutType, numPerfect)
                else:
                    lowQualCount += 1
            else:
                refCount +=1
                bases += readBase
                baseQuals += readBaseQual
        
        if (anIsDebug):
            logging.debug("pysam bases=%s", bases)
            logging.debug("pysam quals=%s", baseQuals)
            logging.debug("final at pos %s:%s, mutSS=%s, mutType=%s, refCount=%s, mutCountReads=%s, mutCountQualReads=%s", aChromList, aPosList, aMutSS, aMutType, refCount, mutCountReads, mutCountQualReads)
            logging.debug("final at pos %s:%s, mutSS=%s, mutType=%s, lowQualCount=%s, improperPairs=%s, ins=%s, dels=%s, neighborBQ=%s, maxMuts=%s, maxSoftClips=%s, numPerfect=%s", aChromList, aPosList, aMutSS, aMutType, lowQualCount, readsWithImproperPairs, readsWithIns, readsWithDels, readsWithNeighborBaseQuals, readsWithMaxMuts, readsWithMaxSoftClips, numPerfect)
        
        # keep track of all filters
        filters = []
        
        mcMultiMappingPct = 0.0
        # testing multiple mapping to the mutCountReads
        if (mutCountReads > 0) and (mutCountReads >= aParamsDict["minMultiMapDepth"]):
            mcMultiMappingPct = round(mutCountReadsWithMultiMap/float(mutCountReads),2)
            if (mcMultiMappingPct >= aParamsDict["maxMultiMapPct"]):
                filters.append("mxmmp")
                if (anIsDebug):
                    logging.debug("checkfilter multiMap applied for %s:%s, mutSS=%s, mutType=%s, perfectReadsWithMultiMap=%s (%s)", aChromList, aPosList, aMutSS, aMutType, mutCountReadsWithMultiMap, mcMultiMappingPct)
        elif (anIsDebug):
            logging.debug("checkfilter multiMap no minDepth for %s:%s, mutSS=%s, mutType=%s, perfectReadsWithMultiMap=%s", aChromList, aPosList, aMutSS, aMutType, mutCountReadsWithMultiMap)
        
        # only apply strand bias to the perfect reads with mutations if we have enough reads
        if (numPerfect > 0) and (numPerfect >= aParamsDict["minStrandBiasDepth"]):
            sbias = round(perfectForStrand/float(numPerfect),2)
            #if sbias < 0.1 or sbias > 0.9:
            if (sbias > (aParamsDict["maxStrandBias"]) or sbias < (1.0 - aParamsDict["maxStrandBias"])):
                filters.append("perfsbias")
            if (anIsDebug):
                logging.debug("checkfilter sbias for %s:%s, numPerfect=%s, filters=%s, forstrand=%s, revstrand=%s, sbias=%s", aChromList, aPosList, numPerfect, filters, perfectForStrand, perfectRevStrand, sbias)
        elif (anIsDebug):
                logging.debug("checkfilter sbias no minDepth for %s:%s, numPerfect=%s", aChromList, aPosList, numPerfect)
        
        # only apply positional bias to the reads with mutations if we have enough reads
        if ((mutCountQualReads > 0) and (mutCountQualReads >= aParamsDict["minPositionBiasDepth"])):
            #if (anIsDebug):
            #    logging.debug("mutCountQualReads=%s, pbiasStarts=%s, pbiasEnds=%s", mutCountQualReads, round(starts/float(mutCountQualReads),2), round(ends/float(mutCountQualReads),2))
            if (round(starts/float(mutCountQualReads),2) >= aParamsDict["maxPositionBias"]):
                #if (anIsDebug):
                #    logging.debug("pbias from starts starts=%s, middles=%s, ends=%s, mutCountQualReads=%s", starts, middles, ends, mutCountQualReads)
                filters.append("pbias")
            elif (round(ends/float(mutCountQualReads),2) >= aParamsDict["maxPositionBias"]):
                #if (anIsDebug):
                #    logging.debug("pbias from ends starts=%s, middles=%s, ends=%s, mutCountQualReads=%s", starts, middles, ends, mutCountQualReads)
                filters.append("pbias")
            if (anIsDebug):
                logging.debug("checkfilter pbias for %s:%s, mutCountQualReads=%s, filters=%s, starts=%s (%s), middles=%s (%s), ends=%s (%s)", aChromList, aPosList, mutCountQualReads, filters, starts, starts/float(mutCountQualReads), middles, middles/float(mutCountQualReads), ends, ends/float(mutCountQualReads))
        elif (anIsDebug):
            logging.debug("checkfilter pbias no minDepth for %s:%s, mutCountQualReads=%s", aChromList, aPosList, mutCountQualReads)
        
        # only apply positional bias to the perfect reads with mutations if we have enough reads
        if ((numPerfect > 0) and (numPerfect >= aParamsDict["minPositionBiasDepth"])):   
            #if (anIsDebug):
            #    logging.debug("numPerfect=%s, perfectpbiasStarts=%s, perfectpbiasEnds=%s", numPerfect, round(perfectStarts/float(numPerfect),2), round(perfectEnds/float(numPerfect),2))
            if (round(perfectStarts/float(numPerfect),2) >= aParamsDict["maxPositionBias"]):
                #if (anIsDebug):
                #    logging.debug("perfectpbias from starts perfectStarts=%s, perfectMiddles=%s, perfectEnds=%s, numPerfect=%s", perfectStarts, perfectMiddles, perfectEnds, numPerfect)
                filters.append("perfpbias")
            elif (round(perfectEnds/float(numPerfect),2) >= aParamsDict["maxPositionBias"]):
                #if (anIsDebug):
                #    logging.debug("perfectpbias from ends perfectStarts=%s, perfectMiddles=%s, perfectEnds=%s, numPerfect=%s", perfectStarts, perfectMiddles, perfectEnds, numPerfect)
                filters.append("perfpbias")
            if (anIsDebug):
                logging.debug("checkfilter perfectpbias for %s:%s, numPerfect=%s, filters=%s, perfectStarts=%s (%s), perfectMiddles=%s (%s), perfectEnds=%s (%s)", aChromList, aPosList, numPerfect, filters, perfectStarts, perfectStarts/float(numPerfect), perfectMiddles, perfectMiddles/float(numPerfect), perfectEnds, perfectEnds/float(numPerfect))
        elif (anIsDebug):
            logging.debug("checkfilter perfectpbias no minDepth for %s:%s, numPerfect=%s", aChromList, aPosList, numPerfect)
        
        # if we don't have enough perfect reads
        if numPerfect < aParamsDict["minPerfectReads"]:
            filters.append("perfmnad")
        # check if the percent of perfect reads is high enough
        if (mutCountReads + refCount > 0):
            perfectPct = round(numPerfect/float(mutCountReads + refCount), 2)
            if perfectPct < aParamsDict["minPerfectReadsPct"]:
                filters.append("perfmnaf")
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
    cmdLineParser.add_option("-s", "--scoreAllVCFCalls", action="store_false", default=True, dest="scorePassingVCFCallsOnly", help="by default the score will only be calculated for all passing VCF calls, include this argument if the score should be calculated for all VCF calls")
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
    cmdLineParser.add_option("", "--maxMultiMapPct", type="float", default=float(0.90), dest="maxMultiMapPct", metavar="MAX_MULTI_MAP_PCT", help="the maximum percentage of secondary mapping reads that support the ALT, %default by default")
    cmdLineParser.add_option("", "--minMultiMapDepth", type="int", default=int(4), dest="minMultiMapDepth", metavar="MIN_MULTI_MAP_DP", help="the minimum total depth needed for the multi map filter to be applied, %default by default")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(1, 46, 1)
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
    i_scorePassingVCFCallsOnly = cmdLineOptions.scorePassingVCFCallsOnly
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
    i_maxMultiMapPct = cmdLineOptions.maxMultiMapPct
    i_minMultiMapDepth = cmdLineOptions.minMultiMapDepth
    
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
        logging.debug("scorePassingVCFCallsOnly=%s", i_scorePassingVCFCallsOnly)
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
    headerLines = "##FILTER=<ID=perfmnad,Description=\"The number of perfect ALT reads is less than the minimum\">\n"
    headerLines += "##FILTER=<ID=perfmnaf,Description=\"The percentage of perfect ALT reads is less than the minimum\">\n"
    headerLines += "##FILTER=<ID=perfsbias,Description=\"A strand bias exists on the perfect ALT reads\">\n"
    headerLines += "##FILTER=<ID=perfpbias,Description=\"A positional bias exists on the perfect ALT reads\">\n"
    headerLines += "##FILTER=<ID=mxmmp,Description=\"The percentage of ALT multi-mapping reads is greater than the maximum\">\n"
    headerLines += "##FILTER=<ID=pbias,Description=\"A positional bias exists\">\n"
    
    # initialize the factorial list
    factLogList = collections.defaultdict(float)
    maxSize = 1000
    factLogList = club.init_factorial(factLogList, 0, maxSize, i_debug)
    
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
            if (curr_data.info["VT"] == "SNP" and (curr_data.info["SS"] == "2" or curr_data.info["SS"] == "4")):
    
                # When merging multiple callers, non-RADIA calls won't have the ORIGIN flag
                if (curr_data.info["ORIGIN"] == None):
                    curr_data.filter = club.filter_by_read_support([curr_data.chrom], [curr_data.pos], [None], curr_data.ref, curr_data.alt, curr_data.info["SS"], curr_data.info["MT"], "DNA", params, i_debug)
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
                                    dnaFilter += club.filter_by_read_support([curr_data.chrom], [curr_data.pos], [None], curr_data.ref, curr_data.alt, curr_data.info["SS"], modType, "DNA", params, i_debug)
                        # if we already passed using the DNA, then don't bother checking the RNA
                        elif ("PASS" not in dnaFilter):
                            if ((i_passedVCFCallsOnlyFlag and "PASS" in curr_data.filter) or (not i_passedVCFCallsOnlyFlag)):
                                # for RNA editing events, a call can have both normal and tumor editing, so loop through them both
                                modTypes = curr_data.info["MT"].split(",")
                                for modType in modTypes:
                                    try:
                                        # if the transcript name, coordinate, and strand should be used instead of the genomic chrom and coordinate
                                        if ((i_transcriptNameTag != None and i_transcriptCoordinateTag != None and i_transcriptStrandTag != None) and 
                                            (curr_data.info[i_transcriptNameTag] != None and curr_data.info[i_transcriptCoordinateTag] != None and curr_data.info[i_transcriptStrandTag] != None)):
                                            rnaFilter += club.filter_by_read_support(curr_data.info[i_transcriptNameTag].split(","), curr_data.info[i_transcriptCoordinateTag].split(","), curr_data.info[i_transcriptStrandTag].split(","), curr_data.ref, curr_data.alt, curr_data.info["SS"], modType, "RNA", params, i_debug)
                                        else:
                                            rnaFilter += club.filter_by_read_support([curr_data.chrom], [curr_data.pos], [None], curr_data.ref, curr_data.alt, curr_data.info["SS"], modType, "RNA", params, i_debug)
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
                    
                    if (i_debug):
                        logging.debug("dnaFilter=%s, rnaFilter=%s, curr_data.filter=%s", dnaFilter, rnaFilter, curr_data.filter)
            
            # calculating the score is computationally expensive, so only do it for all passing calls by default
            pvalue = 0.98
            phred = 0
            if ((i_scorePassingVCFCallsOnly and "PASS" in curr_data.filter) or (not i_scorePassingVCFCallsOnly)):
                # a dict with either
                #    - the format item as key (e.g. DP) and single items as values or
                #    - the format item as key (e.g. AD) and a new dict with the alleles as keys and then values
                # e.g. rnaTumorDict["DP"] = 100
                # e.g. rnaTumorDict["AD"]["A"] = 50
                # e.g. rnaTumorDict["AD"]["T"] = 50
                dnaNormalDict = None
                rnaNormalDict = None
                dnaTumorDict = None
                rnaTumorDict = None
                try:
                    # if we have a 9th column, figure out which dataset it is
                    if (len(currVCF.headers) > 9):
                        if (currVCF.headers[9] == "DNA_NORMAL"):
                            dnaNormalDict = club.parse_sample_data(curr_data.genotype[0], curr_data.genotype[1], [curr_data.ref] + curr_data.alt, i_debug)
                        elif (currVCF.headers[9] == "RNA_NORMAL"):
                            rnaNormalDict = club.parse_sample_data(curr_data.genotype[0], curr_data.genotype[1], [curr_data.ref] + curr_data.alt, i_debug)
                        elif (currVCF.headers[9] == "DNA_TUMOR"):
                            dnaTumorDict = club.parse_sample_data(curr_data.genotype[0], curr_data.genotype[1], [curr_data.ref] + curr_data.alt, i_debug)
                        elif (currVCF.headers[9] == "RNA_TUMOR"):
                            rnaTumorDict = club.parse_sample_data(curr_data.genotype[0], curr_data.genotype[1], [curr_data.ref] + curr_data.alt, i_debug)
                    # if we have a 10th column, figure out which dataset it is
                    if (len(currVCF.headers) > 10):
                        if (currVCF.headers[10] == "RNA_NORMAL"):
                            rnaNormalDict = club.parse_sample_data(curr_data.genotype[0], curr_data.genotype[2], [curr_data.ref] + curr_data.alt, i_debug)
                        elif (currVCF.headers[10] == "DNA_TUMOR"):
                            dnaTumorDict = club.parse_sample_data(curr_data.genotype[0], curr_data.genotype[2], [curr_data.ref] + curr_data.alt, i_debug)
                        elif (currVCF.headers[10] == "RNA_TUMOR"):
                            rnaTumorDict = club.parse_sample_data(curr_data.genotype[0], curr_data.genotype[2], [curr_data.ref] + curr_data.alt, i_debug)
                    # if we have a 11th column, figure out which dataset it is
                    if (len(currVCF.headers) > 11):
                        if (currVCF.headers[11] == "DNA_TUMOR"):
                            dnaTumorDict = club.parse_sample_data(curr_data.genotype[0], curr_data.genotype[3], [curr_data.ref] + curr_data.alt, i_debug)
                        elif (currVCF.headers[11] == "RNA_TUMOR"):
                            rnaTumorDict = club.parse_sample_data(curr_data.genotype[0], curr_data.genotype[3], [curr_data.ref] + curr_data.alt, i_debug)
                    # if we have a 12th column, figure out which dataset it is
                    if (len(currVCF.headers) > 12):
                        if (currVCF.headers[12] == "RNA_TUMOR"):
                            rnaTumorDict = club.parse_sample_data(curr_data.genotype[0], curr_data.genotype[4], [curr_data.ref] + curr_data.alt, i_debug)
                    
                    # parse the info dict
                    infoDict = club.parse_info_field(dataAsList[7], i_debug)
                except:
                    logging.error("Problem with this line: %s", line)
                    raise
                
                if ("GERM" in infoDict["MT"]):
                    for (modType, modChange) in izip(infoDict["MT"], infoDict["MC"]):
                        if (i_debug):
                            logging.debug("modType=%s, modChange=%s", modType, modChange)
                        
                        if (modType == "GERM"):
                            ref, alt = modChange.split(">")
                            totalRefReads = 0
                            totalAltReads = 0
                            if (dnaNormalDict != None):
                                totalRefReads += int(dnaNormalDict["AD"][ref])
                                totalAltReads += int(dnaNormalDict["AD"][alt])
                            if (rnaNormalDict != None):
                                totalRefReads += int(rnaNormalDict["AD"][ref])
                                totalAltReads += int(rnaNormalDict["AD"][alt])
                            if (dnaTumorDict != None):
                                totalRefReads += int(dnaTumorDict["AD"][ref])
                                totalAltReads += int(dnaTumorDict["AD"][alt])
                            if (rnaTumorDict != None):
                                totalRefReads += int(rnaTumorDict["AD"][ref])
                                totalAltReads += int(rnaTumorDict["AD"][alt])
                            
                            totalCoverage = totalRefReads + totalAltReads
                            
                            newMaxSize = totalCoverage * 2
                            if (newMaxSize > maxSize):
                                factLogList = club.init_factorial(factLogList, maxSize, newMaxSize, i_debug)
                                maxSize = newMaxSize
                            pvalue, phred = club.get_score(totalCoverage, 0, totalRefReads, totalAltReads, maxSize, factLogList, i_debug)
                            if (i_debug):
                                logging.debug("pval=%s, phred=%s", pvalue, phred)
                            
                            curr_data.qual = str(phred)
                            curr_data.info["SSC"] = "0"
                            curr_data.info["PVAL"] = str(pvalue)
                
                elif ("SOM" in infoDict["MT"]):
                    for (modType, modChange) in izip(infoDict["MT"], infoDict["MC"]):
                        logging.debug("modType=%s, modChange=%s", modType, modChange)
                        if (modType == "SOM"):
                            ref, alt = modChange.split(">")
                            normalRefReads = 0
                            normalAltReads = 0
                            tumorRefReads = 0
                            tumorAltReads = 0
                            
                            if (dnaNormalDict != None):
                                normalRefReads += int(dnaNormalDict["AD"][ref])
                                normalAltReads += int(dnaNormalDict["AD"][alt])
                            if (rnaNormalDict != None):
                                normalRefReads += int(rnaNormalDict["AD"][ref])
                                normalAltReads += int(rnaNormalDict["AD"][alt])
                            if (dnaTumorDict != None):
                                tumorRefReads += int(dnaTumorDict["AD"][ref])
                                tumorAltReads += int(dnaTumorDict["AD"][alt])
                            if (rnaTumorDict != None):
                                tumorRefReads += int(rnaTumorDict["AD"][ref])
                                tumorAltReads += int(rnaTumorDict["AD"][alt])
                            
                            totalCoverage = normalRefReads + normalAltReads + tumorRefReads + tumorAltReads
                            
                            newMaxSize = totalCoverage
                            if (newMaxSize > maxSize):
                                factLogList = club.init_factorial(factLogList, maxSize, newMaxSize, i_debug)
                                maxSize = newMaxSize
                            pvalue, phred = club.get_score(normalRefReads, normalAltReads, tumorRefReads, tumorAltReads, maxSize, factLogList, i_debug)
                            if (i_debug):
                                logging.debug("pval=%s, phred=%s", pvalue, phred)
                            
                            curr_data.qual = str(phred)
                            curr_data.info["SSC"] = str(phred)
                            curr_data.info["PVAL"] = str(pvalue)
                            
                elif ("TUM_EDIT" in infoDict["MT"]):
                    for (modType, modChange) in izip(infoDict["MT"], infoDict["MC"]):
                        logging.debug("modType=%s, modChange=%s", modType, modChange)
                        if (modType == "TUM_EDIT"):
                            ref, alt = modChange.split(">")
                            totalRefReads = 0
                            totalAltReads = 0
                            editingRefReads = 0
                            editingAltReads = 0
                            if (dnaNormalDict != None):
                                totalRefReads += int(dnaNormalDict["AD"][ref])
                                totalAltReads += int(dnaNormalDict["AD"][alt])
                            if (rnaNormalDict != None):
                                totalRefReads += int(rnaNormalDict["AD"][ref])
                                totalAltReads += int(rnaNormalDict["AD"][alt])
                                editingRefReads += int(rnaNormalDict["AD"][ref])
                                editingAltReads += int(rnaNormalDict["AD"][alt])
                            if (dnaTumorDict != None):
                                totalRefReads += int(dnaTumorDict["AD"][ref])
                                totalAltReads += int(dnaTumorDict["AD"][alt])
                            if (rnaTumorDict != None):
                                totalRefReads += int(rnaTumorDict["AD"][ref])
                                totalAltReads += int(rnaTumorDict["AD"][alt])
                                editingRefReads += int(rnaTumorDict["AD"][ref])
                                editingAltReads += int(rnaTumorDict["AD"][alt])
                            
                            totalCoverage = totalRefReads + totalAltReads
                            
                            newMaxSize = totalCoverage + editingRefReads + editingAltReads
                            if (newMaxSize > maxSize):
                                factLogList = club.init_factorial(factLogList, maxSize, newMaxSize, i_debug)
                                maxSize = newMaxSize
                            pvalue, phred = club.get_score(totalCoverage, 0, editingRefReads, editingAltReads, maxSize, factLogList, i_debug)
                            if (i_debug):
                                logging.debug("pval=%s, phred=%s", pvalue, phred)
                            
                            curr_data.qual = str(phred)
                            curr_data.info["SSC"] = "0"
                            curr_data.info["PVAL"] = str(pvalue)
                else:
                    # these are RNA_TUM_VAR lines that have little or no DNA
                    curr_data.qual = str(phred)
                    curr_data.info["SSC"] = str(phred)
                    curr_data.info["PVAL"] = str(pvalue)
            else:
                # calculating the score is computationally expensive, so only do it for all passing calls by default
                # these are non-passing lines
                curr_data.qual = str(phred)
                curr_data.info["SSC"] = str(phred)
                curr_data.info["PVAL"] = str(pvalue)
            
            if (i_debug):
                logging.debug("pos=%s, qual=%s, pval=%s, phred=%s", curr_data.pos, curr_data.qual, pvalue, phred)
            
            # output the final line
            i_outputFileHandler.write(str(curr_data) + "\n")
        else:
            logging.warning("The number of VCF columns (%s) doesn't equal the number of VCF header columns (%s).", len(dataAsList), len(currVCF.headers))
            logging.warning("Here are the VCF header columns: %s", currVCF.headers)
            logging.warning("Here is the VCF line: %s", line.strip())
    
