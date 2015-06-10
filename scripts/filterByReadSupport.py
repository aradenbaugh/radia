#!/usr/bin/env python

__requires__=['pysam>=0.8.1']
import pkg_resources
import pysam
import sys
import os
from optparse import OptionParser
import myvcf
import gzip

MINNQS = 10
MINBQ = 10
MINMQ = 10

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


def mode(nums):
    numdict = {}
    for num in nums:
        if num not in numdict:
            numdict[num] = 0
        numdict[num] += 1
    maxcount = 0
    maxnum = -1
    
    for num, count in numdict.items():
        if count > maxcount:
            maxcount = count
            maxnum = num
    return maxnum, maxcount


def parsevcf(filename):
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
    
    germlinedict = {}
    mutations = {}
    filterDict = {}
    
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
                
        elif line.startswith("#CHROM"):
            currVCF.set_headers(line[1:].strip().split("\t"))
            break
    
    #loop through and categorize all calls
    for line in vcf:
        dataAsList = line.strip().split("\t")
        if len(dataAsList) == len(currVCF.headers):
            curr_data = currVCF.make_data(dataAsList)
            
            if (curr_data.chrom not in germlinedict):
                germlinedict[curr_data.chrom] = {}
                mutations[curr_data.chrom] = {}
                filterDict[curr_data.chrom] = {}
            
            if curr_data.info["VT"] == "SNP" and (curr_data.info["SS"] == "Germline" or curr_data.info["SS"] == "LOH") and curr_data.filter == ["PASS"]:
                germlinedict[curr_data.chrom][curr_data.pos] = curr_data.alt
            elif curr_data.info["VT"] == "SNP" and curr_data.info["SS"] == "Somatic" and curr_data.filter == ["PASS"]:
                mutations[curr_data.chrom][curr_data.pos] = curr_data
            elif curr_data.info["VT"] == "SNP" and curr_data.info["SS"] == "Somatic" and curr_data.filter != ["PASS"]:
                filterDict[curr_data.chrom][curr_data.pos] = curr_data
        else:
            #print line
            continue
        
    return dnaNormalBam, rnaNormalBam, rnaNormalFasta, rnaNormalChrPrefix, dnaTumorBam, dnaTumorFasta, dnaTumorChrPrefix, rnaTumorBam, rnaTumorFasta, rnaTumorChrPrefix, germlinedict, mutations, filterDict
    

def checkread(pileupread):
    if pileupread.alignment.mapq < MINMQ:
        return False
    if pileupread.is_del:
        return False
    
    if ord(pileupread.alignment.qual[pileupread.query_position]) - ord('!') < MINBQ:
        return False 

    #calculate indices to check, 5 before, and 5 after
    start = pileupread.query_position - 5
    if start < 2:
        start = 2
    stop = pileupread.query_position + 5
    
    if stop > pileupread.alignment.rlen - 2:
        stop = pileupread.alignment.rlen - 2
    
    #if the region examined is somehow negative
    if start > stop:
        return False

    #if anything in neighborhood has too low of base quality
    for q in pileupread.alignment.qual[start:stop]:
        if q < MINNQS:
            return False

    return True


def softclip(alignedread):
    leftcig, leftlen = alignedread.cigar[0]
    rightcig, rightlen = alignedread.cigar[-1]
    if cigardict[leftcig] == "softclipped":
        leftclip = alignedread.pos
    else:
        leftclip = -1
    if cigardict[rightcig] == "softclipped":
        rightclip = alignedread.pos + alignedread.qlen
    else:
        rightclip = -1
    return leftclip, rightclip


def ismut(pileupread, chrom, pos, fastafile, alt):
    #wait, what happens at insertions or deletions?
    if pileupread.indel != 0:
        return True
    if pileupread.alignment.seq[pileupread.query_position] == fastafile.fetch(chrom, pos, pos+1):
	return False
    if pileupread.alignment.seq[pileupread.query_position] in alt:
        return True
    return False


def nummut(alignedread, germline, samfile, fastafile):
    currentpos = alignedread.pos
    currentind = alignedread.qstart
    i = 0
    d = 0
    m = 0
    g = 0
    
    match = 0
    
    for cig in alignedread.cigar:
        if cigardict[cig[0]] == "match":
            #this just means it aligns
            chrom = samfile.getrname(alignedread.tid)
            for offset in range(cig[1]):
                pos = currentpos + offset
                base = alignedread.seq[currentind + offset]
                if base == fastafile.fetch(chrom, pos, pos+1):
                    match += 1
                elif chrom in germline and pos in germline[chrom] and base in germline[chrom][pos]:
                    g += 1
                else:
                    m += 1

            currentpos += cig[1]
            currentind += cig[1]
            continue
        elif cigardict[cig[0]] == "insertion":
            currentind += cig[1]
            i += cig[1]
        elif cigardict[cig[0]] == "deletion":
            currentpos += cig[1]
            d += cig[1]
        elif cigardict[cig[0]] == "skipped":
            #1/0
            continue
        elif cigardict[cig[0]] == "softclipped":
            continue
        elif cigardict[cig[0]] == "hardclipped":
            #1/0
            continue
        elif cigardict[cig[0]] == "padded":
            continue
        elif cigardict[cig[0]] == "seqmatch":
            #match! everything is okay, advance
            currentpos += cig[1]
            currentind += cig[1]
            match += cig[1]
            continue
        elif cigardict[cig[0]] == "seqmismatch":
            #mismatch, check to see if it's in the germline
            chrom = samfile.getrname(alignedread.tid)
            #loop through all mismatches
            for offset in range(cig[1]):
                pos = currentpos + offset
                base = alignedread.seq[currentind + offset]
                if chrom in germline and pos in germline[chrom] and base in germline[chrom][pos]:
                    g += 1
                else:
                    m += 1

            currentpos += cig[1]
            currentind += cig[1]
            continue
        else:
            1/0
    return (i, d, m, g)


class Club():
    def __init__(self, vcffilename):
        self.normalbamfilename, self.rnanormalbamfilename, self.rnanormalfastafilename, self.rnanormalchrprefix, self.tumorbamfilename, self.tumorfastafilename, self.tumorchrprefix, self.rnatumorbamfilename, self.rnatumorfastafilename, self.rnatumorchrprefix, self.germlinedict, self.mutations, self.filter = parsevcf(vcffilename)
        
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
            

    def checkfilter(self, chrom, pos, alt, mutSS, mutType, bamOrigin):
        
        if (mutSS == "Somatic" or mutSS == "2"):
            if (bamOrigin == "RNA"):
                bamfile = self.rnatumorbamfile
                fasta = self.rnatumorfastafile
                if (self.rnatumorchrprefix == "True"):
                    chrom = "chr" + chrom
            else:
                bamfile = self.tumorbamfile
                fasta = self.fastafile
                if (self.tumorchrprefix == "True"):
                    chrom = "chr" + chrom
        elif (mutSS == "5"):
            if (mutType == "NOR_EDIT"):
                bamfile = self.rnanormalbamfile
                fasta = self.rnanormalfastafile
                if (self.rnanormalchrprefix == "True"):
                    chrom = "chr" + chrom
            elif (mutType == "TUM_EDIT"):
                bamfile = self.rnatumorbamfile
                fasta = self.rnatumorfastafile
                if (self.rnatumorchrprefix == "True"):
                    chrom = "chr" + chrom
        
        #print >> sys.stderr, chrom, pos, mutSS, mutType, bamOrigin, bamfile
        
        #expects a 0 based pos
        refcount = 0
        mutcount = 0
        muts = [0] * 10
        dels = 0
        ins = 0
        left = 0
        right = 0
        improperpair = 0
        numperfect = 0
        forstrand = 0
        revstrand = 0
        
        leftclip = []
        rightclip = []
        mutpileupreads = []
        #get the pileup
        for pileupcolumn in bamfile.pileup(chrom, pos, pos+1):
        #move through pileup until at the correct position
            if pileupcolumn.pos != pos:
                continue
            for pileupread in pileupcolumn.pileups:
                alignedread = pileupread.alignment
                if ismut(pileupread, chrom, pos, fasta, alt) and checkread(pileupread):
                    #are orphans secondary?
                    l, r = softclip(alignedread)
                    if l != -1:
                        leftclip.append(l)
                    if r != -1:
                        rightclip.append(r)
                    mutcount +=1
                    mutpileupreads.append(pileupread)
                else:
                    refcount +=1
        
        maxleft, maxleftcount = mode(leftclip)
        maxright, maxrightcount = mode(rightclip)
        
        for pileupread in mutpileupreads:
            perfect = True
            
            alignedread = pileupread.alignment
            l, r = softclip(alignedread)
            if maxleftcount > 1 and l == maxleft:
                left += 1
                perfect = False
            if maxrightcount > 1 and r == maxright:
                right += 1
                perfect = False
            i, d, m, g = nummut(alignedread, self.germlinedict, bamfile, fasta)

	    if m > 9:
                m = 9
            muts[m] += 1
            if i:
                ins += 1
                perfect = False
            if d:
                dels += 1
                perfect = False
            # MapSplice doesn't set the proper pair flag for RNA-Seq reads, so only do this for DNA reads
            if (bamOrigin == "DNA" and not alignedread.is_proper_pair):
                improperpair += 1
                perfect = False
            # RNA editing events with mutSS=5 occur in clusters, so only apply this filter when mutSS!=5
            if mutSS != "5" and m >= 4:
                perfect = False
            if perfect:
                if alignedread.is_reverse:
                    revstrand+=1
                else:
                    forstrand+=1
                numperfect += 1
                
            #print >> sys.stderr, chrom, pos, mutSS, mutType, i, d, m, g, alignedread.is_proper_pair, maxleftcount > 1 and l == maxleft, maxrightcount > 1 and r == maxright, perfect, numperfect
        
        if (forstrand + revstrand) > 4:
            sbias = float(forstrand) / (forstrand + revstrand)
            if sbias < 0.1 or sbias > 0.9:
                sbiasflag = True
            else:
                sbiasflag = False
            #print >> sys.stderr, chrom, pos, mutSS, mutType, forstrand, revstrand, sbias, sbiasflag
        else:
            sbiasflag = False

        if numperfect < 5 and numperfect < mutcount:
            return ["perfect5"]
        elif sbiasflag and numperfect < mutcount:
            return ["perfectsbias"]
        else:
            return ["PASS"]


usage = "usage: python %prog vcfFile [Options]"
cmdLineParser = OptionParser(usage=usage)
cmdLineParser.add_option("-c", "--allVCFCalls", action="store_false", default=True, dest="passedVCFCallsOnly", help="by default only the VCF calls that have passed all filters thus far are processed, include this argument if all of the VCF calls should be processed")

# range(inclusiveFrom, exclusiveTo, by)
i_possibleArgLengths = range(1,5,1)
i_argLength = len(sys.argv)

# check if this is one of the possible correct commands
if (i_argLength not in i_possibleArgLengths):
    cmdLineParser.print_help()
    sys.exit(1)
        
(cmdLineOptions, cmdLineArgs) = cmdLineParser.parse_args()

filename = cmdLineArgs[0]

i_passedVCFCallsOnlyFlag = cmdLineOptions.passedVCFCallsOnly

vcf = get_read_fileHandler(filename)
currVCF = myvcf.VCF()

club = Club(filename)

for line in vcf:
    if line.startswith("#CHROM"):
        print "##FILTER=<ID=perfect5,Description=\"There were not enough perfect reads for this position.\">"
        print "##FILTER=<ID=perfectsbias,Description=\"A strand bias exists on the perfect reads. \">"
        print line[:-1]
        currVCF.set_headers(line[1:].strip().split("\t"))
        break

    print line[:-1]

#loop through and categorize all calls
for line in vcf:
    dataAsList = line.strip().split("\t")
    if len(dataAsList) == len(currVCF.headers):
        curr_data = currVCF.make_data(dataAsList)
        
        # if this is a somatic mutation or an RNA editing event
        if curr_data.info["VT"] == "SNP" and curr_data.filter == ["PASS"] and (curr_data.info["SS"] == "Somatic" or curr_data.info["SS"] == "2" or curr_data.info["SS"] == "5"):

            # Pebbles doesn't have the ORIGIN flag
            if (curr_data.info["ORIGIN"] == None):
                curr_data.filter = club.checkfilter(curr_data.chrom, curr_data.pos-1, curr_data.alt, curr_data.info["SS"], curr_data.info["MT"], "DNA")
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
                                dnaFilter += club.checkfilter(curr_data.chrom, curr_data.pos-1, curr_data.alt, curr_data.info["SS"], modType, "DNA")
                                    
                    # if we already passed using the DNA, then don't bother checking the RNA
                    elif ("PASS" not in dnaFilter):
                        if ((i_passedVCFCallsOnlyFlag and "PASS" in curr_data.filter) or (not i_passedVCFCallsOnlyFlag)):
                            # for RNA editing events, a call can have both normal and tumor editing, so loop through them both
                            modTypes = curr_data.info["MT"].split(",")
                            for modType in modTypes:
                                rnaFilter += club.checkfilter(curr_data.chrom, curr_data.pos-1, curr_data.alt, curr_data.info["SS"], modType, "RNA")
                                #print >> sys.stderr, curr_data.chrom, curr_data.pos-1, modType, rnaFilter

                if ("PASS" in dnaFilter or "PASS" in rnaFilter):
                    curr_data.filter = ["PASS"]
                else:
                    if (len(dnaFilter) > 0):
                        curr_data.filter = dnaFilter
                    elif (len(rnaFilter) > 0):
                        curr_data.filter = rnaFilter
        print curr_data
    else:
        #print line
        continue
    

