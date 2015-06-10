#!/usr/bin/env python

from optparse import OptionParser   # used for parsing command line arguments
import radiaUtil                    # utility functions for rna editing
import sys                          # system module
import collections
import logging
from itertools import izip
import time
import subprocess
import re
import os
from math import floor
import gzip

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


# this regular expression is used to remove insertions and deletions from raw reads
# a read could look like:  "T$TT+3AGGGT+2AG+2AG.-2AGGG..-1A"
# insertions start with a "+", deletions with a "-"
# in theory, there could be multiple digits
i_numOfIndelsRegEx = re.compile("[+-]{1}(\\d+)")

# this regular expression will match any number of valid cDNA strings
i_cDNARegEx = re.compile("[ACGTNacgtn]+")


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


def fix_genotypes(aChrom, aRefList, anAltList, anAlleleDepthsList, aGTMinDepth, aGTMinPct):
    '''
    ' This method assigns the genotype.
    '
    ' aChrom: The chrom
    ' aRefList: The list of ref alleles
    ' anAltList: The list of alt alleles
    ' anAlleleDepthsList: The list of depths for all alleles
    ' aGTMinDepth: The minimum depth needed for the genotyping
    ' aGTMinPct: The minimum percent needed for the genotyping
    '''
    singleGenotypeChroms = ["chrY", "Y"]
    mChroms = ["chrM", "chrMT", "M", "MT"]
    refAltList = aRefList + anAltList
    
    # if it is a single chrom, then we can only assign one allele for the genotype
    # if one of the alts has a depth and percent above the mins, then use it, otherwise use the ref
    if (aChrom in singleGenotypeChroms):
        if (len(anAlleleDepthsList) > 0):
            
            # get the total depth
            totalDepth = sum(anAlleleDepthsList)
        
            # create a dict of just the alt alleles
            altCountsDict = {}
            for alt in anAltList:
                altIndex = refAltList.index(alt)
                altCountsDict[alt] = anAlleleDepthsList[altIndex]
                
            # find the max alt allele    
            (maxAltBase, maxAltDepth) = max(altCountsDict.iteritems(), key=lambda x:x[1])
            maxAltPct = round(maxAltDepth/float(totalDepth), 2)
            
            # if the max alt depth is large enough
            if (maxAltDepth >= aGTMinDepth and maxAltPct >= aGTMinPct):
                # find the index for the max depth on the original list
                maxAltIndex = refAltList.index(maxAltBase)
            else:
                # it wasn't large enough, so just use the ref
                maxAltIndex = 0
            
            # set the single genotype    
            genotypes = [maxAltIndex]
            
        else:
            # we don't have any bases, so just set it to the ref
            genotypes = [0]
    
    # if it is an M chrom, then we can assign as many alleles as we want for the genotype
    # for all bases with a depth and percent above the mins, set the genotype
    elif (aChrom in mChroms):
        if (len(anAlleleDepthsList) > 0):
            
            # get the total depth
            totalDepth = sum(anAlleleDepthsList)
            
            tmpGenotypes = []
            # for each base in the ref and alt
            for index in range(0, len(anAlleleDepthsList)):
                base = refAltList[index]
                depth = anAlleleDepthsList[index]
                # calculate the percent
                percent = round(depth/float(totalDepth), 2)
                # if the max alt depth and percent are large enough
                if (depth >= aGTMinDepth and percent >= aGTMinPct):
                    # add the index to the list
                    index = refAltList.index(base)
                    tmpGenotypes.append(index)
                    
            # if nothing passed the mins, then just take the ref
            if (len(tmpGenotypes) == 0):
                tmpGenotypes = [0]
                
            genotypes = sorted(tmpGenotypes)
        else:
            # we don't have any bases, so just set it to the ref
            genotypes = [0]
    
    # if it is a diploid chrom, then assign the 2 max counts above the min cutoffs
    else:
        # get the total depth
        totalDepth = sum(anAlleleDepthsList)

        # make a copy of the list to manipulate
        depthsTmpList = list(anAlleleDepthsList)
    
        # get the max depth
        max1Depth = max(depthsTmpList)
        
        # find the index for the max depth on the original list
        max1DepthIndex = anAlleleDepthsList.index(max1Depth)

        # remove the max from the tmp list
        depthsTmpList.remove(max1Depth)
    
        # if we still have some depths, find the 2nd max
        if (len(depthsTmpList) > 0):
        
            # get the max depth
            max2Depth = max(depthsTmpList)
            # have to do the rounding here, otherwise the ones that are just below 10% that get rounded up won't have the 0/1 genotype
            max2Pct = round(max2Depth/float(totalDepth), 2)
            
            # if the max depth is large enough
            if (max2Depth >= aGTMinDepth and max2Pct >= aGTMinPct):
                # if the two maxes are the same depth, then return the second index
                if (max1Depth == max2Depth):
                    # get all indices for this depth
                    indices = [i for i, x in enumerate(anAlleleDepthsList) if x == max2Depth]
                    # take the second index, since we already took the first index above
                    max2DepthIndex = indices[1]
                else:
                    # find the index for the max depth on the original list
                    max2DepthIndex = anAlleleDepthsList.index(max2Depth)
            else:
                # it wasn't large enough, so just use previous max
                max2DepthIndex = max1DepthIndex
        else:
            # otherwise it's the same as the first
            max2DepthIndex = max1DepthIndex
            
        genotypes = sorted([max1DepthIndex, max2DepthIndex])
        
    return genotypes


def pre_filter_mod_types(aLine, aRefPlusAltList, anInfoDict, aDNANormalDepths, anRNANormalDepths, aDNATumorDepths, anRNATumorDepths, aModMinDepth, aModMinPct, anLohMaxDepth, anLohMaxPct):
    
    aModTypeList = anInfoDict["MT"]
    aModChangeList = anInfoDict["MC"]
    
    # make copies of the lists to manipulate
    modTypesList = list(aModTypeList)
    modChangesList = list(aModChangeList)
    
    filterSet = set()
    
    try:
        # for each modification type and change
        modTypeIndex = 0

        for (modType, modChange) in izip(aModTypeList, aModChangeList):
            # get the source and target alleles
            (source, target) = modChange.split(">")              
            
            # for every germline or somatic call
            if (modType == "GERM" or modType == "SOM"):
                # get the indices
                sourceIndex = aRefPlusAltList.index(source)
                targetIndex = aRefPlusAltList.index(target)
                
                # get the DNA normal total depth
                totalNormalDepth = sum(aDNANormalDepths)

                # if we don't have any normal DNA, then remove the call and add a filter
                if (totalNormalDepth == 0):
                    modTypesList.remove(modType)
                    modChangesList.remove(modChange)
                   
                    if (modType == "SOM"):
                        # add the filters
                        filterSet.add("dnmntb")
                    
                else:
                    # get the normal source depth
                    sourceDepth = aDNANormalDepths[sourceIndex]
                    sourcePct = round(sourceDepth/float(totalNormalDepth), 2)
                    
                    # if the depth doesn't reach minimum, then remove it and add filters
                    if (sourceDepth < aModMinDepth or sourcePct < aModMinPct):
                        modTypesList.remove(modType)
                        modChangesList.remove(modChange)
                       
                        if (modType == "SOM"):
                            # add the filters
                            if (sourceDepth < aModMinDepth):
                                filterSet.add("dnmnrb")
                            if (sourcePct < aModMinPct):
                                filterSet.add("dnmnrp")
                            continue;
                
                    # if this is labeled as Germline, but the dna normal variant depth is not sufficient, maybe it is a somatic one
                    normalTargetDepth = aDNANormalDepths[targetIndex]
                    if (totalNormalDepth > 0):
                        normalTargetPct = round(normalTargetDepth/float(totalNormalDepth), 2)
                    else:
                        normalTargetPct = 0.0
                    
                    if (modType == "GERM" and (normalTargetDepth < aModMinDepth or normalTargetPct < aModMinPct)):
                
                        # get the tumor depth
                        totalTumorDepth = sum(aDNATumorDepths)
                        
                        # this is a hack for ones that were mis-classified as germline by radia, but should've been classified as somatic
                        if (totalTumorDepth > 0):
                            targetDepth = aDNATumorDepths[targetIndex]
                            targetPct = round(targetDepth/float(totalTumorDepth), 2)
                             
                            # if the tumor depth is above the minimum, then add it (these were mis-classified by original radia script)
                            if (targetDepth >= aModMinDepth or targetPct >= aModMinPct):
                                modTypesList.append("SOM")
                                modChangesList.append(modChange)
                        
            elif (modType == "NOR_EDIT"):
                # get the indices
                sourceIndex = aRefPlusAltList.index(source)
                targetIndex = aRefPlusAltList.index(target)
                
                # if this is classified as a normal edit, but the rna normal variant depth is not sufficient, maybe it is a tumor edit
                normalTargetDepth = anRNANormalDepths[targetIndex]
                if (normalTargetDepth < aModMinDepth):
            
                    # get the tumor depth
                    totalTumorDepth = sum(anRNATumorDepths)
                    
                    # this is a hack for ones that were mis-classified as normal edits by radia, but should've been classified as tumor edits
                    if (totalTumorDepth > 0):
                        targetDepth = anRNATumorDepths[targetIndex]
                        targetPct = round(targetDepth/float(totalTumorDepth), 2)
                         
                        # if the tumor depth is above the minimum, then add it (these were mis-classified by original radia script)
                        if (targetDepth >= aModMinDepth or targetPct >= aModMinPct):
                            modTypesList.append("TUM_EDIT")
                            modChangesList.append(modChange)
            
            elif (modType == "TUM_EDIT"):
                # get the indices
                sourceIndex = aRefPlusAltList.index(source)
                targetIndex = aRefPlusAltList.index(target)
                
                # if this is labeled as a tumor edit, but there are variant reads in the DNA, maybe it is a somatic one                
                # get the tumor depth
                totalTumorDepth = sum(aDNATumorDepths)
                
                # this is a hack for ones that were possibly mis-classified as tumor edits by radia, but could've been classified as somatic
                if (totalTumorDepth > 0):
                    targetDepth = aDNATumorDepths[targetIndex]
                    targetPct = round(targetDepth/float(totalTumorDepth), 2)
                     
                    # if the tumor depth is above the minimum, then add it (these were mis-classified by original radia script)
                    if (targetDepth >= aModMinDepth or targetPct >= aModMinPct):
                        modTypesList.append("SOM")
                        modChangesList.append(modChange)

            elif (modType == "LOH"):
                
                # get the total depths
                totalNormalDepths = sum(aDNANormalDepths)
                totalTumorDepths = sum(aDNATumorDepths)
                
                validSources = []
                
                # for each of the possible sources
                for sourceNormal in source:
                    # if the source is the target, then just skip it, because it hasn't been lost
                    if (source == target):
                        continue
                    
                    # get the source in the normal and tumor
                    sourceIndex = aRefPlusAltList.index(sourceNormal)
                    sourceNormalDepth = aDNANormalDepths[sourceIndex]
                    sourceTumorDepth = aDNATumorDepths[sourceIndex]
                    sourceNormalPct = round(sourceNormalDepth/float(totalNormalDepths), 2)
                    sourceTumorPct = round(sourceTumorDepth/float(totalTumorDepths), 2)
                    
                    # if the normal depth is above the minimum, and the tumor depth is below the maximum, then it's valid
                    if (sourceNormalDepth >= aModMinDepth and sourceNormalPct >= aModMinPct and 
                        sourceTumorDepth <= anLohMaxDepth and sourceTumorPct <= anLohMaxPct):
                        validSources.append(sourceNormal)
                    #else:
                    elif (modType == "SOM"):
                        if (sourceNormalDepth < aModMinDepth):
                            filterSet.add("dnmnrb")
                        if (sourceTumorDepth < anLohMaxDepth):
                            filterSet.add("dtmnrb")
                        if (sourceNormalPct < aModMinPct):
                            filterSet.add("dnmnrp")
                        if (sourceTumorPct < anLohMaxPct):
                            filterSet.add("dtmnrp")
                  
                # if there is a valid source, then change the modType      
                if (len(validSources) > 0): 
                    # remove them
                    modTypesList.remove(modType)
                    modChangesList.remove(modChange)
                    
                    # create a new one with just the valid sources
                    newModChange = "".join(validSources) + ">" + target
                                    
                    # add them back                  
                    modTypesList.append(modType)
                    modChangesList.append(newModChange)
                else:
                    # remove them
                    modTypesList.remove(modType)
                    modChangesList.remove(modChange)
            modTypeIndex += 1
        
    except:
        logging.error("Error: aModTypeList=%s, aModChangeList=%s, aRefPlusAltList=%s, aDNANormalDepths=%s, anRNANormalDepths=%s, aDNATumorDepths=%s, anRNATumorDepths=%s", str(aModTypeList), str(aModChangeList), str(aRefPlusAltList), str(aDNANormalDepths), str(anRNANormalDepths), str(aDNATumorDepths), str(anRNATumorDepths))
        logging.error("UnexpectedError: %s" + sys.exc_info())

    # if there are still some valid mod types, then return them
    if (len(modTypesList) > 0):
        anInfoDict["MT"] = modTypesList
        anInfoDict["MC"] = modChangesList
        return (anInfoDict, set())
    else:
        # otherwise just return the originals, they will get filtered later anyway
        return (anInfoDict, filterSet)


def get_final_mod_type(aRefPlusAltList, anInfoDict, aDNANormalAFs, anRNANormalAFs, aDNATumorAFs, anRNATumorAFs):
    
    aModTypeList = anInfoDict["MT"]
    aModChangeList = anInfoDict["MC"]
    
    # get the final mod type        
    finalModTypeList = list()
    finalModChangeList = list()
    # if there is only one, then use it
    if len(aModTypeList) == 1:
        finalModTypeList = aModTypeList
        finalModChangeList = aModChangeList
    else:
        # at this point, there are multiple events that pass all the filters
        # in this case, pick the passing event in the following order:  GERM, NOR_EDIT, SOM, TUM_EDIT, LOH
        if ("GERM" in aModTypeList):
            finalModTypeList = ["GERM"]
        elif ("NOR_EDIT" in aModTypeList):
            finalModTypeList = ["NOR_EDIT"]
            if ("TUM_EDIT" in aModTypeList):
                finalModTypeList.append("TUM_EDIT")
        elif ("SOM" in aModTypeList):
            finalModTypeList = ["SOM"]
        elif ("TUM_EDIT" in aModTypeList):
            finalModTypeList = ["TUM_EDIT"]
        elif ("LOH" in aModTypeList):
            finalModTypeList = ["LOH"]
        
        for modType in finalModTypeList:
            modIndex = aModTypeList.index(modType)
            modChange = aModChangeList[modIndex]
            finalModChangeList.append(modChange)
        
    finalModType = finalModTypeList[0]
        
    # set the final call
    anInfoDict["SS"] = []
    anInfoDict["SOMATIC"] = []
    if (finalModType == "GERM"):
        anInfoDict["SS"].append("1")
    elif (finalModType == "SOM"):
        anInfoDict["SS"].append("2")
        anInfoDict["SOMATIC"].append("True")
    elif (finalModType.find("EDIT") != -1):
        anInfoDict["SS"].append("5")
    # other
    else:
        anInfoDict["SS"].append("4")
    
    anInfoDict["MT"] = finalModTypeList
    anInfoDict["MC"] = finalModChangeList
    
    return (anInfoDict)


def execute_samtools_cmd(aBamFile, aFastaFile, aBaseQuality, aMappingQuality, aChrom, aUseChrPrefix, aStartCoordinate, aStopCoordinate, anIsDebug):
    '''
    ' This function executes an external command.  The command is the "samtools mpileup" command which returns all 
    ' the information about the sequencing reads for specific coordinates.  There are two things to be careful about
    ' when using the samtools mpileup command.  Some .bam files use the 'chr' prefix when specifying the region to 
    ' select with the -r argument.  If the 'chr' prefix is required, then specify the --useChrPrefix argument and also
    ' make sure that the fasta file that is specified also has the 'chr' prefix.  Here are some examples of the commands
    ' that can be copy/pasted to the command line to view the output:
    '
    ' samtools mpileup -f hg19.fa -Q 10 -q 10 -r 1:10000-20000 myBamfile.bam
    '
    ' aBamFile:            A .bam file to be read from
    ' aFastaFile:          The FASTA file which is needed for the reference base.
    ' aBaseQuality:        The base quality score for the samtools command
    ' aMappingQuality:     The mapping quality score for the samtools command
    ' aChrom:              The chromosome that we are selecting from
    ' aStartCoordinate:    The start coordinate of the selection
    ' aStopCoordinate:     The stop coordinate of the selection
    ' aUseChrPrefix:       Whether the 'chr' should be used in the samtools command
    '''
    
    if (aUseChrPrefix):
        # create the samtools command
        samtoolsSelectStatement = "samtools mpileup -E -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r chr" + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
        #print >> sys.stderr, samtoolsSelectStatement
    else:
        # create the samtools command
        samtoolsSelectStatement = "samtools mpileup -E -f " + aFastaFile + " -Q " + str(aBaseQuality) + " -q " + str(aMappingQuality) + " -r " + aChrom + ":" + str(aStartCoordinate) + "-" + str(aStopCoordinate) + " " + aBamFile
   
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
            logging.debug("%s indels found in %s", str(indelCount), aStringOfRawReads)
            
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

    
def format_bam_output(aChrom, aRefList, anAltList, anOverallAltCountsDict, aStringReads, aStringQualScores, aNumBases, aBaseCountsDict, aQualitySumsOfBasesDict, aPlusStrandCountsDict, anIsDebug):
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
    ' aBaseCountsDict:  A dictionary with the number of bases of each type
    ' aQualitySumsOfBasesDict:  A dictionary with the sum of all quality scores for each type of base
    ' aPlusStrandCountsDict:  The number of bases that occurred on the plus strand
    '''
    
    # initialize the return variables
    sumAltReadSupport = 0
    depths = list()
    readSupports = list()
    baseQuals = list()
    strandBias = list()
    
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
                anOverallAltCountsDict[base] += count
                
            # calculate the allele specific avg quality and plus strand scores
            if (count > 0):
                avgBaseQuality = round(aQualitySumsOfBasesDict[base]/float(count),2)
                avgPlusStrandBias = round(aPlusStrandCountsDict[base]/float(count),2)
            else:
                avgBaseQuality = 0.0
                avgPlusStrandBias = 0.0
                
            baseQuals.append(avgBaseQuality)
            strandBias.append(avgPlusStrandBias)
           
    return (depths, readSupports, baseQuals, strandBias, anOverallAltCountsDict, sumAltReadSupport)

    
def fix_base_qualities(aChrom, aStopCoordinate, aVCFGeneratorPrefix, aVCFGeneratorParamsDict, aRefList, anAltList, anOverallAltCountsDict, aGTMinDepth, aGTMinPct, anIsDebug):
    
    # get the info for executing the samtools mpileup command
    bamFile = aVCFGeneratorParamsDict[aVCFGeneratorPrefix + "Filename"]
    fastaFile = aVCFGeneratorParamsDict[aVCFGeneratorPrefix + "FastaFilename"]
    #fastaFile = "/inside/depot/fa/hg19.fasta"
    baseQualityString = aVCFGeneratorParamsDict[aVCFGeneratorPrefix + "MinBaseQuality"]
    baseQualityInt = int(baseQualityString)
    mappingQuality = aVCFGeneratorParamsDict[aVCFGeneratorPrefix + "MinMappingQuality"]
    useChrPrefixString = aVCFGeneratorParamsDict[aVCFGeneratorPrefix + "UseChrPrefix"]
  
    if (useChrPrefixString == "True"):
        useChrPrefix = True
    else:
        useChrPrefix = False

    if (not os.path.isfile(bamFile)):
        logging.critical("The BAM file specified in the header does not exist: %s", bamFile)
        sys.exit(1)

    if (not os.path.isfile(fastaFile)):
        logging.critical("The FASTA file specified in the header does not exist: %s", fastaFile)
        sys.exit(1)

    pileups = execute_samtools_cmd(bamFile, fastaFile, baseQualityString, mappingQuality, aChrom, useChrPrefix, aStopCoordinate, aStopCoordinate, anIsDebug)
    splitPileups = pileups.split("\n")
    numPileups = len(splitPileups)

    # if there was data for these coordinates
    if (numPileups > 0):
        
        if (anIsDebug):        
            logging.debug("samtools number of lines selected from %s = %s", aStopCoordinate, numPileups)
        
        # for each line representing one coordinate
        for line in splitPileups:
            # if the samtools select statement returns no reads which can happen when the batch size is
            # small and the selection is done in an area with no reads, then a warning message will be
            # returned that starts with "[mpileup]".  We can ignore the message and move on to the next
            # select statement.
            if (line.isspace() or line.startswith("[mpileup]") or line.startswith("<mpileup>")):
                continue;

            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")

            # split the line on the tab
            splitLine = line.split("\t")
        
            if (anIsDebug):    
                logging.debug("Original BAM pileup: %s", line)

            if (len(splitLine) > 1):
                # the coordinate is the second element
                #chr = splitLine[0]
                #coordinate = int(splitLine[1])
                reference = splitLine[2].upper()
                #numOfReads = int(splitLine[3])
                rawReads = splitLine[4]
                rawQualScores = splitLine[5]
            else:
                continue
            
            # convert the raw reads into human-readable reads
            (convertedReads, convertedBaseQuals, numBases, starts, stops, indels) = convert_raw_reads(rawReads, rawQualScores, reference, anIsDebug)
    
            if (anIsDebug):    
                logging.debug("After convert_raw_reads(): %s %s %s %s %s %s %s %s %s", aChrom, aStopCoordinate, reference, numBases, convertedReads, convertedBaseQuals, starts, stops, indels)
            
            # filter out the bases that don't meet the minimum base quality scores
            (convertedReads, convertedBaseQuals, numBases, baseCountsDict, qualitySumsOfBasesDict, plusStrandCountsDict) = filter_by_base_quality(convertedReads, convertedBaseQuals, baseQualityInt, anIsDebug)  
            
            if (anIsDebug):
                logging.debug("After filter_by_base_quality(): %s %s %s %s %s %s %s %s %s %s %s %s", aChrom, aStopCoordinate, reference, numBases, convertedReads, convertedBaseQuals, starts, stops, indels, baseCountsDict, qualitySumsOfBasesDict, plusStrandCountsDict)
            
            # format the bam output         
            (depths, readSupports, baseQuals, strandBias, anOverallAltCountsDict, sumAltReadSupport) = format_bam_output(aChrom, aRefList, anAltList, anOverallAltCountsDict, convertedReads, convertedBaseQuals, numBases, baseCountsDict, qualitySumsOfBasesDict, plusStrandCountsDict, anIsDebug)
            
            if (anIsDebug):
                logging.debug("After format_bam_output(): %s %s %s %s %s %s %s %s %s %s", aChrom, aStopCoordinate, reference, numBases, depths, readSupports, baseQuals, strandBias, anOverallAltCountsDict, sumAltReadSupport)
            
            genotypes = fix_genotypes(aChrom, aRefList, anAltList, depths, aGTMinDepth, aGTMinPct)
            
            outputList = ("/".join(map(str, genotypes)), str(numBases), str(indels), str(starts), str(stops), ",".join(map(str, depths)), ",".join(map(str, readSupports)), ",".join(map(str, baseQuals)), ",".join(map(str, strandBias)))
            returnString = ":".join(outputList)
            
            if (anIsDebug):
                logging.debug("After fix_genotypes(): %s %s", genotypes, returnString)
             
    return (returnString, numBases, sum(qualitySumsOfBasesDict.itervalues()), sum(plusStrandCountsDict.itervalues()), sumAltReadSupport, anOverallAltCountsDict)


def filterByStrandBias(aParamsDict, aSampleDict, aSourceIndex, aTargetIndex):

    isStrandBiased = False
        
    sourceDepth = int(aSampleDict["AD"][aSourceIndex])
    targetDepth = int(aSampleDict["AD"][aTargetIndex])
    sourceStrbias = float(aSampleDict["SB"][aSourceIndex])
    targetStrbias = float(aSampleDict["SB"][aTargetIndex])
    
    # if we have enough depth for both the source and target
    if (sourceDepth >= aParamsDict["MinStrBiasDP"] and targetDepth >= aParamsDict["MinStrBiasDP"]):
        # allow 100% strand bias as long as both the source and target strand bias are 100%
        if ((targetStrbias == 0.0 or targetStrbias == 1.0) and (sourceStrbias == 0.0 or sourceStrbias == 1.0)):
            isStrandBiased = False
        # otherwise, see if the target has a strand bias
        elif (targetStrbias > (aParamsDict["MaxStrandBias"]) or targetStrbias < (1.0 - aParamsDict["MaxStrandBias"])):
            isStrandBiased = True
    # see if the target has a strand bias
    elif (targetDepth >= aParamsDict["MinStrBiasDP"]):
        if (targetStrbias > (aParamsDict["MaxStrandBias"]) or targetStrbias < (1.0 - aParamsDict["MaxStrandBias"])):
            isStrandBiased = True

    return isStrandBiased


def filterByMaxError(aRefPlusAltList, aParamsDict, aSampleDict, aSourceIndex, aTargetIndex, anIncludeTargetAlleles, anIsDebug):
    isMaxError = False
    errorCount = 0
    
    # we only allow a maximum of "other" alleles
    # for most cases, we want to make sure there isn't an over-abundance of a 3rd allele at this position
    # for a somatic mutation, we check to make sure that the normal sample doesn't contain too much of the 
    # target allele - use the anIncludeTargetAlleles flag    
    if (anIncludeTargetAlleles):    
        for allele in aRefPlusAltList:
            alleleIndex = aRefPlusAltList.index(allele)
            # if the allele is not the source, then count it
            if (alleleIndex != aSourceIndex):
                errorCount += int(aSampleDict["AD"][alleleIndex])

    elif (len(aRefPlusAltList) > 2):
        for allele in aRefPlusAltList:
            alleleIndex = aRefPlusAltList.index(allele)
            # if the allele is not the source nor the target, then count it
            if (alleleIndex != aSourceIndex and alleleIndex != aTargetIndex):
                errorCount += int(aSampleDict["AD"][alleleIndex])

    totalDepth = int(aSampleDict["DP"][0])
    
    if ((floor(errorCount/float(totalDepth)*100)/100) > float(aParamsDict["MaxErrPct"])):
        isMaxError = True
    
    if (anIsDebug):
        logging.debug("MaxErrPct: errorCount=%s, totalDepth=%s, errorPct=%s, errorPctParam=%s", str(errorCount), str(totalDepth), str(floor(errorCount/float(totalDepth)*100)/100), str(float(aParamsDict["MaxErrPct"])))
        
    return isMaxError


def get_sample_columns(aFilename, anIsDebug):
    
    # get the file
    i_fileHandler = get_read_fileHandler(aFilename)
    
    for line in i_fileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug):
            logging.debug("Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # if we find the column headers
        elif ("#CHROM" in line):
            columnsLine = line.lstrip("#")
            # the tabs get removed somewhere along the pipeline, so just split on whitespace
            columnsLineSplit = columnsLine.split()
            columnsList = columnsLineSplit[9:len(columnsLineSplit)]
            return columnsList
        
    return


def extract_read_support(aTCGAId, aChrom, aVCFFilename, aHeaderFilename, anOutputFilename, aFilterUsingRNAFlag, anAddOriginFlag, aMinAltAvgBaseQual, aStatsDir, aCmdLineParams, aDnaNormParamsDict, aDnaTumParamsDict, anRnaNormParamsDict, anRnaTumParamsDict, aGTMinDepth, aGTMinPct, aModMinDepth, aModMinPct, anLohMaxDepth, anLohMaxPct, anIsDebug):
    '''
    ' This function filters based on the mpileup read support.
    '
    ' aTCGAId: The TCGA Id for this sample
    ' aChrom: The chromosome being filtered
    ' aVCFFilename: The filename to be filtered
    ' aHeaderFilename: The filename with full header info
    ' anOutputFilename: The output filename that will include the filters
    ' aFilterUsingRNAFlag: If the calls should be filtered by the RNA as well
    ' anAddOriginFlag: If the origin (DNA or RNA) of the call should be added to the INFO tag
    ' aMinAltAvgBaseQual: The minimum alternative allele average base quality
    ' aStatsDir: A directory to output stats on the calls
    ' aCmdLineParams: All the parameters specified by the user
    ' aDnaNormParamsDict: The parameters for the normal DNA
    ' aDnaTumParamsDict: The parameters for the tumor DNA
    ' anRnaNormParamsDict: The parameters for the normal RNA
    ' anRnaTumParamsDict: The parameters for the tumor RNA
    ' aGTMinDepth: The minimum depth needed for the genotype
    ' aGTMinPct: The minimum percent needed for the genotype
    ' aModMinDepth: The minimum depth needed to make a call
    ' aModMinPct: The minimum percent needed to make a call
    ' anLohMaxDepth: 
    ' anLohMaxPct:
    ' aMaxAltDepth:
    ' aMaxAltPct:
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    # initialize some variables
    totalEvents = 0
    includedEvents = 0
    somEventsPassing = 0
    somEventsWithTumorRna = 0
    somEventsWithTumorAltRna = 0
            
    hasAddedHeader = False
    headerLines = "##FILTER=<ID=blat,Description=\"The call did not pass the BLAT filter\">\n"
    headerLines += "##FILTER=<ID=pbias,Description=\"The call did not pass the positional bias filter\">\n"
    headerLines += "##FILTER=<ID=multi,Description=\"There are multiple ALT alleles across all samples at this position\">\n"
    headerLines += "##FILTER=<ID=dnmntb,Description=\"DNA Normal minimum total bases is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=dtmntb,Description=\"DNA Tumor minimum total bases is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rnmntb,Description=\"RNA Normal minimum total bases is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rtmntb,Description=\"RNA Tumor minimum total bases is less than user-specified cut-off\">\n"
    
    headerLines += "##FILTER=<ID=dnmxtb,Description=\"DNA Normal maximum total bases is greater than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=dtmxtb,Description=\"DNA Tumor maximum total bases is greater than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rnmxtb,Description=\"RNA Normal maximum total bases is greater than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rtmxtb,Description=\"RNA Tumor maximum total bases is greater than user-specified cut-off\">\n"
    
    headerLines += "##FILTER=<ID=dnmnab,Description=\"DNA Normal minimum ALT bases is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=dtmnab,Description=\"DNA Tumor minimum ALT bases is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rnmnab,Description=\"RNA Normal minimum ALT bases is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rtmnab,Description=\"RNA Tumor minimum ALT bases is less than user-specified cut-off\">\n"
    
    headerLines += "##FILTER=<ID=dnmnap,Description=\"DNA Normal minimum ALT percentage is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=dtmnap,Description=\"DNA Tumor minimum ALT percentage is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rnmnap,Description=\"RNA Normal minimum ALT percentage is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rtmnap,Description=\"RNA Tumor minimum ALT percentage is less than user-specified cut-off\">\n"
    
    headerLines += "##FILTER=<ID=dnmnbq,Description=\"DNA Normal minimum average ALT base quality is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=dtmnbq,Description=\"DNA Tumor minimum average ALT base quality is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rnmnbq,Description=\"RNA Normal minimum average ALT base quality is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rtmnbq,Description=\"RNA Tumor minimum average ALT base quality is less than user-specified cut-off\">\n"
    
    headerLines += "##FILTER=<ID=dnmnrb,Description=\"DNA Normal minimum REF bases is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=dnmnrp,Description=\"DNA Normal minimum REF percentage is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=dtmnrb,Description=\"DNA Tumor minimum REF bases is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=dtmnrp,Description=\"DNA Tumor minimum REF percentage is less than user-specified cut-off\">\n"
    
    headerLines += "##FILTER=<ID=dnmxerr,Description=\"DNA Normal maximum total ALT percentage attributed to error (sequencing, contamination, etc.) is greater than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=dtmxerr,Description=\"DNA Tumor maximum total ALT percentage attributed to error (sequencing, contamination, etc.) is greater than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rnmxerr,Description=\"RNA Normal maximum total ALT percentage attributed to error (sequencing, contamination, etc.) is greater than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rtmxerr,Description=\"RNA Tumor maximum total ALT percentage attributed to error (sequencing, contamination, etc.) is greater than user-specified cut-off\">\n"
    
    headerLines += "##FILTER=<ID=dnsbias,Description=\"Strand bias, majority of reads supporting ALT are on forward OR reverse strand\">\n"
    headerLines += "##FILTER=<ID=dtsbias,Description=\"Strand bias, majority of reads supporting ALT are on forward OR reverse strand\">\n"
    headerLines += "##FILTER=<ID=rnsbias,Description=\"Strand bias, majority of reads supporting ALT are on forward OR reverse strand\">\n"
    headerLines += "##FILTER=<ID=rtsbias,Description=\"Strand bias, majority of reads supporting ALT are on forward OR reverse strand\">\n"
    
    if (anAddOriginFlag):
        headerLines += "##INFO=<ID=ORIGIN,Number=.,Type=String,Description=\"Where the call originated from, the tumor DNA, RNA, or both\">\n"
    
    # sometimes the header lines are stripped from the file, so get them from a header file if it is specified
    if (aHeaderFilename != None):
        columnsList = get_sample_columns(aHeaderFilename, anIsDebug)
    else:
        columnsList = get_sample_columns(aVCFFilename, anIsDebug)
    
    # get the file
    i_vcfFileHandler = get_read_fileHandler(aVCFFilename)
        
    i_outputFileHandler = None
    if (anOutputFilename):
        i_outputFileHandler = get_write_fileHandler(anOutputFilename)
        
    # for each event in the vcf file 
    for line in i_vcfFileHandler:
        # here are some examples of .vcf lines:
        # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31.4;DP=53;FA=0.04;INDEL=0;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP      
        # GT:DP:INDEL:START:STOP:AD:AF:BQ:SB      0/0:2:0:0:0:2:1.0:35.5:0.0      0/0:1:0:0:0:1:1.0:39.0:1.0      0/1:50:0:2:0:48,2:0.96,0.04:31.65,17.5:0.75,0.5
        
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug):
            logging.debug("VCF Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # if we find the FILTER section, then add the filters from here
        elif ((not hasAddedHeader) and line.startswith("##FILTER")):
            hasAddedHeader = True
            if (i_outputFileHandler != None):        
                i_outputFileHandler.write(headerLines)
                i_outputFileHandler.write(line + "\n")
            else:
                print >> sys.stdout, headerLines.rstrip()
                print >> sys.stdout, line
        
        # if we find the column headers
        elif ("#CHROM" in line):
            if (i_outputFileHandler != None):
                i_outputFileHandler.write(line + "\n")
            else:
                print >> sys.stdout, line
        
        # if we find the vcfGenerator line, then update the params from the user
        elif ("vcfGenerator" in line):
            generatorLine = line[0:(len(line)-1)]
            logging.debug("generatorLine: %s", generatorLine)
            generatorLine = generatorLine[16:len(generatorLine)]
            logging.debug("generatorLine: %s", generatorLine)
            generatorParamsList = generatorLine.split(",")
            generatorParamsDict = {}
            
            # create a dictionary of existing params
            for param in generatorParamsList:
                (key, value) = param.split("=")
                value = value.rstrip(">")
                value = value.lstrip("<")
                generatorParamsDict[key] = value
                            
            # for each new param
            for (paramName, paramValue) in aCmdLineParams.iteritems():
                # don't output the defaults for files that aren't specified
                if (paramName.startswith("dnaNormal") and "DNA_NORMAL" not in columnsList):
                    continue;
                elif (paramName.startswith("rnaNormal") and "RNA_NORMAL" not in columnsList):
                    continue;
                elif (paramName.startswith("dnaTumor") and "DNA_TUMOR" not in columnsList):
                    continue;
                elif (paramName.startswith("rnaTumor") and "RNA_TUMOR" not in columnsList):
                    continue;
                # add new params and overwrite the old params with the new ones 
                else:
                    generatorParamsDict[paramName] = paramValue
                    
            generatorOutput = "##vcfGenerator=<"
            # make it pretty by sorting on the keys
            for (paramName) in sorted(generatorParamsDict.iterkeys()):
                paramValue = generatorParamsDict[paramName]
                if (paramValue != None):
                    generatorOutput += paramName + "=<" + str(paramValue) + ">,"
            generatorOutput = generatorOutput.rstrip(",")
            generatorOutput += ">\n"
    
            if (i_outputFileHandler != None):        
                i_outputFileHandler.write(generatorOutput)
            else:
                print >> sys.stdout, generatorOutput
                
        # these lines are from previous scripts in the pipeline, so output them    
        elif (line.startswith("#")):
            if (i_outputFileHandler != None):
                i_outputFileHandler.write(line + "\n")
            else:
                print >> sys.stdout, line
        # now we are to the data
        else:    
            
            # count the total events    
            totalEvents += 1
            somEventWithTumorRna = False
            somEventWithTumorAltRna = False
            
            # split the line on the tab
            splitLine = line.split("\t")
    
            # sample VCF line
            # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31.4;DP=53;FA=0.04;INDEL=0;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP      
            # GT:DP:INDEL:START:STOP:AD:AF:BQ:SB      0/0:2:0:0:0:2:1.0:35.5:0.0      0/0:1:0:0:0:1:1.0:39.0:1.0      0/1:50:0:2:0:48,2:0.96,0.04:31.65,17.5:0.75,0.5
            
            # the coordinate is the second element
            event_chr = splitLine[0]
            event_stopCoordinate = int(splitLine[1])
            event_idList = splitLine[2].split(";")
            event_refList = splitLine[3].split(",")
            event_altList = splitLine[4].split(",")
            event_score = float(splitLine[5])
            
            # if there are no filters so far, then clear the list
            event_filterSet = set(splitLine[6].split(";"))
            if (len(event_filterSet) == 1 and "PASS" in event_filterSet):
                event_filterSet = set()
            
            # parse the info column and create a dict
            event_infoList = splitLine[7].split(";")
            event_infoDict = collections.defaultdict(list)
            for info in event_infoList:
                keyValueList = info.split("=")
                # some keys are just singular without a value (e.g. DB, SOMATIC, etc.)
                if (len(keyValueList) == 1):
                    event_infoDict[keyValueList[0]] = ["True"]
                else:
                    # the value can be a comma separated list
                    event_infoDict[keyValueList[0]] = keyValueList[1].split(",")
            
            # if we should add the origin to the info column
            if (anAddOriginFlag):
                if (aFilterUsingRNAFlag):
                    origin = "RNA"
                else:
                    origin = "DNA"
                            
                if ("ORIGIN" in event_infoDict):
                    originList = event_infoDict["ORIGIN"]
                    if (origin not in originList):
                        originList.append(origin)
                else:
                    event_infoDict["ORIGIN"] = [origin]
                    
            # get the event format list
            event_formatList = splitLine[8].split(":")
            
            # initialize the optional columns to none
            event_dnaNormalList = None
            event_dnaTumorList = None
            event_rnaNormalList = None
            event_rnaTumorList = None

            # if we have a 9th column, figure out which dataset it is
            if (len(splitLine) > 9):
                if (columnsList[0] == "DNA_NORMAL"):
                    event_dnaNormalList = splitLine[9].split(":")
                elif (columnsList[0] == "RNA_NORMAL"):
                    event_rnaNormalList = splitLine[9].split(":")
                elif (columnsList[0] == "DNA_TUMOR"):
                    event_dnaTumorList = splitLine[9].split(":")
                elif (columnsList[0] == "RNA_TUMOR"):
                    event_rnaTumorList = splitLine[9].split(":")
            # if we have a 10th column, figure out which dataset it is
            if (len(splitLine) > 10):
                if (columnsList[1] == "RNA_NORMAL"):
                    event_rnaNormalList = splitLine[10].split(":")
                elif (columnsList[1] == "DNA_TUMOR"):
                    event_dnaTumorList = splitLine[10].split(":")
                elif (columnsList[1] == "RNA_TUMOR"):
                    event_rnaTumorList = splitLine[10].split(":") 
            # if we have a 11th column, figure out which dataset it is
            if (len(splitLine) > 11):
                if (columnsList[2] == "DNA_TUMOR"):
                    event_dnaTumorList = splitLine[11].split(":")
                elif (columnsList[2] == "RNA_TUMOR"):
                    event_rnaTumorList = splitLine[11].split(":")
            # if we have a 12th column, figure out which dataset it is
            if (len(splitLine) > 12):
                if (columnsList[3] == "RNA_TUMOR"):
                    event_rnaTumorList = splitLine[12].split(":")
            
            haveDnaNormData = True
            haveRnaNormData = True
            haveDnaTumData = True
            haveRnaTumData = True
            
            # if there is no data, then set the flag
            if (event_dnaNormalList == None or event_dnaNormalList[0] == "."):
                haveDnaNormData = False
            # if there is no data, then set the flag
            if (event_rnaNormalList == None or event_rnaNormalList[0] == "."):
                haveRnaNormData = False
            # if there is no data, then set the flag
            if (event_dnaTumorList == None or event_dnaTumorList[0] == "."):
                haveDnaTumData = False
            # if there is no data, then set the flag
            if (event_rnaTumorList == None or event_rnaTumorList[0] == "."):
                haveRnaTumData = False
                         
            # parse the dna and rna columns and create dicts for each
            event_dnaNormalDict = collections.defaultdict(list)
            event_dnaTumorDict = collections.defaultdict(list)
            event_rnaNormalDict = collections.defaultdict(list)
            event_rnaTumorDict = collections.defaultdict(list)
            
            index = 0
            for formatItem in event_formatList:
                if (formatItem == "GT"):
                    sep = "/"
                else:
                    sep = ","
                    
                if (haveDnaNormData):
                    dnaNormalItem = event_dnaNormalList[index]
                    event_dnaNormalDict[formatItem] = dnaNormalItem.split(sep)
                if (haveRnaNormData):
                    rnaNormalItem = event_rnaNormalList[index]
                    event_rnaNormalDict[formatItem] = rnaNormalItem.split(sep)
                if (haveDnaTumData):
                    dnaTumorItem = event_dnaTumorList[index]
                    event_dnaTumorDict[formatItem] = dnaTumorItem.split(sep)
                if (haveRnaTumData):
                    rnaTumorItem = event_rnaTumorList[index]
                    event_rnaTumorDict[formatItem] = rnaTumorItem.split(sep)
                index += 1
            
            # fix the original genotypes
            genotypeIndex = event_formatList.index("GT")    
            if (haveDnaNormData):
                event_dnaNormalDict["GT"] = fix_genotypes(event_chr, event_refList, event_altList, map(int, event_dnaNormalDict["AD"]), aGTMinDepth, aGTMinPct)
                event_dnaNormalList[genotypeIndex] = "/".join(map(str, event_dnaNormalDict["GT"]))
            if (haveRnaNormData):
                event_rnaNormalDict["GT"] = fix_genotypes(event_chr, event_refList, event_altList, map(int, event_rnaNormalDict["AD"]), aGTMinDepth, aGTMinPct)
                event_rnaNormalList[genotypeIndex] = "/".join(map(str, event_rnaNormalDict["GT"]))
            if (haveDnaTumData):
                event_dnaTumorDict["GT"] = fix_genotypes(event_chr, event_refList, event_altList, map(int, event_dnaTumorDict["AD"]), aGTMinDepth, aGTMinPct)
                event_dnaTumorList[genotypeIndex] = "/".join(map(str, event_dnaTumorDict["GT"]))
            if (haveRnaTumData):
                event_rnaTumorDict["GT"] = fix_genotypes(event_chr, event_refList, event_altList, map(int, event_rnaTumorDict["AD"]), aGTMinDepth, aGTMinPct)
                event_rnaTumorList[genotypeIndex] = "/".join(map(str, event_rnaTumorDict["GT"]))
            
            # combine the refs and alts in one list
            refPlusAltList = event_refList + event_altList
            
            # get rid of bad mod types that don't meet the minimum requirements
            (event_infoDict, filterSet) = pre_filter_mod_types(line, refPlusAltList, event_infoDict, map(int, event_dnaNormalDict["AD"]), map(int, event_rnaNormalDict["AD"]), map(int, event_dnaTumorDict["AD"]), map(int, event_rnaTumorDict["AD"]), aModMinDepth, aModMinPct, anLohMaxDepth, anLohMaxPct)
            
            # make copies of the lists to manipulate
            modTypesList = list(event_infoDict["MT"])
            modChangesList = list(event_infoDict["MC"]) 
    
            # for each modification type and change
            for (modType, modChange) in izip(event_infoDict["MT"], event_infoDict["MC"]):
                isValidMod = True
                
                # get the source and target alleles
                (source, target) = modChange.split(">")   
                
                if (modType == "GERM"):
                    sourceIndex = refPlusAltList.index(source)
                    targetIndex = refPlusAltList.index(target)
                    
                    # check to make sure the normal DNA sample is between the min and the max of total bases
                    if (int(event_dnaNormalDict["DP"][0]) < aDnaNormParamsDict["MinTotalNumBases"]):
                        isValidMod = False
                        filterSet.add("dnmntb")
                            
                    elif (int(event_dnaNormalDict["DP"][0]) > aDnaNormParamsDict["MaxTotalNumBases"]):
                        isValidMod = False
                        filterSet.add("dnmxtb")
                        
                    # check to make sure the normal DNA sample number of ALT bases is above the min
                    if (int(event_dnaNormalDict["AD"][targetIndex]) < aDnaNormParamsDict["MinAltNumBases"]):
                        isValidMod = False
                        filterSet.add("dnmnab")
                        
                    # check to make sure the normal DNA sample percentage of ALT bases is above the min
                    if (float(event_dnaNormalDict["AF"][targetIndex]) < aDnaNormParamsDict["MinAltPct"]):
                        isValidMod = False
                        filterSet.add("dnmnap")
                        
                    # check to make sure the normal DNA sample average base quality for ALT bases is above the min
                    if (float(event_dnaNormalDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                        isValidMod = False
                        filterSet.add("dnmnbq")
                       
                    # check to make sure the normal variant reads don't have a strand bias
                    if (filterByStrandBias(aDnaNormParamsDict, event_dnaNormalDict, sourceIndex, targetIndex)):   
                        isValidMod = False
                        filterSet.add("dnsbias")
                    
                    # we want to make sure the normal DNA sample error percentage is below the max
                    # we want to make sure that the percentage of other ALTs in this sample is below the max error
                    # this is the same as making sure that the percentage of the source and target alleles is above one minus the max error percentage
                    if (filterByMaxError(refPlusAltList, aDnaNormParamsDict, event_dnaNormalDict, sourceIndex, targetIndex, False, anIsDebug)):
                        isValidMod = False
                        filterSet.add("dnmxerr")
                    
                    # if we are also filtering using the RNA
                    if (aFilterUsingRNAFlag):
                        # check to make sure the normal RNA sample has data and is between the min and the max of total bases
                        if (haveRnaNormData):
                            if (int(event_rnaNormalDict["DP"][0]) < anRnaNormParamsDict["MinTotalNumBases"]):
                                isValidMod = False
                                filterSet.add("rnmntb")
                            elif (int(event_rnaNormalDict["DP"][0]) > anRnaNormParamsDict["MaxTotalNumBases"]):
                                isValidMod = False
                                filterSet.add("rnmxtb")
    
                            # check to make sure the normal RNA sample number of ALT bases is above the min
                            if (int(event_rnaNormalDict["AD"][targetIndex]) < anRnaNormParamsDict["MinAltNumBases"]):
                                isValidMod = False
                                filterSet.add("rnmnab")
                                
                            # check to make sure the normal RNA sample percentage of ALT bases is above the min
                            if (float(event_rnaNormalDict["AF"][targetIndex]) < anRnaNormParamsDict["MinAltPct"]):
                                isValidMod = False
                                filterSet.add("rnmnap")
                                
                            # check to make sure the normal RNA sample average base quality for ALT bases is above the min
                            if (float(event_rnaNormalDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                                isValidMod = False
                                filterSet.add("rnmnbq")
                            
                            # check to make sure the normal variant reads don't have a strand bias
                            if (filterByStrandBias(anRnaNormParamsDict, event_rnaNormalDict, sourceIndex, targetIndex)):   
                                isValidMod = False
                                filterSet.add("rnsbias")
                            
                            # we want to make sure the normal DNA sample error percentage is below the max
                            # we want to make sure that the percentage of other ALTs in this sample is below the max error
                            if (filterByMaxError(refPlusAltList, anRnaNormParamsDict, event_rnaNormalDict, sourceIndex, targetIndex, False, anIsDebug)):
                                isValidMod = False
                                filterSet.add("rnmxerr")
                                
                        # else if a minimum amount of total bases were required, but none were found, then set the filter
                        elif (anRnaNormParamsDict["MinTotalNumBases"] > 0):
                            isValidMod = False
                            filterSet.add("rnmntb")
                            
                elif (modType.find("NOR_EDIT") != -1 or modType.find("RNA_NOR_VAR") != -1):
                    sourceIndex = refPlusAltList.index(source)
                    targetIndex = refPlusAltList.index(target)
                    
                    # if we are also filtering using the RNA
                    if (aFilterUsingRNAFlag):
                        # if this is a normal edit, then we need to check the DNA
                        # if this is an RNA normal variant, then there isn't any DNA to check
                        if (modType.find("NOR_EDIT") != -1):
                            # check to make sure the normal DNA sample has data and is between the min and the max of total bases
                            if (haveDnaNormData):
                                if (int(event_dnaNormalDict["DP"][0]) < aDnaNormParamsDict["MinTotalNumBases"]):
                                    isValidMod = False
                                    filterSet.add("dnmntb")
                                elif (int(event_dnaNormalDict["DP"][0]) > aDnaNormParamsDict["MaxTotalNumBases"]):
                                    isValidMod = False
                                    filterSet.add("dnmxtb")
                                # we want to make sure the normal DNA sample error percentage is below the max
                                # we want to make sure that the percentage of other ALTs in this sample is below the max error
                                if (filterByMaxError(refPlusAltList, aDnaNormParamsDict, event_dnaNormalDict, sourceIndex, targetIndex, True, anIsDebug)):
                                    isValidMod = False
                                    filterSet.add("dnmxerr")
                            # else if a minimum amount of total bases were required, but none were found, then set the filter
                            elif (aDnaNormParamsDict["MinTotalNumBases"] > 0):
                                isValidMod = False
                                filterSet.add("dnmntb")
                                
                        # check to make sure the normal RNA sample is between the min and the max of total bases
                        if (int(event_rnaNormalDict["DP"][0]) < anRnaNormParamsDict["MinTotalNumBases"]):
                            isValidMod = False
                            filterSet.add("rnmntb")
                        elif (int(event_rnaNormalDict["DP"][0]) > anRnaNormParamsDict["MaxTotalNumBases"]):
                            isValidMod = False
                            filterSet.add("rnmxtb")
                            
                        # check to make sure the tumor RNA sample number of ALT bases is above the min
                        if (int(event_rnaNormalDict["AD"][targetIndex]) < anRnaNormParamsDict["MinAltNumBases"]):
                            isValidMod = False
                            filterSet.add("rnmnab")
                        
                        # check to make sure the tumor DNA sample percentage of ALT bases is above the min
                        if (float(event_rnaNormalDict["AF"][targetIndex]) < anRnaNormParamsDict["MinAltPct"]):
                            isValidMod = False
                            filterSet.add("rnmnap")
                            
                        # check to make sure the tumor DNA sample average base quality for ALT bases is above the min
                        if (float(event_rnaNormalDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                            isValidMod = False
                            filterSet.add("rnmnbq")
                        
                        # check to make sure the tumor variant reads don't have a strand bias
                        if (filterByStrandBias(anRnaNormParamsDict, event_rnaNormalDict, sourceIndex, targetIndex)):   
                            isValidMod = False
                            filterSet.add("rnsbias")
                        
                        # we want to make sure the tumor RNA sample error percentage is below the max
                        # we want to make sure that the percentage of other ALTs in this sample is below the max error
                        if (filterByMaxError(refPlusAltList, anRnaNormParamsDict, event_rnaNormalDict, sourceIndex, targetIndex, False, anIsDebug)):
                            isValidMod = False
                            filterSet.add("rnmxerr")
                    else:
                        isValidMod = False
                        filterSet.add("rnmntb")
                        
                elif (modType == "SOM"):
                    sourceIndex = refPlusAltList.index(source)
                    targetIndex = refPlusAltList.index(target)
                            
                    # check to make sure the tumor DNA sample is between the min and the max of total bases
                    if (int(event_dnaTumorDict["DP"][0]) < aDnaTumParamsDict["MinTotalNumBases"]):
                        isValidMod = False
                        filterSet.add("dtmntb")
                    elif (int(event_dnaTumorDict["DP"][0]) > aDnaTumParamsDict["MaxTotalNumBases"]):
                        isValidMod = False
                        filterSet.add("dtmxtb")
            
                    # check to make sure the tumor DNA sample number of ALT bases is above the min
                    if (int(event_dnaTumorDict["AD"][targetIndex]) < aDnaTumParamsDict["MinAltNumBases"]):
                        isValidMod = False
                        filterSet.add("dtmnab")
                    
                    # check to make sure the tumor DNA sample percentage of ALT bases is above the min
                    if (float(event_dnaTumorDict["AF"][targetIndex]) < aDnaTumParamsDict["MinAltPct"]):
                        isValidMod = False
                        filterSet.add("dtmnap")
                        
                    # check to make sure the tumor DNA sample average base quality for ALT bases is above the min
                    if (float(event_dnaTumorDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                        isValidMod = False
                        filterSet.add("dtmnbq")
                    
                    # check to make sure the tumor variant reads don't have a strand bias
                    if (filterByStrandBias(aDnaTumParamsDict, event_dnaTumorDict, sourceIndex, targetIndex)):   
                        isValidMod = False
                        filterSet.add("dtsbias")
                    
                    # we want to make sure the tumor DNA sample error percentage is below the max
                    # we want to make sure that the percentage of other ALTs in this sample is below the max error
                    if (filterByMaxError(refPlusAltList, aDnaTumParamsDict, event_dnaTumorDict, sourceIndex, targetIndex, False, anIsDebug)):
                        isValidMod = False
                        filterSet.add("dtmxerr")
                    
                    # check to make sure the normal DNA sample has data and is between the min and the max of total bases
                    if (haveDnaNormData):
                        if (int(event_dnaNormalDict["DP"][0]) < aDnaNormParamsDict["MinTotalNumBases"]):
                            isValidMod = False
                            filterSet.add("dnmntb")
                        elif (int(event_dnaNormalDict["DP"][0]) > aDnaNormParamsDict["MaxTotalNumBases"]):
                            isValidMod = False
                            filterSet.add("dnmxtb")
                        # we want to make sure the normal DNA sample error percentage is below the max
                        # we want to make sure that the percentage of other ALTs in this sample is below the max error
                        if (filterByMaxError(refPlusAltList, aDnaNormParamsDict, event_dnaNormalDict, sourceIndex, targetIndex, True, anIsDebug)):
                            isValidMod = False
                            filterSet.add("dnmxerr")        
                    elif (aDnaNormParamsDict["MinTotalNumBases"] > 0):
                        isValidMod = False
                        filterSet.add("dnmntb")
                                                            
                    # set some flags
                    if (haveRnaTumData):
                        # check if there is any RNA
                        if (int(event_rnaTumorDict["DP"][0]) > 1):
                            somEventWithTumorRna = True
                        # check if there are any Alt RNA
                        if (int(event_rnaTumorDict["AD"][targetIndex]) > 1):
                            somEventWithTumorAltRna = True
                            
                    # if we are also filtering using the RNA
                    if (aFilterUsingRNAFlag):
                        # check to make sure the tumor RNA sample has data and is between the min and the max of total bases
                        if (haveRnaTumData):
                            if (int(event_rnaTumorDict["DP"][0]) < anRnaTumParamsDict["MinTotalNumBases"]):
                                isValidMod = False
                                filterSet.add("rtmntb")
                            elif (int(event_rnaTumorDict["DP"][0]) > anRnaTumParamsDict["MaxTotalNumBases"]):
                                isValidMod = False
                                filterSet.add("rtmxtb")
                        
                            # check to make sure the tumor RNA sample number of ALT bases is above the min
                            if (int(event_rnaTumorDict["AD"][targetIndex]) < anRnaTumParamsDict["MinAltNumBases"]):
                                isValidMod = False
                                filterSet.add("rtmnab")
                               
                            # check to make sure the tumor RNA sample percentage of ALT bases is above the min
                            if (float(event_rnaTumorDict["AF"][targetIndex]) < anRnaTumParamsDict["MinAltPct"]):
                                isValidMod = False
                                filterSet.add("rtmnap")
                                
                            # check to make sure the tumor RNA sample average base quality for ALT bases is above the min
                            if (float(event_rnaTumorDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                                isValidMod = False
                                filterSet.add("rtmnbq")
                    
                            # check to make sure the tumor variant reads don't have a strand bias
                            if (filterByStrandBias(anRnaTumParamsDict, event_rnaTumorDict, sourceIndex, targetIndex)):   
                                isValidMod = False
                                filterSet.add("rtsbias")
                            
                            # we want to make sure the tumor RNA sample error percentage is below the max
                            # we want to make sure that the percentage of other ALTs in this sample is below the max error
                            if (filterByMaxError(refPlusAltList, anRnaTumParamsDict, event_rnaTumorDict, sourceIndex, targetIndex, False, anIsDebug)):
                                isValidMod = False
                                filterSet.add("rtmxerr")

                        # else if a minimum amount of total bases were required, but none were found, then set the filter                            
                        elif (anRnaTumParamsDict["MinTotalNumBases"] > 0):
                            isValidMod = False
                            filterSet.add("rtmntb")
                    
                elif (modType.find("TUM_EDIT") != -1 or modType.find("RNA_TUM_VAR") != -1):
                    sourceIndex = refPlusAltList.index(source)
                    targetIndex = refPlusAltList.index(target)
                    
                    # if we are also filtering using the RNA
                    if (aFilterUsingRNAFlag):
                        # if this is a tumor edit, then we need to check the DNA
                        # if this is an RNA tumor variant, then there isn't any DNA to check
                        if (modType.find("TUM_EDIT") != -1):
                            # check to make sure the normal DNA sample has data and is between the min and the max of total bases
                            if (haveDnaNormData):
                                if (int(event_dnaNormalDict["DP"][0]) < aDnaNormParamsDict["MinTotalNumBases"]):
                                    isValidMod = False
                                    filterSet.add("dnmntb")
                                elif (int(event_dnaNormalDict["DP"][0]) > aDnaNormParamsDict["MaxTotalNumBases"]):
                                    isValidMod = False
                                    filterSet.add("dnmxtb")
                                # we want to make sure the normal DNA sample error percentage is below the max
                                # we want to make sure that the percentage of other ALTs in this sample is below the max error
                                if (filterByMaxError(refPlusAltList, aDnaNormParamsDict, event_dnaNormalDict, sourceIndex, targetIndex, True, anIsDebug)):
                                    isValidMod = False
                                    filterSet.add("dnmxerr")
                            # else if a minimum amount of total bases were required, but none were found, then set the filter
                            elif (aDnaNormParamsDict["MinTotalNumBases"] > 0):
                                isValidMod = False
                                filterSet.add("dnmntb")
                                
                            # check to make sure the tumor DNA sample is between the min and the max of total bases
                            if (haveDnaTumData):
                                if (int(event_dnaTumorDict["DP"][0]) < aDnaTumParamsDict["MinTotalNumBases"]):
                                    isValidMod = False
                                    filterSet.add("dtmntb")
                                elif (int(event_dnaTumorDict["DP"][0]) > aDnaTumParamsDict["MaxTotalNumBases"]):
                                    isValidMod = False
                                    filterSet.add("dtmxtb")
                                # we want to make sure the tumor DNA sample error percentage is below the max
                                # we want to make sure that the percentage of other ALTs in this sample is below the max error
                                if (filterByMaxError(refPlusAltList, aDnaTumParamsDict, event_dnaTumorDict, sourceIndex, targetIndex, True, anIsDebug)):
                                    isValidMod = False
                                    filterSet.add("dtmxerr")
                            # else if a minimum amount of total bases were required, but none were found, then set the filter
                            elif (aDnaTumParamsDict["MinTotalNumBases"] > 0):
                                isValidMod = False
                                filterSet.add("dtmntb")
                                
                        # check to make sure the tumor RNA sample is between the min and the max of total bases
                        if (int(event_rnaTumorDict["DP"][0]) < anRnaTumParamsDict["MinTotalNumBases"]):
                            isValidMod = False
                            filterSet.add("rtmntb")
                        elif (int(event_rnaTumorDict["DP"][0]) > anRnaTumParamsDict["MaxTotalNumBases"]):
                            isValidMod = False
                            filterSet.add("rtmxtb")
                            
                        # check to make sure the tumor RNA sample number of ALT bases is above the min
                        if (int(event_rnaTumorDict["AD"][targetIndex]) < anRnaTumParamsDict["MinAltNumBases"]):
                            isValidMod = False
                            filterSet.add("rtmnab")
                        
                        # check to make sure the tumor DNA sample percentage of ALT bases is above the min
                        if (float(event_rnaTumorDict["AF"][targetIndex]) < anRnaTumParamsDict["MinAltPct"]):
                            isValidMod = False
                            filterSet.add("rtmnap")
                            
                        # check to make sure the tumor DNA sample average base quality for ALT bases is above the min
                        if (float(event_rnaTumorDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                            isValidMod = False
                            filterSet.add("rtmnbq")
                        
                        # check to make sure the tumor variant reads don't have a strand bias
                        if (filterByStrandBias(anRnaTumParamsDict, event_rnaTumorDict, sourceIndex, targetIndex)):   
                            isValidMod = False
                            filterSet.add("rtsbias")
                        
                        # we want to make sure the tumor RNA sample error percentage is below the max
                        # we want to make sure that the percentage of other ALTs in this sample is below the max error
                        if (filterByMaxError(refPlusAltList, anRnaTumParamsDict, event_rnaTumorDict, sourceIndex, targetIndex, False, anIsDebug)):
                            isValidMod = False
                            filterSet.add("rtmxerr")
                    else:
                        isValidMod = False
                        filterSet.add("rtmntb")
                
                # if this one is not valid and we have more, 
                # then remove this one and try the next one
                if (not isValidMod):
                    modIndices = range(0, len(modTypesList))
                    for (removeModType, removeModChange, modIndex) in izip(modTypesList, modChangesList, modIndices):
                        if (modType == removeModType and modChange == removeModChange):
                            del modTypesList[modIndex]
                            del modChangesList[modIndex]
                            #modTypesList.remove(modType)
                            #modChangesList.remove(modChange)
                            break;
                     
            # after looping through all of them:  if there are still some valid mod types, then set them in the infoDict and ignore the other filtered calls
            if (len(modTypesList) > 0):
                event_infoDict["MT"] = modTypesList
                event_infoDict["MC"] = modChangesList
            # otherwise add the appropriate filters
            else:
                event_filterSet = event_filterSet.union(filterSet)
                
            # get the final mod type
            event_infoDict = get_final_mod_type(refPlusAltList, event_infoDict, map(float, event_dnaNormalDict["AF"]), map(float, event_rnaNormalDict["AF"]), map(float, event_dnaTumorDict["AF"]), map(float, event_rnaTumorDict["AF"]))              
                
            fixBaseQualities = False
            
            # if all the filters are due to the base qualities, then set the flag to fix the base qualities
            # if there is a filter that is not due to the base qualities, then don't fix the base qualities and move on to the output
            if (len(event_filterSet) == 0):
                fixBaseQualities = False
            else:
                fixBaseQualities = True
                for filterItem in event_filterSet:
                    if (not filterItem.endswith("mnbq")):
                        fixBaseQualities = False
                        break
            
            # don't need this for spikeins, ucec_triplets_mapsplice, and everything going further
            fixBaseQualities = False
            
            # if this is a passing event so far, then let's try to fix the quality scores
            if (anAddOriginFlag and fixBaseQualities and "SOM" in event_infoDict["MT"]):  
                # fix the average base quality scores...which unfortunately affects not only the sample output, but the info tag as well...
                overallAltCountsDict = collections.defaultdict(int)
                totalReadDepth = 0
                totalSumBaseQual = 0
                totalSumStrandBias = 0
                totalAltReadDepth = 0  
                if (haveDnaNormData):
                    (dnaNormalOutput, numBases, totalBaseQual, totalStrandBias, totalAltReadSupport, overallAltCountsDict) = fix_base_qualities(event_chr, event_stopCoordinate, "dnaNormal", generatorParamsDict, event_refList, event_altList, overallAltCountsDict, aGTMinDepth, aGTMinPct, anIsDebug)
                    event_dnaNormalList = dnaNormalOutput.split(":")
                    if (numBases > 0):
                        totalReadDepth += numBases
                        totalSumBaseQual += totalBaseQual
                        totalSumStrandBias += totalStrandBias
                        totalAltReadDepth += totalAltReadSupport
                if (haveRnaNormData):
                    (rnaNormalOutput, numBases, totalBaseQual, totalStrandBias, totalAltReadSupport, overallAltCountsDict) = fix_base_qualities(event_chr, event_stopCoordinate, "rnaNormal", generatorParamsDict, event_refList, event_altList, overallAltCountsDict, aGTMinDepth, aGTMinPct, anIsDebug)
                    event_rnaNormalList = rnaNormalOutput.split(":")
                    if (numBases > 0):
                        totalReadDepth += numBases
                        totalSumBaseQual += totalBaseQual
                        totalSumStrandBias += totalStrandBias
                        totalAltReadDepth += totalAltReadSupport
                if (haveDnaTumData):
                    (dnaTumorOutput, numBases, totalBaseQual, totalStrandBias, totalAltReadSupport, overallAltCountsDict) = fix_base_qualities(event_chr, event_stopCoordinate, "dnaTumor", generatorParamsDict, event_refList, event_altList, overallAltCountsDict, aGTMinDepth, aGTMinPct, anIsDebug)
                    event_dnaTumorList = dnaTumorOutput.split(":")
                    if (numBases > 0):
                        totalReadDepth += numBases
                        totalSumBaseQual += totalBaseQual
                        totalSumStrandBias += totalStrandBias
                        totalAltReadDepth += totalAltReadSupport
                if (haveRnaTumData):
                    (rnaTumorOutput, numBases, totalBaseQual, totalStrandBias, totalAltReadSupport, overallAltCountsDict) = fix_base_qualities(event_chr, event_stopCoordinate, "rnaTumor", generatorParamsDict, event_refList, event_altList, overallAltCountsDict, aGTMinDepth, aGTMinPct, anIsDebug)
                    event_rnaTumorList = rnaTumorOutput.split(":")
                    if (numBases > 0):
                        totalReadDepth += numBases
                        totalSumBaseQual += totalBaseQual
                        totalSumStrandBias += totalStrandBias
                        totalAltReadDepth += totalAltReadSupport
                
                # parse the dna and rna columns and create dicts for each
                event_dnaNormalDict = collections.defaultdict(list)
                event_dnaTumorDict = collections.defaultdict(list)
                event_rnaNormalDict = collections.defaultdict(list)
                event_rnaTumorDict = collections.defaultdict(list)
                
                index = 0
                for formatItem in event_formatList:
                    if (formatItem == "GT"):
                        sep = "/"
                    else:
                        sep = ","
                        
                    if (haveDnaNormData):
                        dnaNormalItem = event_dnaNormalList[index]
                        event_dnaNormalDict[formatItem] = dnaNormalItem.split(sep)
                    if (haveRnaNormData):
                        rnaNormalItem = event_rnaNormalList[index]
                        event_rnaNormalDict[formatItem] = rnaNormalItem.split(sep)
                    if (haveDnaTumData):
                        dnaTumorItem = event_dnaTumorList[index]
                        event_dnaTumorDict[formatItem] = dnaTumorItem.split(sep)
                    if (haveRnaTumData):
                        rnaTumorItem = event_rnaTumorList[index]
                        event_rnaTumorDict[formatItem] = rnaTumorItem.split(sep)
                    index += 1
                 
                # clear the old totals
                event_infoDict["AC"] = list()
                event_infoDict["AF"] = list()
                event_infoDict["DP"] = list()
                event_infoDict["BQ"] = list()
                event_infoDict["SB"] = list()
                event_infoDict["FA"] = list()
                
                # re-calculate with new totals    
                # add the alt counts and frequencies in the same order as the alt list 
                for base in event_altList:
                    event_infoDict["AC"].append(str(overallAltCountsDict[base]))
                    event_infoDict["AF"].append(str(round(overallAltCountsDict[base]/float(totalReadDepth),2)))
                
                # re-calculate with new totals
                event_infoDict["DP"].append(str(totalReadDepth))
                if (totalReadDepth > 0):
                    event_infoDict["BQ"].append(str(round(totalSumBaseQual/float(totalReadDepth),2)))
                    event_infoDict["SB"].append(str(round(totalSumStrandBias/float(totalReadDepth),2)))
                    event_infoDict["FA"].append(str(round(totalAltReadDepth/float(totalReadDepth),2)))
                    
                if (anIsDebug):
                    logging.debug("Info Tags Before = %s, After: AC=%s, AF=%s, DP=%s, BQ=%s, SB=%s, FA=%s", event_infoList, event_infoDict["AC"], event_infoDict["AF"], event_infoDict["DP"], event_infoDict["BQ"], event_infoDict["SB"], event_infoDict["FA"])
            
                # now that we have the correct base qualities, try to filter again
                event_filterSet = set()
                
                # for each modification type and change
                for (modType, modChange) in izip(event_infoDict["MT"], event_infoDict["MC"]):
                    
                    # get the source and target alleles
                    (source, target) = modChange.split(">")   
                    
                    if (modType == "GERM"):
                        sourceIndex = refPlusAltList.index(source)
                        targetIndex = refPlusAltList.index(target)
                        
                        # check to make sure the normal DNA sample average base quality for ALT bases is above the min
                        if (float(event_dnaNormalDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                            event_filterSet.add("dnmnbq")
                            
                        # if we are also filtering using the RNA
                        if (aFilterUsingRNAFlag):
                            # check to make sure the normal RNA sample has data and is between the min and the max of total bases
                            if (haveRnaNormData):
                                
                                # check to make sure the normal RNA sample average base quality for ALT bases is above the min
                                if (float(event_rnaNormalDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                                    event_filterSet.add("rnmnbq")
                                    
                    elif (modType == "SOM"):
                        sourceIndex = refPlusAltList.index(source)
                        targetIndex = refPlusAltList.index(target)
                                
                        # check to make sure the tumor DNA sample average base quality for ALT bases is above the min
                        if (float(event_dnaTumorDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                            event_filterSet.add("dtmnbq")
                                                                            
                        # if we are also filtering using the RNA
                        if (aFilterUsingRNAFlag):
                            # check to make sure the tumor RNA sample has data and is between the min and the max of total bases
                            if (haveRnaTumData):
                                # check to make sure the tumor RNA sample average base quality for ALT bases is above the min
                                if (float(event_rnaTumorDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                                    event_filterSet.add("rtmnbq")
                                
                    elif (modType.find("NOR_EDIT") != -1):
                        sourceIndex = refPlusAltList.index(source)
                        targetIndex = refPlusAltList.index(target)
                        
                        # check to make sure the normal RNA sample average base quality for ALT bases is above the min
                        if (float(event_rnaNormalDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                            event_filterSet.add("rnmnbq")
                        
                    elif (modType.find("TUM_EDIT") != -1 or modType.find("RNA_TUM_VAR") != -1):
                        sourceIndex = refPlusAltList.index(source)
                        targetIndex = refPlusAltList.index(target)
                                  
                        # check to make sure the tumor RNA sample average base quality for ALT bases is above the min
                        if (float(event_rnaTumorDict["BQ"][targetIndex]) < aMinAltAvgBaseQual):
                            event_filterSet.add("rtmnbq")
                                
            # create the output list                        
            vcfOutputList = [event_chr, str(event_stopCoordinate)]
                        
            # add the ref, alt, and score
            vcfOutputList.append(";".join(event_idList))
            vcfOutputList.append(",".join(event_refList))
            vcfOutputList.append(",".join(event_altList))
            vcfOutputList.append(str(event_score))  
              
            # if there are no filters thus far, then pass it    
            if (len(event_filterSet) == 0):
                event_filterSet.add("PASS")
                includedEvents += 1
                
                # check to see if this is a passing somatic event
                if ("SOM" in event_infoDict["MT"]):
                    somEventsPassing += 1
                    
                    # check if there is any RNA
                    if (somEventWithTumorRna):
                        somEventsWithTumorRna += 1
                    # check if there is any Alt RNA
                    if (somEventWithTumorAltRna):
                        somEventsWithTumorAltRna += 1
                
            vcfOutputList.append(";".join(event_filterSet))
            
            # add the modified info dict
            infoField = ""
            for key in sorted(event_infoDict.iterkeys()):
                if (len(event_infoDict[key]) == 0):
                    continue
                elif ("True" in event_infoDict[key]):
                    infoField += key + ";"
                else:    
                    infoField += key + "=" + ",".join(event_infoDict[key]) + ";"
            
            vcfOutputList.append(infoField.rstrip(";"))
                
            #vcfOutputList.append(";".join(event_infoList))
            vcfOutputList.append(":".join(event_formatList))
            if (event_dnaNormalList != None):
                vcfOutputList.append(":".join(event_dnaNormalList))
            if (event_rnaNormalList != None):
                vcfOutputList.append(":".join(event_rnaNormalList))
            if (event_dnaTumorList != None):
                vcfOutputList.append(":".join(event_dnaTumorList))
            if (event_rnaTumorList != None):
                vcfOutputList.append(":".join(event_rnaTumorList))

            if (i_outputFileHandler != None):
                i_outputFileHandler.write("\t".join(vcfOutputList) + "\n")
            else:
                print >> sys.stdout, "\t".join(vcfOutputList)
                
    logging.info("Chrom %s and Id %s: %s events passed out of %s total events", aChrom, aTCGAId, includedEvents, totalEvents)
    logging.info("\t".join([aTCGAId, aChrom, str(somEventsPassing), str(somEventsWithTumorRna), str(somEventsWithTumorAltRna)]))
    
    # close the files 
    if (anOutputFilename != None):
        i_outputFileHandler.close()
      
    i_vcfFileHandler.close()    
    return


def main():
    
    # python filterByReadSupportVCF.py TCGA-AB-2995 12 ../data/test/TCGA-AB-2995.vcf --log=DEBUG
    
    # create the usage statement
    usage = "usage: python %prog id chrom vcfFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    # add the optional params    
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, STDOUT by default")
    i_cmdLineParser.add_option("-n", "--headerFilename", dest="headerFilename", metavar="HEADER_FILE", help="the name of the header file")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    i_cmdLineParser.add_option("-s", "--statsDir", dest="statsDir", metavar="STATS_DIR", help="a stats directory where some basic stats can be output")
    i_cmdLineParser.add_option("-r", "--filterUsingRNA", action="store_true", default=False, dest="filterUsingRNA", help="include this argument if the germline and somatic calls should be filtered by the RNA")
    i_cmdLineParser.add_option("-d", "--filterUsingDNA", action="store_true", default=False, dest="filterUsingDNA", help="include this argument if the germline and somatic calls should be filtered by the DNA")
    i_cmdLineParser.add_option("-a", "--addOrigin", action="store_true", default=False, dest="addOrigin", help="include this argument if the origin of the call should be specified in the INFO tags")
        
    i_cmdLineParser.add_option("", "--genotypeMinDepth", type="int", default=int(4), dest="genotypeMinDepth", metavar="GT_MIN_DP", help="the minimum number of bases required for the genotype, %default by default")
    i_cmdLineParser.add_option("", "--genotypeMinPct", type="float", default=float(0.10), dest="genotypeMinPct", metavar="GT_MIN_PCT", help="the minimum percentage of reads required for the genotype, %default by default")
    i_cmdLineParser.add_option("", "--modMinDepth", type="int", default=int(4), dest="modMinDepth", metavar="MOD_MIN_DP", help="the minimum number of bases required for a modification, %default by default")
    i_cmdLineParser.add_option("", "--modMinPct", type="float", default=float(0.10), dest="modMinPct", metavar="MOD_MIN_PCT", help="the minimum percentage of reads required for a modification, %default by default")
    i_cmdLineParser.add_option("", "--lohMaxDepth", type="int", default=int(2), dest="lohMaxDepth", metavar="LOH_MAX_DP", help="the maximum number of bases allowed in the tumor DNA for an LOH, %default by default")
    i_cmdLineParser.add_option("", "--lohMaxPct", type="float", default=float(0.02), dest="lohMaxPct", metavar="LOH_MAX_PCT", help="the maximum percentage of reads in the tumor DNA for an LOH, %default by default")
    #i_cmdLineParser.add_option("", "--maxAltDepth", type="int", default=int(2), dest="maxAltDepth", metavar="MAX_ALT_DP", help="the maximum number of alterative bases allowed in the normal DNA when a call is somatic, %default by default")
    #i_cmdLineParser.add_option("", "--maxAltPct", type="float", default=float(0.01), dest="maxAltPct", metavar="MAX_ALT_PCT", help="the maximum percentage of reads allowed in the normal DNA when a call is somatic, %default by default")
    i_cmdLineParser.add_option("", "--minAltAvgBaseQual", type="float", default=float(20.0), dest="minAltAvgBaseQual", metavar="MIN_ALT_AVG_BQ", help="the minimum average base quality for the ALT reads, %default by default")
    
    i_cmdLineParser.add_option("", "--dnaNormalMinTotalBases", type="int", default=int(10), dest="dnaNormalMinTotalNumBases", metavar="DNA_NOR_MIN_TOTAL_BASES", help="the minimum number of overall normal DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMaxTotalBases", type="int", default=int(20000), dest="dnaNormalMaxTotalNumBases", metavar="DNA_NOR_MAX_TOTAL_BASES", help="the maximum number of overall normal DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinAltBases", type="int", default=int(4), dest="dnaNormalMinAltNumBases", metavar="DNA_NOR_MIN_ALT_BASES", help="the minimum number of alternative normal DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinAltPct", type="float", default=float(0.10), dest="dnaNormalMinAltPct", metavar="DNA_NOR_MIN_ALT_PCT", help="the minimum percentage of alternative normal DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMaxErrPct", type="float", default=float(0.01), dest="dnaNormalMaxErrPct", metavar="DNA_NOR_MAX_ERR_PCT", help="the maximum percentage of alternative normal DNA reads allowed that support a variant, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMaxStrandBias", type="float", default=float(0.99), dest="dnaNormalMaxStrandBias", metavar="DNA_NOR_MAX_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinStrandBiasDepth", type="int", default=float(4), dest="dnaNormalMinStrandBiasDepth", metavar="DNA_NOR_MIN_STRAND_BIAS_DP", help="the minimum total depth needed for the strand bias filter to be applied, %default by default")
    
    i_cmdLineParser.add_option("", "--dnaTumorMinTotalBases", type="int", default=int(10), dest="dnaTumorMinTotalNumBases", metavar="DNA_TUM_MIN_TOTAL_BASES", help="the minimum number of overall tumor DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMaxTotalBases", type="int", default=int(20000), dest="dnaTumorMaxTotalNumBases", metavar="DNA_TUM_MAX_TOTAL_BASES", help="the maximum number of overall tumor DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinAltBases", type="int", default=int(4), dest="dnaTumorMinAltNumBases", metavar="DNA_TUM_MIN_ALT_BASES", help="the minimum number of alternative tumor DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinAltPct", type="float", default=float(0.10), dest="dnaTumorMinAltPct", metavar="DNA_TUM_MIN_ALT_PCT", help="the minimum percentage of alternative tumor DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMaxErrPct", type="float", default=float(0.01), dest="dnaTumorMaxErrPct", metavar="DNA_TUM_MAX_ERR_PCT", help="the maximum percentage of alternative tumor DNA reads allowed that support a variant in the tumor RNA, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMaxStrandBias", type="float", default=float(0.99), dest="dnaTumorMaxStrandBias", metavar="DNA_TUM_MAX_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinStrandBiasDepth", type="int", default=float(4), dest="dnaTumorMinStrandBiasDepth", metavar="DNA_TUM_MIN_STRAND_BIAS_DP", help="the minimum total depth needed for the strand bias filter to be applied, %default by default")
    
    i_cmdLineParser.add_option("", "--rnaNormalMinTotalBases", type="int", default=int(10), dest="rnaNormalMinTotalNumBases", metavar="RNA_NOR_MIN_TOTAL_BASES", help="the minimum number of overall normal RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMaxTotalBases", type="int", default=int(20000), dest="rnaNormalMaxTotalNumBases", metavar="RNA_NOR_MAX_TOTAL_BASES", help="the maximum number of overall normal RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinAltBases", type="int", default=int(4), dest="rnaNormalMinAltNumBases", metavar="RNA_NOR_MIN_ALT_BASES", help="the minimum number of alternative normal RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinAltPct", type="float", default=float(0.10), dest="rnaNormalMinAltPct", metavar="RNA_NOR_MIN_ALT_PCT", help="the minimum percentage of alternative normal RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMaxErrPct", type="float", default=float(0.01), dest="rnaNormalMaxErrPct", metavar="RNA_NOR_MAX_ERR_PCT", help="the maximum percentage of alternative normal RNA reads allowed that support a variant in the tumor, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMaxStrandBias", type="float", default=float(0.99), dest="rnaNormalMaxStrandBias", metavar="RNA_NOR_MAX_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinStrandBiasDepth", type="int", default=float(4), dest="rnaNormalMinStrandBiasDepth", metavar="RNA_NOR_MIN_STRAND_BIAS_DP", help="the minimum total depth needed for the strand bias filter to be applied, %default by default")
    
    i_cmdLineParser.add_option("", "--rnaTumorMinTotalBases", type="int", default=int(10), dest="rnaTumorMinTotalNumBases", metavar="RNA_TUM_MIN_TOTAL_BASES", help="the minimum number of overall tumor RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMaxTotalBases", type="int", default=int(20000), dest="rnaTumorMaxTotalNumBases", metavar="RNA_TUM_MAX_TOTAL_BASES", help="the maximum number of overall tumor RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinAltBases", type="int", default=int(4), dest="rnaTumorMinAltNumBases", metavar="RNA_TUM_MIN_ALT_BASES", help="the minimum number of alternative tumor RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinAltPct", type="float", default=float(0.10), dest="rnaTumorMinAltPct", metavar="RNA_TUM_MIN_ALT_PCT", help="the minimum percentage of alternative tumor RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMaxErrPct", type="float", default=float(0.01), dest="rnaTumorMaxErrPct", metavar="RNA_TUM_MAX_ERR_PCT", help="the maximum percentage of alternative tumor RNA reads allowed that support a variant in a future data type, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMaxStrandBias", type="float", default=float(0.99), dest="rnaTumorMaxStrandBias", metavar="RNA_TUM_MAX_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinStrandBiasDepth", type="int", default=float(4), dest="rnaTumorMinStrandBiasDepth", metavar="RNA_TUM_MIN_STRAND_BIAS_DP", help="the minimum total depth needed for the strand bias filter to be applied, %default by default")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(3,58,1)
    i_argLength = len(sys.argv)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)
        
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_cmdLineOptionsDict = vars(i_cmdLineOptions)
    i_id = str(i_cmdLineArgs[0])
    i_chrom = str(i_cmdLineArgs[1])
    i_vcfFilename = str(i_cmdLineArgs[2])
    
    # get the optional params with default values
    i_logLevel = i_cmdLineOptions.logLevel
    i_addOrigin = i_cmdLineOptions.addOrigin
    i_filterUsingRNA = i_cmdLineOptions.filterUsingRNA
    i_filterUsingDNA = i_cmdLineOptions.filterUsingDNA
    i_genotypeMinDepth = i_cmdLineOptions.genotypeMinDepth
    i_genotypeMinPct = i_cmdLineOptions.genotypeMinPct
    i_modMinDepth = i_cmdLineOptions.modMinDepth
    i_modMinPct = i_cmdLineOptions.modMinPct
    #i_maxAltDepth = i_cmdLineOptions.maxAltDepth
    #i_maxAltPct = i_cmdLineOptions.maxAltPct
    i_lohMaxDepth = i_cmdLineOptions.lohMaxDepth
    i_lohMaxPct = i_cmdLineOptions.lohMaxPct
    i_minAltAvgBaseQual = i_cmdLineOptions.minAltAvgBaseQual
    
    i_dnaNormParams = {}
    i_dnaNormParams["MinTotalNumBases"] = i_cmdLineOptions.dnaNormalMinTotalNumBases
    i_dnaNormParams["MaxTotalNumBases"] = i_cmdLineOptions.dnaNormalMaxTotalNumBases
    i_dnaNormParams["MinAltNumBases"] = i_cmdLineOptions.dnaNormalMinAltNumBases
    i_dnaNormParams["MinAltPct"] = i_cmdLineOptions.dnaNormalMinAltPct
    i_dnaNormParams["MaxErrPct"] = i_cmdLineOptions.dnaNormalMaxErrPct
    i_dnaNormParams["MaxStrandBias"] = i_cmdLineOptions.dnaNormalMaxStrandBias
    i_dnaNormParams["MinStrBiasDP"] = i_cmdLineOptions.dnaNormalMinStrandBiasDepth
    
    i_dnaTumParams = {}
    i_dnaTumParams["MinTotalNumBases"] = i_cmdLineOptions.dnaTumorMinTotalNumBases
    i_dnaTumParams["MaxTotalNumBases"] = i_cmdLineOptions.dnaTumorMaxTotalNumBases
    i_dnaTumParams["MinAltNumBases"] = i_cmdLineOptions.dnaTumorMinAltNumBases
    i_dnaTumParams["MinAltPct"] = i_cmdLineOptions.dnaTumorMinAltPct
    i_dnaTumParams["MaxErrPct"] = i_cmdLineOptions.dnaTumorMaxErrPct
    i_dnaTumParams["MaxStrandBias"] = i_cmdLineOptions.dnaTumorMaxStrandBias
    i_dnaTumParams["MinStrBiasDP"] = i_cmdLineOptions.dnaTumorMinStrandBiasDepth
    
    i_rnaNormParams = {}
    i_rnaNormParams["MinTotalNumBases"] = i_cmdLineOptions.rnaNormalMinTotalNumBases
    i_rnaNormParams["MaxTotalNumBases"] = i_cmdLineOptions.rnaNormalMaxTotalNumBases
    i_rnaNormParams["MinAltNumBases"] = i_cmdLineOptions.rnaNormalMinAltNumBases
    i_rnaNormParams["MinAltPct"] = i_cmdLineOptions.rnaNormalMinAltPct
    i_rnaNormParams["MaxErrPct"] = i_cmdLineOptions.rnaNormalMaxErrPct
    i_rnaNormParams["MaxStrandBias"] = i_cmdLineOptions.rnaNormalMaxStrandBias
    i_rnaNormParams["MinStrBiasDP"] = i_cmdLineOptions.rnaNormalMinStrandBiasDepth
    
    i_rnaTumParams = {}
    i_rnaTumParams["MinTotalNumBases"] = i_cmdLineOptions.rnaTumorMinTotalNumBases
    i_rnaTumParams["MaxTotalNumBases"] = i_cmdLineOptions.rnaTumorMaxTotalNumBases
    i_rnaTumParams["MinAltNumBases"] = i_cmdLineOptions.rnaTumorMinAltNumBases
    i_rnaTumParams["MinAltPct"] = i_cmdLineOptions.rnaTumorMinAltPct
    i_rnaTumParams["MaxErrPct"] = i_cmdLineOptions.rnaTumorMaxErrPct
    i_rnaTumParams["MaxStrandBias"] = i_cmdLineOptions.rnaTumorMaxStrandBias
    i_rnaTumParams["MinStrBiasDP"] = i_cmdLineOptions.rnaTumorMinStrandBiasDepth
    
    # try to get any optional parameters with no defaults    
    i_readFilenameList = [i_vcfFilename]
    i_writeFilenameList = []
    i_dirList = []
    
    i_outputFilename = None
    i_headerFilename = None
    i_logFilename = None
    i_statsDir = None
    if (i_cmdLineOptions.outputFilename != None):
        i_outputFilename = str(i_cmdLineOptions.outputFilename)
        i_writeFilenameList += [i_outputFilename]
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
        i_writeFilenameList += [i_logFilename]
    if (i_cmdLineOptions.headerFilename != None):
        i_headerFilename = str(i_cmdLineOptions.headerFilename)
        i_writeFilenameList += [i_headerFilename]
    if (i_cmdLineOptions.statsDir != None):
        i_statsDir = str(i_cmdLineOptions.statsDir)
        i_dirList += [i_statsDir]
            
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
        logging.debug("headerFilename=%s" % i_headerFilename)
        logging.debug("logLevel=%s" % i_logLevel)
        logging.debug("logFilename=%s" % i_logFilename)
        logging.debug("statsDir=%s" % i_statsDir)
        logging.debug("filterUsingRNA=%s" % i_filterUsingRNA)
        logging.debug("filterUsingDNA=%s" % i_filterUsingDNA)
        logging.debug("addOrigin=%s" % i_addOrigin)
        
        logging.debug("genotypeMinDepth=%s" % i_genotypeMinDepth)
        logging.debug("genotypeMinPct=%s" % i_genotypeMinPct)
        logging.debug("modMinDepth=%s" % i_modMinDepth)
        logging.debug("modMinPct=%s" % i_modMinPct)
        logging.debug("lohMaxDepth=%s" % i_lohMaxDepth)
        logging.debug("lohMaxPct=%s" % i_lohMaxPct)
        logging.debug("minAltAvgBaseQual=%s" % i_minAltAvgBaseQual)
        
        logging.debug("dna normal minTotalBases: %s" % i_dnaNormParams["MinTotalNumBases"])
        logging.debug("dna normal maxTotalBases: %s" % i_dnaNormParams["MaxTotalNumBases"])
        logging.debug("dna normal minAltBases: %s" % i_dnaNormParams["MinAltNumBases"])
        logging.debug("dna normal minAltPct: %s" % i_dnaNormParams["MinAltPct"])
        logging.debug("dna normal maxErrPct: %s" % i_dnaNormParams["MaxErrPct"])
        logging.debug("dna normal strandBias: %s" % i_dnaNormParams["MaxStrandBias"])
        logging.debug("dna normal strandBias minDP: %s" % i_dnaNormParams["MinStrBiasDP"])
        
        logging.debug("dna tumor minTotalBases: %s" % i_dnaTumParams["MinTotalNumBases"])
        logging.debug("dna tumor maxTotalBases: %s" % i_dnaTumParams["MaxTotalNumBases"])
        logging.debug("dna tumor minAltBases: %s" % i_dnaTumParams["MinAltNumBases"])
        logging.debug("dna tumor minAltPct: %s" % i_dnaTumParams["MinAltPct"])
        logging.debug("dna tumor maxErrPct: %s" % i_dnaTumParams["MaxErrPct"])
        logging.debug("dna tumor strandBias: %s" % i_dnaTumParams["MaxStrandBias"])
        logging.debug("dna tumor strandBias minDP: %s" % i_dnaTumParams["MinStrBiasDP"])
        
        logging.debug("rna normal minTotalBases: %s" % i_rnaNormParams["MinTotalNumBases"])
        logging.debug("rna normal maxTotalBases: %s" % i_rnaNormParams["MaxTotalNumBases"])
        logging.debug("rna normal minAltBases: %s" % i_rnaNormParams["MinAltNumBases"])
        logging.debug("rna normal minAltPct: %s" % i_rnaNormParams["MinAltPct"])
        logging.debug("rna normal maxErrPct: %s" % i_rnaNormParams["MaxErrPct"])
        logging.debug("rna normal strandBias: %s" % i_rnaNormParams["MaxStrandBias"])
        logging.debug("rna normal strandBias minDP: %s" % i_rnaNormParams["MinStrBiasDP"])
        
        logging.debug("rna tumor minTotalBases: %s" % i_rnaTumParams["MinTotalNumBases"])
        logging.debug("rna tumor maxTotalBases: %s" % i_rnaTumParams["MaxTotalNumBases"])
        logging.debug("rna tumor minAltBases: %s" % i_rnaTumParams["MinAltNumBases"])
        logging.debug("rna tumor minAltPct: %s" % i_rnaTumParams["MinAltPct"])
        logging.debug("rna tumor maxErrPct: %s" % i_rnaTumParams["MaxErrPct"])
        logging.debug("rna tumor strandBias: %s" % i_rnaTumParams["MaxStrandBias"])
        logging.debug("rna tumor strandBias minDP: %s" % i_rnaTumParams["MinStrBiasDP"])
                    
    # check for any errors
    if (not radiaUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)
    
    extract_read_support(i_id, i_chrom, i_vcfFilename, i_headerFilename, i_outputFilename, i_filterUsingRNA, i_addOrigin, i_minAltAvgBaseQual, i_statsDir, i_cmdLineOptionsDict, i_dnaNormParams, i_dnaTumParams, i_rnaNormParams, i_rnaTumParams, i_genotypeMinDepth, i_genotypeMinPct, i_modMinDepth, i_modMinPct, i_lohMaxDepth, i_lohMaxPct, i_debug)
    return

main()
sys.exit(0)
