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
from math import floor
import gzip

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

# this regular expression is used to extract the ID tag from the INFO, FORMAT, and FILTER fields
i_headerIDRegEx = re.compile("ID=(\\w)*,")
# this regular expression is used to extract the Type tag from the INFO and FORMAT fields
i_headerTypeRegEx = re.compile("Type=(\\w)*,")


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


def pre_filter_mod_types(aRefPlusAltList, anAllFiltersSet, anInfoDict, aDNANormalDepths, anRNANormalDepths, aDNATumorDepths, anRNATumorDepths, aModMinDepth, aModMinPct, anLohMaxDepth, anLohMaxPct):
    
    aModTypeList = anInfoDict["MT"]
    aModChangeList = anInfoDict["MC"]
    
    # make copies of the lists to manipulate
    modTypesList = list(aModTypeList)
    modChangesList = list(aModChangeList)
    
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
                        anAllFiltersSet.add("dnmntb")
                    
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
                                anAllFiltersSet.add("dnmnrb")
                            if (sourcePct < aModMinPct):
                                anAllFiltersSet.add("dnmnrp")
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
                        
            elif (modType == "NOR_EDIT" or modType == "RNA_NOR_VAR"):
                # get the indices
                sourceIndex = aRefPlusAltList.index(source)
                targetIndex = aRefPlusAltList.index(target)
                
                # if this is labeled as a normal edit, but there are variant reads in the normal DNA, maybe it is a germline                
                # get the normal depth
                totalDNANormalDepth = sum(aDNANormalDepths)
                
                # this is a hack for ones that were possibly mis-classified as normal edits by radia, but could've been classified as germline
                if (totalDNANormalDepth > 0):
                    targetDNANormalDepth = aDNANormalDepths[targetIndex]
                    targetDNANormalPct = round(targetDNANormalDepth/float(totalDNANormalDepth), 2)
                     
                    # if the normal depth is above the minimum, then add it (these were mis-classified by original radia script)
                    if (targetDNANormalDepth >= aModMinDepth or targetDNANormalPct >= aModMinPct):
                        modTypesList.append("GERM")
                        modChangesList.append(modChange)
                    # these are calls with normal DNA reads, but not enough variant reads to be considered as germline, maybe they're edits
                    elif (modType == "RNA_NOR_VAR"):
                        modTypesList.append("NOR_EDIT")
                        modChangesList.append(modChange)
                
                # if this is classified as a normal edit, but the rna normal variant depth is not sufficient, maybe it is a tumor edit
                targetRNANormalDepth = anRNANormalDepths[targetIndex]
                if (targetRNANormalDepth < aModMinDepth):
            
                    # get the tumor depth
                    totalRNATumorDepth = sum(anRNATumorDepths)
                    
                    # this is a hack for ones that were mis-classified as normal edits by radia, but should've been classified as tumor edits
                    if (totalRNATumorDepth > 0):
                        targetRNATumorDepth = anRNATumorDepths[targetIndex]
                        targetRNATumorPct = round(targetRNATumorDepth/float(totalRNATumorDepth), 2)
                         
                        # if the tumor depth is above the minimum, then add it (these were mis-classified by original radia script)
                        if (targetRNATumorDepth >= aModMinDepth or targetRNATumorPct >= aModMinPct):
                            modTypesList.append("TUM_EDIT")
                            modChangesList.append(modChange)
                        
            elif (modType == "TUM_EDIT" or modType == "RNA_TUM_VAR"):
                # get the indices
                sourceIndex = aRefPlusAltList.index(source)
                targetIndex = aRefPlusAltList.index(target)
                
                # if this is labeled as a tumor edit, but there are variant reads in the DNA, maybe it is a somatic one                
                # get the tumor depth
                totalDNATumorDepth = sum(aDNATumorDepths)
                
                # this is a hack for ones that were possibly mis-classified as tumor edits by radia, but could've been classified as somatic
                if (totalDNATumorDepth > 0):
                    targetDNATumorDepth = aDNATumorDepths[targetIndex]
                    targetDNATumorPct = round(targetDNATumorDepth/float(totalDNATumorDepth), 2)
                     
                    # if the tumor depth is above the minimum, then add it (these were mis-classified by original radia script)
                    if (targetDNATumorDepth >= aModMinDepth or targetDNATumorPct >= aModMinPct):
                        modTypesList.append("SOM")
                        modChangesList.append(modChange)
                    # these are calls with DNA, but not enough DNA (by default 1 read) to be considered as somatic, maybe they're edits
                    elif (modType == "RNA_TUM_VAR" and targetDNATumorDepth == 0):
                        modTypesList.append("TUM_EDIT")
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
                            anAllFiltersSet.add("dnmnrb")
                        if (sourceTumorDepth < anLohMaxDepth):
                            anAllFiltersSet.add("dtmnrb")
                        if (sourceNormalPct < aModMinPct):
                            anAllFiltersSet.add("dnmnrp")
                        if (sourceTumorPct < anLohMaxPct):
                            anAllFiltersSet.add("dtmnrp")
                  
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
        logging.error("Filtering Error in pre_filter_mod_types(): aModTypeList=%s, aModChangeList=%s, aRefPlusAltList=%s, aDNANormalDepths=%s, anRNANormalDepths=%s, aDNATumorDepths=%s, anRNATumorDepths=%s", str(aModTypeList), str(aModChangeList), str(aRefPlusAltList), str(aDNANormalDepths), str(anRNANormalDepths), str(aDNATumorDepths), str(anRNATumorDepths))
        raise

    # if there are still some valid mod types, then return them
    if (len(modTypesList) > 0):
        anInfoDict["MT"] = modTypesList
        anInfoDict["MC"] = modChangesList
        return (anInfoDict, set())
    else:
        # otherwise just return the originals, they will get filtered later anyway
        return (anInfoDict, anAllFiltersSet)


def get_final_mod_type(anInfoDict, anIsDebug):
    
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
        # in this case, pick the passing event in the following order:  GERM, NOR_EDIT, SOM, TUM_EDIT, RNA_TUM_VAR, LOH
        if ("GERM" in aModTypeList):
            finalModTypeList = ["GERM"]
        elif ("NOR_EDIT" in aModTypeList):
            finalModTypeList = ["NOR_EDIT"]
            if ("TUM_EDIT" in aModTypeList):
                finalModTypeList.append("TUM_EDIT")
        elif ("RNA_NOR_VAR" in aModTypeList):
            finalModTypeList = ["RNA_NOR_VAR"]
            if ("RNA_TUM_VAR" in aModTypeList):
                finalModTypeList.append("RNA_TUM_VAR")
        elif ("SOM" in aModTypeList):
            finalModTypeList = ["SOM"]
        elif ("TUM_EDIT" in aModTypeList):
            finalModTypeList = ["TUM_EDIT"]
        elif ("RNA_TUM_VAR" in aModTypeList):
            finalModTypeList = ["RNA_TUM_VAR"]
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
    # germline
    if (finalModType == "GERM"):
        anInfoDict["SS"].append("1")
    # somatic
    elif (finalModType == "SOM"):
        anInfoDict["SS"].append("2")
        anInfoDict["SOMATIC"].append("True")
    # rna-editing
    elif (finalModType.find("EDIT") != -1):
        anInfoDict["SS"].append("4")
    # unknown
    elif (finalModType.find("RNA_NOR_VAR") != -1 or finalModType.find("RNA_TUM_VAR") != 1):
        anInfoDict["SS"].append("5")
    # unknown
    else:
        anInfoDict["SS"].append("5")
    
    anInfoDict["MT"] = finalModTypeList
    anInfoDict["MC"] = finalModChangeList
    
    return (anInfoDict)


def filterByMapQualZero(aParamsDict, aSampleDict, aTargetIndex):
    
    if "MQ0" in aSampleDict:
        # check to make sure the sample does not exceed the maximum percentage of MQ0 reads supporting the ALT
        mapQualZeroReads = int(aSampleDict["MQ0"][aTargetIndex])
        totalAltReads = int(aSampleDict["AD"][aTargetIndex])
        
        if (totalAltReads > 0):
            if ((floor(mapQualZeroReads/float(totalAltReads)*100)/100) > float(aParamsDict["MaxAltMapQualZeroPct"])):
                return True
    
    return False


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
    
    # if the number of "other" alleles is above the minimum count and
    # the "other" allele percentage is above the max percent
    if ((errorCount >= int(aParamsDict["MinErrPctDP"])) and
        ((floor(errorCount/float(totalDepth)*100)/100) > float(aParamsDict["MaxErrPct"]))):
        isMaxError = True
    
    if (anIsDebug):
        logging.debug("MaxErrPct: isMaxError=%s, errorCount=%s, totalDepth=%s, errorPct=%s, errorPctParam=%s", str(isMaxError), str(errorCount), str(totalDepth), str(floor(errorCount/float(totalDepth)*100)/100), str(float(aParamsDict["MaxErrPct"])))
        
    return isMaxError


def get_sample_columns(aFilename, aHeaderDict, anIsDebug):
    
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
            aHeaderDict["chrom"] = line + "\n"
            columnsLine = line.lstrip("#")
            # the tabs get removed somewhere along the pipeline, so just split on whitespace
            columnsLineSplit = columnsLine.split()
            columnsList = columnsLineSplit[9:len(columnsLineSplit)]
            return columnsList, aHeaderDict
        
    return


def get_meta_id(aLine, anIsDebug):
    metaId = aLine.split("=")[0]
    return metaId


def get_id(aLine, anIsDebug):
    matchObj = i_headerIDRegEx.search(aLine)
    matchId = matchObj.group()
    matchId = matchId.rstrip(",")
    (tag, matchId) = matchId.split("=")
    #matchId = matchId.lstrip("ID=")
    return matchId


def add_header_data(aHeaderDict, aLine, anIsDebug):
    
    if (not aLine.startswith("#")):
        return aHeaderDict
    
    # extract the metadata
    if (aLine.startswith("##FORMAT")):
        formatId = get_id(aLine, anIsDebug)
        aHeaderDict["format"][formatId] = aLine
        
    elif (aLine.startswith("##INFO")):
        infoId = get_id(aLine, anIsDebug)
        aHeaderDict["info"][infoId] = aLine
        
    elif (aLine.startswith("##FILTER")):
        filterId = get_id(aLine, anIsDebug)
        aHeaderDict["filter"][filterId] = aLine
        
    elif (aLine.startswith("##SAMPLE")):
        sampleId = get_id(aLine, anIsDebug)
        aHeaderDict["sample"][sampleId] = aLine
        
    elif (aLine.startswith("##")):
        metadataId = get_meta_id(aLine, anIsDebug)
        aHeaderDict["metadata"][metadataId] = aLine
        
    elif (aLine.startswith("#CHROM")):
        aHeaderDict["chrom"] = aLine
        
    return aHeaderDict


def get_mpileup_header(anAddOriginFlag):
    
    headerDict = dict()
    headerDict["metadata"] = collections.defaultdict(str)
    headerDict["sample"] = collections.defaultdict(str)
    headerDict["format"] = collections.defaultdict(str)
    headerDict["info"] = collections.defaultdict(str)
    headerDict["filter"] = collections.defaultdict(str)
    
    # if we're adding the origin tag, then add the INFO tag
    if (anAddOriginFlag):
        headerDict["info"]["ORIGIN"] = "##INFO=<ID=ORIGIN,Number=.,Type=String,Description=\"Where the call originated from, the tumor DNA, RNA, or both\">\n"
        
    # add all the filters
    headerDict["filter"]["blat"] = "##FILTER=<ID=blat,Description=\"The call did not pass the BLAT filter\">\n"
    headerDict["filter"]["pbias"] = "##FILTER=<ID=pbias,Description=\"A positional bias exists\">\n"
    headerDict["filter"]["multi"] = "##FILTER=<ID=multi,Description=\"There are multiple ALT alleles across all samples at this position\">\n"
    headerDict["filter"]["rnacall"] = "##FILTER=<ID=rnacall,Description=\"This is a dummy filter for a call that originated in the RNA being filtered by the DNA\">\n"
    headerDict["filter"]["dnacall"] = "##FILTER=<ID=dnacall,Description=\"This is a dummy filter for a call that originated in the DNA being filtered by the RNA\">\n"
    
    headerDict["filter"]["dnmntb"] = "##FILTER=<ID=dnmntb,Description=\"DNA Normal total bases is less than the minimum\">\n"
    headerDict["filter"]["dtmntb"] = "##FILTER=<ID=dtmntb,Description=\"DNA Tumor total bases is less than the minimum\">\n"
    headerDict["filter"]["rnmntb"] = "##FILTER=<ID=rnmntb,Description=\"RNA Normal total bases is less than the minimum\">\n"
    headerDict["filter"]["rtmntb"] = "##FILTER=<ID=rtmntb,Description=\"RNA Tumor total bases is less than the minimum\">\n"
    
    headerDict["filter"]["dnmxtb"] = "##FILTER=<ID=dnmxtb,Description=\"DNA Normal total bases is greater than the maximum\">\n"
    headerDict["filter"]["dtmxtb"] = "##FILTER=<ID=dtmxtb,Description=\"DNA Tumor total bases is greater than the maximum\">\n"
    headerDict["filter"]["rnmxtb"] = "##FILTER=<ID=rnmxtb,Description=\"RNA Normal total bases is greater than the maximum\">\n"
    headerDict["filter"]["rtmxtb"] = "##FILTER=<ID=rtmxtb,Description=\"RNA Tumor total bases is greater than the maximum\">\n"
    
    headerDict["filter"]["dnmnab"] = "##FILTER=<ID=dnmnab,Description=\"DNA Normal ALT bases is less than the minimum\">\n"
    headerDict["filter"]["dtmnab"] = "##FILTER=<ID=dtmnab,Description=\"DNA Tumor ALT bases is less than the minimum\">\n"
    headerDict["filter"]["rnmnab"] = "##FILTER=<ID=rnmnab,Description=\"RNA Normal ALT bases is less than the minimum\">\n"
    headerDict["filter"]["rtmnab"] = "##FILTER=<ID=rtmnab,Description=\"RNA Tumor ALT bases is less than the minimum\">\n"
    
    headerDict["filter"]["dnmnap"] = "##FILTER=<ID=dnmnap,Description=\"DNA Normal ALT percentage is less than the minimum\">\n"
    headerDict["filter"]["dtmnap"] = "##FILTER=<ID=dtmnap,Description=\"DNA Tumor ALT percentage is less than the minimum\">\n"
    headerDict["filter"]["rnmnap"] = "##FILTER=<ID=rnmnap,Description=\"RNA Normal ALT percentage is less than the minimum\">\n"
    headerDict["filter"]["rtmnap"] = "##FILTER=<ID=rtmnap,Description=\"RNA Tumor ALT percentage is less than the minimum\">\n"
    
    headerDict["filter"]["dnmnbq"] = "##FILTER=<ID=dnmnbq,Description=\"DNA Normal average ALT base quality is less than the minimum\">\n"
    headerDict["filter"]["dtmnbq"] = "##FILTER=<ID=dtmnbq,Description=\"DNA Tumor average ALT base quality is less than the minimum\">\n"
    headerDict["filter"]["rnmnbq"] = "##FILTER=<ID=rnmnbq,Description=\"RNA Normal average ALT base quality is less than the minimum\">\n"
    headerDict["filter"]["rtmnbq"] = "##FILTER=<ID=rtmnbq,Description=\"RNA Tumor average ALT base quality is less than the minimum\">\n"
    
    headerDict["filter"]["dnmnmqa"] = "##FILTER=<ID=dnmnmqa,Description=\"DNA Normal average ALT mapping quality is less than the minimum\">\n"
    headerDict["filter"]["dtmnmqa"] = "##FILTER=<ID=dtmnmqa,Description=\"DNA Tumor average ALT mapping quality is less than the minimum\">\n"
    headerDict["filter"]["rnmnmqa"] = "##FILTER=<ID=rnmnmqa,Description=\"RNA Normal average ALT mapping quality is less than the minimum\">\n"
    headerDict["filter"]["rtmnmqa"] = "##FILTER=<ID=rtmnmqa,Description=\"RNA Tumor average ALT mapping quality is less than the minimum\">\n"
    
    headerDict["filter"]["dnmnmq"] = "##FILTER=<ID=dnmnmq,Description=\"DNA Normal has no ALT reads with the minimum mapping quality\">\n"
    headerDict["filter"]["dtmnmq"] = "##FILTER=<ID=dtmnmq,Description=\"DNA Tumor has no ALT reads with the minimum mapping quality\">\n"
    headerDict["filter"]["rnmnmq"] = "##FILTER=<ID=rnmnmq,Description=\"RNA Normal has no ALT reads with the minimum mapping quality\">\n"
    headerDict["filter"]["rtmnmq"] = "##FILTER=<ID=rtmnmq,Description=\"RNA Tumor has no ALT reads with the minimum mapping quality\">\n"
    
    headerDict["filter"]["dnmxmq0"] = "##FILTER=<ID=dnmxmq0,Description=\"DNA Normal percentage of mapping quality zero reads for the ALT is above the maximum\">\n"
    headerDict["filter"]["dtmxmq0"] = "##FILTER=<ID=dtmxmq0,Description=\"DNA Tumor percentage of mapping quality zero reads for the ALT is above the maximum\">\n"
    headerDict["filter"]["rnmxmq0"] = "##FILTER=<ID=rnmxmq0,Description=\"RNA Normal percentage of mapping quality zero reads for the ALT is above the maximum\">\n"
    headerDict["filter"]["rtmxmq0"] = "##FILTER=<ID=rtmxmq0,Description=\"RNA Tumor percentage of mapping quality zero reads for the ALT is above the maximum\">\n"
    
    headerDict["filter"]["dnmnrb"] = "##FILTER=<ID=dnmnrb,Description=\"DNA Normal REF bases is less than the minimum\">\n"
    headerDict["filter"]["dnmnrp"] = "##FILTER=<ID=dnmnrp,Description=\"DNA Normal REF percentage is less than the minimum\">\n"
    headerDict["filter"]["dtmnrb"] = "##FILTER=<ID=dtmnrb,Description=\"DNA Tumor REF bases is less than the minimum\">\n"
    headerDict["filter"]["dtmnrp"] = "##FILTER=<ID=dtmnrp,Description=\"DNA Tumor REF percentage is less than the minimum\">\n"
    
    headerDict["filter"]["dnmxerr"] = "##FILTER=<ID=dnmxerr,Description=\"DNA Normal total ALT percentage attributed to error (sequencing, contamination, etc.) is greater than the maximum\">\n"
    headerDict["filter"]["dtmxerr"] = "##FILTER=<ID=dtmxerr,Description=\"DNA Tumor total ALT percentage attributed to error (sequencing, contamination, etc.) is greater than the maximum\">\n"
    headerDict["filter"]["rnmxerr"] = "##FILTER=<ID=rnmxerr,Description=\"RNA Normal total ALT percentage attributed to error (sequencing, contamination, etc.) is greater than the maximum\">\n"
    headerDict["filter"]["rtmxerr"] = "##FILTER=<ID=rtmxerr,Description=\"RNA Tumor total ALT percentage attributed to error (sequencing, contamination, etc.) is greater than the maximum\">\n"
    
    headerDict["filter"]["dnsbias"] = "##FILTER=<ID=dnsbias,Description=\"Strand bias, majority of reads supporting ALT are on forward OR reverse strand\">\n"
    headerDict["filter"]["dtsbias"] = "##FILTER=<ID=dtsbias,Description=\"Strand bias, majority of reads supporting ALT are on forward OR reverse strand\">\n"
    headerDict["filter"]["rnsbias"] = "##FILTER=<ID=rnsbias,Description=\"Strand bias, majority of reads supporting ALT are on forward OR reverse strand\">\n"
    headerDict["filter"]["rtsbias"] = "##FILTER=<ID=rtsbias,Description=\"Strand bias, majority of reads supporting ALT are on forward OR reverse strand\">\n"

    return headerDict


def get_vcf_header(aHeaderDict, aFilename, aCmdLineParams, aColumnsList, anIsDebug):
    
    # open the file
    vcfFileHandler = get_read_fileHandler(aFilename)
    
    for line in vcfFileHandler:
        
        # strip the carriage return and newline characters
        #line = line.rstrip("\r\n")

        '''
        if (anIsDebug):
            logging.debug("read vcfLine: %s", line)
        '''
        
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # get the header info
        if (line.startswith("##FORMAT")):
            formatId = get_id(line, anIsDebug)
            aHeaderDict["format"][formatId] = line
            
        elif (line.startswith("##INFO")):
            infoId = get_id(line, anIsDebug)
            aHeaderDict["info"][infoId] = line
            
        elif (line.startswith("##FILTER")):
            filterId = get_id(line, anIsDebug)
            aHeaderDict["filter"][filterId] = line
            
        elif (line.startswith("##SAMPLE")):
            sampleId = get_id(line, anIsDebug)
            aHeaderDict["sample"][sampleId] = line
            
        # if we find the vcfGenerator line, then update the params from the user
        elif ("vcfGenerator" in line):
            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")
            
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
                if (paramName.startswith("dnaNormal") and "DNA_NORMAL" not in aColumnsList):
                    continue;
                elif (paramName.startswith("rnaNormal") and "RNA_NORMAL" not in aColumnsList):
                    continue;
                elif (paramName.startswith("dnaTumor") and "DNA_TUMOR" not in aColumnsList):
                    continue;
                elif (paramName.startswith("rnaTumor") and "RNA_TUMOR" not in aColumnsList):
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
            aHeaderDict["metadata"]["vcfGenerator"] = generatorOutput
            
        elif (line.startswith("##")):
            metadataId = get_meta_id(line, anIsDebug)
            aHeaderDict["metadata"][metadataId] = line
            
        elif (line.startswith("#CHROM")):
            aHeaderDict["chrom"] = line
            
        if (not line.startswith("#")):
            break
            
    # close the file
    vcfFileHandler.close()

    return aHeaderDict


def output_header(aHeaderDict, aSortFlag, anOutputFileHandler):
    
    # output the header
    keys = aHeaderDict.keys()
    if aSortFlag:
        keys.sort(key=str)
    for key in keys:
        anOutputFileHandler.write(aHeaderDict[key])



def filter_by_mpileup_support(anId, aChrom, aVCFFilename, aHeaderFilename, anOutputFilename, aFilterUsingRNAFlag, anAddOriginFlag, aCmdLineParams, aDnaNormParamsDict, aDnaTumParamsDict, anRnaNormParamsDict, anRnaTumParamsDict, aGTMinDepth, aGTMinPct, aModMinDepth, aModMinPct, anLohMaxDepth, anLohMaxPct, anIsDebug):
    '''
    ' This function filters based on the mpileup read support.
    '
    ' anId:                    The Id for this sample
    ' aChrom:                  The chromosome being filtered
    ' aVCFFilename:            The filename to be filtered
    ' aHeaderFilename:         The filename with full header info
    ' anOutputFilename:        The output filename that will include the filters
    ' aFilterUsingRNAFlag:     If the calls should be filtered by the RNA as well
    ' anAddOriginFlag:         If the origin (DNA or RNA) of the call should be added to the INFO tag
    ' aCmdLineParams:          All the parameters specified by the user
    ' aDnaNormParamsDict:      The parameters for the normal DNA
    ' aDnaTumParamsDict:       The parameters for the tumor DNA
    ' anRnaNormParamsDict:     The parameters for the normal RNA
    ' anRnaTumParamsDict:      The parameters for the tumor RNA
    ' aGTMinDepth:             The minimum depth needed for the genotype
    ' aGTMinPct:               The minimum percent needed for the genotype
    ' aModMinDepth:            The minimum depth needed to make a call
    ' aModMinPct:              The minimum percent needed to make a call
    ' anIsDebug:               A flag for outputting debug messages to STDERR
    '''
    
    # initialize some variables
    totalEvents = 0
    includedEvents = 0
    somEventsPassing = 0
    somEventsWithTumorRna = 0
    somEventsWithTumorAltRna = 0
    
    # get the output file handler
    i_outputFileHandler = None
    if (anOutputFilename):
        i_outputFileHandler = get_write_fileHandler(anOutputFilename)
    
    # get the mpileup filter header lines
    headerDict = get_mpileup_header(anAddOriginFlag)
    
    # sometimes the header lines are stripped from the file, so get the necessary columnsList from the #CHROM line from a header file if it is specified
    if (aHeaderFilename != None):
        columnsList, headerDict = get_sample_columns(aHeaderFilename, headerDict, anIsDebug)
        headerDict = get_vcf_header(headerDict, aHeaderFilename, aCmdLineParams, columnsList, anIsDebug)
    else:
        columnsList, headerDict = get_sample_columns(aVCFFilename, headerDict, anIsDebug)
        headerDict = get_vcf_header(headerDict, aVCFFilename, aCmdLineParams, columnsList, anIsDebug)
    
    # output the header information
    output_header(headerDict["metadata"], False, i_outputFileHandler)
    output_header(headerDict["sample"], False, i_outputFileHandler)
    output_header(headerDict["filter"], True, i_outputFileHandler)
    output_header(headerDict["info"], True, i_outputFileHandler)
    output_header(headerDict["format"], True, i_outputFileHandler)
    i_outputFileHandler.write(headerDict["chrom"])
    
    # get the file
    i_vcfFileHandler = get_read_fileHandler(aVCFFilename)
        
    # for each event in the vcf file 
    for line in i_vcfFileHandler:
        # here are some examples of .vcf lines:
        # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31;DP=53;FA=0.04;INS=0;DEL=0;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP
        # GT:DP:AD:AF:INS:DEL:START:STOP:BQ:SB      0/0:2:2,0:1.0,0.0:0:0:0:0::36,0:0.0,0.0      0/0:1:1,0:1.0,0.0:0:0:0:0:39,0:1.0,0.0      0/1:50:48,2:0.96,0.04:0:0:2:0:32,18:0.75,0.5
        
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug):
            logging.debug("VCF Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # skip the header lines that are taken care of above
        if (line.startswith("#")):
            continue
        
        # now we are to the data
        # count the total events
        totalEvents += 1
        somEventWithTumorRna = False
        somEventWithTumorAltRna = False
        
        # split the line on the tab
        splitLine = line.split("\t")

        # sample VCF line
        # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31;DP=53;FA=0.04;INS=0;DEL=0;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP
        # GT:DP:AD:AF:INS:DEL:START:STOP:BQ:SB      0/0:2:2,0:1.0,0.0:0:0:0:0::36,0:0.0,0.0      0/0:1:1,0:1.0,0.0:0:0:0:0:39,0:1.0,0.0      0/1:50:48,2:0.96,0.04:0:0:2:0:32,18:0.75,0.5
        
        # the coordinate is the second element
        event_chr = splitLine[0]
        event_stopCoordinate = int(splitLine[1])
        event_idList = splitLine[2].split(";")
        event_refList = splitLine[3].split(",")
        event_altList = splitLine[4].split(",")
        event_score = splitLine[5]
        
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
        if (event_dnaNormalList == None or event_dnaNormalList[0] == "." or event_dnaNormalList[0] == "./."):
            haveDnaNormData = False
        # if there is no data, then set the flag
        if (event_rnaNormalList == None or event_rnaNormalList[0] == "." or event_rnaNormalList[0] == "./."):
            haveRnaNormData = False
        # if there is no data, then set the flag
        if (event_dnaTumorList == None or event_dnaTumorList[0] == "." or event_dnaTumorList[0] == "./."):
            haveDnaTumData = False
        # if there is no data, then set the flag
        if (event_rnaTumorList == None or event_rnaTumorList[0] == "." or event_rnaTumorList[0] == "./."):
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
        
        allFiltersSet = set()
        
        # get rid of bad mod types that don't meet the minimum requirements
        (event_infoDict, allFiltersSet) = pre_filter_mod_types(refPlusAltList, allFiltersSet, event_infoDict, map(int, event_dnaNormalDict["AD"]), map(int, event_rnaNormalDict["AD"]), map(int, event_dnaTumorDict["AD"]), map(int, event_rnaTumorDict["AD"]), aModMinDepth, aModMinPct, anLohMaxDepth, anLohMaxPct)
        if (anIsDebug):
            logging.debug("after pre_filter_mod_types(): modTypes=%s, modChanges=%s", list(event_infoDict["MT"]), list(event_infoDict["MC"]))
            
        # make copies of the lists to manipulate
        modTypesList = list(event_infoDict["MT"])
        modChangesList = list(event_infoDict["MC"])
        
        # keep track of filters for each mod to add to INFO
        modFilterTypes = []
        modFilters = []
        
        # for each modification type and change
        for (modType, modChange) in izip(event_infoDict["MT"], event_infoDict["MC"]):
            isValidMod = True
            filterSet = set()
            
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
                if (int(event_dnaNormalDict["BQ"][targetIndex]) < aDnaNormParamsDict["MinAltAvgBaseQual"]):
                    isValidMod = False
                    filterSet.add("dnmnbq")
                
                # check to make sure the normal DNA sample average mapping quality for ALT reads is above the min
                if (int(event_dnaNormalDict["MQA"][targetIndex]) < aDnaNormParamsDict["MinAltAvgMapQual"]):
                    isValidMod = False
                    filterSet.add("dnmnmqa")
                
                # check to make sure the normal DNA sample has at least 1 ALT read with a mapping quality above the min
                if (int(event_dnaNormalDict["MMQ"][targetIndex]) < aDnaNormParamsDict["MinAltMapQual"]):
                    isValidMod = False
                    filterSet.add("dnmnmq")
                
                # check to make sure the normal DNA sample has a maximum percentage of MQ0 reads supporting the ALT
                if (filterByMapQualZero(aDnaNormParamsDict, event_dnaNormalDict, targetIndex)):
                    isValidMod = False
                    filterSet.add("dnmxmq0")
                
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
                        if (int(event_rnaNormalDict["BQ"][targetIndex]) < anRnaNormParamsDict["MinAltAvgBaseQual"]):
                            isValidMod = False
                            filterSet.add("rnmnbq")
                        
                        # check to make sure the normal RNA sample average mapping quality for ALT reads is above the min
                        if (int(event_rnaNormalDict["MQA"][targetIndex]) < anRnaNormParamsDict["MinAltAvgMapQual"]):
                            isValidMod = False
                            filterSet.add("rnmnmqa")
                        
                        # check to make sure the normal RNA sample has at least 1 ALT read with a mapping quality above the min
                        if (int(event_rnaNormalDict["MMQ"][targetIndex]) < anRnaNormParamsDict["MinAltMapQual"]):
                            isValidMod = False
                            filterSet.add("rnmnmq")
                        
                        # check to make sure the normal RNA sample has a maximum percentage of MQ0 reads supporting the ALT
                        if (filterByMapQualZero(anRnaNormParamsDict, event_rnaNormalDict, targetIndex)):
                            isValidMod = False
                            filterSet.add("rnmxmq0")
                        
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
                        filterSet.add("dnacall")
                        
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
                    
                    # check to make sure the normal RNA sample number of ALT bases is above the min
                    if (int(event_rnaNormalDict["AD"][targetIndex]) < anRnaNormParamsDict["MinAltNumBases"]):
                        isValidMod = False
                        filterSet.add("rnmnab")
                    
                    # check to make sure the normal RNA sample percentage of ALT bases is above the min
                    if (float(event_rnaNormalDict["AF"][targetIndex]) < anRnaNormParamsDict["MinAltPct"]):
                        isValidMod = False
                        filterSet.add("rnmnap")
                    
                    # check to make sure the normal RNA sample average base quality for ALT bases is above the min
                    if (int(event_rnaNormalDict["BQ"][targetIndex]) < anRnaNormParamsDict["MinAltAvgBaseQual"]):
                        isValidMod = False
                        filterSet.add("rnmnbq")
                    
                    # check to make sure the normal RNA sample average map quality for ALT reads is above the min
                    if (int(event_rnaNormalDict["MQA"][targetIndex]) < anRnaNormParamsDict["MinAltAvgMapQual"]):
                        isValidMod = False
                        filterSet.add("rnmnmqa")
                    
                    # check to make sure the normal RNA sample has at least 1 ALT read with a mapping quality above the min
                    if (int(event_rnaNormalDict["MMQ"][targetIndex]) < anRnaNormParamsDict["MinAltMapQual"]):
                        isValidMod = False
                        filterSet.add("rnmnmq")
                    
                    # check to make sure the normal RNA sample has a maximum percentage of MQ0 reads supporting the ALT
                    if (filterByMapQualZero(anRnaNormParamsDict, event_rnaNormalDict, targetIndex)):
                        isValidMod = False
                        filterSet.add("rnmxmq0")
                        
                    # check to make sure the normal variant reads don't have a strand bias
                    if (filterByStrandBias(anRnaNormParamsDict, event_rnaNormalDict, sourceIndex, targetIndex)):
                        isValidMod = False
                        filterSet.add("rnsbias")
                    
                    # we want to make sure the normal RNA sample error percentage is below the max
                    # we want to make sure that the percentage of other ALTs in this sample is below the max error
                    if (filterByMaxError(refPlusAltList, anRnaNormParamsDict, event_rnaNormalDict, sourceIndex, targetIndex, False, anIsDebug)):
                        isValidMod = False
                        filterSet.add("rnmxerr")
                # we are filtering via the DNA, so put in a dummy filter so that they don't pass
                else:
                    isValidMod = False
                    filterSet.add("rnacall")
                    
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
                if (int(event_dnaTumorDict["BQ"][targetIndex]) < aDnaTumParamsDict["MinAltAvgBaseQual"]):
                    isValidMod = False
                    filterSet.add("dtmnbq")
                
                # check to make sure the tumor DNA sample average map quality for ALT reads is above the min
                if (int(event_dnaTumorDict["MQA"][targetIndex]) < aDnaTumParamsDict["MinAltAvgMapQual"]):
                    isValidMod = False
                    filterSet.add("dtmnmqa")
                    
                # check to make sure the tumor DNA sample has at least 1 ALT read with a mapping quality above the min
                if (int(event_dnaTumorDict["MMQ"][targetIndex]) < aDnaTumParamsDict["MinAltMapQual"]):
                    isValidMod = False
                    filterSet.add("dtmnmq")
                
                # check to make sure the tumor DNA sample has a maximum percentage of MQ0 reads supporting the ALT
                if (filterByMapQualZero(aDnaTumParamsDict, event_dnaTumorDict, targetIndex)):
                    isValidMod = False
                    filterSet.add("dtmxmq0")
                
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
                        if (int(event_rnaTumorDict["BQ"][targetIndex]) < anRnaTumParamsDict["MinAltAvgBaseQual"]):
                            isValidMod = False
                            filterSet.add("rtmnbq")
                        
                        # check to make sure the tumor RNA sample average map quality for ALT reads is above the min
                        if (int(event_rnaTumorDict["MQA"][targetIndex]) < anRnaTumParamsDict["MinAltAvgMapQual"]):
                            isValidMod = False
                            filterSet.add("rtmnmqa")
                        
                        # check to make sure the tumor RNA sample has at least 1 ALT read with a mapping quality above the min
                        if (int(event_rnaTumorDict["MMQ"][targetIndex]) < anRnaTumParamsDict["MinAltMapQual"]):
                            isValidMod = False
                            filterSet.add("rtmnmq")
                        
                        # check to make sure the tumor RNA sample has a maximum percentage of MQ0 reads supporting the ALT
                        if (filterByMapQualZero(anRnaTumParamsDict, event_rnaTumorDict, targetIndex)):
                            isValidMod = False
                            filterSet.add("rtmxmq0")
                        
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
                    # if this is an RNA tumor variant, then don't filter on the DNA
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
                    
                    # check to make sure the tumor RNA sample percentage of ALT bases is above the min
                    if (float(event_rnaTumorDict["AF"][targetIndex]) < anRnaTumParamsDict["MinAltPct"]):
                        isValidMod = False
                        filterSet.add("rtmnap")
                    
                    # check to make sure the tumor RNA sample average base quality for ALT bases is above the min
                    if (int(event_rnaTumorDict["BQ"][targetIndex]) < anRnaTumParamsDict["MinAltAvgBaseQual"]):
                        isValidMod = False
                        filterSet.add("rtmnbq")
                    
                    # check to make sure the tumor RNA sample average map quality for ALT reads is above the min
                    if (int(event_rnaTumorDict["MQA"][targetIndex]) < anRnaTumParamsDict["MinAltAvgMapQual"]):
                        isValidMod = False
                        filterSet.add("rtmnmqa")
                        
                    # check to make sure the tumor RNA sample has at least 1 ALT read with a mapping quality above the min
                    if (int(event_rnaTumorDict["MMQ"][targetIndex]) < anRnaTumParamsDict["MinAltMapQual"]):
                        isValidMod = False
                        filterSet.add("rtmnmq")
                        
                    # check to make sure the tumor RNA sample has a maximum percentage of MQ0 reads supporting the ALT
                    if (filterByMapQualZero(anRnaTumParamsDict, event_rnaTumorDict, targetIndex)):
                        isValidMod = False
                        filterSet.add("rtmxmq0")
                    
                    # check to make sure the tumor variant reads don't have a strand bias
                    if (filterByStrandBias(anRnaTumParamsDict, event_rnaTumorDict, sourceIndex, targetIndex)):
                        isValidMod = False
                        filterSet.add("rtsbias")
                    
                    # we want to make sure the tumor RNA sample error percentage is below the max
                    # we want to make sure that the percentage of other ALTs in this sample is below the max error
                    if (filterByMaxError(refPlusAltList, anRnaTumParamsDict, event_rnaTumorDict, sourceIndex, targetIndex, False, anIsDebug)):
                        isValidMod = False
                        filterSet.add("rtmxerr")
                # we are filtering via the DNA, so put in a dummy filter so that they don't pass
                else:
                    isValidMod = False
                    filterSet.add("rnacall")
                
            if (anIsDebug):
                logging.debug("modType=%s, modChange=%s, isValidMod=%s, filters=%s", modType, modChange, isValidMod, filterSet)
                
            # if this one is not valid and we have more,
            # then set the filters for this one, remove it,
            # and try the next one
            if (not isValidMod):
                allFiltersSet = allFiltersSet.union(filterSet)
                
                # find the origin
                origin = "DNA"
                if (aFilterUsingRNAFlag):
                    origin = "RNA"

                modFilterTypes.append("_".join([origin, modType, modChange]))
                modFilters.append("_".join(filterSet))
                
                # remove it and try the next one
                modIndices = range(0, len(modTypesList))
                for (removeModType, removeModChange, modIndex) in izip(modTypesList, modChangesList, modIndices):
                    if (modType == removeModType and modChange == removeModChange):
                        del modTypesList[modIndex]
                        del modChangesList[modIndex]
                        break;
                 
        # after looping through all of them:  if there are still some valid mod types, then set them in the infoDict and ignore the other filtered calls
        if (len(modTypesList) > 0):
            event_infoDict["MT"] = modTypesList
            event_infoDict["MC"] = modChangesList
            
            # if an event passed, get the final mod type
            event_infoDict = get_final_mod_type(event_infoDict, anIsDebug)
        
        # otherwise add the appropriate filters
        else:
            event_infoDict["MFT"] = modFilterTypes
            event_infoDict["MF"] = modFilters
            event_filterSet = event_filterSet.union(allFiltersSet)
        
        # create the output list
        vcfOutputList = [event_chr, str(event_stopCoordinate)]
        
        # add the ref, alt, and score
        vcfOutputList.append(";".join(event_idList))
        vcfOutputList.append(",".join(event_refList))
        vcfOutputList.append(",".join(event_altList))
        vcfOutputList.append(event_score)
        
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
                
    logging.info("Chrom %s and Id %s: %s events passed out of %s total events", aChrom, anId, includedEvents, totalEvents)
    logging.info("\t".join([anId, aChrom, str(somEventsPassing), str(somEventsWithTumorRna), str(somEventsWithTumorAltRna)]))
    
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
    
    i_cmdLineParser.add_option("", "--dnaNormalMinTotalBases", type="int", default=int(10), dest="dnaNormalMinTotalNumBases", metavar="DNA_NOR_MIN_TOTAL_BASES", help="the minimum number of overall normal DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMaxTotalBases", type="int", default=int(20000), dest="dnaNormalMaxTotalNumBases", metavar="DNA_NOR_MAX_TOTAL_BASES", help="the maximum number of overall normal DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinAltBases", type="int", default=int(4), dest="dnaNormalMinAltNumBases", metavar="DNA_NOR_MIN_ALT_BASES", help="the minimum number of alternative normal DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinAltPct", type="float", default=float(0.10), dest="dnaNormalMinAltPct", metavar="DNA_NOR_MIN_ALT_PCT", help="the minimum percentage of alternative normal DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMaxErrPct", type="float", default=float(0.01), dest="dnaNormalMaxErrPct", metavar="DNA_NOR_MAX_ERR_PCT", help="the maximum percentage of alternative normal DNA reads allowed that support a variant, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinErrPctDepth", type="int", default=float(2), dest="dnaNormalMinErrPctDepth", metavar="DNA_NOR_MIN_ERR_PCT_DEPTH", help="the minimum error count depth needed for the max error percent filter to be applied, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMaxStrandBias", type="float", default=float(0.99), dest="dnaNormalMaxStrandBias", metavar="DNA_NOR_MAX_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinStrandBiasDepth", type="int", default=float(4), dest="dnaNormalMinStrandBiasDepth", metavar="DNA_NOR_MIN_STRAND_BIAS_DP", help="the minimum total depth needed for the strand bias filter to be applied, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinAltAvgBaseQual", type="int", default=int(20), dest="dnaNormalMinAltAvgBaseQual", metavar="DNA_NOR_MIN_ALT_AVG_BQ", help="the minimum average base quality for the ALT reads, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinAltAvgMapQual", type="int", default=int(15), dest="dnaNormalMinAltAvgMapQual", metavar="DNA_NOR_MIN_ALT_AVG_MQ", help="the minimum average mapping quality for the ALT reads, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinAltMapQual", type="int", default=int(20), dest="dnaNormalMinAltMapQual", metavar="DNA_NOR_MIN_ALT_MQ", help="at least 1 ALT read needs this minimum mapping quality, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMaxAltMapQualZeroPct", type="float", default=float(0.50), dest="dnaNormalMaxAltMapQualZeroPct", metavar="DNA_NOR_MAX_ALT_MQ0_PCT", help="the maximum percentage of mapping quality zero reads for the ALT reads, %default by default")
    
    i_cmdLineParser.add_option("", "--dnaTumorMinTotalBases", type="int", default=int(10), dest="dnaTumorMinTotalNumBases", metavar="DNA_TUM_MIN_TOTAL_BASES", help="the minimum number of overall tumor DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMaxTotalBases", type="int", default=int(20000), dest="dnaTumorMaxTotalNumBases", metavar="DNA_TUM_MAX_TOTAL_BASES", help="the maximum number of overall tumor DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinAltBases", type="int", default=int(4), dest="dnaTumorMinAltNumBases", metavar="DNA_TUM_MIN_ALT_BASES", help="the minimum number of alternative tumor DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinAltPct", type="float", default=float(0.10), dest="dnaTumorMinAltPct", metavar="DNA_TUM_MIN_ALT_PCT", help="the minimum percentage of alternative tumor DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMaxErrPct", type="float", default=float(0.01), dest="dnaTumorMaxErrPct", metavar="DNA_TUM_MAX_ERR_PCT", help="the maximum percentage of alternative tumor DNA reads allowed that support a variant in the tumor RNA, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinErrPctDepth", type="int", default=float(2), dest="dnaTumorMinErrPctDepth", metavar="DNA_TUM_MIN_ERR_PCT_DEPTH", help="the minimum error count depth needed for the max error percent filter to be applied, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMaxStrandBias", type="float", default=float(0.99), dest="dnaTumorMaxStrandBias", metavar="DNA_TUM_MAX_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinStrandBiasDepth", type="int", default=float(4), dest="dnaTumorMinStrandBiasDepth", metavar="DNA_TUM_MIN_STRAND_BIAS_DP", help="the minimum total depth needed for the strand bias filter to be applied, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinAltAvgBaseQual", type="int", default=int(20), dest="dnaTumorMinAltAvgBaseQual", metavar="DNA_TUM_MIN_ALT_AVG_BQ", help="the minimum average base quality for the ALT reads, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinAltAvgMapQual", type="int", default=int(15), dest="dnaTumorMinAltAvgMapQual", metavar="DNA_TUM_MIN_ALT_AVG_MQ", help="the minimum average mapping quality for the ALT reads, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinAltMapQual", type="int", default=int(20), dest="dnaTumorMinAltMapQual", metavar="DNA_TUM_MIN_ALT_MQ", help="at least 1 ALT read needs this minimum mapping quality, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMaxAltMapQualZeroPct", type="float", default=float(0.50), dest="dnaTumorMaxAltMapQualZeroPct", metavar="DNA_TUM_MAX_ALT_MQ0_PCT", help="the maximum percentage of mapping quality zero reads for the ALT reads, %default by default")
    
    i_cmdLineParser.add_option("", "--rnaNormalMinTotalBases", type="int", default=int(10), dest="rnaNormalMinTotalNumBases", metavar="RNA_NOR_MIN_TOTAL_BASES", help="the minimum number of overall normal RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMaxTotalBases", type="int", default=int(20000), dest="rnaNormalMaxTotalNumBases", metavar="RNA_NOR_MAX_TOTAL_BASES", help="the maximum number of overall normal RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinAltBases", type="int", default=int(4), dest="rnaNormalMinAltNumBases", metavar="RNA_NOR_MIN_ALT_BASES", help="the minimum number of alternative normal RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinAltPct", type="float", default=float(0.10), dest="rnaNormalMinAltPct", metavar="RNA_NOR_MIN_ALT_PCT", help="the minimum percentage of alternative normal RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMaxErrPct", type="float", default=float(0.01), dest="rnaNormalMaxErrPct", metavar="RNA_NOR_MAX_ERR_PCT", help="the maximum percentage of alternative normal RNA reads allowed that support a variant in the tumor, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinErrPctDepth", type="int", default=float(2), dest="rnaNormalMinErrPctDepth", metavar="RNA_NOR_MIN_ERR_PCT_DEPTH", help="the minimum error count depth needed for the max error percent filter to be applied, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMaxStrandBias", type="float", default=float(0.99), dest="rnaNormalMaxStrandBias", metavar="RNA_NOR_MAX_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinStrandBiasDepth", type="int", default=float(4), dest="rnaNormalMinStrandBiasDepth", metavar="RNA_NOR_MIN_STRAND_BIAS_DP", help="the minimum total depth needed for the strand bias filter to be applied, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinAltAvgBaseQual", type="int", default=int(20), dest="rnaNormalMinAltAvgBaseQual", metavar="RNA_NOR_MIN_ALT_AVG_BQ", help="the minimum average base quality for the ALT reads, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinAltAvgMapQual", type="int", default=int(15), dest="rnaNormalMinAltAvgMapQual", metavar="RNA_NOR_MIN_ALT_AVG_MQ", help="the minimum average mapping quality for the ALT reads, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinAltMapQual", type="int", default=int(20), dest="rnaNormalMinAltMapQual", metavar="RNA_NOR_MIN_ALT_MQ", help="at least 1 ALT read needs this minimum mapping quality, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMaxAltMapQualZeroPct", type="float", default=float(0.50), dest="rnaNormalMaxAltMapQualZeroPct", metavar="RNA_NOR_MAX_ALT_MQ0_PCT", help="the maximum percentage of mapping quality zero reads for the ALT reads, %default by default")
    
    i_cmdLineParser.add_option("", "--rnaTumorMinTotalBases", type="int", default=int(10), dest="rnaTumorMinTotalNumBases", metavar="RNA_TUM_MIN_TOTAL_BASES", help="the minimum number of overall tumor RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMaxTotalBases", type="int", default=int(20000), dest="rnaTumorMaxTotalNumBases", metavar="RNA_TUM_MAX_TOTAL_BASES", help="the maximum number of overall tumor RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinAltBases", type="int", default=int(4), dest="rnaTumorMinAltNumBases", metavar="RNA_TUM_MIN_ALT_BASES", help="the minimum number of alternative tumor RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinAltPct", type="float", default=float(0.10), dest="rnaTumorMinAltPct", metavar="RNA_TUM_MIN_ALT_PCT", help="the minimum percentage of alternative tumor RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMaxErrPct", type="float", default=float(0.01), dest="rnaTumorMaxErrPct", metavar="RNA_TUM_MAX_ERR_PCT", help="the maximum percentage of alternative tumor RNA reads allowed that support a variant in a future data type, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinErrPctDepth", type="int", default=float(2), dest="rnaTumorMinErrPctDepth", metavar="RNA_TUM_MIN_ERR_PCT_DEPTH", help="the minimum error count depth needed for the max error percent filter to be applied, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMaxStrandBias", type="float", default=float(0.99), dest="rnaTumorMaxStrandBias", metavar="RNA_TUM_MAX_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinStrandBiasDepth", type="int", default=float(4), dest="rnaTumorMinStrandBiasDepth", metavar="RNA_TUM_MIN_STRAND_BIAS_DP", help="the minimum total depth needed for the strand bias filter to be applied, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinAltAvgBaseQual", type="int", default=int(20), dest="rnaTumorMinAltAvgBaseQual", metavar="RNA_TUM_MIN_ALT_AVG_BQ", help="the minimum average base quality for the ALT reads, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinAltAvgMapQual", type="int", default=int(15), dest="rnaTumorMinAltAvgMapQual", metavar="RNA_TUM_MIN_ALT_AVG_MQ", help="the minimum average mapping quality for the ALT reads, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinAltMapQual", type="int", default=int(20), dest="rnaTumorMinAltMapQual", metavar="RNA_TUM_MIN_ALT_MQ", help="at least 1 ALT read needs this minimum mapping quality, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMaxAltMapQualZeroPct", type="float", default=float(0.50), dest="rnaTumorMaxAltMapQualZeroPct", metavar="RNA_TUM_MAX_ALT_MQ0_PCT", help="the maximum percentage of mapping quality zero reads for the ALT reads, %default by default")
    
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
    
    i_dnaNormParams = {}
    i_dnaNormParams["MinTotalNumBases"] = i_cmdLineOptions.dnaNormalMinTotalNumBases
    i_dnaNormParams["MaxTotalNumBases"] = i_cmdLineOptions.dnaNormalMaxTotalNumBases
    i_dnaNormParams["MinAltNumBases"] = i_cmdLineOptions.dnaNormalMinAltNumBases
    i_dnaNormParams["MinAltPct"] = i_cmdLineOptions.dnaNormalMinAltPct
    i_dnaNormParams["MaxErrPct"] = i_cmdLineOptions.dnaNormalMaxErrPct
    i_dnaNormParams["MinErrPctDP"] = i_cmdLineOptions.dnaNormalMinErrPctDepth
    i_dnaNormParams["MaxStrandBias"] = i_cmdLineOptions.dnaNormalMaxStrandBias
    i_dnaNormParams["MinStrBiasDP"] = i_cmdLineOptions.dnaNormalMinStrandBiasDepth
    i_dnaNormParams["MinAltAvgBaseQual"] = i_cmdLineOptions.dnaNormalMinAltAvgBaseQual
    i_dnaNormParams["MinAltAvgMapQual"] = i_cmdLineOptions.dnaNormalMinAltAvgMapQual
    i_dnaNormParams["MinAltMapQual"] = i_cmdLineOptions.dnaNormalMinAltMapQual
    i_dnaNormParams["MaxAltMapQualZeroPct"] = i_cmdLineOptions.dnaNormalMaxAltMapQualZeroPct
    
    i_dnaTumParams = {}
    i_dnaTumParams["MinTotalNumBases"] = i_cmdLineOptions.dnaTumorMinTotalNumBases
    i_dnaTumParams["MaxTotalNumBases"] = i_cmdLineOptions.dnaTumorMaxTotalNumBases
    i_dnaTumParams["MinAltNumBases"] = i_cmdLineOptions.dnaTumorMinAltNumBases
    i_dnaTumParams["MinAltPct"] = i_cmdLineOptions.dnaTumorMinAltPct
    i_dnaTumParams["MaxErrPct"] = i_cmdLineOptions.dnaTumorMaxErrPct
    i_dnaTumParams["MinErrPctDP"] = i_cmdLineOptions.dnaTumorMinErrPctDepth
    i_dnaTumParams["MaxStrandBias"] = i_cmdLineOptions.dnaTumorMaxStrandBias
    i_dnaTumParams["MinStrBiasDP"] = i_cmdLineOptions.dnaTumorMinStrandBiasDepth
    i_dnaTumParams["MinAltAvgBaseQual"] = i_cmdLineOptions.dnaTumorMinAltAvgBaseQual
    i_dnaTumParams["MinAltAvgMapQual"] = i_cmdLineOptions.dnaTumorMinAltAvgMapQual
    i_dnaTumParams["MinAltMapQual"] = i_cmdLineOptions.dnaTumorMinAltMapQual
    i_dnaTumParams["MaxAltMapQualZeroPct"] = i_cmdLineOptions.dnaTumorMaxAltMapQualZeroPct
    
    i_rnaNormParams = {}
    i_rnaNormParams["MinTotalNumBases"] = i_cmdLineOptions.rnaNormalMinTotalNumBases
    i_rnaNormParams["MaxTotalNumBases"] = i_cmdLineOptions.rnaNormalMaxTotalNumBases
    i_rnaNormParams["MinAltNumBases"] = i_cmdLineOptions.rnaNormalMinAltNumBases
    i_rnaNormParams["MinAltPct"] = i_cmdLineOptions.rnaNormalMinAltPct
    i_rnaNormParams["MaxErrPct"] = i_cmdLineOptions.rnaNormalMaxErrPct
    i_rnaNormParams["MinErrPctDP"] = i_cmdLineOptions.rnaNormalMinErrPctDepth
    i_rnaNormParams["MaxStrandBias"] = i_cmdLineOptions.rnaNormalMaxStrandBias
    i_rnaNormParams["MinStrBiasDP"] = i_cmdLineOptions.rnaNormalMinStrandBiasDepth
    i_rnaNormParams["MinAltAvgBaseQual"] = i_cmdLineOptions.rnaNormalMinAltAvgBaseQual
    i_rnaNormParams["MinAltAvgMapQual"] = i_cmdLineOptions.rnaNormalMinAltAvgMapQual
    i_rnaNormParams["MinAltMapQual"] = i_cmdLineOptions.rnaNormalMinAltMapQual
    i_rnaNormParams["MaxAltMapQualZeroPct"] = i_cmdLineOptions.rnaNormalMaxAltMapQualZeroPct
    
    i_rnaTumParams = {}
    i_rnaTumParams["MinTotalNumBases"] = i_cmdLineOptions.rnaTumorMinTotalNumBases
    i_rnaTumParams["MaxTotalNumBases"] = i_cmdLineOptions.rnaTumorMaxTotalNumBases
    i_rnaTumParams["MinAltNumBases"] = i_cmdLineOptions.rnaTumorMinAltNumBases
    i_rnaTumParams["MinAltPct"] = i_cmdLineOptions.rnaTumorMinAltPct
    i_rnaTumParams["MaxErrPct"] = i_cmdLineOptions.rnaTumorMaxErrPct
    i_rnaTumParams["MinErrPctDP"] = i_cmdLineOptions.rnaTumorMinErrPctDepth
    i_rnaTumParams["MaxStrandBias"] = i_cmdLineOptions.rnaTumorMaxStrandBias
    i_rnaTumParams["MinStrBiasDP"] = i_cmdLineOptions.rnaTumorMinStrandBiasDepth
    i_rnaTumParams["MinAltAvgBaseQual"] = i_cmdLineOptions.rnaTumorMinAltAvgBaseQual
    i_rnaTumParams["MinAltAvgMapQual"] = i_cmdLineOptions.rnaTumorMinAltAvgMapQual
    i_rnaTumParams["MinAltMapQual"] = i_cmdLineOptions.rnaTumorMinAltMapQual
    i_rnaTumParams["MaxAltMapQualZeroPct"] = i_cmdLineOptions.rnaTumorMaxAltMapQualZeroPct
    
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
    i_debug = (i_numericLogLevel == logging.DEBUG)
        
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
        
        logging.debug("dna normal minTotalBases: %s" % i_dnaNormParams["MinTotalNumBases"])
        logging.debug("dna normal maxTotalBases: %s" % i_dnaNormParams["MaxTotalNumBases"])
        logging.debug("dna normal minAltBases: %s" % i_dnaNormParams["MinAltNumBases"])
        logging.debug("dna normal minAltPct: %s" % i_dnaNormParams["MinAltPct"])
        logging.debug("dna normal maxErrPct: %s" % i_dnaNormParams["MaxErrPct"])
        logging.debug("dna normal minErrPctDP: %s" % i_dnaNormParams["MinErrPctDP"])
        logging.debug("dna normal strandBias: %s" % i_dnaNormParams["MaxStrandBias"])
        logging.debug("dna normal strandBias minDP: %s" % i_dnaNormParams["MinStrBiasDP"])
        logging.debug("dna normal minAltAvgBaseQual: %s" % i_dnaNormParams["MinAltAvgBaseQual"])
        logging.debug("dna normal minAltAvgMapQual: %s" % i_dnaNormParams["MinAltAvgMapQual"])
        logging.debug("dna normal minAltMapQual: %s" % i_dnaNormParams["MinAltMapQual"])
        logging.debug("dna normal maxAltMapQualZeroPct: %s" % i_dnaNormParams["MaxAltMapQualZeroPct"])
        
        logging.debug("dna tumor minTotalBases: %s" % i_dnaTumParams["MinTotalNumBases"])
        logging.debug("dna tumor maxTotalBases: %s" % i_dnaTumParams["MaxTotalNumBases"])
        logging.debug("dna tumor minAltBases: %s" % i_dnaTumParams["MinAltNumBases"])
        logging.debug("dna tumor minAltPct: %s" % i_dnaTumParams["MinAltPct"])
        logging.debug("dna tumor maxErrPct: %s" % i_dnaTumParams["MaxErrPct"])
        logging.debug("dna tumor minErrPctDP: %s" % i_dnaTumParams["MinErrPctDP"])
        logging.debug("dna tumor strandBias: %s" % i_dnaTumParams["MaxStrandBias"])
        logging.debug("dna tumor strandBias minDP: %s" % i_dnaTumParams["MinStrBiasDP"])
        logging.debug("dna tumor minAltAvgBaseQual: %s" % i_dnaTumParams["MinAltAvgBaseQual"])
        logging.debug("dna tumor minAltAvgMapQual: %s" % i_dnaTumParams["MinAltAvgMapQual"])
        logging.debug("dna tumor minAltMapQual: %s" % i_dnaTumParams["MinAltMapQual"])
        logging.debug("dna tumor maxAltMapQualZeroPct: %s" % i_dnaTumParams["MaxAltMapQualZeroPct"])
        
        logging.debug("rna normal minTotalBases: %s" % i_rnaNormParams["MinTotalNumBases"])
        logging.debug("rna normal maxTotalBases: %s" % i_rnaNormParams["MaxTotalNumBases"])
        logging.debug("rna normal minAltBases: %s" % i_rnaNormParams["MinAltNumBases"])
        logging.debug("rna normal minAltPct: %s" % i_rnaNormParams["MinAltPct"])
        logging.debug("rna normal maxErrPct: %s" % i_rnaNormParams["MaxErrPct"])
        logging.debug("rna normal minErrPctDP: %s" % i_rnaNormParams["MinErrPctDP"])
        logging.debug("rna normal strandBias: %s" % i_rnaNormParams["MaxStrandBias"])
        logging.debug("rna normal strandBias minDP: %s" % i_rnaNormParams["MinStrBiasDP"])
        logging.debug("rna normal minAltAvgBaseQual: %s" % i_rnaNormParams["MinAltAvgBaseQual"])
        logging.debug("rna normal minAltAvgMapQual: %s" % i_rnaNormParams["MinAltAvgMapQual"])
        logging.debug("rna normal minAltMapQual: %s" % i_rnaNormParams["MinAltMapQual"])
        logging.debug("rna normal maxAltMapQualZeroPct: %s" % i_rnaNormParams["MaxAltMapQualZeroPct"])
        
        logging.debug("rna tumor minTotalBases: %s" % i_rnaTumParams["MinTotalNumBases"])
        logging.debug("rna tumor maxTotalBases: %s" % i_rnaTumParams["MaxTotalNumBases"])
        logging.debug("rna tumor minAltBases: %s" % i_rnaTumParams["MinAltNumBases"])
        logging.debug("rna tumor minAltPct: %s" % i_rnaTumParams["MinAltPct"])
        logging.debug("rna tumor maxErrPct: %s" % i_rnaTumParams["MaxErrPct"])
        logging.debug("rna tumor minErrPctDP: %s" % i_rnaTumParams["MinErrPctDP"])
        logging.debug("rna tumor strandBias: %s" % i_rnaTumParams["MaxStrandBias"])
        logging.debug("rna tumor strandBias minDP: %s" % i_rnaTumParams["MinStrBiasDP"])
        logging.debug("rna tumor minAltAvgBaseQual: %s" % i_rnaTumParams["MinAltAvgBaseQual"])
        logging.debug("rna tumor minAltAvgMapQual: %s" % i_rnaTumParams["MinAltAvgMapQual"])
        logging.debug("rna tumor minAltMapQual: %s" % i_rnaTumParams["MinAltMapQual"])
        logging.debug("rna tumor maxAltMapQualZeroPct: %s" % i_rnaTumParams["MaxAltMapQualZeroPct"])
                    
    # check for any errors
    if (not radiaUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)
    
    filter_by_mpileup_support(i_id, i_chrom, i_vcfFilename, i_headerFilename, i_outputFilename, i_filterUsingRNA, i_addOrigin, i_cmdLineOptionsDict, i_dnaNormParams, i_dnaTumParams, i_rnaNormParams, i_rnaTumParams, i_genotypeMinDepth, i_genotypeMinPct, i_modMinDepth, i_modMinPct, i_lohMaxDepth, i_lohMaxPct, i_debug)
    return

main()
sys.exit(0)
