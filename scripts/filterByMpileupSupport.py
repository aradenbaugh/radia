#!/usr/bin/env python

from optparse import OptionParser   # used for parsing command line arguments
import radiaUtil                    # utility functions for rna editing
import sys                          # system module
import collections
import logging
from itertools import izip
import re
from math import floor

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

# this regular expression is used to extract the
# ID tag from the INFO, FORMAT, and FILTER fields
i_headerIDRegEx = re.compile("ID=(\\w)*,")
# this regular expression is used to extract the
# Type tag from the INFO and FORMAT fields
i_headerTypeRegEx = re.compile("Type=(\\w)*,")
# this regular expression is used to extract the
# Derived tag from the PEDIGREE fields
i_headerDerivedRegEx = re.compile("Derived=(\\w)*,")


def fix_genotypes(aChrom, aRefList, anAltList,
                  anAlleleDepthsList, aParamsDict):
    '''
    ' This method assigns the genotype.
    '
    ' aChrom:             The chrom
    ' aRefList:           The list of ref alleles
    ' anAltList:          The list of alt alleles
    ' anAlleleDepthsList: The list of depths for all alleles
    ' aParamsDict:        Contains the minGTdp and minGTpct
    '''
    singleGenotypeChroms = ["chrY", "Y"]
    mChroms = ["chrM", "chrMT", "M", "MT"]
    refAltList = aRefList + anAltList
    gtMinDepth = aParamsDict["MinGenotypeDepth"]
    gtMinPct = aParamsDict["MinGenotypePct"]

    # if it is a single chrom, then we can only
    # assign one allele for the genotype
    # if one of the alts has a depth and percent above
    # the mins, then use it, otherwise use the ref
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
            (maxAltBase, maxAltDepth) = max(altCountsDict.iteritems(),
                                            key=lambda x:x[1])
            maxAltPct = round(maxAltDepth/float(totalDepth), 2)

            # if the max alt depth is large enough
            if (maxAltDepth >= gtMinDepth and maxAltPct >= gtMinPct):
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

    # if it is an M chrom, then we can assign
    # as many alleles as we want for the genotype
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
                if (depth >= gtMinDepth and percent >= gtMinPct):
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

    # if it is a diploid chrom, then assign the
    # 2 max counts above the min cutoffs
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
            # have to do the rounding here, otherwise the ones that are
            # just below 10% that get rounded up won't have the 0/1 genotype
            max2Pct = round(max2Depth/float(totalDepth), 2)

            # if the max depth is large enough
            if (max2Depth >= gtMinDepth and max2Pct >= gtMinPct):
                # if the two maxes are the same depth,
                # then return the second index
                if (max1Depth == max2Depth):
                    # get all indices for this depth
                    indices = [i for i,
                               x in enumerate(anAlleleDepthsList)
                               if x == max2Depth]
                    # take the second index, since we
                    # already took the first index above
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


def pre_filter_mod_types(aRefPlusAltList, anAllFiltersSet, anInfoDict,
                         aDNANormalDepths, anRNANormalDepths, aDNATumorDepths,
                         anRNATumorDepths, aParamsDict, anIsDebug):

    aModTypeList = anInfoDict["MT"]
    aModChangeList = anInfoDict["MC"]

    # make copies of the lists to manipulate
    modTypesList = list(aModTypeList)
    modChangesList = list(aModChangeList)

    modMinDepth = aParamsDict["MinModDepth"]
    modMinPct = aParamsDict["MinModPct"]
    lohMaxDepth = aParamsDict["MaxLohDepth"]
    lohMaxPct = aParamsDict["MaxLohPct"]

    try:
        # for each modification type and change
        for (modType, modChange) in izip(aModTypeList, aModChangeList):
            # get the source and target alleles
            (source, target) = modChange.split(">")

            # if (anIsDebug):
            #    logging.debug("modType=%s, modChange=%s", modType, modChange)

            # for every germline
            if (modType == "GERM"):
                # get the indices
                sIndex = aRefPlusAltList.index(source)
                tIndex = aRefPlusAltList.index(target)

                # get the DNA normal total depth
                dnDepth = sum(aDNANormalDepths)

                # if this is labeled as Germline, but the dna normal variant
                # depth is not sufficient, maybe it is a somatic one
                dnTargetDepth = aDNANormalDepths[tIndex]
                if (dnDepth > 0):
                    dnTargetPct = round(dnTargetDepth/float(dnDepth), 2)
                else:
                    dnTargetPct = 0.0

                if (dnTargetDepth < modMinDepth or
                    dnTargetPct < modMinPct):
                    modTypesList.remove(modType)
                    modChangesList.remove(modChange)

                    # add the filters
                    if (dnTargetDepth < modMinDepth):
                        anAllFiltersSet.add("dnmnad")
                    if (dnTargetPct < modMinPct):
                        anAllFiltersSet.add("dnmnaf")

                    # get the tumor depth
                    dtDepth = sum(aDNATumorDepths)

                    # this is a hack for ones that were mis-classified
                    # as germline by radia, but should've been
                    # classified as somatic
                    if (dtDepth > 0):
                        dtTargetDepth = aDNATumorDepths[tIndex]
                        dtTargetPct = round(dtTargetDepth/float(dtDepth), 2)

                        # if the tumor depth is above the minimum,
                        # then add it (these were mis-classified by
                        # the original radia.py script)
                        if (dtTargetDepth >= modMinDepth or
                            dtTargetPct >= modMinPct):
                            '''
                            if (anIsDebug):
                                logging.debug("GERM call with not enough " +
                                              "normal DNA " +
                                              "(altDepth=%s < minDepth=%s " +
                                              "or altPct=%s < minPct=%s) is " +
                                              "being tested as SOM",
                                              dnTargetDepth, modMinDepth,
                                              dnTargetPct, modMinPct)
                            '''
                            modTypesList.append("SOM")
                            modChangesList.append(modChange)

            # for every somatic call
            elif (modType == "SOM"):
                # get the indices
                sIndex = aRefPlusAltList.index(source)
                tIndex = aRefPlusAltList.index(target)

                # get the DNA normal total depth
                dnDepth = sum(aDNANormalDepths)

                # if we don't have any normal DNA,
                # then remove the call and add a filter
                if (dnDepth == 0):
                    '''
                    if (anIsDebug):
                        logging.debug("SOM call with 0 total normal " +
                                      "DNA reads is being removed")
                    '''
                    modTypesList.remove(modType)
                    modChangesList.remove(modChange)
                    anAllFiltersSet.add("dnmndp")

                else:
                    # get the normal source depth
                    dnSourceDepth = aDNANormalDepths[sIndex]
                    dnSourcePct = round(dnSourceDepth/float(dnDepth), 2)

                    # if the normal ref depth doesn't reach the minimum,
                    # then remove the somatic call and add filters
                    if (dnSourceDepth < modMinDepth or
                        dnSourcePct < modMinPct):
                        modTypesList.remove(modType)
                        modChangesList.remove(modChange)
                        '''
                        if (anIsDebug):
                            logging.debug("SOM call with not enough normal " +
                                          "DNA (refDepth=%s < minDepth=%s " +
                                          "or refPct=%s < minPct=%s) is " +
                                          "being removed",
                                          dnSourceDepth, modMinDepth,
                                          dnSourcePct, modMinPct)
                            logging.debug("normalDepths=%s, tumorDepths=%s, " +
                                          "rnaDepths=%s", aDNANormalDepths,
                                          aDNATumorDepths, anRNATumorDepths)
                        '''
                        # add the filters
                        if (dnSourceDepth < modMinDepth):
                            anAllFiltersSet.add("dnmnrefad")
                        if (dnSourcePct < modMinPct):
                            anAllFiltersSet.add("dnmnrefaf")

            elif (modType == "NOR_EDIT" or modType == "RNA_NOR_VAR"):
                # get the indices
                sIndex = aRefPlusAltList.index(source)
                tIndex = aRefPlusAltList.index(target)

                # if this is labeled as a normal edit, but there are
                # variant reads in the normal DNA, maybe it is a germline
                # get the normal depth
                dnDepth = sum(aDNANormalDepths)

                # this is a hack for ones that were possibly
                # mis-classified as normal edits by radia,
                # but could've been classified as germline
                if (dnDepth > 0):
                    dnTargetDepth = aDNANormalDepths[tIndex]
                    dnTargetPct = round(dnTargetDepth/float(dnDepth), 2)

                    # if the normal depth is above the minimum, then add it
                    # (these were mis-classified by original radia script)
                    if (dnTargetDepth >= modMinDepth or
                        dnTargetPct >= modMinPct):
                        '''
                        if (anIsDebug):
                            logging.debug("NOR_EDIT or RNA_NOR_VAR call " +
                                          "with variant reads in the normal " +
                                          "DNA (altDepth=%s < minDepth=%s " +
                                          "or altPct=%s < minPct=%s) is " +
                                          "being tested as GERM",
                                          dnTargetDepth, modMinDepth,
                                          dnTargetPct, modMinPct)
                        '''
                        modTypesList.append("GERM")
                        modChangesList.append(modChange)
                    # these are calls with normal DNA reads,
                    # but not enough variant reads to be considered
                    # as germline, maybe they're edits
                    elif (modType == "RNA_NOR_VAR"):
                        '''
                        if (anIsDebug):
                            logging.debug("RNA_NOR_VAR call without variant " +
                                          "reads in the normal DNA " +
                                          "(altDepth=%s < minDepth=%s or " +
                                          "altPct=%s < minPct=%s) is being " +
                                          "tested as NOR_EDIT",
                                          dnTargetDepth, modMinDepth,
                                          dnTargetPct, modMinPct)
                        '''
                        modTypesList.append("NOR_EDIT")
                        modChangesList.append(modChange)

                # if this is classified as a normal edit, but the rna normal
                # variant depth is not sufficient, maybe it is a tumor edit
                rnTargetDepth = anRNANormalDepths[tIndex]
                if (rnTargetDepth < modMinDepth):

                    # get the tumor depth
                    rtDepth = sum(anRNATumorDepths)

                    # this is a hack for ones that were
                    # mis-classified as normal edits by radia,
                    # but should've been classified as tumor edits
                    if (rtDepth > 0):
                        rtTargetDepth = anRNATumorDepths[tIndex]
                        rtTargetPct = round(rtTargetDepth/float(rtDepth), 2)

                        # if the tumor depth is above the minimum, then add it
                        # (these were mis-classified by original radia script)
                        if (rtTargetDepth >= modMinDepth or
                            rtTargetPct >= modMinPct):
                            '''
                            if (anIsDebug):
                                logging.debug("NOR_EDIT or RNA_NOR_VAR call " +
                                              "with not enough variant " +
                                              "reads in the normal RNA " +
                                              "(altDepth=%s < minDepth=%s " +
                                              "or altPct=%s < minPct=%s) " +
                                              "is being tested as TUM_EDIT",
                                              rtTargetDepth, modMinDepth,
                                              rtTargetPct, modMinPct)
                            '''
                            modTypesList.append("TUM_EDIT")
                            modChangesList.append(modChange)

            elif (modType == "TUM_EDIT" or modType == "RNA_TUM_VAR"):
                # get the indices
                sIndex = aRefPlusAltList.index(source)
                tIndex = aRefPlusAltList.index(target)

                # if this is labeled as a tumor edit, but there are
                # variant reads in the DNA, maybe it is a somatic one
                # get the tumor depth
                dtDepth = sum(aDNATumorDepths)

                # this is a hack for ones that were possibly
                # mis-classified as tumor edits by radia,
                # but could've been classified as somatic
                if (dtDepth > 0):
                    dtTargetDepth = aDNATumorDepths[tIndex]
                    dtTargetPct = round(dtTargetDepth/float(dtDepth), 2)

                    # if the tumor depth is above the minimum, then add it
                    # (these were mis-classified by original radia script)
                    if (dtTargetDepth >= modMinDepth or
                        dtTargetPct >= modMinPct):
                        '''
                        if (anIsDebug):
                            logging.debug("TUM_EDIT or RNA_TUM_VAR call " +
                                          "with variant reads in the tumor " +
                                          "DNA (altDepth=%s < minDepth=%s " +
                                          "or altPct=%s < minPct=%s) is " +
                                          "being tested as SOM",
                                          dtTargetDepth, modMinDepth,
                                          dtTargetPct, modMinPct)
                        '''
                        modTypesList.append("SOM")
                        modChangesList.append(modChange)
                    # these are calls with DNA, but not enough DNA (by default
                    # 1 read) to be considered as somatic, maybe they're edits
                    elif (modType == "RNA_TUM_VAR" and
                          dtTargetDepth == 0):
                        '''
                        if (anIsDebug):
                            logging.debug("RNA_TUM_VAR call without variant " +
                                          "reads in the tumor DNA " +
                                          "(altDepth=%s < minDepth=%s or " +
                                          "altPct=%s < minPct=%s) is being " +
                                          "tested as TUM_EDIT",
                                          dtTargetDepth, modMinDepth,
                                          dtTargetPct, modMinPct)
                            logging.debug("normalDepths=%s, tumorDepths=%s, " +
                                          "rnaDepths=%s", aDNANormalDepths,
                                          aDNATumorDepths, anRNATumorDepths)
                        '''
                        modTypesList.append("TUM_EDIT")
                        modChangesList.append(modChange)

            elif (modType == "LOH"):

                # get the total depths
                dnDepths = sum(aDNANormalDepths)
                dtDepths = sum(aDNATumorDepths)

                validSources = []

                # for each of the possible sources
                for sourceNormal in source:
                    # if the source is the target, then just skip it,
                    # because it hasn't been lost
                    if (source == target):
                        continue

                    # get the source in the normal and tumor
                    sIndex = aRefPlusAltList.index(sourceNormal)
                    dnSourceDepth = aDNANormalDepths[sIndex]
                    dtSourceDepth = aDNATumorDepths[sIndex]
                    dnSourcePct = round(dnSourceDepth/float(dnDepths), 2)
                    dtSourcePct = round(dtSourceDepth/float(dtDepths), 2)

                    # if the normal depth is above the minimum, and the
                    # tumor depth is below the maximum, then it's valid
                    if (dnSourceDepth >= modMinDepth and
                        dnSourcePct >= modMinPct and
                        dtSourceDepth <= lohMaxDepth and
                        dtSourcePct <= lohMaxPct):
                        validSources.append(sourceNormal)
                    elif (modType == "SOM"):
                        if (dnSourceDepth < modMinDepth):
                            anAllFiltersSet.add("dnmnrefad")
                        if (dtSourceDepth < lohMaxDepth):
                            anAllFiltersSet.add("dtmnrefad")
                        if (dnSourcePct < modMinPct):
                            anAllFiltersSet.add("dnmnrefaf")
                        if (dtSourcePct < lohMaxPct):
                            anAllFiltersSet.add("dtmnrefaf")

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

            '''
            if (anIsDebug):
                logging.debug("final modTypes=%s, modChanges=%s",
                              modTypesList, modChangesList)
            '''
    except:
        logging.error("Error in pre_filter_mod_types(): aModTypeList=%s, " +
                      "aModChangeList=%s, aRefPlusAltList=%s, " +
                      "aDNANormalDepths=%s, anRNANormalDepths=%s, " +
                      "aDNATumorDepths=%s, anRNATumorDepths=%s",
                      str(aModTypeList), str(aModChangeList),
                      str(aRefPlusAltList), str(aDNANormalDepths),
                      str(anRNANormalDepths), str(aDNATumorDepths),
                      str(anRNATumorDepths))
        raise

    # if there are still some valid mod types, then return them
    if (len(modTypesList) > 0):
        anInfoDict["MT"] = modTypesList
        anInfoDict["MC"] = modChangesList
        return (anInfoDict, set())
    else:
        # otherwise just return the originals,
        # they will get filtered later anyway
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
        # in this case, pick the passing event in the following order:
        # GERM, NOR_EDIT, SOM, TUM_EDIT, RNA_TUM_VAR, LOH
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
    elif (finalModType.find("RNA_NOR_VAR") != -1 or
          finalModType.find("RNA_TUM_VAR") != 1):
        anInfoDict["SS"].append("5")
    # unknown
    else:
        anInfoDict["SS"].append("5")

    anInfoDict["MT"] = finalModTypeList
    anInfoDict["MC"] = finalModChangeList

    return (anInfoDict)


def filterByMapQualZero(aParamsDict, aSampleDict, aTargetIndex):

    if "MQ0" in aSampleDict:
        # check to make sure the sample does not exceed the maximum
        # percentage of MQ0 reads supporting the ALT
        mapQualZeroReads = int(aSampleDict["MQ0"][aTargetIndex])
        totalAltReads = int(aSampleDict["AD"][aTargetIndex])

        if (totalAltReads > 0):
            mq0Pct = (floor(mapQualZeroReads/float(totalAltReads)*100)/100)
            if (mq0Pct > float(aParamsDict["MaxAltMQ0Pct"])):
                return True

    return False


def filterByStrandBias(aParamsDict, aSampleDict, aSourceIndex, aTargetIndex):

    isStrandBiased = False

    sourceDepth = int(aSampleDict["AD"][aSourceIndex])
    targetDepth = int(aSampleDict["AD"][aTargetIndex])
    sourceStrbias = float(aSampleDict["SB"][aSourceIndex])
    targetStrbias = float(aSampleDict["SB"][aTargetIndex])

    # if we have enough depth for both the source and target
    if (sourceDepth >= aParamsDict["MinStrBiasDP"] and
        targetDepth >= aParamsDict["MinStrBiasDP"]):
        # allow 100% strand bias as long as both the
        # source and target strand bias are 100%
        if ((targetStrbias == 0.0 or targetStrbias == 1.0) and
            (sourceStrbias == 0.0 or sourceStrbias == 1.0)):
            isStrandBiased = False
        # otherwise, see if the target has a strand bias
        elif (targetStrbias > (aParamsDict["MaxStrandBias"]) or
              targetStrbias < (1.0 - aParamsDict["MaxStrandBias"])):
            isStrandBiased = True
    # see if the target has a strand bias
    elif (targetDepth >= aParamsDict["MinStrBiasDP"]):
        if (targetStrbias > (aParamsDict["MaxStrandBias"]) or
            targetStrbias < (1.0 - aParamsDict["MaxStrandBias"])):
            isStrandBiased = True

    return isStrandBiased


def filterByMaxError(aRefPlusAltList, aParamsDict, aSampleDict, aSourceIndex,
                     aTargetIndex, anIncludeTargetAlleles, anIsDebug):
    isMaxError = False
    errorCount = 0

    # we only allow a certain maximum of "other" alleles
    # for most cases, we want to make sure there isn't
    # an over-abundance of a 3rd allele at this position
    # for a somatic mutation, we check to make sure that
    # the normal sample doesn't contain too much of the
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

    maxErrPct = 0.0
    totalDepth = int(aSampleDict["DP"][0])
    if (totalDepth > 0):
        maxErrPct = (floor(errorCount/float(totalDepth)*100)/100)

    # if the number of "other" alleles is above the minimum count and
    # the "other" allele percentage is above the max percent
    if ((errorCount >= int(aParamsDict["MinErrPctDP"])) and
        (maxErrPct > float(aParamsDict["MaxErrPct"]))):
        isMaxError = True

    '''
    if (anIsDebug):
        logging.debug("MaxErrPct: isMaxError=%s, errorCount=%s, " +
                      "totalDepth=%s, errorPct=%s, errorPctParam=%s",
                      str(isMaxError), str(errorCount), str(totalDepth),
                      str(floor(errorCount/float(totalDepth)*100)/100),
                      str(float(aParamsDict["MaxErrPct"])))
    '''
    return isMaxError


def get_sample_columns(aFilename, aHeaderDict, anIsDebug):

    # get the file
    i_fileHandler = radiaUtil.get_read_fileHandler(aFilename)

    for line in i_fileHandler:

        # if it is an empty line, then just continue
        if (line.isspace()):
            continue

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        # if (anIsDebug):
        #    logging.debug("Line: %s", line)

        # if we find the column headers
        if ("#CHROM" in line):
            aHeaderDict["chrom"] = line + "\n"
            columnsLine = line.lstrip("#")
            # the tabs get removed somewhere along the
            # pipeline, so just split on whitespace
            columnsLineSplit = columnsLine.split()
            columnsList = columnsLineSplit[9:len(columnsLineSplit)]
            return columnsList, aHeaderDict

    return


def get_meta_id(aLine):
    metaId = aLine.split("=")[0]
    return metaId


def get_derived_sample(aLine):
    matchObj = i_headerDerivedRegEx.search(aLine)
    matchId = matchObj.group()
    matchId = matchId.rstrip(",")
    (tag, matchId) = matchId.split("=")
    # matchId = matchId.lstrip("Derived=")
    return matchId


def get_id(aLine):
    matchObj = i_headerIDRegEx.search(aLine)
    matchId = matchObj.group()
    matchId = matchId.rstrip(",")
    (tag, matchId) = matchId.split("=")
    # matchId = matchId.lstrip("ID=")
    return matchId


def add_header_data(aHeaderDict, aLine, anIsDebug):

    if (not aLine.startswith("#")):
        return aHeaderDict

    # extract the metadata
    if (aLine.startswith("##FORMAT")):
        formatId = get_id(aLine)
        aHeaderDict["format"][formatId] = aLine

    elif (aLine.startswith("##INFO")):
        infoId = get_id(aLine)
        aHeaderDict["info"][infoId] = aLine

    elif (aLine.startswith("##FILTER")):
        filterId = get_id(aLine)
        aHeaderDict["filter"][filterId] = aLine

    elif (aLine.startswith("##SAMPLE")):
        sampleId = get_id(aLine)
        aHeaderDict["sample"][sampleId] = aLine

    elif (aLine.startswith("##PEDIGREE")):
            derivedSample = get_derived_sample(aLine)
            aHeaderDict["pedigree"][derivedSample] = aLine

    elif (aLine.startswith("##")):
        metadataId = get_meta_id(aLine)
        aHeaderDict["metadata"][metadataId] = aLine

    elif (aLine.startswith("#CHROM")):
        aHeaderDict["chrom"] = aLine

    return aHeaderDict


def get_mpileup_header(anAddOriginFlag):

    headerDict = dict()
    headerDict["metadata"] = collections.defaultdict(str)
    headerDict["sample"] = collections.defaultdict(str)
    headerDict["pedigree"] = collections.defaultdict(str)
    headerDict["format"] = collections.defaultdict(str)
    headerDict["info"] = collections.defaultdict(str)
    headerDict["filter"] = collections.defaultdict(str)

    # if we're adding the origin tag, then add the INFO tag
    if (anAddOriginFlag):
        headerDict["info"]["ORIGIN"] = (
            "##INFO=<ID=ORIGIN,Number=.,Type=String,Description=\"Where " +
            "the call originated from, the tumor DNA, RNA, or both\">\n")

    # add all the filters
    headerDict["filter"]["blat"] = (
        "##FILTER=<ID=blat,Description=\"The call " +
        "did not pass the BLAT filter\">\n")
    headerDict["filter"]["indel"] = (
        "##FILTER=<ID=indel,Description=\"The number of INDELS " +
        "across all samples is above the maximum\">\n")
    headerDict["filter"]["multi"] = (
        "##FILTER=<ID=multi,Description=\"There are multiple " +
        "ALT alleles across all samples\">\n")
    headerDict["filter"]["rnacall"] = (
        "##FILTER=<ID=rnacall,Description=\"This is a dummy filter for a " +
        "call that originated in the RNA being filtered by the DNA\">\n")
    headerDict["filter"]["dnacall"] = (
        "##FILTER=<ID=dnacall,Description=\"This is a dummy filter for a " +
        "call that originated in the DNA being filtered by the RNA\">\n")

    headerDict["filter"]["dnmndp"] = (
        "##FILTER=<ID=dnmndp,Description=\"DNA Normal total depth " +
        "is less than the minimum\">\n")
    headerDict["filter"]["dtmndp"] = (
        "##FILTER=<ID=dtmndp,Description=\"DNA Tumor total depth " +
        "is less than the minimum\">\n")
    headerDict["filter"]["rnmndp"] = (
        "##FILTER=<ID=rnmndp,Description=\"RNA Normal total depth " +
        "is less than the minimum\">\n")
    headerDict["filter"]["rtmndp"] = (
        "##FILTER=<ID=rtmndp,Description=\"RNA Tumor total depth " +
        "is less than the minimum\">\n")

    headerDict["filter"]["dnmxdp"] = (
        "##FILTER=<ID=dnmxdp,Description=\"DNA Normal total depth " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["dtmxdp"] = (
        "##FILTER=<ID=dtmxdp,Description=\"DNA Tumor total depth " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["rnmxdp"] = (
        "##FILTER=<ID=rnmxdp,Description=\"RNA Normal total depth " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["rtmxdp"] = (
        "##FILTER=<ID=rtmxdp,Description=\"RNA Tumor total depth " +
        "is greater than the maximum\">\n")

    headerDict["filter"]["dnmnad"] = (
        "##FILTER=<ID=dnmnad,Description=\"DNA Normal ALT depth " +
        "is less than the minimum\">\n")
    headerDict["filter"]["dtmnad"] = (
        "##FILTER=<ID=dtmnad,Description=\"DNA Tumor ALT depth " +
        "is less than the minimum\">\n")
    headerDict["filter"]["rnmnad"] = (
        "##FILTER=<ID=rnmnad,Description=\"RNA Normal ALT depth " +
        "is less than the minimum\">\n")
    headerDict["filter"]["rtmnad"] = (
        "##FILTER=<ID=rtmnad,Description=\"RNA Tumor ALT depth " +
        "is less than the minimum\">\n")

    headerDict["filter"]["dnmnaf"] = (
        "##FILTER=<ID=dnmnaf,Description=\"DNA Normal ALT frequency " +
        "is less than the minimum\">\n")
    headerDict["filter"]["dtmnaf"] = (
        "##FILTER=<ID=dtmnaf,Description=\"DNA Tumor ALT frequency " +
        "is less than the minimum\">\n")
    headerDict["filter"]["rnmnaf"] = (
        "##FILTER=<ID=rnmnaf,Description=\"RNA Normal ALT frequency " +
        "is less than the minimum\">\n")
    headerDict["filter"]["rtmnaf"] = (
        "##FILTER=<ID=rtmnaf,Description=\"RNA Tumor ALT frequency " +
        "is less than the minimum\">\n")

    headerDict["filter"]["dnmnbq"] = (
        "##FILTER=<ID=dnmnbq,Description=\"DNA Normal average ALT " +
        "base quality is less than the minimum\">\n")
    headerDict["filter"]["dtmnbq"] = (
        "##FILTER=<ID=dtmnbq,Description=\"DNA Tumor average ALT " +
        "base quality is less than the minimum\">\n")
    headerDict["filter"]["rnmnbq"] = (
        "##FILTER=<ID=rnmnbq,Description=\"RNA Normal average ALT " +
        "base quality is less than the minimum\">\n")
    headerDict["filter"]["rtmnbq"] = (
        "##FILTER=<ID=rtmnbq,Description=\"RNA Tumor average ALT " +
        "base quality is less than the minimum\">\n")

    headerDict["filter"]["dnmnmqa"] = (
        "##FILTER=<ID=dnmnmqa,Description=\"DNA Normal average ALT " +
        "mapping quality is less than the minimum\">\n")
    headerDict["filter"]["dtmnmqa"] = (
        "##FILTER=<ID=dtmnmqa,Description=\"DNA Tumor average ALT " +
        "mapping quality is less than the minimum\">\n")
    headerDict["filter"]["rnmnmqa"] = (
        "##FILTER=<ID=rnmnmqa,Description=\"RNA Normal average ALT " +
        "mapping quality is less than the minimum\">\n")
    headerDict["filter"]["rtmnmqa"] = (
        "##FILTER=<ID=rtmnmqa,Description=\"RNA Tumor average ALT " +
        "mapping quality is less than the minimum\">\n")

    headerDict["filter"]["dnmnmq"] = (
        "##FILTER=<ID=dnmnmq,Description=\"DNA Normal has no " +
        "ALT reads with the minimum mapping quality\">\n")
    headerDict["filter"]["dtmnmq"] = (
        "##FILTER=<ID=dtmnmq,Description=\"DNA Tumor has no " +
        "ALT reads with the minimum mapping quality\">\n")
    headerDict["filter"]["rnmnmq"] = (
        "##FILTER=<ID=rnmnmq,Description=\"RNA Normal has no " +
        "ALT reads with the minimum mapping quality\">\n")
    headerDict["filter"]["rtmnmq"] = (
        "##FILTER=<ID=rtmnmq,Description=\"RNA Tumor has no " +
        "ALT reads with the minimum mapping quality\">\n")

    headerDict["filter"]["dnmxmq0"] = (
        "##FILTER=<ID=dnmxmq0,Description=\"DNA Normal " +
        "percentage of mapping quality zero reads for the ALT " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["dtmxmq0"] = (
        "##FILTER=<ID=dtmxmq0,Description=\"DNA Tumor " +
        "percentage of mapping quality zero reads for the ALT " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["rnmxmq0"] = (
        "##FILTER=<ID=rnmxmq0,Description=\"RNA Normal " +
        "percentage of mapping quality zero reads for the ALT " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["rtmxmq0"] = (
        "##FILTER=<ID=rtmxmq0,Description=\"RNA Tumor " +
        "percentage of mapping quality zero reads for the ALT " +
        "is greater than the maximum\">\n")

    headerDict["filter"]["dnmnrefad"] = (
        "##FILTER=<ID=dnmnrefad,Description=\"DNA Normal REF depth " +
        "is less than the minimum\">\n")
    headerDict["filter"]["dnmnrefaf"] = (
        "##FILTER=<ID=dnmnrefaf,Description=\"DNA Normal REF frequency " +
        "is less than the minimum\">\n")
    headerDict["filter"]["dtmnrefad"] = (
        "##FILTER=<ID=dtmnrefad,Description=\"DNA Tumor REF depth " +
        "is less than the minimum\">\n")
    headerDict["filter"]["dtmnrefaf"] = (
        "##FILTER=<ID=dtmnrefaf,Description=\"DNA Tumor REF frequency " +
        "is less than the minimum\">\n")

    headerDict["filter"]["dnmxerr"] = (
        "##FILTER=<ID=dnmxerr,Description=\"DNA Normal total ALT percentage " +
        "attributed to error (sequencing, contamination, etc.) " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["dtmxerr"] = (
        "##FILTER=<ID=dtmxerr,Description=\"DNA Tumor total ALT percentage " +
        "attributed to error (sequencing, contamination, etc.) " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["rnmxerr"] = (
        "##FILTER=<ID=rnmxerr,Description=\"RNA Normal total ALT percentage " +
        "attributed to error (sequencing, contamination, etc.) " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["rtmxerr"] = (
        "##FILTER=<ID=rtmxerr,Description=\"RNA Tumor total ALT percentage " +
        "attributed to error (sequencing, contamination, etc.) " +
        "is greater than the maximum\">\n")

    headerDict["filter"]["dnsbias"] = (
        "##FILTER=<ID=dnsbias,Description=\"DNA Normal strand bias, " +
        "majority of reads supporting ALT " +
        "are on forward OR reverse strand\">\n")
    headerDict["filter"]["dtsbias"] = (
        "##FILTER=<ID=dtsbias,Description=\"DNA Tumor strand bias, " +
        "majority of reads supporting ALT " +
        "are on forward OR reverse strand\">\n")
    headerDict["filter"]["rnsbias"] = (
        "##FILTER=<ID=rnsbias,Description=\"RNA Normal strand bias, " +
        "majority of reads supporting ALT " +
        "are on forward OR reverse strand\">\n")
    headerDict["filter"]["rtsbias"] = (
        "##FILTER=<ID=rtsbias,Description=\"RNA Tumor strand bias, " +
        "majority of reads supporting ALT " +
        "are on forward OR reverse strand\">\n")

    headerDict["filter"]["dnmxindel"] = (
        "##FILTER=<ID=dnmxindel,Description=\"DNA Normal INDEL count " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["dtmxindel"] = (
        "##FILTER=<ID=dtmxindel,Description=\"DNA Tumor INDEL count " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["rnmxindel"] = (
        "##FILTER=<ID=rnmxindel,Description=\"RNA Normal INDEL count " +
        "is greater than the maximum\">\n")
    headerDict["filter"]["rtmxindel"] = (
        "##FILTER=<ID=rtmxindel,Description=\"RNA Tumor INDEL count " +
        "is greater than the maximum\">\n")

    return headerDict


def get_vcf_header(aHeaderDict, aFilename, aCmdLineParams,
                   aColumnsList, anIsDebug):

    # open the file
    vcfFileHandler = radiaUtil.get_read_fileHandler(aFilename)

    for line in vcfFileHandler:

        # if (anIsDebug):
        #    logging.debug("read vcfLine: %s", line)

        # if it is an empty line, then just continue
        if (line.isspace()):
            continue

        # get the header info
        if (line.startswith("##FORMAT")):
            formatId = get_id(line)
            aHeaderDict["format"][formatId] = line

        elif (line.startswith("##INFO")):
            infoId = get_id(line)
            aHeaderDict["info"][infoId] = line

        elif (line.startswith("##FILTER")):
            filterId = get_id(line)
            aHeaderDict["filter"][filterId] = line

        elif (line.startswith("##SAMPLE")):
            sampleId = get_id(line)
            aHeaderDict["sample"][sampleId] = line

        elif (line.startswith("##PEDIGREE")):
            derivedSample = get_derived_sample(line)
            aHeaderDict["pedigree"][derivedSample] = line

        # if we find the vcfGenerator line,
        # then update the params from the user
        elif ("vcfGenerator" in line):
            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")

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

            # for each new param
            for (paramName, paramValue) in aCmdLineParams.iteritems():
                # don't output the defaults for files that aren't specified
                if (paramName.startswith("dnaNormal") and
                    "DNA_NORMAL" not in aColumnsList):
                    continue
                elif (paramName.startswith("rnaNormal") and
                      "RNA_NORMAL" not in aColumnsList):
                    continue
                elif (paramName.startswith("dnaTumor") and
                      "DNA_TUMOR" not in aColumnsList):
                    continue
                elif (paramName.startswith("rnaTumor") and
                      "RNA_TUMOR" not in aColumnsList):
                    continue
                # add new params and overwrite the old params with the new ones
                else:
                    generatorParamsDict[paramName] = paramValue

            generatorOut = "##vcfGenerator=<"
            # make it pretty by sorting on the keys
            for (paramName) in sorted(generatorParamsDict.iterkeys()):
                paramValue = generatorParamsDict[paramName]
                if (paramValue is not None):
                    generatorOut += paramName + "=<" + str(paramValue) + ">,"
            generatorOut = generatorOut.rstrip(",")
            generatorOut += ">\n"
            aHeaderDict["metadata"]["vcfGenerator"] = generatorOut

        elif (line.startswith("##")):
            metadataId = get_meta_id(line)
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


def filter_by_mpileup_support(anId, aChrom, aVCFFilename, aHeaderFilename,
                              anOutputFilename, aFilterUsingRNAFlag,
                              anAddOriginFlag, aCmdLineParams,
                              aDNPs, aDTPs, anRNPs, anRTPs,
                              aParamsDict, anIsDebug):
    '''
    ' This function filters based on the mpileup read support.
    '
    ' anId:
    '    The Id for this sample
    ' aChrom:
    '    The chromosome being filtered
    ' aVCFFilename:
    '    The filename to be filtered
    ' aHeaderFilename:
    '    The filename with full header info
    ' anOutputFilename:
    '    The output filename
    ' aFilterUsingRNAFlag:
    '    If the calls should be filtered by the RNA as well
    ' anAddOriginFlag:
    '    If the origin (DNA or RNA) of the call should be added to the INFO tag
    ' aCmdLineParams:
    '    All the parameters specified by the user
    ' aDNPs:
    '    The parameters for the normal DNA
    ' aDTPs:
    '    The parameters for the tumor DNA
    ' anRNPs:
    '    The parameters for the normal RNA
    ' anRTPs:
    '    The parameters for the tumor RNA
    ' aParamsDict:
    '    General parameters
    ' anIsDebug:
    '    A flag for outputting debug messages to STDERR
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
        i_outputFileHandler = radiaUtil.get_write_fileHandler(anOutputFilename)

    # get the mpileup filter header lines
    headerDict = get_mpileup_header(anAddOriginFlag)

    # sometimes the header lines are stripped from the file,
    # so get the necessary columnsList from the #CHROM line
    # from a header file if it is specified
    if (aHeaderFilename is not None):
        columnsList, headerDict = get_sample_columns(aHeaderFilename,
                                                     headerDict,
                                                     anIsDebug)
        headerDict = get_vcf_header(headerDict,
                                    aHeaderFilename,
                                    aCmdLineParams,
                                    columnsList,
                                    anIsDebug)
    else:
        columnsList, headerDict = get_sample_columns(aVCFFilename,
                                                     headerDict,
                                                     anIsDebug)
        headerDict = get_vcf_header(headerDict,
                                    aVCFFilename,
                                    aCmdLineParams,
                                    columnsList,
                                    anIsDebug)

    # output the header information
    output_header(headerDict["metadata"], False, i_outputFileHandler)
    output_header(headerDict["sample"], False, i_outputFileHandler)
    output_header(headerDict["pedigree"], False, i_outputFileHandler)
    output_header(headerDict["filter"], True, i_outputFileHandler)
    output_header(headerDict["info"], True, i_outputFileHandler)
    output_header(headerDict["format"], True, i_outputFileHandler)
    i_outputFileHandler.write(headerDict["chrom"])

    # get the file
    i_vcfFileHandler = radiaUtil.get_read_fileHandler(aVCFFilename)

    # for each event in the vcf file
    for line in i_vcfFileHandler:

        # if it is an empty line or header line, then just continue
        if (line.isspace() or line.startswith("#")):
            continue

        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        # if (anIsDebug):
        #    logging.debug("VCF Line: %s", line)

        # now we are to the data
        # count the total events
        totalEvents += 1
        somEventWithTumorRna = False
        somEventWithTumorAltRna = False

        # split the line on the tab
        splitLine = line.split("\t")

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

        # if we should add the origin to the info column
        if (anAddOriginFlag):
            if (aFilterUsingRNAFlag):
                origin = "RNA"
            else:
                origin = "DNA"

            if ("ORIGIN" in infoDict):
                originList = infoDict["ORIGIN"]
                if (origin not in originList):
                    originList.append(origin)
            else:
                infoDict["ORIGIN"] = [origin]

        # get the event format list
        formatList = splitLine[8].split(":")

        # initialize the optional columns to none
        dnaNormalList = None
        dnaTumorList = None
        rnaNormalList = None
        rnaTumorList = None

        # if we have a 9th column, figure out which dataset it is
        if (len(splitLine) > 9):
            if (columnsList[0] == "DNA_NORMAL"):
                dnaNormalList = splitLine[9].split(":")
            elif (columnsList[0] == "RNA_NORMAL"):
                rnaNormalList = splitLine[9].split(":")
            elif (columnsList[0] == "DNA_TUMOR"):
                dnaTumorList = splitLine[9].split(":")
            elif (columnsList[0] == "RNA_TUMOR"):
                rnaTumorList = splitLine[9].split(":")
        # if we have a 10th column, figure out which dataset it is
        if (len(splitLine) > 10):
            if (columnsList[1] == "RNA_NORMAL"):
                rnaNormalList = splitLine[10].split(":")
            elif (columnsList[1] == "DNA_TUMOR"):
                dnaTumorList = splitLine[10].split(":")
            elif (columnsList[1] == "RNA_TUMOR"):
                rnaTumorList = splitLine[10].split(":")
        # if we have a 11th column, figure out which dataset it is
        if (len(splitLine) > 11):
            if (columnsList[2] == "DNA_TUMOR"):
                dnaTumorList = splitLine[11].split(":")
            elif (columnsList[2] == "RNA_TUMOR"):
                rnaTumorList = splitLine[11].split(":")
        # if we have a 12th column, figure out which dataset it is
        if (len(splitLine) > 12):
            if (columnsList[3] == "RNA_TUMOR"):
                rnaTumorList = splitLine[12].split(":")

        haveDnaNormData = True
        haveRnaNormData = True
        haveDnaTumData = True
        haveRnaTumData = True

        # if there is no data, then set the flag
        if (dnaNormalList is None or
            dnaNormalList[0] == "." or
            dnaNormalList[0] == "./."):
            haveDnaNormData = False
        # if there is no data, then set the flag
        if (rnaNormalList is None or
            rnaNormalList[0] == "." or
            rnaNormalList[0] == "./."):
            haveRnaNormData = False
        # if there is no data, then set the flag
        if (dnaTumorList is None or
            dnaTumorList[0] == "." or
            dnaTumorList[0] == "./."):
            haveDnaTumData = False
        # if there is no data, then set the flag
        if (rnaTumorList is None or
            rnaTumorList[0] == "." or
            rnaTumorList[0] == "./."):
            haveRnaTumData = False

        # parse the dna and rna columns and create dicts for each
        dnDict = collections.defaultdict(list)
        dtDict = collections.defaultdict(list)
        rnDict = collections.defaultdict(list)
        rtDict = collections.defaultdict(list)

        index = 0
        for formatItem in formatList:
            if (formatItem == "GT"):
                sep = "/"
            else:
                sep = ","

            if (haveDnaNormData):
                dnaNormalItem = dnaNormalList[index]
                dnDict[formatItem] = dnaNormalItem.split(sep)
            if (haveRnaNormData):
                rnaNormalItem = rnaNormalList[index]
                rnDict[formatItem] = rnaNormalItem.split(sep)
            if (haveDnaTumData):
                dnaTumorItem = dnaTumorList[index]
                dtDict[formatItem] = dnaTumorItem.split(sep)
            if (haveRnaTumData):
                rnaTumorItem = rnaTumorList[index]
                rtDict[formatItem] = rnaTumorItem.split(sep)
            index += 1

        # fix the original genotypes
        gtIndex = formatList.index("GT")
        if (haveDnaNormData):
            dnDict["GT"] = fix_genotypes(event_chr,
                                         event_refList,
                                         event_altList,
                                         map(int, dnDict["AD"]),
                                         aParamsDict)
            dnaNormalList[gtIndex] = "/".join(map(str, dnDict["GT"]))
        if (haveRnaNormData):
            rnDict["GT"] = fix_genotypes(event_chr,
                                         event_refList,
                                         event_altList,
                                         map(int, rnDict["AD"]),
                                         aParamsDict)
            rnaNormalList[gtIndex] = "/".join(map(str, rnDict["GT"]))
        if (haveDnaTumData):
            dtDict["GT"] = fix_genotypes(event_chr,
                                         event_refList,
                                         event_altList,
                                         map(int, dtDict["AD"]),
                                         aParamsDict)
            dnaTumorList[gtIndex] = "/".join(map(str, dtDict["GT"]))
        if (haveRnaTumData):
            rtDict["GT"] = fix_genotypes(event_chr,
                                         event_refList,
                                         event_altList,
                                         map(int, rtDict["AD"]),
                                         aParamsDict)
            rnaTumorList[gtIndex] = "/".join(map(str, rtDict["GT"]))

        # combine the refs and alts in one list
        refPlusAltList = event_refList + event_altList

        allFiltersSet = set()

        # get rid of bad mod types that don't meet the minimum requirements
        '''
        if (anIsDebug):
            logging.debug("before pre_filter_mod_types(): " +
                          "modTypes=%s, modChanges=%s",
                          list(infoDict["MT"]),
                          list(infoDict["MC"]))
        '''
        (infoDict,
         allFiltersSet) = pre_filter_mod_types(refPlusAltList,
                                               allFiltersSet,
                                               infoDict,
                                               map(int, dnDict["AD"]),
                                               map(int, rnDict["AD"]),
                                               map(int, dtDict["AD"]),
                                               map(int, rtDict["AD"]),
                                               aParamsDict,
                                               anIsDebug)
        '''
        if (anIsDebug):
            logging.debug("after pre_filter_mod_types(): " +
                          "modTypes=%s, modChanges=%s",
                          list(infoDict["MT"]),
                          list(infoDict["MC"]))
        '''

        # make copies of the lists to manipulate
        modTypesList = list(infoDict["MT"])
        modChangesList = list(infoDict["MC"])

        # keep track of filters for each mod to add to INFO
        modFilterTypes = []
        modFilters = []

        # for each modification type and change
        for (modType, modChange) in izip(infoDict["MT"],
                                         infoDict["MC"]):
            isValidMod = True
            filterSet = set()

            # get the source and target alleles
            (source, target) = modChange.split(">")

            if (modType == "GERM"):
                sIndex = refPlusAltList.index(source)
                tIndex = refPlusAltList.index(target)

                # check to make sure the normal DNA sample
                # is between the min and the max of total bases
                dnDepth = int(dnDict["DP"][0])
                if (dnDepth < aDNPs["MinDepth"]):
                    isValidMod = False
                    filterSet.add("dnmndp")

                elif (dnDepth >= aDNPs["MaxDepth"]):
                    isValidMod = False
                    filterSet.add("dnmxdp")

                # check to make sure the normal DNA sample
                # number of ALT bases is above the min
                if (int(dnDict["AD"][tIndex]) < aDNPs["MinAltDepth"]):
                    isValidMod = False
                    filterSet.add("dnmnad")

                # check to make sure the normal DNA sample
                # percentage of ALT bases is above the min
                # adjust the minAltPct for purity
                dnMinAltPct = round(aDNPs["MinAltPct"] / aDNPs["Purity"], 2)
                if (float(dnDict["AF"][tIndex]) < dnMinAltPct):
                    isValidMod = False
                    filterSet.add("dnmnaf")

                # check to make sure the normal DNA sample
                # average base quality for ALT bases is above the min
                if (int(dnDict["BQ"][tIndex]) < aDNPs["MinAltAvgBQ"]):
                    isValidMod = False
                    filterSet.add("dnmnbq")

                # check to make sure the normal DNA sample
                # average mapping quality for ALT reads is above the min
                if (int(dnDict["MQA"][tIndex]) < aDNPs["MinAltMQA"]):
                    isValidMod = False
                    filterSet.add("dnmnmqa")

                # check to make sure the normal DNA sample has at least
                # 1 ALT read with a mapping quality above the min
                if (int(dnDict["MMQ"][tIndex]) < aDNPs["MinAltMMQ"]):
                    isValidMod = False
                    filterSet.add("dnmnmq")

                # check to make sure the normal DNA sample has a
                # maximum percentage of MQ0 reads supporting the ALT
                if (filterByMapQualZero(aDNPs,
                                        dnDict,
                                        tIndex)):
                    isValidMod = False
                    filterSet.add("dnmxmq0")

                # check to make sure the number of normal DNA
                # sample INDELs is below the maximum
                indels = int(dnDict["INS"][0]) + int(dnDict["DEL"][0])
                if (indels >= aDNPs["MaxIndels"]):
                    isValidMod = False
                    filterSet.add("dnmxindel")

                # check to make sure the normal variant
                # reads don't have a strand bias
                if (filterByStrandBias(aDNPs,
                                       dnDict,
                                       sIndex,
                                       tIndex)):
                    isValidMod = False
                    filterSet.add("dnsbias")

                # we want to make sure the normal DNA sample
                # error percentage is below the max
                # we want to make sure that the percentage of
                # other ALTs in this sample is below the max error
                # this is the same as making sure that the percentage
                # of the source and target alleles is above one minus
                # the max error percentage
                if (filterByMaxError(refPlusAltList,
                                     aDNPs,
                                     dnDict,
                                     sIndex,
                                     tIndex,
                                     False,
                                     anIsDebug)):
                    isValidMod = False
                    filterSet.add("dnmxerr")

                # if we are also filtering using the RNA
                if (aFilterUsingRNAFlag):
                    # check to make sure the normal RNA sample has data
                    if (haveRnaNormData):
                        # check to make sure the normal RNA sample
                        # total depth is between the min and the max
                        if (int(rnDict["DP"][0]) < anRNPs["MinDepth"]):
                            isValidMod = False
                            filterSet.add("rnmndp")
                        elif (int(rnDict["DP"][0]) >= anRNPs["MaxDepth"]):
                            isValidMod = False
                            filterSet.add("rnmxdp")

                        # check to make sure the normal RNA sample
                        # number of ALT bases is above the min
                        if (int(rnDict["AD"][tIndex]) < anRNPs["MinAltDepth"]):
                            isValidMod = False
                            filterSet.add("rnmnad")

                        # check to make sure the normal RNA sample
                        # percentage of ALT bases is above the min
                        if (float(rnDict["AF"][tIndex]) < anRNPs["MinAltPct"]):
                            isValidMod = False
                            filterSet.add("rnmnaf")

                        # check to make sure the normal RNA sample average
                        # base quality for ALT bases is above the min
                        if (int(rnDict["BQ"][tIndex]) < anRNPs["MinAltAvgBQ"]):
                            isValidMod = False
                            filterSet.add("rnmnbq")

                        # check to make sure the normal RNA sample average
                        # mapping quality for ALT reads is above the min
                        if (int(rnDict["MQA"][tIndex]) < anRNPs["MinAltMQA"]):
                            isValidMod = False
                            filterSet.add("rnmnmqa")

                        # check to make sure the normal RNA sample has at
                        # least 1 ALT read with a mapping quality above the min
                        if (int(rnDict["MMQ"][tIndex]) < anRNPs["MinAltMMQ"]):
                            isValidMod = False
                            filterSet.add("rnmnmq")

                        # check to make sure the normal RNA sample has a
                        # maximum percentage of MQ0 reads supporting the ALT
                        if (filterByMapQualZero(anRNPs, rnDict, tIndex)):
                            isValidMod = False
                            filterSet.add("rnmxmq0")

                        # check to make sure the number of normal RNA
                        # sample INDELs is below the maximum
                        indels = int(rnDict["INS"][0]) + int(rnDict["DEL"][0])
                        if (indels >= anRNPs["MaxIndels"]):
                            isValidMod = False
                            filterSet.add("rnmxindel")

                        # check to make sure the normal variant
                        # reads don't have a strand bias
                        if (filterByStrandBias(anRNPs,
                                               rnDict,
                                               sIndex,
                                               tIndex)):
                            isValidMod = False
                            filterSet.add("rnsbias")

                        # we want to make sure the normal DNA sample
                        # error percentage is below the max
                        # we want to make sure that the percentage of
                        # other ALTs in this sample is below the max error
                        if (filterByMaxError(refPlusAltList,
                                             anRNPs,
                                             rnDict,
                                             sIndex,
                                             tIndex,
                                             False,
                                             anIsDebug)):
                            isValidMod = False
                            filterSet.add("rnmxerr")

                    # else if a minimum amount of total bases were required,
                    # but none were found, then set the filter
                    elif (anRNPs["MinDepth"] > 0):
                        isValidMod = False
                        filterSet.add("dnacall")

            elif (modType.find("NOR_EDIT") != -1 or
                  modType.find("RNA_NOR_VAR") != -1):
                sIndex = refPlusAltList.index(source)
                tIndex = refPlusAltList.index(target)

                # if we are also filtering using the RNA
                if (aFilterUsingRNAFlag):
                    # if this is a normal edit, then we need to check the DNA
                    # if this is an RNA normal variant,
                    # then there isn't any DNA to check
                    if (modType.find("NOR_EDIT") != -1):
                        # check to make sure the normal DNA sample has data
                        # and is between the min and the max of total bases
                        if (haveDnaNormData):
                            if (int(dnDict["DP"][0]) < aDNPs["MinDepth"]):
                                isValidMod = False
                                filterSet.add("dnmndp")
                            elif (int(dnDict["DP"][0]) >= aDNPs["MaxDepth"]):
                                isValidMod = False
                                filterSet.add("dnmxdp")
                            # we want to make sure the normal DNA sample
                            # error percentage is below the max
                            # we want to make sure that the percentage of
                            # other ALTs in this sample is below the max error
                            if (filterByMaxError(refPlusAltList,
                                                 aDNPs,
                                                 dnDict,
                                                 sIndex,
                                                 tIndex,
                                                 True,
                                                 anIsDebug)):
                                isValidMod = False
                                filterSet.add("dnmxerr")
                        # else if a minimum amount of total bases were required
                        # but none were found, then set the filter
                        elif (aDNPs["MinDepth"] > 0):
                            isValidMod = False
                            filterSet.add("dnmndp")

                    # check to make sure the normal RNA sample is
                    # between the min and the max of total bases
                    if (int(rnDict["DP"][0]) < anRNPs["MinDepth"]):
                        isValidMod = False
                        filterSet.add("rnmndp")
                    elif (int(rnDict["DP"][0]) >= anRNPs["MaxDepth"]):
                        isValidMod = False
                        filterSet.add("rnmxdp")

                    # check to make sure the normal RNA sample
                    # number of ALT bases is above the min
                    if (int(rnDict["AD"][tIndex]) < anRNPs["MinAltDepth"]):
                        isValidMod = False
                        filterSet.add("rnmnad")

                    # check to make sure the normal RNA sample
                    # percentage of ALT bases is above the min
                    if (float(rnDict["AF"][tIndex]) < anRNPs["MinAltPct"]):
                        isValidMod = False
                        filterSet.add("rnmnaf")

                    # check to make sure the normal RNA sample average
                    # base quality for ALT bases is above the min
                    if (int(rnDict["BQ"][tIndex]) < anRNPs["MinAltAvgBQ"]):
                        isValidMod = False
                        filterSet.add("rnmnbq")

                    # check to make sure the normal RNA sample average
                    # map quality for ALT reads is above the min
                    if (int(rnDict["MQA"][tIndex]) < anRNPs["MinAltMQA"]):
                        isValidMod = False
                        filterSet.add("rnmnmqa")

                    # check to make sure the normal RNA sample has at least
                    # 1 ALT read with a mapping quality above the min
                    if (int(rnDict["MMQ"][tIndex]) < anRNPs["MinAltMMQ"]):
                        isValidMod = False
                        filterSet.add("rnmnmq")

                    # check to make sure the normal RNA sample has a maximum
                    # percentage of MQ0 reads supporting the ALT
                    if (filterByMapQualZero(anRNPs,
                                            rnDict,
                                            tIndex)):
                        isValidMod = False
                        filterSet.add("rnmxmq0")

                    # check to make sure the number of normal RNA sample
                    # INDELs is below the maximum
                    indels = int(rnDict["INS"][0]) + int(rnDict["DEL"][0])
                    if (indels >= anRNPs["MaxIndels"]):
                        isValidMod = False
                        filterSet.add("rnmxindel")

                    # check to make sure the normal variant
                    # reads don't have a strand bias
                    if (filterByStrandBias(anRNPs,
                                           rnDict,
                                           sIndex,
                                           tIndex)):
                        isValidMod = False
                        filterSet.add("rnsbias")

                    # we want to make sure the normal RNA sample
                    # error percentage is below the max
                    # we want to make sure that the percentage of
                    # other ALTs in this sample is below the max error
                    if (filterByMaxError(refPlusAltList,
                                         anRNPs,
                                         rnDict,
                                         sIndex,
                                         tIndex,
                                         False,
                                         anIsDebug)):
                        isValidMod = False
                        filterSet.add("rnmxerr")
                # we are filtering via the DNA, so put in a
                # dummy filter so that they don't pass
                else:
                    isValidMod = False
                    filterSet.add("rnacall")

            elif (modType == "SOM"):
                sIndex = refPlusAltList.index(source)
                tIndex = refPlusAltList.index(target)

                # check to make sure the tumor DNA sample is
                # between the min and the max of total bases
                if (int(dtDict["DP"][0]) < aDTPs["MinDepth"]):
                    isValidMod = False
                    filterSet.add("dtmndp")
                elif (int(dtDict["DP"][0]) >= aDTPs["MaxDepth"]):
                    isValidMod = False
                    filterSet.add("dtmxdp")

                # check to make sure the tumor DNA sample
                # number of ALT bases is above the min
                if (int(dtDict["AD"][tIndex]) < aDTPs["MinAltDepth"]):
                    isValidMod = False
                    filterSet.add("dtmnad")
                # check to make sure the tumor DNA sample
                # percentage of ALT bases is above the min
                # adjust the minAltPct for purity
                dtMinAltPct = round(aDTPs["MinAltPct"] * aDTPs["Purity"], 2)
                if (float(dtDict["AF"][tIndex]) < dtMinAltPct):
                    isValidMod = False
                    filterSet.add("dtmnaf")

                # check to make sure the tumor DNA sample average
                # base quality for ALT bases is above the min
                if (int(dtDict["BQ"][tIndex]) < aDTPs["MinAltAvgBQ"]):
                    isValidMod = False
                    filterSet.add("dtmnbq")

                # check to make sure the tumor DNA sample average
                # map quality for ALT reads is above the min
                if (int(dtDict["MQA"][tIndex]) < aDTPs["MinAltMQA"]):
                    isValidMod = False
                    filterSet.add("dtmnmqa")

                # check to make sure the tumor DNA sample has at least
                # 1 ALT read with a mapping quality above the min
                if (int(dtDict["MMQ"][tIndex]) < aDTPs["MinAltMMQ"]):
                    isValidMod = False
                    filterSet.add("dtmnmq")

                # check to make sure the tumor DNA sample has a
                # maximum percentage of MQ0 reads supporting the ALT
                if (filterByMapQualZero(aDTPs,
                                        dtDict,
                                        tIndex)):
                    isValidMod = False
                    filterSet.add("dtmxmq0")

                # check to make sure the number of tumor DNA sample
                # INDELs is below the maximum
                indels = int(dtDict["INS"][0]) + int(dtDict["DEL"][0])
                if (indels >= aDTPs["MaxIndels"]):
                    isValidMod = False
                    filterSet.add("dtmxindel")

                # check to make sure the tumor variant reads
                # don't have a strand bias
                if (filterByStrandBias(aDTPs,
                                       dtDict,
                                       sIndex,
                                       tIndex)):
                    isValidMod = False
                    filterSet.add("dtsbias")

                # we want to make sure the tumor DNA sample
                # error percentage is below the max
                # we want to make sure that the percentage of
                # other ALTs in this sample is below the max error
                if (filterByMaxError(refPlusAltList,
                                     aDTPs,
                                     dtDict,
                                     sIndex,
                                     tIndex,
                                     False,
                                     anIsDebug)):
                    isValidMod = False
                    filterSet.add("dtmxerr")

                # check to make sure the normal DNA sample has data
                # and is between the min and the max of total bases
                if (haveDnaNormData):
                    if (int(dnDict["DP"][0]) < aDNPs["MinDepth"]):
                        isValidMod = False
                        filterSet.add("dnmndp")
                    elif (int(dnDict["DP"][0]) >= aDNPs["MaxDepth"]):
                        isValidMod = False
                        filterSet.add("dnmxdp")
                    # we want to make sure the normal DNA sample
                    # error percentage is below the max
                    # we want to make sure that the percentage of
                    # other ALTs in this sample is below the max error
                    if (filterByMaxError(refPlusAltList,
                                         aDNPs,
                                         dnDict,
                                         sIndex,
                                         tIndex,
                                         True,
                                         anIsDebug)):
                        isValidMod = False
                        filterSet.add("dnmxerr")
                elif (aDNPs["MinDepth"] > 0):
                    isValidMod = False
                    filterSet.add("dnmndp")

                # set some flags
                if (haveRnaTumData):
                    # check if there is any RNA
                    if (int(rtDict["DP"][0]) > 1):
                        somEventWithTumorRna = True
                    # check if there are any Alt RNA
                    if (int(rtDict["AD"][tIndex]) > 1):
                        somEventWithTumorAltRna = True

                # if we are also filtering using the RNA
                if (aFilterUsingRNAFlag):
                    # check to make sure the tumor RNA sample has data
                    # and is between the min and the max of total bases
                    if (haveRnaTumData):
                        if (int(rtDict["DP"][0]) < anRTPs["MinDepth"]):
                            isValidMod = False
                            filterSet.add("rtmndp")
                        elif (int(rtDict["DP"][0]) >= anRTPs["MaxDepth"]):
                            isValidMod = False
                            filterSet.add("rtmxdp")

                        # check to make sure the tumor RNA sample
                        # number of ALT bases is above the min
                        if (int(rtDict["AD"][tIndex]) < anRTPs["MinAltDepth"]):
                            isValidMod = False
                            filterSet.add("rtmnad")

                        # check to make sure the tumor RNA sample
                        # percentage of ALT bases is above the min
                        if (float(rtDict["AF"][tIndex]) < anRTPs["MinAltPct"]):
                            isValidMod = False
                            filterSet.add("rtmnaf")

                        # check to make sure the tumor RNA sample average
                        # base quality for ALT bases is above the min
                        if (int(rtDict["BQ"][tIndex]) < anRTPs["MinAltAvgBQ"]):
                            isValidMod = False
                            filterSet.add("rtmnbq")

                        # check to make sure the tumor RNA sample average
                        # map quality for ALT reads is above the min
                        if (int(rtDict["MQA"][tIndex]) < anRTPs["MinAltMQA"]):
                            isValidMod = False
                            filterSet.add("rtmnmqa")

                        # check to make sure the tumor RNA sample has at
                        # least 1 ALT read with a mapping quality above the min
                        if (int(rtDict["MMQ"][tIndex]) < anRTPs["MinAltMMQ"]):
                            isValidMod = False
                            filterSet.add("rtmnmq")

                        # check to make sure the tumor RNA sample has a
                        # maximum percentage of MQ0 reads supporting the ALT
                        if (filterByMapQualZero(anRTPs,
                                                rtDict,
                                                tIndex)):
                            isValidMod = False
                            filterSet.add("rtmxmq0")

                        # check to make sure the number of tumor RNA
                        # sample INDELs is below the maximum
                        indels = int(rtDict["INS"][0]) + int(rtDict["DEL"][0])
                        if (indels >= anRTPs["MaxIndels"]):
                            isValidMod = False
                            filterSet.add("rtmxindel")

                        # check to make sure the tumor variant reads
                        # don't have a strand bias
                        if (filterByStrandBias(anRTPs,
                                               rtDict,
                                               sIndex,
                                               tIndex)):
                            isValidMod = False
                            filterSet.add("rtsbias")

                        # we want to make sure the tumor RNA sample
                        # error percentage is below the max
                        # we want to make sure that the percentage of
                        # other ALTs in this sample is below the max error
                        if (filterByMaxError(refPlusAltList,
                                             anRTPs,
                                             rtDict,
                                             sIndex,
                                             tIndex,
                                             False,
                                             anIsDebug)):
                            isValidMod = False
                            filterSet.add("rtmxerr")

                    # else if a minimum amount of total bases were required,
                    # but none were found, then set the filter
                    elif (anRTPs["MinDepth"] > 0):
                        isValidMod = False
                        filterSet.add("rtmndp")

            elif (modType.find("TUM_EDIT") != -1 or
                  modType.find("RNA_TUM_VAR") != -1):

                sIndex = refPlusAltList.index(source)
                tIndex = refPlusAltList.index(target)

                # if we are also filtering using the RNA
                if (aFilterUsingRNAFlag):
                    # if this is a tumor edit, then we need to check the DNA
                    # if this is an RNA tumor variant,
                    # then don't filter on the DNA
                    if (modType.find("TUM_EDIT") != -1):
                        # check to make sure the normal DNA sample has data
                        # and is between the min and the max of total bases
                        if (haveDnaNormData):
                            if (int(dnDict["DP"][0]) < aDNPs["MinDepth"]):
                                isValidMod = False
                                filterSet.add("dnmndp")
                            elif (int(dnDict["DP"][0]) >= aDNPs["MaxDepth"]):
                                isValidMod = False
                                filterSet.add("dnmxdp")
                            # we want to make sure the normal DNA sample
                            # error percentage is below the max
                            # we want to make sure that the percentage of
                            # other ALTs in this sample is below the max error
                            if (filterByMaxError(refPlusAltList,
                                                 aDNPs,
                                                 dnDict,
                                                 sIndex,
                                                 tIndex,
                                                 True,
                                                 anIsDebug)):
                                isValidMod = False
                                filterSet.add("dnmxerr")
                        # else if a minimum amount of total bases were
                        # required, but none were found, then set the filter
                        elif (aDNPs["MinDepth"] > 0):
                            isValidMod = False
                            filterSet.add("dnmndp")

                        # check to make sure the tumor DNA sample is
                        # between the min and the max of total bases
                        if (haveDnaTumData):
                            if (int(dtDict["DP"][0]) < aDTPs["MinDepth"]):
                                isValidMod = False
                                filterSet.add("dtmndp")
                            elif (int(dtDict["DP"][0]) >= aDTPs["MaxDepth"]):
                                isValidMod = False
                                filterSet.add("dtmxdp")
                            # we want to make sure the tumor DNA sample
                            # error percentage is below the max
                            # we want to make sure that the percentage of
                            # other ALTs in this sample is below the max error
                            if (filterByMaxError(refPlusAltList,
                                                 aDTPs,
                                                 dtDict,
                                                 sIndex,
                                                 tIndex,
                                                 True,
                                                 anIsDebug)):
                                isValidMod = False
                                filterSet.add("dtmxerr")
                        # else if a minimum amount of total bases were
                        # required, but none were found, then set the filter
                        elif (aDTPs["MinDepth"] > 0):
                            isValidMod = False
                            filterSet.add("dtmndp")

                    # check to make sure the tumor RNA sample is
                    # between the min and the max of total bases
                    if (int(rtDict["DP"][0]) < anRTPs["MinDepth"]):
                        isValidMod = False
                        filterSet.add("rtmndp")
                    elif (int(rtDict["DP"][0]) >= anRTPs["MaxDepth"]):
                        isValidMod = False
                        filterSet.add("rtmxdp")

                    # check to make sure the tumor RNA sample
                    # number of ALT bases is above the min
                    if (int(rtDict["AD"][tIndex]) < anRTPs["MinAltDepth"]):
                        isValidMod = False
                        filterSet.add("rtmnad")

                    # check to make sure the tumor RNA sample
                    # percentage of ALT bases is above the min
                    if (float(rtDict["AF"][tIndex]) < anRTPs["MinAltPct"]):
                        isValidMod = False
                        filterSet.add("rtmnaf")

                    # check to make sure the tumor RNA sample average
                    # base quality for ALT bases is above the min
                    if (int(rtDict["BQ"][tIndex]) < anRTPs["MinAltAvgBQ"]):
                        isValidMod = False
                        filterSet.add("rtmnbq")

                    # check to make sure the tumor RNA sample average
                    # map quality for ALT reads is above the min
                    if (int(rtDict["MQA"][tIndex]) < anRTPs["MinAltMQA"]):
                        isValidMod = False
                        filterSet.add("rtmnmqa")

                    # check to make sure the tumor RNA sample has at least
                    # 1 ALT read with a mapping quality above the min
                    if (int(rtDict["MMQ"][tIndex]) < anRTPs["MinAltMMQ"]):
                        isValidMod = False
                        filterSet.add("rtmnmq")

                    # check to make sure the tumor RNA sample has a maximum
                    # percentage of MQ0 reads supporting the ALT
                    if (filterByMapQualZero(anRTPs,
                                            rtDict,
                                            tIndex)):
                        isValidMod = False
                        filterSet.add("rtmxmq0")

                    # check to make sure the number of tumor RNA sample
                    # INDELs is below the maximum
                    indels = int(rtDict["INS"][0]) + int(rtDict["DEL"][0])
                    if (indels >= anRTPs["MaxIndels"]):
                        isValidMod = False
                        filterSet.add("rtmxindel")

                    # check to make sure the tumor variant reads
                    # don't have a strand bias
                    if (filterByStrandBias(anRTPs,
                                           rtDict,
                                           sIndex,
                                           tIndex)):
                        isValidMod = False
                        filterSet.add("rtsbias")

                    # we want to make sure the tumor RNA sample
                    # error percentage is below the max
                    # we want to make sure that the percentage of
                    # other ALTs in this sample is below the max error
                    if (filterByMaxError(refPlusAltList,
                                         anRTPs,
                                         rtDict,
                                         sIndex,
                                         tIndex,
                                         False,
                                         anIsDebug)):
                        isValidMod = False
                        filterSet.add("rtmxerr")
                # we are filtering via the DNA, so put in a
                # dummy filter so that they don't pass
                else:
                    isValidMod = False
                    filterSet.add("rnacall")

            # filters for all mod types start here

            # check to make sure the number of INDELs
            # across all samples is below the maximum
            indels = int(infoDict["INS"][0]) + int(infoDict["DEL"][0])
            if (indels >= aParamsDict["MaxIndels"]):
                isValidMod = False
                filterSet.add("indel")

            if (anIsDebug):
                logging.debug("modType=%s, modChange=%s, " +
                              "isValidMod=%s, filters=%s",
                              modType, modChange, isValidMod, filterSet)

            # if this one is not valid and we have more,
            # then set the filters for this one, remove it,
            # and try the next one
            if (not isValidMod):
                allFiltersSet = allFiltersSet.union(filterSet)

                # remove it and try the next one
                modIndices = range(0, len(modTypesList))
                for (rmModType, rmModChange, modIndex) in izip(modTypesList,
                                                               modChangesList,
                                                               modIndices):
                    if (modType == rmModType and modChange == rmModChange):
                        del modTypesList[modIndex]
                        del modChangesList[modIndex]
                        break

            # keep track of the mod filters for all calls
            # whether they are valid or not. this is needed
            # to include filters upstream from this script in the final calls
            origin = "DNA"
            if (aFilterUsingRNAFlag):
                origin = "RNA"

            modFilterTypes.append("_".join([origin, modType, modChange]))
            modFilters.append("_".join(event_filterSet.union(filterSet)))

        # after looping through all of them:  if there are
        # still some valid mod types, then set them in the
        # infoDict and ignore the other filtered calls
        if (len(modTypesList) > 0):
            infoDict["MT"] = modTypesList
            infoDict["MC"] = modChangesList

            # if an event passed, get the final mod type
            infoDict = get_final_mod_type(infoDict, anIsDebug)

            # if this call passes this script but did not pass previous scripts
            # then set the mod filters from previous scripts
            if (len(event_filterSet) != 0):
                infoDict["MFT"] = modFilterTypes
                infoDict["MF"] = modFilters
        # otherwise add the appropriate filters
        else:
            infoDict["MFT"] = modFilterTypes
            infoDict["MF"] = modFilters
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
            if ("SOM" in infoDict["MT"]):
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
        for key in sorted(infoDict.iterkeys()):
            if (len(infoDict[key]) == 0):
                continue
            elif ("True" in infoDict[key]):
                infoField += key + ";"
            else:
                infoField += key + "=" + ",".join(infoDict[key]) + ";"

        vcfOutputList.append(infoField.rstrip(";"))

        vcfOutputList.append(":".join(formatList))
        if (dnaNormalList is not None):
            vcfOutputList.append(":".join(dnaNormalList))
        if (rnaNormalList is not None):
            vcfOutputList.append(":".join(rnaNormalList))
        if (dnaTumorList is not None):
            vcfOutputList.append(":".join(dnaTumorList))
        if (rnaTumorList is not None):
            vcfOutputList.append(":".join(rnaTumorList))

        if (i_outputFileHandler is not None):
            i_outputFileHandler.write("\t".join(vcfOutputList) + "\n")
        else:
            print >> sys.stdout, "\t".join(vcfOutputList)

    logging.info("filterByMpileupSupport.py Chrom %s and Id %s: " +
                 "%s events passed out of %s total events",
                 aChrom, anId, includedEvents, totalEvents)
    '''
    logging.info("\t".join([anId, aChrom,
                            str(somEventsPassing),
                            str(somEventsWithTumorRna),
                            str(somEventsWithTumorAltRna)]))
    '''

    # close the files
    if (anOutputFilename is not None):
        i_outputFileHandler.close()

    i_vcfFileHandler.close()
    return


def main():

    # python filterByReadSupport.py TCGA-AB-2995 12
    # ../data/test/TCGA-AB-2995.vcf

    # create the usage statement
    usage = "usage: python %prog id chrom vcfFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)

    # add the optional params
    i_cmdLineParser.add_option(
        "-o", "--outputFilename",
        dest="outputFilename", metavar="OUTPUT_FILE",
        help="the name of the output file, STDOUT by default")
    i_cmdLineParser.add_option(
        "-n", "--headerFilename",
        dest="headerFilename", metavar="HEADER_FILE",
        help="the name of the header file")
    i_cmdLineParser.add_option(
        "-l", "--log",
        dest="logLevel", default="WARNING", metavar="LOG",
        help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), " +
             "%default by default")
    i_cmdLineParser.add_option(
        "-g", "--logFilename",
        dest="logFilename", metavar="LOG_FILE",
        help="the name of the log file, STDOUT by default")
    i_cmdLineParser.add_option(
        "-s", "--statsDir",
        dest="statsDir", metavar="STATS_DIR",
        help="a stats directory where some basic stats can be output")
    i_cmdLineParser.add_option(
        "-r", "--filterUsingRNA",
        dest="filterUsingRNA", action="store_true", default=False,
        help="include this argument if the germline and somatic calls " +
             "should be filtered by the RNA")
    i_cmdLineParser.add_option(
        "-d", "--filterUsingDNA",
        dest="filterUsingDNA", action="store_true", default=False,
        help="include this argument if the germline and somatic calls " +
             "should be filtered by the DNA")
    i_cmdLineParser.add_option(
        "-a", "--addOrigin",
        dest="addOrigin", action="store_true", default=False,
        help="include this argument if the origin of the call should be " +
             "specified in the INFO tags")

    i_cmdLineParser.add_option(
        "", "--genotypeMinDepth",
        type="int", default=int(4),
        dest="genotypeMinDepth", metavar="GT_MIN_DP",
        help="the minimum number of bases required for the genotype, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--genotypeMinPct",
        type="float", default=float(0.10),
        dest="genotypeMinPct", metavar="GT_MIN_PCT",
        help="the minimum percentage of reads required for the genotype, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--modMinDepth",
        type="int", default=int(4),
        dest="modMinDepth", metavar="MOD_MIN_DP",
        help="the minimum number of bases required for a modification, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--modMinPct",
        type="float", default=float(0.10),
        dest="modMinPct", metavar="MOD_MIN_PCT",
        help="the minimum percentage of reads required for a modification, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--lohMaxDepth",
        type="int", default=int(2),
        dest="lohMaxDepth", metavar="LOH_MAX_DP",
        help="the maximum number of bases allowed in the tumor DNA for " +
             "an LOH, %default by default")
    i_cmdLineParser.add_option(
        "", "--lohMaxPct",
        type="float", default=float(0.02),
        dest="lohMaxPct", metavar="LOH_MAX_PCT",
        help="the maximum percentage of reads in the tumor DNA for " +
             "an LOH, %default by default")
    i_cmdLineParser.add_option(
        "", "--maxIndels",
        type="int", default=int(3),
        dest="maxIndels", metavar="MAX_INDELS",
        help="the maximum number of INDELS allowed at this position " +
             "across all samples, %default by default")

    i_cmdLineParser.add_option(
        "", "--dnaNormalMinTotalBases",
        type="int", default=int(8),
        dest="dnaNormalMinTotalNumBases", metavar="DNA_NOR_MIN_TOTAL_BASES",
        help="the minimum number of overall normal DNA reads covering " +
             "a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMaxTotalBases",
        type="int", default=int(8000),
        dest="dnaNormalMaxTotalNumBases", metavar="DNA_NOR_MAX_TOTAL_BASES",
        help="the maximum number of overall normal DNA reads covering " +
             "a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMinAltBases",
        type="int", default=int(4),
        dest="dnaNormalMinAltNumBases", metavar="DNA_NOR_MIN_ALT_BASES",
        help="the minimum number of alternative normal DNA reads " +
             "supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMinAltPct",
        type="float", default=float(0.03),
        dest="dnaNormalMinAltPct", metavar="DNA_NOR_MIN_ALT_PCT",
        help="the minimum percentage of alternative normal DNA reads " +
             "supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMaxErrPct",
        type="float", default=float(0.01),
        dest="dnaNormalMaxErrPct", metavar="DNA_NOR_MAX_ERR_PCT",
        help="the maximum percentage of alternative normal DNA reads " +
             "allowed that support a variant, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMinErrPctDepth",
        type="int", default=float(2),
        dest="dnaNormalMinErrPctDepth", metavar="DNA_NOR_MIN_ERR_PCT_DEPTH",
        help="the minimum error count depth needed for the max error " +
             "percent filter to be applied, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMaxStrandBias",
        type="float", default=float(0.99),
        dest="dnaNormalMaxStrandBias", metavar="DNA_NOR_MAX_STRAND_BIAS",
        help="the maximum percentage of strand bias on reads that " +
             "support the ALT, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMinStrandBiasDepth",
        type="int", default=float(4),
        dest="dnaNormalMinStrandBiasDepth",
        metavar="DNA_NOR_MIN_STRAND_BIAS_DP",
        help="the minimum total depth needed for the strand bias filter " +
             "to be applied, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMinAltAvgBaseQual",
        type="int", default=int(20),
        dest="dnaNormalMinAltAvgBaseQual", metavar="DNA_NOR_MIN_ALT_AVG_BQ",
        help="the minimum average base quality for the ALT reads, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMinAltAvgMapQual",
        type="int", default=int(15),
        dest="dnaNormalMinAltAvgMapQual", metavar="DNA_NOR_MIN_ALT_AVG_MQ",
        help="the minimum average mapping quality for the ALT reads, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMinAltMapQual",
        type="int", default=int(20),
        dest="dnaNormalMinAltMapQual", metavar="DNA_NOR_MIN_ALT_MQ",
        help="at least 1 ALT read needs this minimum mapping quality, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMaxAltMapQualZeroPct",
        type="float", default=float(0.50),
        dest="dnaNormalMaxAltMapQualZeroPct",
        metavar="DNA_NOR_MAX_ALT_MQ0_PCT",
        help="the maximum percentage of mapping quality zero reads " +
             "for the ALT reads, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaNormalMaxIndels",
        type="int", default=int(3),
        dest="dnaNormalMaxIndels", metavar="DNA_NOR_MAX_INDELS",
        help="the maximum number of INDELS allowed at a position, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--danNormalPurity",
        type="float", default=float(1.0),
        dest="dnaNormalPurity", metavar="DNA_NOR_PURITY",
        help="estimated purity (non-tumor content) of normal DNA sample, " +
             "%default by default")

    i_cmdLineParser.add_option(
        "", "--dnaTumorMinTotalBases",
        type="int", default=int(8),
        dest="dnaTumorMinTotalNumBases", metavar="DNA_TUM_MIN_TOTAL_BASES",
        help="the minimum number of overall tumor DNA reads covering " +
             "a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMaxTotalBases",
        type="int", default=int(8000),
        dest="dnaTumorMaxTotalNumBases", metavar="DNA_TUM_MAX_TOTAL_BASES",
        help="the maximum number of overall tumor DNA reads covering " +
             "a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMinAltBases",
        type="int", default=int(4),
        dest="dnaTumorMinAltNumBases", metavar="DNA_TUM_MIN_ALT_BASES",
        help="the minimum number of alternative tumor DNA reads " +
             "supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMinAltPct",
        type="float", default=float(0.10),
        dest="dnaTumorMinAltPct", metavar="DNA_TUM_MIN_ALT_PCT",
        help="the minimum percentage of alternative tumor DNA reads " +
             "supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMaxErrPct",
        type="float", default=float(0.01),
        dest="dnaTumorMaxErrPct", metavar="DNA_TUM_MAX_ERR_PCT",
        help="the maximum percentage of alternative tumor DNA reads allowed " +
             "that support a variant in the tumor RNA, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMinErrPctDepth",
        type="int", default=float(2),
        dest="dnaTumorMinErrPctDepth", metavar="DNA_TUM_MIN_ERR_PCT_DEPTH",
        help="the minimum error count depth needed for the max error " +
             "percent filter to be applied, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMaxStrandBias",
        type="float", default=float(0.99),
        dest="dnaTumorMaxStrandBias", metavar="DNA_TUM_MAX_STRAND_BIAS",
        help="the maximum percentage of strand bias on reads that support " +
             "the ALT, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMinStrandBiasDepth",
        type="int", default=float(4),
        dest="dnaTumorMinStrandBiasDepth",
        metavar="DNA_TUM_MIN_STRAND_BIAS_DP",
        help="the minimum total depth needed for the strand bias filter " +
             "to be applied, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMinAltAvgBaseQual",
        type="int", default=int(20),
        dest="dnaTumorMinAltAvgBaseQual", metavar="DNA_TUM_MIN_ALT_AVG_BQ",
        help="the minimum average base quality for the ALT reads, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMinAltAvgMapQual",
        type="int", default=int(15),
        dest="dnaTumorMinAltAvgMapQual", metavar="DNA_TUM_MIN_ALT_AVG_MQ",
        help="the minimum average mapping quality for the ALT reads, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMinAltMapQual",
        type="int", default=int(20),
        dest="dnaTumorMinAltMapQual", metavar="DNA_TUM_MIN_ALT_MQ",
        help="at least 1 ALT read needs this minimum mapping quality, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMaxAltMapQualZeroPct",
        type="float", default=float(0.50),
        dest="dnaTumorMaxAltMapQualZeroPct", metavar="DNA_TUM_MAX_ALT_MQ0_PCT",
        help="the maximum percentage of mapping quality zero reads " +
             "for the ALT reads, %default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorMaxIndels",
        type="int", default=int(3),
        dest="dnaTumorMaxIndels", metavar="DNA_TUM_MAX_INDELS",
        help="the maximum number of INDELS allowed at a position, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--dnaTumorPurity",
        type="float", default=float(1.0),
        dest="dnaTumorPurity", metavar="DNA_TUM_PURITY",
        help="estimated purity (tumor content) of tumor DNA sample, " +
             "%default by default")

    i_cmdLineParser.add_option(
        "", "--rnaNormalMinTotalBases",
        type="int", default=int(8),
        dest="rnaNormalMinTotalNumBases", metavar="RNA_NOR_MIN_TOTAL_BASES",
        help="the minimum number of overall normal RNA-Seq reads covering " +
             "a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMaxTotalBases",
        type="int", default=int(8000),
        dest="rnaNormalMaxTotalNumBases", metavar="RNA_NOR_MAX_TOTAL_BASES",
        help="the maximum number of overall normal RNA-Seq reads covering " +
             "a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMinAltBases",
        type="int", default=int(4),
        dest="rnaNormalMinAltNumBases", metavar="RNA_NOR_MIN_ALT_BASES",
        help="the minimum number of alternative normal RNA-Seq reads " +
             "supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMinAltPct",
        type="float", default=float(0.10),
        dest="rnaNormalMinAltPct", metavar="RNA_NOR_MIN_ALT_PCT",
        help="the minimum percentage of alternative normal RNA-Seq reads " +
             "supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMaxErrPct",
        type="float", default=float(0.01),
        dest="rnaNormalMaxErrPct", metavar="RNA_NOR_MAX_ERR_PCT",
        help="the maximum percentage of alternative normal RNA reads " +
             "allowed that support a variant in the tumor, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMinErrPctDepth",
        type="int", default=float(2),
        dest="rnaNormalMinErrPctDepth", metavar="RNA_NOR_MIN_ERR_PCT_DEPTH",
        help="the minimum error count depth needed for the max error " +
             "percent filter to be applied, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMaxStrandBias",
        type="float", default=float(0.99),
        dest="rnaNormalMaxStrandBias", metavar="RNA_NOR_MAX_STRAND_BIAS",
        help="the maximum percentage of strand bias on reads that support " +
             "the ALT, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMinStrandBiasDepth",
        type="int", default=float(4),
        dest="rnaNormalMinStrandBiasDepth",
        metavar="RNA_NOR_MIN_STRAND_BIAS_DP",
        help="the minimum total depth needed for the strand bias filter " +
             "to be applied, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMinAltAvgBaseQual",
        type="int", default=int(20),
        dest="rnaNormalMinAltAvgBaseQual", metavar="RNA_NOR_MIN_ALT_AVG_BQ",
        help="the minimum average base quality for the ALT reads, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMinAltAvgMapQual",
        type="int", default=int(15),
        dest="rnaNormalMinAltAvgMapQual", metavar="RNA_NOR_MIN_ALT_AVG_MQ",
        help="the minimum average mapping quality for the ALT reads, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMinAltMapQual",
        type="int", default=int(20),
        dest="rnaNormalMinAltMapQual", metavar="RNA_NOR_MIN_ALT_MQ",
        help="at least 1 ALT read needs this minimum mapping quality, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMaxAltMapQualZeroPct",
        type="float", default=float(0.50),
        dest="rnaNormalMaxAltMapQualZeroPct",
        metavar="RNA_NOR_MAX_ALT_MQ0_PCT",
        help="the maximum percentage of mapping quality zero reads for the " +
             "ALT reads, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaNormalMaxIndels",
        type="int", default=int(3),
        dest="rnaNormalMaxIndels", metavar="RNA_NOR_MAX_INDELS",
        help="the maximum number of INDELS allowed at a position, " +
             "%default by default")

    i_cmdLineParser.add_option(
        "", "--rnaTumorMinTotalBases",
        type="int", default=int(8),
        dest="rnaTumorMinTotalNumBases", metavar="RNA_TUM_MIN_TOTAL_BASES",
        help="the minimum number of overall tumor RNA-Seq reads covering " +
             "a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMaxTotalBases",
        type="int", default=int(8000),
        dest="rnaTumorMaxTotalNumBases", metavar="RNA_TUM_MAX_TOTAL_BASES",
        help="the maximum number of overall tumor RNA-Seq reads covering " +
             "a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMinAltBases",
        type="int", default=int(4),
        dest="rnaTumorMinAltNumBases", metavar="RNA_TUM_MIN_ALT_BASES",
        help="the minimum number of alternative tumor RNA-Seq reads " +
             "supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMinAltPct",
        type="float", default=float(0.10),
        dest="rnaTumorMinAltPct", metavar="RNA_TUM_MIN_ALT_PCT",
        help="the minimum percentage of alternative tumor RNA-Seq reads " +
             "supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMaxErrPct",
        type="float", default=float(0.01),
        dest="rnaTumorMaxErrPct", metavar="RNA_TUM_MAX_ERR_PCT",
        help="the maximum percentage of alternative tumor RNA reads " +
             "allowed that support a variant in a future data type, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMinErrPctDepth",
        type="int", default=float(2),
        dest="rnaTumorMinErrPctDepth", metavar="RNA_TUM_MIN_ERR_PCT_DEPTH",
        help="the minimum error count depth needed for the max error " +
             "percent filter to be applied, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMaxStrandBias",
        type="float", default=float(0.99),
        dest="rnaTumorMaxStrandBias", metavar="RNA_TUM_MAX_STRAND_BIAS",
        help="the maximum percentage of strand bias on reads that support " +
             "the ALT, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMinStrandBiasDepth",
        type="int", default=float(4),
        dest="rnaTumorMinStrandBiasDepth",
        metavar="RNA_TUM_MIN_STRAND_BIAS_DP",
        help="the minimum total depth needed for the strand bias filter " +
             "to be applied, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMinAltAvgBaseQual",
        type="int", default=int(20),
        dest="rnaTumorMinAltAvgBaseQual", metavar="RNA_TUM_MIN_ALT_AVG_BQ",
        help="the minimum average base quality for the ALT reads, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMinAltAvgMapQual",
        type="int", default=int(15),
        dest="rnaTumorMinAltAvgMapQual", metavar="RNA_TUM_MIN_ALT_AVG_MQ",
        help="the minimum average mapping quality for the ALT reads, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMinAltMapQual",
        type="int", default=int(20),
        dest="rnaTumorMinAltMapQual", metavar="RNA_TUM_MIN_ALT_MQ",
        help="at least 1 ALT read needs this minimum mapping quality, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMaxAltMapQualZeroPct",
        type="float", default=float(0.50),
        dest="rnaTumorMaxAltMapQualZeroPct", metavar="RNA_TUM_MAX_ALT_MQ0_PCT",
        help="the maximum percentage of mapping quality zero reads " +
             "for the ALT reads, %default by default")
    i_cmdLineParser.add_option(
        "", "--rnaTumorMaxIndels",
        type="int", default=int(3),
        dest="rnaTumorMaxIndels", metavar="RNA_TUM_MAX_INDELS",
        help="the maximum number of INDELS allowed at a position, " +
             "%default by default")

    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(3, 58, 1)
    i_argLength = len(sys.argv)

    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)

    # get the required parameters
    (cmdLineOpts, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    cmdLineOptsDict = vars(cmdLineOpts)
    i_id = str(i_cmdLineArgs[0])
    i_chrom = str(i_cmdLineArgs[1])
    i_vcfFilename = str(i_cmdLineArgs[2])

    # get the optional params with default values
    i_logLevel = cmdLineOpts.logLevel
    i_addOrigin = cmdLineOpts.addOrigin
    i_filterUsingRNA = cmdLineOpts.filterUsingRNA
    i_filterUsingDNA = cmdLineOpts.filterUsingDNA

    i_paramsDict = {}
    i_paramsDict["MinGenotypeDepth"] = cmdLineOpts.genotypeMinDepth
    i_paramsDict["MinGenotypePct"] = cmdLineOpts.genotypeMinPct
    i_paramsDict["MinModDepth"] = cmdLineOpts.modMinDepth
    i_paramsDict["MinModPct"] = cmdLineOpts.modMinPct
    # i_paramsDict["MaxModDepth"] = cmdLineOpts.modMaxDepth
    # i_paramsDict["MaxModPct"] = cmdLineOpts.modMaxPct
    i_paramsDict["MaxLohDepth"] = cmdLineOpts.lohMaxDepth
    i_paramsDict["MaxLohPct"] = cmdLineOpts.lohMaxPct
    i_paramsDict["MaxIndels"] = cmdLineOpts.maxIndels

    dnParams = {}
    dnParams["MinDepth"] = cmdLineOpts.dnaNormalMinTotalNumBases
    dnParams["MaxDepth"] = cmdLineOpts.dnaNormalMaxTotalNumBases
    dnParams["MinAltDepth"] = cmdLineOpts.dnaNormalMinAltNumBases
    dnParams["MinAltPct"] = cmdLineOpts.dnaNormalMinAltPct
    dnParams["MaxErrPct"] = cmdLineOpts.dnaNormalMaxErrPct
    dnParams["MinErrPctDP"] = cmdLineOpts.dnaNormalMinErrPctDepth
    dnParams["MaxStrandBias"] = cmdLineOpts.dnaNormalMaxStrandBias
    dnParams["MinStrBiasDP"] = cmdLineOpts.dnaNormalMinStrandBiasDepth
    dnParams["MinAltAvgBQ"] = cmdLineOpts.dnaNormalMinAltAvgBaseQual
    dnParams["MinAltMQA"] = cmdLineOpts.dnaNormalMinAltAvgMapQual
    dnParams["MinAltMMQ"] = cmdLineOpts.dnaNormalMinAltMapQual
    dnParams["MaxAltMQ0Pct"] = cmdLineOpts.dnaNormalMaxAltMapQualZeroPct
    dnParams["MaxIndels"] = cmdLineOpts.dnaNormalMaxIndels
    dnParams["Purity"] = cmdLineOpts.dnaNormalPurity

    dtParams = {}
    dtParams["MinDepth"] = cmdLineOpts.dnaTumorMinTotalNumBases
    dtParams["MaxDepth"] = cmdLineOpts.dnaTumorMaxTotalNumBases
    dtParams["MinAltDepth"] = cmdLineOpts.dnaTumorMinAltNumBases
    dtParams["MinAltPct"] = cmdLineOpts.dnaTumorMinAltPct
    dtParams["MaxErrPct"] = cmdLineOpts.dnaTumorMaxErrPct
    dtParams["MinErrPctDP"] = cmdLineOpts.dnaTumorMinErrPctDepth
    dtParams["MaxStrandBias"] = cmdLineOpts.dnaTumorMaxStrandBias
    dtParams["MinStrBiasDP"] = cmdLineOpts.dnaTumorMinStrandBiasDepth
    dtParams["MinAltAvgBQ"] = cmdLineOpts.dnaTumorMinAltAvgBaseQual
    dtParams["MinAltMQA"] = cmdLineOpts.dnaTumorMinAltAvgMapQual
    dtParams["MinAltMMQ"] = cmdLineOpts.dnaTumorMinAltMapQual
    dtParams["MaxAltMQ0Pct"] = cmdLineOpts.dnaTumorMaxAltMapQualZeroPct
    dtParams["MaxIndels"] = cmdLineOpts.dnaTumorMaxIndels
    dtParams["Purity"] = cmdLineOpts.dnaTumorPurity

    rnParams = {}
    rnParams["MinDepth"] = cmdLineOpts.rnaNormalMinTotalNumBases
    rnParams["MaxDepth"] = cmdLineOpts.rnaNormalMaxTotalNumBases
    rnParams["MinAltDepth"] = cmdLineOpts.rnaNormalMinAltNumBases
    rnParams["MinAltPct"] = cmdLineOpts.rnaNormalMinAltPct
    rnParams["MaxErrPct"] = cmdLineOpts.rnaNormalMaxErrPct
    rnParams["MinErrPctDP"] = cmdLineOpts.rnaNormalMinErrPctDepth
    rnParams["MaxStrandBias"] = cmdLineOpts.rnaNormalMaxStrandBias
    rnParams["MinStrBiasDP"] = cmdLineOpts.rnaNormalMinStrandBiasDepth
    rnParams["MinAltAvgBQ"] = cmdLineOpts.rnaNormalMinAltAvgBaseQual
    rnParams["MinAltMQA"] = cmdLineOpts.rnaNormalMinAltAvgMapQual
    rnParams["MinAltMMQ"] = cmdLineOpts.rnaNormalMinAltMapQual
    rnParams["MaxAltMQ0Pct"] = cmdLineOpts.rnaNormalMaxAltMapQualZeroPct
    rnParams["MaxIndels"] = cmdLineOpts.rnaNormalMaxIndels

    rtParams = {}
    rtParams["MinDepth"] = cmdLineOpts.rnaTumorMinTotalNumBases
    rtParams["MaxDepth"] = cmdLineOpts.rnaTumorMaxTotalNumBases
    rtParams["MinAltDepth"] = cmdLineOpts.rnaTumorMinAltNumBases
    rtParams["MinAltPct"] = cmdLineOpts.rnaTumorMinAltPct
    rtParams["MaxErrPct"] = cmdLineOpts.rnaTumorMaxErrPct
    rtParams["MinErrPctDP"] = cmdLineOpts.rnaTumorMinErrPctDepth
    rtParams["MaxStrandBias"] = cmdLineOpts.rnaTumorMaxStrandBias
    rtParams["MinStrBiasDP"] = cmdLineOpts.rnaTumorMinStrandBiasDepth
    rtParams["MinAltAvgBQ"] = cmdLineOpts.rnaTumorMinAltAvgBaseQual
    rtParams["MinAltMQA"] = cmdLineOpts.rnaTumorMinAltAvgMapQual
    rtParams["MinAltMMQ"] = cmdLineOpts.rnaTumorMinAltMapQual
    rtParams["MaxAltMQ0Pct"] = cmdLineOpts.rnaTumorMaxAltMapQualZeroPct
    rtParams["MaxIndels"] = cmdLineOpts.rnaTumorMaxIndels

    # try to get any optional parameters with no defaults
    i_readFilenameList = [i_vcfFilename]
    i_writeFilenameList = []
    i_dirList = []

    i_outputFilename = None
    i_headerFilename = None
    i_logFilename = None
    i_statsDir = None
    if (cmdLineOpts.outputFilename is not None):
        i_outputFilename = str(cmdLineOpts.outputFilename)
        i_writeFilenameList += [i_outputFilename]
    if (cmdLineOpts.logFilename is not None):
        i_logFilename = str(cmdLineOpts.logFilename)
        i_writeFilenameList += [i_logFilename]
    if (cmdLineOpts.headerFilename is not None):
        i_headerFilename = str(cmdLineOpts.headerFilename)
        i_writeFilenameList += [i_headerFilename]
    if (cmdLineOpts.statsDir is not None):
        i_statsDir = str(cmdLineOpts.statsDir)
        i_dirList += [i_statsDir]

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

        logging.debug("genotypeMinDepth=%s" % i_paramsDict["MinGenotypeDepth"])
        logging.debug("genotypeMinPct=%s" % i_paramsDict["MinGenotypePct"])
        logging.debug("modMinDepth=%s" % i_paramsDict["MinModDepth"])
        logging.debug("modMinPct=%s" % i_paramsDict["MinModPct"])
        logging.debug("lohMaxDepth=%s" % i_paramsDict["MaxLohDepth"])
        logging.debug("lohMaxPct=%s" % i_paramsDict["MaxLohPct"])
        logging.debug("maxIndels=%s" % i_paramsDict["MaxIndels"])

        logging.debug("dn minDepth: %s" % dnParams["MinDepth"])
        logging.debug("dn maxDepth: %s" % dnParams["MaxDepth"])
        logging.debug("dn minAltDepth: %s" % dnParams["MinAltDepth"])
        logging.debug("dn minAltPct: %s" % dnParams["MinAltPct"])
        logging.debug("dn maxErrPct: %s" % dnParams["MaxErrPct"])
        logging.debug("dn minErrPctDP: %s" % dnParams["MinErrPctDP"])
        logging.debug("dn strandBias: %s" % dnParams["MaxStrandBias"])
        logging.debug("dn strandBiasMinDP: %s" % dnParams["MinStrBiasDP"])
        logging.debug("dn minAltAvgBQ: %s" % dnParams["MinAltAvgBQ"])
        logging.debug("dn minAltAvgMapQual: %s" % dnParams["MinAltMQA"])
        logging.debug("dn minAltMapQual: %s" % dnParams["MinAltMMQ"])
        logging.debug("dn maxAltMQ0Pct: %s" % dnParams["MaxAltMQ0Pct"])
        logging.debug("dn maxIndels: %s" % dnParams["MaxIndels"])
        logging.debug("dn purity=%s" % dnParams["Purity"])

        logging.debug("dt minDepth: %s" % dtParams["MinDepth"])
        logging.debug("dt maxDepth: %s" % dtParams["MaxDepth"])
        logging.debug("dt minAltDepth: %s" % dtParams["MinAltDepth"])
        logging.debug("dt minAltPct: %s" % dtParams["MinAltPct"])
        logging.debug("dt maxErrPct: %s" % dtParams["MaxErrPct"])
        logging.debug("dt minErrPctDP: %s" % dtParams["MinErrPctDP"])
        logging.debug("dt strandBias: %s" % dtParams["MaxStrandBias"])
        logging.debug("dt strandBiasMinDP: %s" % dtParams["MinStrBiasDP"])
        logging.debug("dt minAltAvgBQ: %s" % dtParams["MinAltAvgBQ"])
        logging.debug("dt minAltAvgMapQual: %s" % dtParams["MinAltMQA"])
        logging.debug("dt minAltMapQual: %s" % dtParams["MinAltMMQ"])
        logging.debug("dt maxAltMQ0Pct: %s" % dtParams["MaxAltMQ0Pct"])
        logging.debug("dt maxIndels: %s" % dtParams["MaxIndels"])
        logging.debug("dt purity=%s" % dtParams["Purity"])

        logging.debug("rn minDepth: %s" % rnParams["MinDepth"])
        logging.debug("rn maxDepth: %s" % rnParams["MaxDepth"])
        logging.debug("rn minAltDepth: %s" % rnParams["MinAltDepth"])
        logging.debug("rn minAltPct: %s" % rnParams["MinAltPct"])
        logging.debug("rn maxErrPct: %s" % rnParams["MaxErrPct"])
        logging.debug("rn minErrPctDP: %s" % rnParams["MinErrPctDP"])
        logging.debug("rn strandBias: %s" % rnParams["MaxStrandBias"])
        logging.debug("rn strandBiasMinDP: %s" % rnParams["MinStrBiasDP"])
        logging.debug("rn minAltAvgBQ: %s" % rnParams["MinAltAvgBQ"])
        logging.debug("rn minAltAvgMapQual: %s" % rnParams["MinAltMQA"])
        logging.debug("rn minAltMapQual: %s" % rnParams["MinAltMMQ"])
        logging.debug("rn maxAltMQ0Pct: %s" % rnParams["MaxAltMQ0Pct"])
        logging.debug("rn maxIndels: %s" % rnParams["MaxIndels"])

        logging.debug("rt minDepth: %s" % rtParams["MinDepth"])
        logging.debug("rt maxDepth: %s" % rtParams["MaxDepth"])
        logging.debug("rt minAltDepth: %s" % rtParams["MinAltDepth"])
        logging.debug("rt minAltPct: %s" % rtParams["MinAltPct"])
        logging.debug("rt maxErrPct: %s" % rtParams["MaxErrPct"])
        logging.debug("rt minErrPctDP: %s" % rtParams["MinErrPctDP"])
        logging.debug("rt strandBias: %s" % rtParams["MaxStrandBias"])
        logging.debug("rt strandBiasMinDP: %s" % rtParams["MinStrBiasDP"])
        logging.debug("rt minAltAvgBQ: %s" % rtParams["MinAltAvgBQ"])
        logging.debug("rt minAltAvgMapQual: %s" % rtParams["MinAltMQA"])
        logging.debug("rt minAltMapQual: %s" % rtParams["MinAltMMQ"])
        logging.debug("rt maxAltMQ0Pct: %s" % rtParams["MaxAltMQ0Pct"])
        logging.debug("rt maxIndels: %s" % rtParams["MaxIndels"])

    # check for any errors
    if (not radiaUtil.check_for_argv_errors(i_dirList,
                                            i_readFilenameList,
                                            i_writeFilenameList)):
        sys.exit(1)

    filter_by_mpileup_support(i_id, i_chrom, i_vcfFilename, i_headerFilename,
                              i_outputFilename, i_filterUsingRNA, i_addOrigin,
                              cmdLineOptsDict, dnParams, dtParams, rnParams,
                              rtParams, i_paramsDict, i_debug)
    return


main()
sys.exit(0)
