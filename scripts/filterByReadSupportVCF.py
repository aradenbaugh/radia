#!/usr/bin/env python2.7

from optparse import OptionParser   # used for parsing command line arguments
import rnaEditingUtil               # utility functions for rna editing
import sys                          # system module
import re
import collections
import logging
from itertools import izip

'''
'   Amie Radenbaugh - 02/17/2012
'   UCSC - RNA Editing  
'   Program name: "filterByReadSupportVCF.py"
'''


def fix_genotypes(aChrom, aGenotypeList, anAlleleDepthsList, aGTMinDepth, aGTMinPct):
    
    singleGenotypeChroms = ["chrX", "chrY", "chrM", "chrMT", "X", "Y", "M", "MT"]
    
    # get the total depth
    totalDepth = sum(anAlleleDepthsList)
    
    # make a copy of the list to manipulate
    depthsTmpList = list(anAlleleDepthsList)
            
    if (aChrom in singleGenotypeChroms):
        # the single genotypes are ok, b/c they are always the max
        return aGenotypeList
    elif (len(depthsTmpList) == 0):
        # there are no alleles
        return aGenotypeList
        
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
        max2Pct = int(max2Depth/float(totalDepth) * 100)
        
        # if the max depth is large enough
        if (max2Depth >= aGTMinDepth and max2Pct >= aGTMinPct):
            # find the index for the max depth on the original list
            max2DepthIndex = anAlleleDepthsList.index(max2Depth)
        else:
            # it wasn't large enough, so just use previous max
            max2DepthIndex = max1DepthIndex
    else:
        # otherwise it's the same as the first
        max2DepthIndex = max1DepthIndex
       
    # return the genotypes
    return sorted([max1DepthIndex, max2DepthIndex])


def pre_filter_mod_types(aRefPlusAltList, anInfoDict, aDNANormalDepths, anRNANormalDepths, aDNATumorDepths, anRNATumorDepths, aModMinDepth, aModMinPct, anLohMaxDepth, anLohMaxPct):

    aModTypeList = anInfoDict["MT"]
    aModChangeList = anInfoDict["MC"]
    
    # if we only have one mod type, then just return it
    #if (len(aModTypeList) == 1):
    #    return (anInfoDict, set())
    
    # make copies of the lists to manipulate
    modTypesList = list(aModTypeList)
    modChangesList = list(aModChangeList)
    
    filterSet = set()
    
    try:
        # for each modification type and change
        for (modType, modChange) in izip(aModTypeList, aModChangeList):
            # get the source and target alleles
            (source, target) = modChange.split(">")              
            
            # for every germline or somatic call
            if (modType == "GERM" or modType == "SOM"):
                # get the indices
                sourceIndex = aRefPlusAltList.index(source)
                
                # get the DNA normal total depth
                totalDepth = sum(aDNANormalDepths)
                
                # if we don't have any normal DNA, then move on
                if (totalDepth == 0):
                    modTypesList.remove(modType)
                    modChangesList.remove(modChange)
                    filterSet.add("dnmntb")
                else:
                    # get the depth
                    sourceDepth = aDNANormalDepths[sourceIndex]
                    sourcePct = int(sourceDepth/float(totalDepth) * 100)
                     
                    # if the depth doesn't reach minimum, then remove it
                    if (sourceDepth < aModMinDepth or sourcePct < aModMinPct):
                        modTypesList.remove(modType)
                        modChangesList.remove(modChange)
                        
                        if (sourceDepth < aModMinDepth):
                            filterSet.add("dnmnrb")
                        if (sourcePct < aModMinPct):
                            filterSet.add("dnmnrp")
    
            elif (modType == "LOH"):
                
                # get the total depths
                totalNormalDepths = sum(aDNANormalDepths)
                totalTumorDepths = sum(aDNATumorDepths)
                
                validSources = []
                
                # for each of the possible sources
                for sourceNormal in source:
                    # get the source in the normal and tumor
                    sourceIndex = aRefPlusAltList.index(sourceNormal)
                    sourceNormalDepth = aDNANormalDepths[sourceIndex]
                    sourceTumorDepth = aDNATumorDepths[sourceIndex]
                    sourceNormalPct = int(sourceNormalDepth/float(totalNormalDepths) * 100)
                    sourceTumorPct = int(sourceTumorDepth/float(totalTumorDepths) * 100)
                    
                    # if the normal depth is above the minimum, and the tumor depth is below the minimum, then it's valid
                    if (sourceNormalDepth >= aModMinDepth and sourceNormalPct >= aModMinPct and 
                        sourceTumorDepth <= anLohMaxDepth and sourceTumorPct <= anLohMaxPct):
                        validSources.append(sourceNormal)
                    else:
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
            
    except:
        print "Error:", aModTypeList, aModChangeList, aRefPlusAltList, aDNANormalDepths, aDNATumorDepths
            
    # if there are still some valid mod types, then return them
    if (len(modTypesList) > 0):
        anInfoDict["MT"] = modTypesList
        anInfoDict["MC"] = modChangesList
    # otherwise just return the originals, they will get filtered later anyway
    else:
        #print >> sys.stderr, "no more valid modTypes", aRefPlusAltList, aDNANormalDepths, aDNATumorDepths, anInfoDict
        return (anInfoDict, filterSet)
        
    #if (len(modTypesList) > 1):
    #    print >> sys.stderr, "still more than one modType", aRefPlusAltList, aDNANormalDepths, aDNATumorDepths, anInfoDict
        
    finalModChange = ""
    if len(modTypesList) == 1:
        finalModChange = modTypesList[0]
    # you can have some strange cases where a call is both somatic and an LOH, 
    # but the TCGA VCF format only allows one call, so let's choose SOM and
    # leave the LOH details in the info tag
    elif "SOM" in modTypesList and "LOH" in modTypesList:
        finalModChange = "SOM"
    else:
        finalModChange = modTypesList[0]
        
    # get the final call
    anInfoDict["SS"] = []
    anInfoDict["SOMATIC"] = []
    if (finalModChange == "GERM"):
        anInfoDict["SS"].append("1")
    elif (finalModChange == "SOM"):
        anInfoDict["SS"].append("2")
        anInfoDict["SOMATIC"].append("True")
    elif (finalModChange.find("EDIT") != -1):
        anInfoDict["SS"].append("5")
    # other
    else:
        anInfoDict["SS"].append("4")
                                      
    return (anInfoDict, set())


def extract_read_support(aTCGAId, aChrom, aVCFFilename, anOutputFilename, aFilterUsingRNAFlag, aMinAltAvgBaseQual, aStatsDir, aCmdLineParams, aDnaNormParamsDict, aDnaTumParamsDict, anRnaNormParamsDict, anRnaTumParamsDict, aGTMinDepth, aGTMinPct, aModMinDepth, aModMinPct, anLohMaxDepth, anLohMaxPct, anIsDebug):
    '''
    ' This function...
    '
    ' aTCGAId: The TCGA Id for this sample
    ' aChrom: The chromosome being filtered
    '
    '
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    # get the files
    i_vcfFileHandler = open(aVCFFilename, "r")
        
    i_outputFileHandler = None
    if (anOutputFilename):
        i_outputFileHandler = open(anOutputFilename, "w")
    
    # initialize some variables
    totalEvents = 0
    includedEvents = 0
    hasAddedHeader = False
    headerLines = "##FILTER=<ID=multi,Description=\"There are multiple ALT alleles across all samples at this position\">\n"
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
    
    headerLines += "##FILTER=<ID=dnmnabq,Description=\"DNA Normal minimum average ALT base quality is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=dtmnabq,Description=\"DNA Tumor minimum average ALT base quality is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rnmnabq,Description=\"RNA Normal minimum average ALT base quality is less than user-specified cut-off\">\n"
    headerLines += "##FILTER=<ID=rtmnabq,Description=\"RNA Tumor minimum average ALT base quality is less than user-specified cut-off\">\n"
    
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
    
    
    sampleList = ["DNA_NORMAL", "DNA_TUMOR", "RNA_TUMOR"]
    #columnsList = ["DNA_NORMAL", "DNA_TUMOR"]
    
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
            columnsLine = line.lstrip("#")
            columnsList = columnsLine.split("\t")
            columnsList = columnsList[9:len(columnsList)]
            
            if (i_outputFileHandler != None):
                i_outputFileHandler.write(line + "\n")
            else:
                print >> sys.stdout, line
        
            '''       
            # unindent this 
            # if we find the radia line, then update the params from the user
            elif ("SAMPLE" in line):
                sampleLine = line.rstrip(">")
                sampleLine = sampleLine.lstrip("##SAMPLE=<")
                sampleParams = sampleLine.split(",")
                for param in sampleParams:
                    if ("ID" in param):
                        (id, idName) = param.split("=")
                        sampleList.append(idName)
                        break;
            '''
        
        # if we find the vcfGenerator line, then update the params from the user
        elif ("vcfGenerator" in line):
            #generatorLine = line.rstrip(">")
            #generatorLine = generatorLine.lstrip("##vcfGenerator=<")
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
                if (paramName.startswith("dnaNormal") and "DNA_NORMAL" not in sampleList):
                    continue;
                elif (paramName.startswith("rnaNormal") and "RNA_NORMAL" not in sampleList):
                    continue;
                elif (paramName.startswith("dnaTumor") and "DNA_TUMOR" not in sampleList):
                    continue;
                elif (paramName.startswith("rnaTumor") and "RNA_TUMOR" not in sampleList):
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
            
            # split the line on the tab
            splitLine = re.split("\t", line)
    
            # sample VCF line
            # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31.4;DP=53;FA=0.04;INDEL=0;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP      
            # GT:DP:INDEL:START:STOP:AD:AF:BQ:SB      0/0:2:0:0:0:2:1.0:35.5:0.0      0/0:1:0:0:0:1:1.0:39.0:1.0      0/1:50:0:2:0:48,2:0.96,0.04:31.65,17.5:0.75,0.5
            
            # the coordinate is the second element
            event_chr = splitLine[0]
            event_stopCoordinate = int(splitLine[1])
            event_idList = re.split(";", splitLine[2])
            event_refList = re.split(",", splitLine[3])
            event_altList = re.split(",", splitLine[4])
            event_score = int(splitLine[5])
            
            # if there are no filters so far, then clear the list
            event_filterSet = set(re.split(";", splitLine[6]))
            if (len(event_filterSet) == 1 and "PASS" in event_filterSet):
                event_filterSet = set()
            
            # parse the info column and create a dict
            event_infoList = re.split(";", splitLine[7])
            event_infoDict = collections.defaultdict(list)
            for info in event_infoList:
                keyValueList = re.split("=", info)
                # some keys are just singular without a value (e.g. DB, SOMATIC, etc.)
                if (len(keyValueList) == 1):
                    event_infoDict[keyValueList[0]] = ["True"]
                else:
                    # the value can be a comma separated list
                    event_infoDict[keyValueList[0]] = re.split(",", keyValueList[1])

            # get the event format list
            event_formatList = re.split(":", splitLine[8])
            
            # initialize the optional columns to none
            event_dnaNormalList = None
            event_dnaTumorList = None
            event_rnaNormalList = None
            event_rnaTumorList = None

            # if we have a 9th column, figure out which dataset it is
            if (len(splitLine) > 9):
                if (columnsList[0] == "DNA_NORMAL"):
                    event_dnaNormalList = re.split(":", splitLine[9])
                elif (columnsList[0] == "RNA_NORMAL"):
                    event_rnaNormalList = re.split(":", splitLine[9])
                elif (columnsList[0] == "DNA_TUMOR"):
                    event_dnaTumorList = re.split(":", splitLine[9])
                elif (columnsList[0] == "RNA_TUMOR"):
                    event_rnaTumorList = re.split(":", splitLine[9])
            # if we have a 10th column, figure out which dataset it is
            if (len(splitLine) > 10):
                if (columnsList[1] == "RNA_NORMAL"):
                    event_rnaNormalList = re.split(":", splitLine[10])
                elif (columnsList[1] == "DNA_TUMOR"):
                    event_dnaTumorList = re.split(":", splitLine[10])
                elif (columnsList[1] == "RNA_TUMOR"):
                    event_rnaTumorList = re.split(":", splitLine[10]) 
            # if we have a 11th column, figure out which dataset it is
            if (len(splitLine) > 11):
                if (columnsList[2] == "DNA_TUMOR"):
                    event_dnaTumorList = re.split(":", splitLine[11])
                elif (columnsList[2] == "RNA_TUMOR"):
                    event_rnaTumorList = re.split(":", splitLine[11])
            # if we have a 12th column, figure out which dataset it is
            if (len(splitLine) > 12):
                if (columnsList[3] == "RNA_TUMOR"):
                    event_rnaTumorList = re.split(":", splitLine[12]) 
                
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
            for format in event_formatList:
                if (format == "GT"):
                    sep = "/"
                else:
                    sep = ","
                    
                if (haveDnaNormData):
                    dnaNormalItem = event_dnaNormalList[index]
                    event_dnaNormalDict[format] = re.split(sep, dnaNormalItem)
                if (haveRnaNormData):
                    rnaNormalItem = event_rnaNormalList[index]
                    event_rnaNormalDict[format] = re.split(sep, rnaNormalItem)
                if (haveDnaTumData):
                    dnaTumorItem = event_dnaTumorList[index]
                    event_dnaTumorDict[format] = re.split(sep, dnaTumorItem)
                if (haveRnaTumData):
                    rnaTumorItem = event_rnaTumorList[index]
                    event_rnaTumorDict[format] = re.split(sep, rnaTumorItem)
                index += 1
            
            '''
            # get the ref and alt counts
            numAlts = len(event_altList)
            
            # if we have more than 1 alt, then add a filter
            if (numAlts > 1):
                event_filterList.append("multi")
            '''
            
            # fix the original genotypes
            genotypeIndex = event_formatList.index("GT")
            if (haveDnaNormData):
                event_dnaNormalDict["GT"] = fix_genotypes(event_chr, event_dnaNormalDict["GT"], map(int, event_dnaNormalDict["AD"]), aGTMinDepth, aGTMinPct)
                event_dnaNormalList[genotypeIndex] = "/".join(map(str, event_dnaNormalDict["GT"]))
            if (haveRnaNormData):
                event_rnaNormalDict["GT"] = fix_genotypes(event_chr, event_rnaNormalDict["GT"], map(int, event_rnaNormalDict["AD"]), aGTMinDepth, aGTMinPct)
                event_rnaNormalList[genotypeIndex] = "/".join(map(str, event_rnaNormalDict["GT"]))
            if (haveDnaTumData):
                event_dnaTumorDict["GT"] = fix_genotypes(event_chr, event_dnaTumorDict["GT"], map(int, event_dnaTumorDict["AD"]), aGTMinDepth, aGTMinPct)
                event_dnaTumorList[genotypeIndex] = "/".join(map(str, event_dnaTumorDict["GT"]))
            if (haveRnaTumData):
                event_rnaTumorDict["GT"] = fix_genotypes(event_chr, event_rnaTumorDict["GT"], map(int, event_rnaTumorDict["AD"]), aGTMinDepth, aGTMinPct)
                event_rnaTumorList[genotypeIndex] = "/".join(map(str, event_rnaTumorDict["GT"]))
            
            # combine the refs and alts in one list
            refPlusAltList = event_refList + event_altList
            
            # get rid of bad mod types that don't meet the minimum requirements
            (event_infoDict, filterSet) = pre_filter_mod_types(refPlusAltList, event_infoDict, map(int, event_dnaNormalDict["AD"]), map(int, event_rnaNormalDict["AD"]), map(int, event_dnaTumorDict["AD"]), map(int, event_rnaTumorDict["AD"]), aModMinDepth, aModMinPct, anLohMaxDepth, anLohMaxPct)
            
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
                        filterSet.add("dnmnabq")
                        
                    sourceDepth = int(event_dnaNormalDict["AD"][0])
                    targetDepth = int(event_dnaNormalDict["AD"][targetIndex])
                    sourceStrbias = float(event_dnaNormalDict["SB"][0])
                    targetStrbias = float(event_dnaNormalDict["SB"][targetIndex])
                    
                    # check to make sure the normal DNA strand bias of ALT bases is in range
                    if (sourceDepth > 0 and targetDepth > 0):
                        if (targetStrbias == 0.0 and sourceStrbias != targetStrbias):
                            isValidMod = False
                            filterSet.add("dnsbias")
                        elif (targetStrbias == 1.0 and sourceStrbias != targetStrbias):
                            isValidMod = False
                            filterSet.add("dnsbias")
                    elif (targetStrbias < (aDnaNormParamsDict["StrandBias"]) or targetStrbias > (1.0 - aDnaNormParamsDict["StrandBias"])):
                        isValidMod = False
                        filterSet.add("dnsbias")
                        
                    # we want to make sure the normal DNA sample error percentage is below the max
                    # we want to make sure that the percentage of other ALTs in this sample is below the max error
                    # this is the same as making sure that the percentage of the source and target alleles is above one minus the max error percentage
                    if ((float(event_dnaNormalDict["AF"][sourceIndex]) + float(event_dnaNormalDict["AF"][targetIndex])) < (1.0-aDnaNormParamsDict["MaxErrPct"])):
                        isValidMod = False
                        filterSet.add("dnmxerr")
                        
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
                                filterSet.add("rtmnabq")
                                
                            # check to make sure the tumor DNA strand bias of ALT bases is in range
                            if (float(event_rnaTumorDict["SB"][targetIndex]) < (anRnaTumParamsDict["StrandBias"]) or float(event_rnaTumorDict["SB"][targetIndex]) > (1.0 - anRnaTumParamsDict["StrandBias"])):
                                isValidMod = False
                                filterSet.add("rtbias")
                                
                            # we want to make sure the sample error percentage is below the max
                            # we want to make sure that the total percentage of ALTs in the sample is below the max error
                            # this is the same as making sure that the percentage of the source is above one minus the max error percentage
                            #if (float(event_rnaTumorDict["AF"][sourceIndex]) < (1.0-anRnaTumParamsDict["MaxErrPct"])):
                            if ((float(event_rnaTumorDict["AF"][sourceIndex]) + float(event_rnaTumorDict["AF"][targetIndex])) < (1.0-anRnaTumParamsDict["MaxErrPct"])):
                                isValidMod = False
                                filterSet.add("rtmxerr")
                        elif (anRnaTumParamsDict["MinTotalNumBases"] > 0):
                            isValidMod = False
                            filterSet.add("rtmntb")
                        
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
                        filterSet.add("dtmnabq")
                    
                    sourceDepth = int(event_dnaTumorDict["AD"][0])
                    targetDepth = int(event_dnaTumorDict["AD"][targetIndex])
                    sourceStrbias = float(event_dnaTumorDict["SB"][0])
                    targetStrbias = float(event_dnaTumorDict["SB"][targetIndex])
                    
                    # check to make sure the normal DNA strand bias of ALT bases is in range
                    if (sourceDepth > 0 and targetDepth > 0):
                        if (targetStrbias == 0.0 and sourceStrbias != targetStrbias):
                            isValidMod = False
                            filterSet.add("dtbias")
                        elif (targetStrbias == 1.0 and sourceStrbias != targetStrbias):
                            isValidMod = False
                            filterSet.add("dtbias")
                    elif (targetStrbias < (aDnaTumParamsDict["StrandBias"]) or targetStrbias > (1.0 - aDnaTumParamsDict["StrandBias"])):
                        isValidMod = False
                        filterSet.add("dtbias")
                            
                    # we want to make sure the sample error percentage is below the max
                    # we want to make sure that the percentage of other ALTs in the tumor sample is below the max error
                    # this is the same as making sure that the percentage of the source and target alleles is above one minus the max error percentage
                    if ((float(event_dnaTumorDict["AF"][sourceIndex]) + float(event_dnaTumorDict["AF"][targetIndex])) < (1.0-aDnaTumParamsDict["MaxErrPct"])):
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
                            
                        # we want to make sure the sample error percentage is below the max
                        # we want to make sure that the total percentage of ALTs in the sample is below the max error
                        # this is the same as making sure that the percentage of the source is above one minus the max error percentage
                        if (float(event_dnaNormalDict["AF"][sourceIndex]) < (1.0-aDnaNormParamsDict["MaxErrPct"])):
                            isValidMod = False
                            filterSet.add("dnmxerr")
                    else:
                        isValidMod = False
                        filterSet.add("dnmntb")
                    
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
                                filterSet.add("rtmnabq")
                                
                            sourceDepth = int(event_rnaTumorDict["AD"][0])
                            targetDepth = int(event_rnaTumorDict["AD"][targetIndex])
                            sourceStrbias = float(event_rnaTumorDict["SB"][0])
                            targetStrbias = float(event_rnaTumorDict["SB"][targetIndex])
                            
                            # check to make sure the normal DNA strand bias of ALT bases is in range
                            if (sourceDepth > 0 and targetDepth > 0):
                                if (targetStrbias == 0.0 and sourceStrbias != targetStrbias):
                                    isValidMod = False
                                    filterSet.add("rtbias")
                                elif (targetStrbias == 1.0 and sourceStrbias != targetStrbias):
                                    isValidMod = False
                                    filterSet.add("rtbias")
                            elif (targetStrbias < (anRnaTumParamsDict["StrandBias"]) or targetStrbias > (1.0 - anRnaTumParamsDict["StrandBias"])):
                                isValidMod = False
                                filterSet.add("rtbias")    
                                
                            # we want to make sure the sample error percentage is below the max
                            # we want to make sure that the total percentage of ALTs in the sample is below the max error
                            # this is the same as making sure that the percentage of the source is above one minus the max error percentage
                            #if (float(event_rnaTumorDict["AF"][sourceIndex]) < (1.0-anRnaTumParamsDict["MaxErrPct"])):
                            if ((float(event_rnaTumorDict["AF"][sourceIndex]) + float(event_rnaTumorDict["AF"][targetIndex])) < (1.0-anRnaTumParamsDict["MaxErrPct"])):
                                isValidMod = False
                                filterSet.add("rtmxerr")
                        elif (anRnaTumParamsDict["MinTotalNumBases"] > 0):
                            isValidMod = False
                            filterSet.add("rtmntb")
        
                elif (modType.find("NOR_EDIT") != -1):
                    sourceIndex = refPlusAltList.index(source)
                    targetIndex = refPlusAltList.index(target)
                    
                elif (modType.find("TUM_EDIT") != -1):
                    sourceIndex = refPlusAltList.index(source)
                    targetIndex = refPlusAltList.index(target)
                    
                    # check to make sure the normal DNA sample has data and is between the min and the max of total bases
                    if (haveDnaNormData):
                        if (int(event_dnaNormalDict["DP"][0]) < aDnaNormParamsDict["MinTotalNumBases"]):
                            isValidMod = False
                            filterSet.add("dnmntb")
                        elif (int(event_dnaNormalDict["DP"][0]) > aDnaNormParamsDict["MaxTotalNumBases"]):
                            isValidMod = False
                            filterSet.add("dnmxtb")
                            
                        # we want to make sure the sample error percentage is below the max
                        # we want to make sure that the total percentage of ALTs in the normal sample is below the max error
                        # this is the same as making sure that the percentage of the source is above one minus the max error percentage
                        if ((len(event_dnaNormalDict["AF"]) > sourceIndex) and float(event_dnaNormalDict["AF"][sourceIndex]) < (1.0-aDnaNormParamsDict["MaxErrPct"])):
                            isValidMod = False
                            filterSet.add("dnmxerr")
                    else:
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
                            
                        # we want to make sure the sample error percentage is below the max
                        # we want to make sure that the total percentage of ALTs in the tumor sample is below the max error
                        # this is the same as making sure that the percentage of the source is above one minus the max error percentage
                        if ((len(event_dnaTumorDict["AF"]) > sourceIndex) and float(event_dnaTumorDict["AF"][sourceIndex]) < (1.0-aDnaTumParamsDict["MaxErrPct"])):
                            isValidMod = False
                            filterSet.add("dtmxerr")
                    else:
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
                        filterSet.add("rtmnabq")
                    
                    sourceDepth = int(event_rnaTumorDict["AD"][0])
                    targetDepth = int(event_rnaTumorDict["AD"][targetIndex])
                    sourceStrbias = float(event_rnaTumorDict["SB"][0])
                    targetStrbias = float(event_rnaTumorDict["SB"][targetIndex])
                    
                    # check to make sure the normal DNA strand bias of ALT bases is in range
                    if (sourceDepth > 0 and targetDepth > 0):
                        if (targetStrbias == 0.0 and sourceStrbias != targetStrbias):
                            isValidMod = False
                            filterSet.add("rtbias")
                        elif (targetStrbias == 1.0 and sourceStrbias != targetStrbias):
                            isValidMod = False
                            filterSet.add("rtbias")
                    elif (targetStrbias < (anRnaTumParamsDict["StrandBias"]) or targetStrbias > (1.0 - anRnaTumParamsDict["StrandBias"])):
                        isValidMod = False
                        filterSet.add("rtbias")
                    
                    # we want to make sure the sample error percentage is below the max
                    # we want to make sure that the percentage of other ALTs in the tumor sample is below the max error
                    # this is the same as making sure that the percentage of the source and target alleles is above one minus the max error percentage
                    if ((float(event_rnaTumorDict["AF"][sourceIndex]) + float(event_rnaTumorDict["AF"][targetIndex])) < (1.0-anRnaTumParamsDict["MaxErrPct"])):
                        isValidMod = False
                        filterSet.add("rtmxerr")
                
                # if this one is not valid and we have more, 
                # then remove this one and try the next one
                if (len(modTypesList) > 1 and not isValidMod):
                    modTypesList.remove(modType)
                    modChangesList.remove(modChange)
                    #print "one not valid, try next", refPlusAltList, event_dnaNormalDict["AD"], event_dnaTumorDict["AD"], event_infoDict
                # otherwise add the appropriate filters
                else:
                    event_filterSet = event_filterSet.union(filterSet)
                    
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
                
            vcfOutputList.append(";".join(event_filterSet))
            vcfOutputList.append(";".join(event_infoList))
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
                
    
    if (aStatsDir != None):
        '''
        #i_readDPFileHandler = open("/inside/grotto/users/aradenba/data/hg18/benchmark/radia/stats/" + i_id + "_" + i_chrom + "readDepth.tab", "w")
        #i_readDPFileHandler = open("../stats/" + i_id + "_" + i_chrom + "_readDepth.tab", "w")
        i_readDPFileHandler = open(aStatsDir + aTCGAId + "_" + aChrom + "_readDepth.tab", "w")
        i_readDPFileHandler.write("ReadDepth\tDNANormDP\tRNANormDP\tDNATumDP\tRNATumDP\n")
        
        dnaNormalMaxDP = 0
        rnaNormalMaxDP = 0
        dnaTumorMaxDP = 0
        rnaTumorMaxDP = 0
        if dnaNormalReadDPDict:
            dnaNormalMaxDP = max(map(int, dnaNormalReadDPDict.iterkeys()))
        if rnaNormalReadDPDict:
            rnaNormalMaxDP = max(map(int, rnaNormalReadDPDict.iterkeys()))
        if dnaTumorReadDPDict:
            dnaTumorMaxDP = max(map(int, dnaTumorReadDPDict.iterkeys()))
        if rnaTumorReadDPDict:
            rnaTumorMaxDP = max(map(int, rnaTumorReadDPDict.iterkeys()))
        
        maxDP = max(dnaNormalMaxDP, dnaTumorMaxDP, rnaNormalMaxDP, rnaTumorMaxDP)
        for index in range(1,(maxDP+1)):
            i_readDPFileHandler.write(str(index) + "\t" + str(dnaNormalReadDPDict[str(index)]) + "\t" + str(rnaNormalReadDPDict[str(index)]) + "\t" + str(dnaTumorReadDPDict[str(index)]) + "\t" + str(rnaTumorReadDPDict[str(index)]) +"\n")
        i_readDPFileHandler.close()
        
        #i_altPerFileHandler = open("/inside/grotto/users/aradenba/data/hg18/benchmark/radia/stats/" + i_id + "_" + i_chrom + "altPercent.tab", "w")
        #i_altPerFileHandler = open("../stats/" + i_id + "_" + i_chrom + "_altPercent.tab", "w")
        i_altPerFileHandler = open(aStatsDir + aTCGAId + "_" + aChrom + "_altPercent.tab", "w")
        i_altPerFileHandler.write("AltPercent\tDNANormAP\tRNANormAP\tDNATumAP\RNATumAP\n")
        
        for index in range(1,101):
            i_altPerFileHandler.write(str(index) + "\t" + str(dnaNormalAltPercentDict[str(index)]) + "\t" + str(rnaNormalAltPercentDict[str(index)]) + "\t" + str(dnaTumorAltPercentDict[str(index)]) + "\t" + str(rnaTumorAltPercentDict[str(index)])+"\n")
        i_altPerFileHandler.close()
        '''
                  
    logging.warning("Chrom %s and Id %s: %s events passed out of %s total events", aChrom, aTCGAId, includedEvents, totalEvents)
    
    # close the files 
    if (anOutputFilename != None):
        i_outputFileHandler.close()
      
    i_vcfFileHandler.close()    
    return


def main():
    
    # python2.7 filterByReadSupportVCF.py TCGA-AB-2995 12 ../data/test/TCGA-AB-2995.vcf --log=DEBUG
    
    # create the usage statement
    usage = "usage: python2.7 %prog id chrom vcfFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    # add the optional params    
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, STDOUT by default")
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    i_cmdLineParser.add_option("-s", "--statsDir", dest="statsDir", metavar="STATS_DIR", help="a stats directory where some basic stats can be output")
    i_cmdLineParser.add_option("-f", "--filterUsingRNA", action="store_true", default=False, dest="filterUsingRNA", help="include this argument if the germline and somatic calls should be filtered by the RNA")
    
    i_cmdLineParser.add_option("", "--genotypeMinDepth", type="int", default=int(4), dest="genotypeMinDepth", metavar="GT_MIN_DP", help="the minimum number of bases required for the genotype, %default by default")
    i_cmdLineParser.add_option("", "--genotypeMinPct", type="float", default=float(.10), dest="genotypeMinPct", metavar="GT_MIN_PCT", help="the minimum percentage of reads required for the genotype, %default by default")
    i_cmdLineParser.add_option("", "--modMinDepth", type="int", default=int(4), dest="modMinDepth", metavar="MOD_MIN_DP", help="the minimum number of bases required for a modification, %default by default")
    i_cmdLineParser.add_option("", "--modMinPct", type="float", default=float(.10), dest="modMinPct", metavar="MOD_MIN_PCT", help="the minimum percentage of reads required for a modification, %default by default")
    i_cmdLineParser.add_option("", "--lohMaxDepth", type="int", default=int(2), dest="lohMaxDepth", metavar="LOH_MAX_DP", help="the maximum number of bases allowed in the tumor DNA for an LOH, %default by default")
    i_cmdLineParser.add_option("", "--lohMaxPct", type="float", default=float(.02), dest="lohMaxPct", metavar="LOH_MAX_PCT", help="the maximum percentage of reads in the tumor DNA for an LOH, %default by default")
    i_cmdLineParser.add_option("", "--minAltAvgBaseQual", type="float", default=float(20.0), dest="minAltAvgBaseQual", metavar="MIN_ALT_AVG_BQ", help="the minimum average base quality for the ALT reads, %default by default")
    
    i_cmdLineParser.add_option("", "--dnaNormalMinTotalBases", type="int", default=int(10), dest="dnaNormalMinTotalNumBases", metavar="DNA_NOR_MIN_TOTAL_BASES", help="the minimum number of overall normal DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMaxTotalBases", type="int", default=int(20000), dest="dnaNormalMaxTotalNumBases", metavar="DNA_NOR_MAX_TOTAL_BASES", help="the maximum number of overall normal DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinAltBases", type="int", default=int(4), dest="dnaNormalMinAltNumBases", metavar="DNA_NOR_MIN_ALT_BASES", help="the minimum number of alternative normal DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMinAltPct", type="float", default=float(0.10), dest="dnaNormalMinAltPct", metavar="DNA_NOR_MIN_ALT_PCT", help="the minimum percentage of alternative normal DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalMaxErrPct", type="float", default=float(0.01), dest="dnaNormalMaxErrPct", metavar="DNA_NOR_MAX_ERR_PCT", help="the maximum percentage of alternative normal DNA reads allowed that support a variant, %default by default")
    i_cmdLineParser.add_option("", "--dnaNormalStrandBias", type="float", default=float(0.01), dest="dnaNormalStrandBias", metavar="DNA_NOR_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    
    i_cmdLineParser.add_option("", "--dnaTumorMinTotalBases", type="int", default=int(10), dest="dnaTumorMinTotalNumBases", metavar="DNA_TUM_MIN_TOTAL_BASES", help="the minimum number of overall tumor DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMaxTotalBases", type="int", default=int(20000), dest="dnaTumorMaxTotalNumBases", metavar="DNA_TUM_MAX_TOTAL_BASES", help="the maximum number of overall tumor DNA reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinAltBases", type="int", default=int(4), dest="dnaTumorMinAltNumBases", metavar="DNA_TUM_MIN_ALT_BASES", help="the minimum number of alternative tumor DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMinAltPct", type="float", default=float(0.10), dest="dnaTumorMinAltPct", metavar="DNA_TUM_MIN_ALT_PCT", help="the minimum percentage of alternative tumor DNA reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorMaxErrPct", type="float", default=float(0.01), dest="dnaTumorMaxErrPct", metavar="DNA_TUM_MAX_ERR_PCT", help="the maximum percentage of alternative tumor DNA reads allowed that support a variant in the tumor RNA, %default by default")
    i_cmdLineParser.add_option("", "--dnaTumorStrandBias", type="float", default=float(0.01), dest="dnaTumorStrandBias", metavar="DNA_TUM_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    
    i_cmdLineParser.add_option("", "--rnaNormalMinTotalBases", type="int", default=int(10), dest="rnaNormalMinTotalNumBases", metavar="RNA_NOR_MIN_TOTAL_BASES", help="the minimum number of overall normal RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMaxTotalBases", type="int", default=int(20000), dest="rnaNormalMaxTotalNumBases", metavar="RNA_NOR_MAX_TOTAL_BASES", help="the maximum number of overall normal RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinAltBases", type="int", default=int(4), dest="rnaNormalMinAltNumBases", metavar="RNA_NOR_MIN_ALT_BASES", help="the minimum number of alternative normal RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMinAltPct", type="float", default=float(0.10), dest="rnaNormalMinAltPct", metavar="RNA_NOR_MIN_ALT_PCT", help="the minimum percentage of alternative normal RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalMaxErrPct", type="float", default=float(0.01), dest="rnaNormalMaxErrPct", metavar="RNA_NOR_MAX_ERR_PCT", help="the maximum percentage of alternative normal RNA reads allowed that support a variant in the tumor, %default by default")
    i_cmdLineParser.add_option("", "--rnaNormalStrandBias", type="float", default=float(0.01), dest="rnaNormalStrandBias", metavar="RNA_NOR_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    
    i_cmdLineParser.add_option("", "--rnaTumorMinTotalBases", type="int", default=int(10), dest="rnaTumorMinTotalNumBases", metavar="RNA_TUM_MIN_TOTAL_BASES", help="the minimum number of overall tumor RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMaxTotalBases", type="int", default=int(20000), dest="rnaTumorMaxTotalNumBases", metavar="RNA_TUM_MAX_TOTAL_BASES", help="the maximum number of overall tumor RNA-Seq reads covering a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinAltBases", type="int", default=int(4), dest="rnaTumorMinAltNumBases", metavar="RNA_TUM_MIN_ALT_BASES", help="the minimum number of alternative tumor RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMinAltPct", type="float", default=float(0.10), dest="rnaTumorMinAltPct", metavar="RNA_TUM_MIN_ALT_PCT", help="the minimum percentage of alternative tumor RNA-Seq reads supporting a variant at a position, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorMaxErrPct", type="float", default=float(0.01), dest="rnaTumorMaxErrPct", metavar="RNA_TUM_MAX_ERR_PCT", help="the maximum percentage of alternative tumor RNA reads allowed that support a variant in a future data type, %default by default")
    i_cmdLineParser.add_option("", "--rnaTumorStrandBias", type="float", default=float(0.01), dest="rnaTumorStrandBias", metavar="RNA_TUM_STRAND_BIAS", help="the maximum percentage of strand bias on reads that support the ALT, %default by default")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(4,51,1)
    i_argLength = len(sys.argv)
    #print >> sys.stderr, "argLen", i_argLength
    
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
    i_filterUsingRNA = i_cmdLineOptions.filterUsingRNA
    i_genotypeMinDepth = i_cmdLineOptions.genotypeMinDepth
    i_genotypeMinPct = i_cmdLineOptions.genotypeMinPct
    i_modMinDepth = i_cmdLineOptions.modMinDepth
    i_modMinPct = i_cmdLineOptions.modMinPct
    i_lohMaxDepth = i_cmdLineOptions.lohMaxDepth
    i_lohMaxPct = i_cmdLineOptions.lohMaxPct
    i_minAltAvgBaseQual = i_cmdLineOptions.minAltAvgBaseQual
    
    i_dnaNormParams = {}
    i_dnaNormParams["MinTotalNumBases"] = i_cmdLineOptions.dnaNormalMinTotalNumBases
    i_dnaNormParams["MaxTotalNumBases"] = i_cmdLineOptions.dnaNormalMaxTotalNumBases
    i_dnaNormParams["MinAltNumBases"] = i_cmdLineOptions.dnaNormalMinAltNumBases
    i_dnaNormParams["MinAltPct"] = i_cmdLineOptions.dnaNormalMinAltPct
    i_dnaNormParams["MaxErrPct"] = i_cmdLineOptions.dnaNormalMaxErrPct
    i_dnaNormParams["StrandBias"] = i_cmdLineOptions.dnaNormalStrandBias
    
    i_dnaTumParams = {}
    i_dnaTumParams["MinTotalNumBases"] = i_cmdLineOptions.dnaTumorMinTotalNumBases
    i_dnaTumParams["MaxTotalNumBases"] = i_cmdLineOptions.dnaTumorMaxTotalNumBases
    i_dnaTumParams["MinAltNumBases"] = i_cmdLineOptions.dnaTumorMinAltNumBases
    i_dnaTumParams["MinAltPct"] = i_cmdLineOptions.dnaTumorMinAltPct
    i_dnaTumParams["MaxErrPct"] = i_cmdLineOptions.dnaTumorMaxErrPct
    i_dnaTumParams["StrandBias"] = i_cmdLineOptions.dnaTumorStrandBias
    
    i_rnaNormParams = {}
    i_rnaNormParams["MinTotalNumBases"] = i_cmdLineOptions.rnaNormalMinTotalNumBases
    i_rnaNormParams["MaxTotalNumBases"] = i_cmdLineOptions.rnaNormalMaxTotalNumBases
    i_rnaNormParams["MinAltNumBases"] = i_cmdLineOptions.rnaNormalMinAltNumBases
    i_rnaNormParams["MinAltPct"] = i_cmdLineOptions.rnaNormalMinAltPct
    i_rnaNormParams["MaxErrPct"] = i_cmdLineOptions.rnaNormalMaxErrPct
    i_rnaNormParams["StrandBias"] = i_cmdLineOptions.rnaNormalStrandBias
    
    i_rnaTumParams = {}
    i_rnaTumParams["MinTotalNumBases"] = i_cmdLineOptions.rnaTumorMinTotalNumBases
    i_rnaTumParams["MaxTotalNumBases"] = i_cmdLineOptions.rnaTumorMaxTotalNumBases
    i_rnaTumParams["MinAltNumBases"] = i_cmdLineOptions.rnaTumorMinAltNumBases
    i_rnaTumParams["MinAltPct"] = i_cmdLineOptions.rnaTumorMinAltPct
    i_rnaTumParams["MaxErrPct"] = i_cmdLineOptions.rnaTumorMaxErrPct
    i_rnaTumParams["StrandBias"] = i_cmdLineOptions.rnaTumorStrandBias
    
    # try to get any optional parameters with no defaults    
    i_readFilenameList = [i_vcfFilename]
    i_writeFilenameList = []
    i_dirList = []
    
    i_outputFilename = None
    i_logFilename = None
    i_statsDir = None
    if (i_cmdLineOptions.outputFilename != None):
        i_outputFilename = str(i_cmdLineOptions.outputFilename)
        i_writeFilenameList += [i_outputFilename]
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
        i_writeFilenameList += [i_logFilename]
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
        logging.debug("logLevel=%s" % i_logLevel)
        logging.debug("logFilename=%s" % i_logFilename)
        logging.debug("i_statsDir=%s" % i_statsDir)
        logging.debug("i_filterUsingRNA=%s" % i_filterUsingRNA)
        
        logging.debug("i_genotypeMinDepth=%s" % i_genotypeMinDepth)
        logging.debug("i_genotypeMinPct=%s" % i_genotypeMinPct)
        logging.debug("i_modMinDepth=%s" % i_modMinDepth)
        logging.debug("i_modMinPct=%s" % i_modMinPct)
        logging.debug("i_lohMaxDepth=%s" % i_lohMaxDepth)
        logging.debug("i_lohMaxPct=%s" % i_lohMaxPct)
        logging.debug("i_minAltAvgBaseQual=%s" % i_minAltAvgBaseQual)
        
        logging.debug("dna normal minTotalBases: %s" % i_dnaNormParams["MinTotalNumBases"])
        logging.debug("dna normal maxTotalBases: %s" % i_dnaNormParams["MaxTotalNumBases"])
        logging.debug("dna normal minAltBases: %s" % i_dnaNormParams["MinAltNumBases"])
        logging.debug("dna normal minAltPct: %s" % i_dnaNormParams["MinAltPct"])
        logging.debug("dna normal maxErrPct: %s" % i_dnaNormParams["MaxErrPct"])
        logging.debug("dna normal strandBias: %s" % i_dnaNormParams["StrandBias"])
        
        logging.debug("dna tumor minTotalBases: %s" % i_dnaTumParams["MinTotalNumBases"])
        logging.debug("dna tumor maxTotalBases: %s" % i_dnaTumParams["MaxTotalNumBases"])
        logging.debug("dna tumor minAltBases: %s" % i_dnaTumParams["MinAltNumBases"])
        logging.debug("dna tumor minAltPct: %s" % i_dnaTumParams["MinAltPct"])
        logging.debug("dna tumor maxErrPct: %s" % i_dnaTumParams["MaxErrPct"])
        logging.debug("dna tumor strandBias: %s" % i_dnaTumParams["StrandBias"])
        
        logging.debug("rna normal minTotalBases: %s" % i_rnaNormParams["MinTotalNumBases"])
        logging.debug("rna normal maxTotalBases: %s" % i_rnaNormParams["MaxTotalNumBases"])
        logging.debug("rna normal minAltBases: %s" % i_rnaNormParams["MinAltNumBases"])
        logging.debug("rna normal minAltPct: %s" % i_rnaNormParams["MinAltPct"])
        logging.debug("rna normal maxErrPct: %s" % i_rnaNormParams["MaxErrPct"])
        logging.debug("rna normal strandBias: %s" % i_rnaNormParams["StrandBias"])
        
        logging.debug("rna tumor minTotalBases: %s" % i_rnaTumParams["MinTotalNumBases"])
        logging.debug("rna tumor maxTotalBases: %s" % i_rnaTumParams["MaxTotalNumBases"])
        logging.debug("rna tumor minAltBases: %s" % i_rnaTumParams["MinAltNumBases"])
        logging.debug("rna tumor minAltPct: %s" % i_rnaTumParams["MinAltPct"])
        logging.debug("rna tumor maxErrPct: %s" % i_rnaTumParams["MaxErrPct"])
        logging.debug("rna tumor strandBias: %s" % i_rnaTumParams["StrandBias"])
                    
    # check for any errors
    if (not rnaEditingUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)
    
    extract_read_support(i_id, i_chrom, i_vcfFilename, i_outputFilename, i_filterUsingRNA, i_minAltAvgBaseQual, i_statsDir, i_cmdLineOptionsDict, i_dnaNormParams, i_dnaTumParams, i_rnaNormParams, i_rnaTumParams, i_genotypeMinDepth, i_genotypeMinPct, i_modMinDepth, i_modMinPct, i_lohMaxDepth, i_lohMaxPct, i_debug)
    return

main()
sys.exit(0)