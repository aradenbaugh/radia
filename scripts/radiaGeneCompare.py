#!/usr/bin/env python

import sys                          # system module
from optparse import OptionParser   # used for parsing command line arguments
import radiaUtil                    # utility functions for rna editing
import logging
import time
import collections
import gzip


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
    

def get_append_fileHandler(aFilename):
    '''
    ' Open aFilename for appending and return
    ' the file handler.  The file can be 
    ' gzipped or not.
    '''
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename,'ab')
    else:
        return open(aFilename,'a')


def get_extracted_data(anInputFilename, aMissenseOnlyFlag, anIsDebug):
    
    inputFileHandler = get_read_fileHandler(anInputFilename)
    geneRescueDict = collections.defaultdict(set)
    geneRescueCountDict = dict()
    geneEditingDict = collections.defaultdict(set)
    geneEditingCountDict = dict()
     
    for line in inputFileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
                
        # these lines are from previous scripts in the pipeline, so output them    
        elif (line.startswith("#")):
            continue;
        
        # now we are to the data
        else:    
            
            # split the line on the tab
            splitLine = line.split("\t")
            
            # get the fields to yield
            # sample chrom coordinate pass oldModType newModType modChange 
            # dnDepth dnRefDp dnAltDp dnAltPct    
            # dtDepth dtRefDp dtAltDp dtAltPct  
            # rtDepth rtRefDp rtAltDp rtAltPct    
            # mmp gene vc cc pc    
            # malawiMutSig tcgaMutSig filters
            sample = splitLine[0]
            chrom = splitLine[1]
            coordinate = splitLine[2]
            passing = splitLine[3]
            modType = splitLine[5]
            modChange = splitLine[6]
            mmp = float(splitLine[19])
            gene = splitLine[20].split("/")[0]
            vc = splitLine[21]
            cc = splitLine[22]
            pc = splitLine[23]
            
            if (anIsDebug):
                logging.debug("sample=%s, gene=%s, pass?=%s, modType=%s, vc=%s", sample, gene, passing, modType, vc)
            
            if (passing == "False"):
                continue
            
            if (mmp >= 0.95):
                continue
                
            if (aMissenseOnlyFlag and (vc != "Missense" and vc != "Nonsense" and vc != "Readthrough" and vc != "SpliceSite")):
                continue
            
            if (modType == "RADIASomRNARescue"):
                geneRescueDict[gene].add(sample)
                
                if (gene not in geneRescueCountDict):
                    geneRescueCountDict[gene] = dict()
                    geneRescueCountDict[gene]["misCount"] = 0
                    geneRescueCountDict[gene]["otherCount"] = 0
                    
                if (vc == "Missense" or vc == "Nonsense" or vc == "Readthrough" or vc == "SpliceSite"):
                    geneRescueCountDict[gene]["misCount"] += 1
                else:
                    geneRescueCountDict[gene]["otherCount"] += 1
                    
            if (modType == "RADIATumEdit"):
                
                geneEditingDict[gene].add(sample)
                
                if (gene not in geneEditingCountDict):
                    geneEditingCountDict[gene] = dict()
                    geneEditingCountDict[gene]["misCount"] = 0
                    geneEditingCountDict[gene]["otherCount"] = 0
                    geneEditingCountDict[gene]["sites"] = dict()
                    
                if (vc == "Missense" or vc == "Nonsense" or vc == "Readthrough" or vc == "SpliceSite"):
                    geneEditingCountDict[gene]["misCount"] += 1
                else:
                    geneEditingCountDict[gene]["otherCount"] += 1
                
                if (chrom + "_" + coordinate) not in geneEditingCountDict[gene]["sites"]:
                    geneEditingCountDict[gene]["sites"][chrom + "_" + coordinate] = dict()
                    geneEditingCountDict[gene]["sites"][chrom + "_" + coordinate]["misCount"] = 0
                    geneEditingCountDict[gene]["sites"][chrom + "_" + coordinate]["otherCount"] = 0
                    
                if (vc == "Missense" or vc == "Nonsense" or vc == "Readthrough" or vc == "SpliceSite"):
                    geneEditingCountDict[gene]["sites"][chrom + "_" + coordinate]["misCount"] += 1
                else:
                    geneEditingCountDict[gene]["sites"][chrom + "_" + coordinate]["otherCount"] += 1
                geneEditingCountDict[gene]["sites"][chrom + "_" + coordinate]["modChange"] = modChange
                                        
    inputFileHandler.close()
        
    return (geneRescueDict, geneRescueCountDict, geneEditingDict, geneEditingCountDict)


def get_maf_data(anInputFilename, aMissenseOnlyFlag, anIsDebug):
    
    inputFileHandler = get_read_fileHandler(anInputFilename)
    outputDict = collections.defaultdict(set)
     
    for line in inputFileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
                
        # these lines are from previous scripts in the pipeline, so output them    
        elif (line.startswith("#")):
            continue;
        
        # now we are to the data
        else:    
            
            # split the line on the tab
            splitLine = line.split("\t")
            
            # get the fields to yield
            sample = splitLine[39]
            gene = splitLine[1]
            variantType = splitLine[10]
            vc = splitLine[9]
            
            if (anIsDebug):
                logging.debug("sample=%s, gene=%s, variantType=%s, vc=%s", sample, gene, variantType, vc)
            
            if (variantType != "SNP"):
                continue
            
            if (aMissenseOnlyFlag and (vc != "Missense_Mutation" and vc != "Nonsense_Mutation" and vc != "Readthrough" and vc != "Splice_Site")):
                continue
            
            outputDict[gene].add(sample)
                                    
    inputFileHandler.close()
        
    return (outputDict)
        
    
def main(): 
        
    # create the usage statement
    usage = "usage: python %prog dnaMutrickMaf rnaRescueFile [Options]"
    i_cmdLineParser = OptionParser(usage=usage)
    
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDERR by default")
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file")
    i_cmdLineParser.add_option("-c", "--geneCountFilename", dest="geneCountFilename", metavar="GENE_COUNT_FILE", help="the name of the gene count file")
    i_cmdLineParser.add_option("-s", "--siteCountFilename", dest="siteCountFilename", metavar="SITE_COUNT_FILE", help="the name of the site count file")
    i_cmdLineParser.add_option("-m", "--missenseOnly", action="store_true", default=False, dest="missenseOnly", help="by default all the calls are processed, include this argument if only the (missense, nonsense, readthrough) calls should be processed")
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(3,11,1)
    i_argLength = len(sys.argv)
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)
    
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_dnaFile = str(i_cmdLineArgs[0])
    i_rnaFile = str(i_cmdLineArgs[1])
    
    # get the optional params with default values   
    i_logLevel = i_cmdLineOptions.logLevel
    i_missenseOnly = i_cmdLineOptions.missenseOnly
    
    # try to get any optional parameters with no defaults   
    # check for any errors
    writeFilenameList = []
    readFilenameList = [i_dnaFile, i_rnaFile]
     
    i_logFilename = None
    i_outputFilename = None
    i_geneCountFilename = None
    i_siteCountFilename = None
    if (i_cmdLineOptions.outputFilename != None):
        i_outputFilename = str(i_cmdLineOptions.outputFilename)
        writeFilenameList += [i_outputFilename]
    if (i_cmdLineOptions.geneCountFilename != None):
        i_geneCountFilename = str(i_cmdLineOptions.geneCountFilename)
        writeFilenameList += [i_geneCountFilename]
    if (i_cmdLineOptions.siteCountFilename != None):
        i_siteCountFilename = str(i_cmdLineOptions.siteCountFilename)
        writeFilenameList += [i_siteCountFilename]
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
        writeFilenameList += [i_logFilename]
        
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
        logging.debug("i_rnaFile=%s", i_rnaFile)
        logging.debug("i_dnaFile=%s" % i_dnaFile)
        logging.debug("i_outputFilename=%s" % i_outputFilename)
        logging.debug("i_geneCountFilename=%s" % i_geneCountFilename)
        logging.debug("i_siteCountFilename=%s" % i_siteCountFilename)
        logging.debug("logLevel=%s", i_logLevel)
        logging.debug("logFile=%s", i_logFilename)
        logging.debug("missenseOnly?=%s", i_missenseOnly)
        
    if (not radiaUtil.check_for_argv_errors(None, readFilenameList, writeFilenameList)):
        sys.exit(1)           
    
    # create the generators for the filter and vcf files
    (i_rescueDict, i_rescueGeneCountDict, i_editingDict, i_editingGeneCountDict) = get_extracted_data(i_rnaFile, i_missenseOnly, i_debug)
    (i_dnaDict) = get_maf_data(i_dnaFile, i_missenseOnly, i_debug)
    
    genes = set(i_rescueDict.keys() + i_editingDict.keys() + i_dnaDict.keys())
    
    geneOutputFileHandler = get_write_fileHandler(i_geneCountFilename)
    siteOutputFileHandler = get_write_fileHandler(i_siteCountFilename)
    for gene in genes:
        if gene in i_rescueGeneCountDict:
            rescueMisCount = i_rescueGeneCountDict[gene]["misCount"]
            rescueOtherCount = i_rescueGeneCountDict[gene]["otherCount"]
        else:
            rescueMisCount = 0
            rescueOtherCount = 0
            
        if gene in i_editingGeneCountDict:
            editingMisCount = i_editingGeneCountDict[gene]["misCount"]
            editingOtherCount = i_editingGeneCountDict[gene]["otherCount"]
        else:
            editingMisCount = 0
            editingOtherCount = 0
            
        geneOutputFileHandler.write("\t".join([gene, str(rescueMisCount), str(rescueOtherCount), str(editingMisCount), str(editingOtherCount)]) + "\n")
        
        if gene in i_editingGeneCountDict and "sites" in i_editingGeneCountDict[gene]:
            for site in i_editingGeneCountDict[gene]["sites"]:
                (chrom, coordinate) = site.split("_")
                siteMisCount = i_editingGeneCountDict[gene]["sites"][site]["misCount"]
                siteOtherCount = i_editingGeneCountDict[gene]["sites"][site]["otherCount"]
                modChange = i_editingGeneCountDict[gene]["sites"][site]["modChange"]
                siteOutputFileHandler.write("\t".join([gene, chrom, coordinate, str(siteMisCount), str(siteOtherCount), str(modChange)]) + "\n")
            
    geneOutputFileHandler.close()
    siteOutputFileHandler.close()
    
    outputFileHandler = get_write_fileHandler(i_outputFilename)
    for gene in genes:
        dnaCount = len(i_dnaDict[gene])
        rescueCount = len(i_rescueDict[gene])
        editingCount = len(i_editingDict[gene])
        outputFileHandler.write("\t".join([gene, str(dnaCount), str(rescueCount), str(editingCount)]) + "\n")    
    outputFileHandler.close()
    
    
    '''
    malawiFileHandler = open("/hottub/research/aradenbaugh/radia/data/common/malawi_escc_top_mutated_genes.list", "r")
    tcgaFileHandler = open("/hottub/research/aradenbaugh/radia/data/common/tcga_esca_top_mutated_genes.list", "r")
    
    sigMutGeneSet = set()
    for line in malawiFileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        sigMutGeneSet.add(line)
    
    for line in tcgaFileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        sigMutGeneSet.add(line)
        
    outputFileHandler = get_write_fileHandler(i_geneCountFilename)
    for gene in sigMutGeneSet:
        if gene in i_rescueGeneCountDict:
            rescueMisCount = i_rescueGeneCountDict[gene]["misCount"]
            rescueOtherCount = i_rescueGeneCountDict[gene]["otherCount"]
        else:
            rescueMisCount = 0
            rescueOtherCount = 0
            
        if gene in i_editingGeneCountDict:
            editingMisCount = i_editingGeneCountDict[gene]["misCount"]
            editingOtherCount = i_editingGeneCountDict[gene]["otherCount"]
        else:
            editingMisCount = 0
            editingOtherCount = 0
            
        outputFileHandler.write("\t".join([gene, str(rescueMisCount), str(rescueOtherCount), str(editingMisCount), str(editingOtherCount)]) + "\n")    
    outputFileHandler.close()
    
    outputFileHandler = get_write_fileHandler(i_outputFilename)
    for gene in sigMutGeneSet:
        dnaCount = len(i_dnaDict[gene])
        rescueCount = len(i_rescueDict[gene].difference(i_dnaDict[gene]))
        editingCount = len(i_editingDict[gene].difference((i_dnaDict[gene].union(i_rescueDict[gene]))))
        outputFileHandler.write("\t".join([gene, str(dnaCount), str(rescueCount), str(editingCount)]) + "\n")    
    outputFileHandler.close()
    '''
    
    
    
    return

main()
sys.exit(0)
