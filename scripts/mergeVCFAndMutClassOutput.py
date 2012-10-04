#!/usr/bin/env python2.7

from optparse import OptionParser   # used for parsing command line arguments
import rnaEditingUtil               # utility functions for rna editing
import sys                          # system module
import re
import collections
import logging
from OrderedSet import OrderedSet

'''
'   Amie Radenbaugh - 02/29/2012
'   UCSC - RNA Editing  
'   Program name: "mergeVCFAndMutClassOutput.py"
'''


def parse_mutClass_file(aMutClassFilename, anIsDebug):
    '''
    ' anInputFileHandler: The input stream for the file
    ' anIsDebug: A flag for outputting debug messages to STDERR
    '''
    
    mutClassFileHandler = open(aMutClassFilename, "r")
    mutClassDict = collections.defaultdict(list)
     
    for line in mutClassFileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (anIsDebug):
            logging.debug("mutClass Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
                
        # these lines are from previous scripts in the pipeline, so ignore them    
        elif (line.startswith("#")):
            continue;
        
        # now we are to the data
        else:    
            
            # split the line on the tab
            splitLine = re.split("\t", line)
            
            # get the fields to yield
            chr = splitLine[0]
            chrNum = chr.lstrip("chr")
            startCoordinate = int(splitLine[1])
            
            # add the coordinate to the output
            mutClassDict[chrNum + "_" + str(startCoordinate+1)].append(line)
                            
    mutClassFileHandler.close()
        
    return (mutClassDict)


def get_mutClass_data(aMutClassLineList, anIsDebug):
    
    # INFO SID = ucsc id
    # INFO GENE = gene name
    # INFO RGN = {5_utr, 3_utr, exon, intron, ncds, sp}
    # RE = known RNA-edits occur
    # FORMAT TE = translational effect {SIL, MIS, NSNS, NSTP, FSH, NA}
    # INFO NSUB = nucleotide substitution "g.Y:2889273A>C"
    # INFO CC = codon change "c.220A>C_c.(220-222)ATA>CTA"
    # INFO AAC = amino acid change "p.I74L"
    geneIdSet = OrderedSet()
    geneNameSet = OrderedSet()
    transEffectSet = OrderedSet()
    regionSet = OrderedSet()
    nucSubSet = OrderedSet()
    codonChangeSet = OrderedSet()
    aaChangeSet = OrderedSet()
    
    for line in aMutClassLineList:
        if (anIsDebug):
            logging.debug("mutClass Line: %s", line)
            
        # split the line on the tab
        splitLine = re.split("\t", line)
            
        # get the fields
        geneNameSet.add(splitLine[3])
        geneIdSet.add(splitLine[4])
        
        transEffect = splitLine[5]
        if (transEffect == "Silent"):
            transEffectSet.add("SIL")
            regionSet.add("exon")
        elif (transEffect == "Missense_Mutation"):
            transEffectSet.add("MIS")
            regionSet.add("exon")
        elif (transEffect == "Nonsense_Mutation" or transEffect == "Nonsense_Mutation(Last_Exon)"):
            transEffectSet.add("NSNS")
            regionSet.add("exon")
        elif (transEffect == "Read_Through"):
            transEffectSet.add("NSTP")
            regionSet.add("exon")
        elif (transEffect == "Frame_Shift_Del" or transEffect == "Frame_Shift_Ins"):
            transEffectSet.add("FSH")
            regionSet.add("exon")
        elif (transEffect == "5'UTR"):
            transEffectSet.add("NA")
            regionSet.add("5_utr")
        elif (transEffect == "3'UTR"):
            transEffectSet.add("NA")
            regionSet.add("3_utr")
        elif (transEffect == "Intron"):
            transEffectSet.add("NA")
            regionSet.add("intron")
        elif (transEffect == "Splice_Site_SNP"):
            transEffectSet.add("NA")
            regionSet.add("sp")
        else:
            transEffectSet.add("NA")
            regionSet.add("NA")
        
        if (transEffect == "Missense_Mutation" or transEffect == "Nonsense_Mutation" or transEffect == "Nonsense_Mutation(Last_Exon)" or transEffect == "Read_Through"):
            nucSubSet.add(splitLine[6])
            codonChangeSet.add(splitLine[7] + "_" + splitLine[8])
            aaChangeSet.add(splitLine[9])
        
    # INFO SID = ucsc id
    # INFO GENE = gene name
    # INFO RGN = {5_utr, 3_utr, exon, intron, ncds, sp}
    # RE = known RNA-edits occur
    # FORMAT TE = translational effect {SIL, MIS, NSNS, NSTP, FSH, NA}
    # INFO NSUB = nucleotide substitution "g.Y:2889273A>C"
    # INFO CC = codon change "c.220A>C_c.(220-222)ATA>CTA"
    # INFO AAC = amino acid change "p.I74L"
    
    # add the coordinate to the output
    mutClassInfo =  "SID=" + ",".join(geneIdSet) + ";"
    mutClassInfo += "GENE=" + ",".join(geneNameSet) + ";" 
    mutClassInfo += "RGN=" + ",".join(regionSet) + ";"
    mutClassInfo += "TE=" + ",".join(transEffectSet) + ";"
    if (len(nucSubSet) > 0):
        mutClassInfo += "NC=" + ",".join(nucSubSet) + ";"
        mutClassInfo += "CC=" + ",".join(codonChangeSet) + ";"
        mutClassInfo += "AAC=" + ",".join(aaChangeSet) + ";"
                                    
    return (mutClassInfo)



def main():
    
    # python2.7 mergeVCFAndMutClassOutput.py TCGA-AB-2995 12 ../data/test/TCGA-AB-2995.vcf ../data/test/TCGA-AB-2995_mutClass_input.bed --log=DEBUG
    
    # create the usage statement
    usage = "usage: python2.7 %prog id chrom vcfInputFile mutClassOutputBedFile"
    i_cmdLineParser = OptionParser(usage=usage)
    
    # add the optional params  
    i_cmdLineParser.add_option("-o", "--outputFilename", dest="outputFilename", metavar="OUTPUT_FILE", help="the name of the output file, STDOUT by default")  
    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")
    
    
    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(4,10,1)
    i_argLength = len(sys.argv)
    #print >> sys.stderr, "argLen", i_argLength
    
    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)
        
    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    i_id = str(i_cmdLineArgs[0])
    i_chrom = str(i_cmdLineArgs[1])
    i_vcfInputFilename = str(i_cmdLineArgs[2])
    i_mutClassFilename = str(i_cmdLineArgs[3])
    
    # get the optional params with default values
    i_logLevel = i_cmdLineOptions.logLevel
    
    # try to get any optional parameters with no defaults 
    i_readFilenameList = [i_vcfInputFilename, i_mutClassFilename]
    i_writeFilenameList = []
    i_dirList = []
    
    i_logFilename = None
    i_outputFilename = None
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)
        i_writeFilenameList += [i_logFilename]
    if (i_cmdLineOptions.outputFilename != None):
        i_outputFilename = str(i_cmdLineOptions.outputFilename)
        i_writeFilenameList += [i_outputFilename]
            
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
        logging.debug("i_vcfInputFilename=%s" % i_vcfInputFilename)
        logging.debug("i_mutClassOutputFilename=%s" % i_mutClassFilename)
        logging.debug("i_vcfOutputFilename=%s" % i_outputFilename)
        logging.debug("logLevel=%s" % i_logLevel)
        logging.debug("logFilename=%s" % i_logFilename)
                        
    # check for any errors
    if (not rnaEditingUtil.check_for_argv_errors(i_dirList, i_readFilenameList, i_writeFilenameList)):
        sys.exit(1)

    mutClassDict = parse_mutClass_file(i_mutClassFilename, i_debug)
    
    # header geneAnno field
    # INFO SID = ucsc id
    # INFO GENE = gene name
    # INFO RGN = {5_utr, 3_utr, exon, intron, ncds, sp}
    # RE = known RNA-edits occur
    # FORMAT TE = translational effect {SIL, MIS, NSNS, NSTP, FSH, NA}
    # INFO NSUB = nucleotide substitution "g.Y:2889273A>C"
    # INFO CC = codon change "c.220A>C_c.(220-222)ATA>CTA"
    # INFO AAC = amino acid change "p.I74L"
        
    mainHeaderLines = "##geneAnno=UCSC Internal mutClass\n"
    
    infoHeaderLines = "##INFO=<ID=SID,Number=.,Type=String,Description=\"Unique identifier (UCSC) from gene annotation source or unknown\">\n"
    infoHeaderLines += "##INFO=<ID=GENE,Number=.,Type=String,Description=\"Gene Name (UCSC) that overlaps variant\">\n"
    infoHeaderLines += "##INFO=<ID=RGN,Number=.,Type=String,Description=\"Region where nucleotide variant occurs in relation to a gene\">\n"
    infoHeaderLines += "##INFO=<ID=RE,Number=0,Type=Flag,Description=\"Position known to have RNA-edits to occur\">\n"
    infoHeaderLines += "##INFO=<ID=TE,Number=.,Type=String,Description=\"Translational effect of the variant in a codon\">\n"
    infoHeaderLines += "##INFO=<ID=NC,Number=.,Type=String,Description=\"Nucleotide substitution specifying the source DNA, chrom, coordinate, and substitution\">\n"
    infoHeaderLines += "##INFO=<ID=CC,Number=.,Type=String,Description=\"Codon change specified by both the short and long-hand notation\">\n"
    infoHeaderLines += "##INFO=<ID=AAC,Number=.,Type=String,Description=\"Amino acid change where '*' represents a stop codon\">\n"
    
    #formatHeaderLines = "##FORMAT=<ID=TE,Number=.,Type=String,Description=\"Translational effect of the variant in a codon\">\n"

    # get the files
    i_vcfInputFileHandler = open(i_vcfInputFilename, "r")
    i_mutClassFileHandler = open(i_mutClassFilename, "r")
    if (i_outputFilename != None):
        i_vcfOutputFileHandler = open(i_outputFilename, "w")
            
    infoHeaderAdded = False
    formatHeaderAdded = False
    
    # for each event in the vcf file 
    for line in i_vcfInputFileHandler:
        # here are some examples of .vcf lines:
        # 20      199696  .       G       T       0       PASS    AC=2;AF=0.04;AN=2;BQ=31.4;DP=53;FA=0.04;INDEL=0;MC=G>T;MT=TUM_EDIT;NS=3;SB=0.72;SS=5;START=2;STOP=0;VT=SNP      
        # GT:DP:INDEL:START:STOP:AD:AF:BQ:SB      0/0:2:0:0:0:2:1.0:35.5:0.0      0/0:1:0:0:0:1:1.0:39.0:1.0      0/1:50:0:2:0:48,2:0.96,0.04:31.65,17.5:0.75,0.5
        
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")
        
        if (i_debug):
            logging.debug("VCF Line: %s", line)
            
        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;
        
        # these lines are from previous scripts in the pipeline, so output them    
        elif (line.startswith("##phasing")):
            if (i_outputFilename != None):
                i_vcfOutputFileHandler.write(line + "\n")
                i_vcfOutputFileHandler.write(mainHeaderLines)
            else:
                print >> sys.stdout, line
                print >> sys.stdout, mainHeaderLines.rstrip()
        
        # if we find the INFO section, then add the info from here
        elif ((not infoHeaderAdded) and line.startswith("##INFO")):
            infoHeaderAdded = True
            if (i_outputFilename != None):        
                i_vcfOutputFileHandler.write(infoHeaderLines)
                i_vcfOutputFileHandler.write(line + "\n")
            else:
                print >> sys.stdout, infoHeaderLines.rstrip()
                print >> sys.stdout, line
                
            '''
            # UNINDENT THIS: if we find the FORMAT section, then add the format from here
            elif ((not formatHeaderAdded) and line.startswith("##FORMAT")):
                formatHeaderAdded = True
                if (i_outputFilename != None):        
                    i_vcfOutputFileHandler.write(formatHeaderLines)
                    i_vcfOutputFileHandler.write(line + "\n")
                else:
                    print >> sys.stdout, formatHeaderLines.rstrip()
                    print >> sys.stdout, line
            '''
        
        # these lines are from previous scripts in the pipeline, so output them    
        elif (line.startswith("#")):
            if (i_outputFilename != None):
                i_vcfOutputFileHandler.write(line + "\n")
            else:
                print >> sys.stdout, line
                
        # now we are to the data
        else:        
            
            # split the line on the tab
            splitLine = re.split("\t", line)
    
            # the coordinate is the second element
            event_chr = splitLine[0]
            event_stopCoordinate = int(splitLine[1])
            
            # the key is the chrom number and the coordinate
            event_key = (event_chr + "_" + str(event_stopCoordinate))
            
            # check if there are any mutClass results for this coordinate
            if (event_key in mutClassDict):
                
                # get the info the format field data
                (mutClassInfo) = get_mutClass_data(mutClassDict[event_key], i_debug)
                event_infoOutput = splitLine[7] + ";" + mutClassInfo    
                
                '''
                start temp code
                
                infoTagList = re.split(";", event_infoList)
                for tag in infoTagList:
                    if (tag.startswith("MT")):
                        mtStartIndex = infoTagList.index(tag)
                        mcStopIndex = mtStartIndex
                    elif (tag.startswith("NS") and not tag.startswith("NSUB")):
                        mtStopIndex = infoTagList.index(tag)
                    elif (tag.startswith("MC")):
                        mcStartIndex = infoTagList.index(tag)
                
                tags = infoTagList[mtStartIndex:mtStopIndex]
                modTypesList = list()
                for tag in tags:
                    if ("=" in tag):
                        crap = re.split("=", tag)
                        modTypesList.append(crap[1])
                    else:
                        modTypesList.append(tag)
   
                tags = infoTagList[mcStartIndex:mcStopIndex]
                modChangesList = list()
                for tag in tags:
                    if ("=" in tag):
                        crap = re.split("=", tag)
                        modChangesList.append(crap[1])
                    else:
                        modChangesList.append(tag)
            
                event_infoOutput = ";".join(infoTagList[0:mcStartIndex])
                event_infoOutput += ";MC=" + ",".join(modChangesList) + ";MT=" + ",".join(modTypesList) + ";"
                event_infoOutput += ";".join(infoTagList[mtStopIndex:len(infoTagList)])
                
                
                stop temp code
                '''
                
                event_infoOutput = event_infoOutput.rstrip(";")
                event_infoList = re.split(";", event_infoOutput)
                finalOutput = ""
                for event in sorted(event_infoList):
                    finalOutput += event + ";"
                finalOutput = finalOutput.rstrip(";")
                
                
                # create the output list                        
                vcfOutputList = splitLine[0:7] + [finalOutput] + splitLine[8:len(splitLine)]
                            
                # output the new info
                if (i_outputFilename != None):
                    i_vcfOutputFileHandler.write("\t".join(vcfOutputList) + "\n")
                else:
                    print >> sys.stdout, "\t".join(vcfOutputList)
            # if there is no mutClass info to add
            else:
                if (i_outputFilename != None):
                    i_vcfOutputFileHandler.write(line + "\n")
                else:
                    print >> sys.stdout, line
        
    i_vcfInputFileHandler.close()
    i_mutClassFileHandler.close()
    if (i_outputFilename != None):
        i_vcfOutputFileHandler.close()
    return

main()
sys.exit(0)