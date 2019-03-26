#!/usr/bin/env python

import sys                          # system module
from optparse import OptionParser   # used for parsing command line arguments
import radiaUtil                    # utility functions for rna editing
import logging
import os
# import time
import subprocess


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


def filter_blacklist(aPythonExecutable, anId, aChromId, anInputFilename,
                     anOutputDir, aPrefix, aBlacklistDir, aScriptsDir,
                     anIgnoreScriptsDirFlag, aJobListFileHandler,
                     aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByPybed.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByPybed.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    filterFilename = os.path.join(aBlacklistDir, "chr" + aChromId + ".bed.gz")
    if (not os.path.isfile(filterFilename)):
        filterFilename = os.path.join(aBlacklistDir, "chr" + aChromId + ".bed")

    if (aGzipFlag):
        baseFilename = aPrefix + "_blacklist_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_blacklist_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + filterFilename
    cmd += " " + anInputFilename
    cmd += " blck --includeFilterName"
    cmd += (" -f \"##FILTER=<ID=blck,Description=\\\"Position overlaps " +
            "1000 Genomes Project blacklist\\\">\" -o " + outputFilename)

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("FilterFilename: %s", filterFilename)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename, filterFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def flag_radar(aPythonExecutable, anId, aChromId, anInputFilename, anOutputDir,
               aPrefix, aDbSnpDir, aScriptsDir, anIgnoreScriptsDirFlag,
               aJobListFileHandler, aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByPybed.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByPybed.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    filterFilename = os.path.join(aDbSnpDir, "chr" + aChromId + ".bed.gz")
    if (not os.path.isfile(filterFilename)):
        filterFilename = os.path.join(aDbSnpDir, "chr" + aChromId + ".bed")

    if (aGzipFlag):
        baseFilename = aPrefix + "_radar_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_radar_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + filterFilename
    cmd += " " + anInputFilename
    cmd += " RADAR"
    cmd += " --includeOverlaps"
    cmd += " --includeFilterName"
    cmd += " --includeIdName"
    cmd += " --idField INFO"
    cmd += " -d INFO"
    cmd += (" -f \"##INFO=<ID=RADAR,Number=.,Type=String," +
            "Description=\\\"Overlaps with Rigorously Annotated Database " +
            "of A-to-I RNA Editing (RADAR)\\\">\" -o " + outputFilename)

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("FilterFilename: %s", filterFilename)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename, filterFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def flag_darned(aPythonExecutable, anId, aChromId, anInputFilename,
                anOutputDir, aPrefix, aDbSnpDir, aScriptsDir,
                anIgnoreScriptsDirFlag, aJobListFileHandler,
                aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByPybed.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByPybed.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    filterFilename = os.path.join(aDbSnpDir, "chr" + aChromId + ".bed.gz")
    if (not os.path.isfile(filterFilename)):
        filterFilename = os.path.join(aDbSnpDir, "chr" + aChromId + ".bed")

    if (aGzipFlag):
        baseFilename = aPrefix + "_darned_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_darned_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + filterFilename
    cmd += " " + anInputFilename
    cmd += " DARNED "
    cmd += " --includeOverlaps"
    cmd += " --includeFilterName"
    cmd += " --includeIdName"
    cmd += " --idField INFO"
    cmd += " -d INFO"
    cmd += (" -f \"##INFO=<ID=DARNED,Number=.,Type=String," +
            "Description=\\\"Overlaps with Database of RNA Editing " +
            "(DARNED)\\\">\" -o " + outputFilename)

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("FilterFilename: %s", filterFilename)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename, filterFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def flag_dbSnp(aPythonExecutable, anId, aChromId, anInputFilename, anOutputDir,
               aPrefix, aDbSnpDir, aScriptsDir, anIgnoreScriptsDirFlag,
               aJobListFileHandler, aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByPybed.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByPybed.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    filterFilename = os.path.join(aDbSnpDir, "chr" + aChromId + ".bed.gz")
    if (not os.path.isfile(filterFilename)):
        filterFilename = os.path.join(aDbSnpDir, "chr" + aChromId + ".bed")

    if (aGzipFlag):
        baseFilename = aPrefix + "_dbsnp_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_dbsnp_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + filterFilename
    cmd += " " + anInputFilename
    cmd += " DB"
    cmd += " --includeOverlaps"
    cmd += " --includeFilterName"
    cmd += " --includeIdName"
    cmd += " -d INFO"
    cmd += (" -f \"##INFO=<ID=DB,Number=0,Type=Flag," +
            "Description=\\\"dbSNP common SNP membership\\\">\" -o " +
            outputFilename)

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("FilterFilename: %s", filterFilename)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename, filterFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def flag_retroGenes(aPythonExecutable, anId, aChromId, anInputFilename,
                    anOutputDir, aPrefix, aRetroGeneDir, aScriptsDir,
                    anIgnoreScriptsDirFlag, aJobListFileHandler,
                    aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByPybed.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByPybed.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    filterFilename = os.path.join(aRetroGeneDir, "chr" + aChromId + ".bed.gz")
    if (not os.path.isfile(filterFilename)):
        filterFilename = os.path.join(aRetroGeneDir, "chr" + aChromId + ".bed")

    if (aGzipFlag):
        baseFilename = aPrefix + "_retroGene_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_retroGene_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + filterFilename
    cmd += " " + anInputFilename
    cmd += " RTPS"
    cmd += " --includeOverlaps"
    cmd += " --includeFilterName"
    cmd += " -d INFO"
    cmd += (" -f \"##INFO=<ID=RTPS,Number=0,Type=Flag," +
            "Description=\\\"Overlaps with retrotransposon or " +
            "pseudogene\\\">\" -o " + outputFilename)

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("FilterFilename: %s", filterFilename)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename, filterFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def flag_pseudoGenes(aPythonExecutable, anId, aChromId, anInputFilename,
                     anOutputDir, aPrefix, aPseudoGeneDir, aScriptsDir,
                     anIgnoreScriptsDirFlag, aJobListFileHandler,
                     aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByPybed.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByPybed.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    filterFile = os.path.join(aPseudoGeneDir, "chr" + aChromId + ".bed.gz")
    if (not os.path.isfile(filterFile)):
        filterFile = os.path.join(aPseudoGeneDir, "chr" + aChromId + ".bed")

    if (aGzipFlag):
        baseFilename = aPrefix + "_pseudoGene_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_pseudoGene_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + filterFile
    cmd += " " + anInputFilename
    cmd += " EGPS"
    cmd += " --includeOverlaps"
    cmd += " --includeFilterName"
    cmd += " -d INFO"
    cmd += (" -f \"##INFO=<ID=EGPS,Number=0,Type=Flag," +
            "Description=\\\"Overlaps with ENCODE/GENCODE " +
            "pseudogenes\\\">\" -o " + outputFilename)

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("FilterFilename: %s", filterFile)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename, filterFile]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def flag_cosmic(aPythonExecutable, anId, aChromId, anInputFilename,
                anOutputDir, aPrefix, aCosmicDir, aScriptsDir,
                anIgnoreScriptsDirFlag, aJobListFileHandler,
                aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByPybed.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByPybed.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    filterFilename = os.path.join(aCosmicDir, "chr" + aChromId + ".bed.gz")
    if (not os.path.isfile(filterFilename)):
        filterFilename = os.path.join(aCosmicDir, "chr" + aChromId + ".bed")

    if (aGzipFlag):
        baseFilename = aPrefix + "_cosmic_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_cosmic_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + filterFilename
    cmd += " " + anInputFilename
    cmd += " COSMIC"
    cmd += " --includeOverlaps"
    cmd += " --includeFilterName"
    cmd += " --includeFilterCount"
    cmd += " -d INFO"
    cmd += (" -f \"##INFO=<ID=COSMIC,Number=1,Type=Integer," +
            "Description=\\\"Number of overlaps with the Catalogue Of " +
            "Somatic Mutations In Cancer (COSMIC)\\\">\" -o " + outputFilename)

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("FilterFilename: %s", filterFilename)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename, filterFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def filter_targets(aPythonExecutable, anId, aChromId, anInputFilename,
                   anOutputDir, aPrefix, aTargetDir, aScriptsDir,
                   anIgnoreScriptsDirFlag, aJobListFileHandler,
                   aGzipFlag, aTargetsInfoFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByPybed.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByPybed.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    filterFilename = os.path.join(aTargetDir, "chr" + aChromId + ".bed.gz")
    if (not os.path.isfile(filterFilename)):
        filterFilename = os.path.join(aTargetDir, "chr" + aChromId + ".bed")

    if (aGzipFlag):
        baseFilename = aPrefix + "_targets_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_targets_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + filterFilename
    cmd += " " + anInputFilename

    # if the target should be added to the INFO instead of the FILTER
    if (aTargetsInfoFlag):
        cmd += " NTR "
        cmd += " --includeFilterName"
        cmd += " -d INFO"
        cmd += (" -f \"##INFO=<ID=NTR,Number=0,Type=Flag," +
                "Description=\\\"Position does not overlap with a " +
                "GENCODE gene region\\\">\" -o " + outputFilename)
    else:
        cmd += " ntr"
        cmd += " --includeFilterName"
        cmd += (" -f \"##FILTER=<ID=ntr,Description=\\\"Position does not " +
                "overlap with a GENCODE gene region\\\">\" -o " +
                outputFilename)

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("FilterFilename: %s", filterFilename)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename, filterFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def filter_mpileupSupport_dna(aPythonExecutable, anId, aChromId,
                              anInputFilename, aHeaderFilename,
                              anOriginFlag, anOutputDir, aPrefix, aScriptsDir,
                              anIgnoreScriptsDirFlag, aJobListFileHandler,
                              aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByMpileupSupport.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByMpileupSupport.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    dnaParamList = ["--genotypeMinPct=0.10"]
    dnaParamList += ["--modMinDepth=4", "--modMinPct=0.10"]
    dnaParamList += ["--dnaNormalMinTotalBases=10", "--dnaNormalMinAltBases=4"]
    dnaParamList += ["--dnaNormalMinAltPct=0.10", "--dnaNormalMaxErrPct=0.01"]
    dnaParamList += ["--dnaTumorMinTotalBases=10", "--dnaTumorMinAltBases=4"]
    dnaParamList += ["--dnaTumorMinAltPct=0.10", "--dnaTumorMaxErrPct=0.01"]
    dnaParamString = " ".join(dnaParamList)

    if (anOriginFlag):
        if (aGzipFlag):
            baseFilename = (aPrefix + "_mpileup_dna_origin_chr" +
                            aChromId + ".vcf.gz")
        else:
            baseFilename = (aPrefix + "_mpileup_dna_origin_chr" +
                            aChromId + ".vcf")
    else:
        if (aGzipFlag):
            baseFilename = aPrefix + "_mpileup_dna_chr" + aChromId + ".vcf.gz"
        else:
            baseFilename = aPrefix + "_mpileup_dna_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + anInputFilename
    cmd += " " + dnaParamString
    cmd += " -o " + outputFilename

    if (anOriginFlag):
        cmd += " --addOrigin"
    if (aHeaderFilename is not None):
        cmd += " -n " + aHeaderFilename

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def filter_mpileupSupport_rna(aPythonExecutable, anId, aChromId,
                              anInputFilename, anOriginFlag, anRnaMinMapQual,
                              anRnaMinAvgMapQual, anOutputDir, aPrefix,
                              aScriptsDir, anIgnoreScriptsDirFlag,
                              aJobListFileHandler, aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByMpileupSupport.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByMpileupSupport.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    rnaParamList = ["--genotypeMinPct=0.0"]
    rnaParamList += ["--modMinDepth=1", "--modMinPct=0.01"]
    rnaParamList += ["--dnaNormalMinTotalBases=10", "--dnaNormalMinAltBases=0"]
    rnaParamList += ["--dnaNormalMinAltPct=0.0", "--dnaNormalMaxErrPct=1.0"]
    rnaParamList += ["--dnaNormalMinAltAvgBaseQual=15"]
    '''
    rnaParamList += ["--dnaNormalMinTotalBases=10", "--dnaNormalMinAltBases=0"]
    rnaParamList += ["--dnaNormalMinAltPct=0.0", "--dnaNormalMaxErrPct=0.01"]
    rnaParamList += ["--dnaNormalMinAltAvgBaseQual=15"]
    '''
    rnaParamList += ["--dnaTumorMinTotalBases=1", "--dnaTumorMinAltBases=1"]
    rnaParamList += ["--dnaTumorMinAltPct=0.01", "--dnaTumorMaxErrPct=1.0"]
    rnaParamList += ["--dnaTumorMinAltAvgBaseQual=15"]
    '''
    rnaParamList += ["--dnaTumorMinTotalBases=1", "--dnaTumorMinAltBases=1"]
    rnaParamList += ["--dnaTumorMinAltPct=0.01", "--dnaTumorMaxErrPct=0.01"]
    rnaParamList += ["--dnaTumorMinAltAvgBaseQual=15"]
    '''
    rnaParamList += ["--rnaNormalMinTotalBases=10", "--rnaNormalMinAltBases=4"]
    rnaParamList += ["--rnaNormalMinAltPct=0.10", "--rnaNormalMaxErrPct=0.01"]
    rnaParamList += ["--rnaNormalMinAltAvgBaseQual=15"]
    rnaParamList += ["--rnaNormalMinAltMapQual=" + str(anRnaMinMapQual)]
    rnaParamList += ["--rnaNormalMinAltAvgMapQual=" + str(anRnaMinAvgMapQual)]
    rnaParamList += ["--rnaTumorMinTotalBases=10", "--rnaTumorMinAltBases=4"]
    rnaParamList += ["--rnaTumorMinAltPct=0.10", "--rnaTumorMaxErrPct=0.01"]
    rnaParamList += ["--rnaTumorMinAltAvgBaseQual=15"]
    rnaParamList += ["--rnaTumorMinAltMapQual=" + str(anRnaMinMapQual)]
    rnaParamList += ["--rnaTumorMinAltAvgMapQual=" + str(anRnaMinAvgMapQual)]
    rnaParamString = " ".join(rnaParamList)

    if (anOriginFlag):
        if (aGzipFlag):
            baseFilename = (aPrefix + "_mpileup_rna_origin_chr" +
                            aChromId + ".vcf.gz")
        else:
            baseFilename = (aPrefix + "_mpileup_rna_origin_chr" +
                            aChromId + ".vcf")
    else:
        if (aGzipFlag):
            baseFilename = aPrefix + "_mpileup_rna_chr" + aChromId + ".vcf.gz"
        else:
            baseFilename = aPrefix + "_mpileup_rna_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + anInputFilename
    cmd += " --filterUsingRNA"
    cmd += " -o " + outputFilename
    cmd += " " + rnaParamString

    if (anOriginFlag):
        cmd += " --addOrigin"

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def radia_compare(aPythonExecutable, anId, aChromId, anRnaFilename,
                  aDnaFilename, anOutputDir, aPrefix, aScriptsDir,
                  anIgnoreScriptsDirFlag, aJobListFileHandler,
                  aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "radiaCompare.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "radiaCompare.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    if (aGzipFlag):
        overlapBaseFilename = aPrefix + "_overlap_chr" + aChromId + ".vcf.gz"
    else:
        overlapBaseFilename = aPrefix + "_overlap_chr" + aChromId + ".vcf"

    if (aGzipFlag):
        nonOverlapBaseFilename = (aPrefix + "_nonoverlap_chr" +
                                  aChromId + ".vcf.gz")
    else:
        nonOverlapBaseFilename = (aPrefix + "_nonoverlap_chr" +
                                  aChromId + ".vcf")

    overlapFilename = os.path.join(anOutputDir, overlapBaseFilename)
    nonOverlapFilename = os.path.join(anOutputDir, nonOverlapBaseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + anRnaFilename
    cmd += " " + aDnaFilename
    cmd += " -c \"SOM=SOM\""
    # cmd += (" -c \"SOM=SOM,TUM_EDIT=TUM_EDIT,NOR_EDIT=NOR_EDIT\"")
    # cmd += (" -c \"SOM=SOM,TUM_EDIT=TUM_EDIT,NOR_EDIT=NOR_EDIT,"
    #         "RNA_TUM_VAR=RNA_TUM_VAR,RNA_NOR_VAR=RNA_NOR_VAR\"")
    cmd += " -o " + overlapFilename
    cmd += " -n " + nonOverlapFilename

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("Input: %s %s", anRnaFilename, aDnaFilename)
        logging.debug("Output: %s %s", overlapFilename, nonOverlapFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anRnaFilename, aDnaFilename]
    writeFilenameList = [overlapFilename, nonOverlapFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return overlapFilename, nonOverlapFilename


def filter_rnaOnly(aChromId, anInputFilename, anOutputDir, aPrefix,
                   aJobListFileHandler, aGzipFlag, anIsDebug):

    # if we pipe the grep output to gzip, the return code and error messages
    # from grep get overwritten by gzip, and we no longer detect when there's
    # a problem with grep

    baseFilename = aPrefix + "_dnaFiltered_chr" + aChromId + ".vcf"
    outputFilename = os.path.join(anOutputDir, baseFilename)

    if (aGzipFlag):
        '''
        cmd = "zcat " + anInputFilename
        cmd += " | grep -v \"dnm\" "
        cmd += " | grep -v \"DB\""
        cmd += " | grep -v \"EGPS\""
        cmd += " | grep -v \"RTPS\""
        cmd += " | grep \"[SOM,EDIT]\" | gzip > " + outputFilename

        cmd = "zcat " + anInputFilename
        cmd += " | grep -v \"dnm\" "
        cmd += " | grep -v \"DB\""
        cmd += " | grep -v \"EGPS\""
        cmd += " | grep -v \"RTPS\""
        cmd += " | grep -v \"GERM\""
        cmd += " | awk '{if ($1 ~ /^#/) {print} else {print $1\"\\t\"$2\""
        cmd += "\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\tPASS\\t\"$8\"\\t\"$9\""
        cmd += "\\t\"$10\"\\t\"$11\"\\t\"$12}}' > " + outputFilename
        '''
        cmd = "zcat " + anInputFilename
        cmd += " | grep -v \"dnm\""
        cmd += " | grep -v \"DB\""
        cmd += " | grep -v \"EGPS\""
        cmd += " | grep -v \"RTPS\""
        cmd += " | grep -v \"GERM\""
        cmd += " | awk '{OFS=\"\\t\"}"
        cmd += " {if ($1 ~ /^#/) {print} else {$7=\"PASS\"; print}}'"
        cmd += " > " + outputFilename
    else:
        '''
        cmd = "grep -v \"dnm\" " + anInputFilename
        cmd += " | grep -v \"DB\""
        cmd += " | grep -v \"EGPS\""
        cmd += " | grep -v \"RTPS\""
        cmd += " | grep \"[SOM,EDIT]\" > " + outputFilename

        cmd = "grep -v \"dnm\" " + anInputFilename
        cmd += " | grep -v \"DB\""
        cmd += " | grep -v \"EGPS\""
        cmd += " | grep -v \"RTPS\""
        cmd += " | grep -v \"GERM\""
        cmd += " | awk '{if ($1 ~ /^#/) {print} else {print $1\"\\t\"$2\"" +
        cmd += "\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\tPASS\\t\"$8\"\\t\"$9\""
        cmd += "\\t\"$10\"\\t\"$11\"\\t\"$12}}' > " + outputFilename
        '''
        cmd = "grep -v \"dnm\" " + anInputFilename
        cmd += " | grep -v \"DB\""
        cmd += " | grep -v \"EGPS\""
        cmd += " | grep -v \"RTPS\""
        cmd += " | grep -v \"GERM\""
        cmd += " | awk '{OFS=\"\\t\"}"
        cmd += " {if ($1 ~ /^#/) {print} else {$7=\"PASS\"; print}}'"
        cmd += " > " + outputFilename

    if (anIsDebug):
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList = [anInputFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def extract_passing(aChromId, anInputFilename, anOutputDir, aPrefix,
                    aJobListFileHandler, aGzipFlag, anIsDebug):

    if (aGzipFlag):
        baseFilename = aPrefix + "_passing_chr" + aChromId + ".vcf"
        outputFilename = os.path.join(anOutputDir, baseFilename)
        cmd = "zcat " + anInputFilename
        cmd += " | grep -e \"PASS\" -e \"^#\" " + " > " + outputFilename
    else:
        baseFilename = aPrefix + "_passing_chr" + aChromId + ".vcf"
        outputFilename = os.path.join(anOutputDir, baseFilename)
        cmd = "grep -e \"PASS\" -e \"^#\" "
        cmd += anInputFilename + " > " + outputFilename

    if (anIsDebug):
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList = [anInputFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def run_snpEff(aChromId, anInputFilename, aSnpEffDir, aSnpEffGenome,
               aSnpEffCanonical, anOutputDir, aPrefix,
               aJobListFileHandler, anIsDebug):

    snpEffJar = os.path.join(aSnpEffDir, "snpEff.jar")
    snpEffConfig = os.path.join(aSnpEffDir, "snpEff.config")

    # by default, snpEff does not gzip the output
    # if we pipe the snpEff output to gzip, the return code and error messages
    # from snpEff get overwritten by gzip, and we no longer detect when
    # there's a problem with snpEff
    baseFilename = aPrefix + "_snpEff_chr" + aChromId + ".vcf"
    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = "java -Xmx4g"
    cmd += " -jar " + snpEffJar + " eff"
    cmd += " -c " + snpEffConfig
    cmd += " -formatEff -cancer"
    cmd += " -no-downstream"
    cmd += " -no-upstream"
    cmd += " -no-intergenic"
    cmd += " -no-intron"
    cmd += " " + aSnpEffGenome
    cmd += " " + anInputFilename
    cmd += " > " + outputFilename

    if (aSnpEffCanonical):
        cmd += " -canon"

    if (anIsDebug):
        logging.debug("Input: %s", anInputFilename)
        logging.debug("SnpEff: %s", snpEffJar)
        logging.debug("SnpEff: %s", snpEffConfig)
        logging.debug("SnpEff: %s", aSnpEffGenome)
        logging.debug("SnpEff: %s", aSnpEffCanonical)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList = [anInputFilename, snpEffJar, snpEffConfig]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def create_blat_input(aPythonExecutable, anId, aChromId, anInputFilename,
                      aTxNameTag, aTxCoordinateTag, aTxStrandTag,
                      anRnaIncludeSecondaryAlignmentsFlag, anOutputDir,
                      aPrefix, aScriptsDir, anIgnoreScriptsDirFlag,
                      aJobListFileHandler, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "createBlatFile.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "createBlatFile.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    # we can't gzip the blat input file
    baseFilename = aPrefix + "_blatInput_chr" + aChromId + ".fa"
    outputFilename = os.path.join(anOutputDir, baseFilename)

    # get the basic cmd
    cmd = script
    cmd += " " + anId
    cmd += " " + anInputFilename
    cmd += " -o " + outputFilename
    cmd += " --blatRnaNormalReads"
    cmd += " --blatRnaTumorReads"

    # if the transcript names, coordinates, and strands should be used
    if (aTxNameTag is not None and aTxCoordinateTag is not None):
        cmd += " --transcriptNameTag " + aTxNameTag
        cmd += " --transcriptCoordinateTag " + aTxCoordinateTag
        cmd += " --transcriptStrandTag " + aTxStrandTag
    # if the RNA secondary alignments should be included
    if (anRnaIncludeSecondaryAlignmentsFlag):
        cmd += " --rnaIncludeSecondaryAlignments"

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def run_blat(aChromId, aBlatInputFilename, aFastaFile, anOutputDir,
             aPrefix, aJobListFileHandler, anIsDebug):

    baseFilename = aPrefix + "_blatOutput_chr" + aChromId + ".blast"
    blatOutputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = "blat -stepSize=5 -repMatch=2253 -t=dna -q=rna"
    # cmd = "blat -t=dna -q=rna"
    # cmd += " -minScore=0 -minIdentity=0"
    cmd += " " + aFastaFile
    cmd += " " + aBlatInputFilename
    cmd += " -out=blast8 " + blatOutputFilename

    '''
    # default output from blat is PSL, but we're using -out=blast8 for now
    cmd = "blat -stepSize=5 -repMatch=2253 -t=dna -q=rna"
    cmd += " " + aFastaFile
    cmd += " " + aBlatInputFilename
    cmd += " " + blatOutputFilename
    '''

    if (anIsDebug):
        logging.debug("Input: %s", aBlatInputFilename)
        logging.debug("Output: %s", blatOutputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList = [aBlatInputFilename, aFastaFile]
    writeFilenameList = [blatOutputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return blatOutputFilename


def filter_blat(aPythonExecutable, anId, aChromId, anInputFilename,
                aHeaderFilename, aBlatInputFilename, aFastaFile,
                aTxNameTag, aTxCoordinateTag, aTxStrandTag,
                anRnaIncludeSecondaryAlignmentsFlag, anOutputDir, aPrefix,
                aScriptsDir, anIgnoreScriptsDirFlag, aJobListFileHandler,
                aGzipFlag, anIsDebug):

    # if no fasta file was specified, try to get it from the header file
    if (aFastaFile is None):
        fileHandler = radiaUtil.get_read_fileHandler(aHeaderFilename)

        for line in fileHandler:

            # strip the carriage return and newline characters
            line = line.rstrip("\r\n")

            # if we find the vcfGenerator line, then create the dict of params
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

                break

        fileHandler.close()

        if (("rnaTumorFastaFilename") in generatorParamsDict):
            aFastaFile = generatorParamsDict["rnaTumorFastaFilename"]

            if (not os.path.isfile(aFastaFile)):
                logging.critical("The FASTA file specified in the header " +
                                 "does not exist: %s", aFastaFile,
                                 ". Specify a FASTA file for the RNA " +
                                 "using the -f option.")
                sys.exit(1)

    blatOutputFilename = run_blat(aChromId,
                                  aBlatInputFilename,
                                  aFastaFile,
                                  anOutputDir,
                                  aPrefix,
                                  aJobListFileHandler,
                                  anIsDebug)

    if (aGzipFlag):
        baseFilename = aPrefix + "_blatFiltered_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_blatFiltered_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByBlat.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByBlat.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    # get the basic cmd
    cmd = script
    cmd += " " + anId
    cmd += " " + anInputFilename
    cmd += " " + blatOutputFilename
    cmd += " -o " + outputFilename
    cmd += " --blatRnaNormalReads"
    cmd += " --blatRnaTumorReads"

    # if the transcript names, coordinates, and strands should be used
    if (aTxNameTag is not None and aTxCoordinateTag is not None):
        cmd += " --transcriptNameTag " + aTxNameTag
        cmd += " --transcriptCoordinateTag " + aTxCoordinateTag
        cmd += " --transcriptStrandTag " + aTxStrandTag
    # if the RNA secondary alignments should be included
    if (anRnaIncludeSecondaryAlignmentsFlag):
        cmd += " --rnaIncludeSecondaryAlignments"

    if (anIsDebug):
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return (blatOutputFilename, outputFilename)


def filter_rnaBlacklist(aPythonExecutable, aChromId, anInputFilename,
                        aGeneBlckFilename, aGeneFamilyBlckFilename,
                        anOutputDir, aPrefix, aScriptsDir,
                        anIgnoreScriptsDirFlag, aJobListFileHandler,
                        aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByRnaBlacklist.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByRnaBlacklist.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    if (aGzipFlag):
        baseFilename = aPrefix + "_rna_genes_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_rna_genes_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anInputFilename
    cmd += " " + aGeneBlckFilename
    cmd += " " + aGeneFamilyBlckFilename
    cmd += " -o " + outputFilename

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Input: %s", aGeneBlckFilename)
        logging.debug("Input: %s", aGeneFamilyBlckFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename,
                         aGeneBlckFilename,
                         aGeneFamilyBlckFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def merge_rnaAndDna(aPythonExecutable, anId, aChromId, aDnaFilename,
                    anRnaFilename, anOverlapsFilname, aNonoverlapsFilename,
                    aDnaHeaderOnlyFlag, anOutputDir, aPrefix, aScriptsDir,
                    anIgnoreScriptsDirFlag, aJobListFileHandler,
                    aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "mergeRnaAndDnaFiles.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "mergeRnaAndDnaFiles.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    if (aGzipFlag):
        baseFilename = aPrefix + "_merged_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_merged_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anId
    cmd += " " + aChromId
    cmd += " " + aDnaFilename
    cmd += " " + anRnaFilename
    cmd += " " + anOverlapsFilname
    cmd += " " + aNonoverlapsFilename
    cmd += " " + outputFilename

    if (aDnaHeaderOnlyFlag):
        cmd += " --dnaHeaderOnly"

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("aDnaFilename: %s", aDnaFilename)
        logging.debug("anOverlapsFilname: %s", anOverlapsFilname)
        logging.debug("aNonoverlapsFilename: %s", aNonoverlapsFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [aDnaFilename, anOverlapsFilname]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def merge_passingAndOriginals(aPythonExecutable, aChromId,
                              aPassingCallsFilename, anOriginalFilname,
                              anOutputDir, aPrefix, aScriptsDir,
                              anIgnoreScriptsDirFlag, aJobListFileHandler,
                              aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "mergePassingAndOriginals.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "mergePassingAndOriginals.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    if (aGzipFlag):
        baseFilename = aPrefix + "_mergedFinal_chr" + aChromId + ".vcf.gz"
    else:
        baseFilename = aPrefix + "_mergedFinal_chr" + aChromId + ".vcf"

    outputFilename = os.path.join(anOutputDir, baseFilename)
    cmd = script
    cmd += " " + aPassingCallsFilename
    cmd += " " + anOriginalFilname
    cmd += " " + outputFilename

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("Input: %s %s", aPassingCallsFilename, anOriginalFilname)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [aPassingCallsFilename, anOriginalFilname]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def filter_readSupport(aPythonExecutable, aChromId, anInputFilename,
                       aTxNameTag, aTxCoordinateTag, aTxStrandTag,
                       anRnaIncludeSecondaryAlignmentsFlag, aMinMapQual,
                       anOutputDir, aPrefix, anOutputFilename, aScriptsDir,
                       anIgnoreScriptsDirFlag, aJobListFileHandler,
                       aGzipFlag, anIsDebug):

    readFilenameList = []
    if (anIgnoreScriptsDirFlag):
        script = "filterByReadSupport.py"
    else:
        pythonScript = os.path.join(aScriptsDir, "filterByReadSupport.py")
        readFilenameList.append(pythonScript)
        script = aPythonExecutable + " " + pythonScript

    if (anOutputFilename is not None):
        outputFilename = anOutputFilename
    else:
        if (aGzipFlag):
            baseFilename = aPrefix + "_chr" + aChromId + ".vcf.gz"
        else:
            baseFilename = aPrefix + "_chr" + aChromId + ".vcf"

        outputFilename = os.path.join(anOutputDir, baseFilename)

    cmd = script
    cmd += " " + anInputFilename
    cmd += " -o " + outputFilename
    cmd += " --minMapQual " + str(aMinMapQual)

    # if the transcript names, coordinates, and strands should be used
    if (aTxNameTag is not None and aTxCoordinateTag is not None):
        cmd += " --transcriptNameTag " + aTxNameTag
        cmd += " --transcriptCoordinateTag " + aTxCoordinateTag
        cmd += " --transcriptStrandTag " + aTxStrandTag
    # if the RNA secondary alignments should be included
    if (anRnaIncludeSecondaryAlignmentsFlag):
        cmd += " --rnaIncludeSecondaryAlignments"

    if (anIsDebug):
        logging.debug("Script: %s", script)
        logging.debug("Input: %s", anInputFilename)
        logging.debug("Output: %s", outputFilename)
        logging.debug("Filter: %s", cmd)

    readFilenameList += [anInputFilename]
    writeFilenameList = [outputFilename]

    run_cmd(cmd, None, readFilenameList, writeFilenameList,
            aJobListFileHandler, anIsDebug)

    return outputFilename


def remove_tmpFiles(aRmTmpFilesList, aJobListFileHandler, anIsDebug):

    finalList = list()
    for tmpFile in aRmTmpFilesList:
        if (os.path.exists(tmpFile)):
            finalList.append("rm " + tmpFile)

    cmd = ";".join(finalList)

    run_cmd(cmd, None, None, None, aJobListFileHandler, anIsDebug)

    return


def run_cmd(aCommand, aDirList, aReadFileList, aWriteFileList,
            aJobListFileHandler, anIsDebug):

    if (not radiaUtil.check_for_argv_errors(aDirList,
                                            aReadFileList,
                                            aWriteFileList)):
        sys.exit(1)

    if (anIsDebug):
        logging.debug("Command: %s", aCommand)

    if (aJobListFileHandler is not None):
        aJobListFileHandler.write(aCommand + "\n")
    else:
        subprocessCall = subprocess.Popen(aCommand,
                                          shell=True,
                                          bufsize=-1,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE,
                                          close_fds=True)
        (stdOut, stdErr) = subprocessCall.communicate()
        if (subprocessCall.returncode != 0):
            logging.error("The return code of '%s' from the following " +
                          "filter command indicates an error.",
                          subprocessCall.returncode)
            logging.error("Error from %s:\n%s", aCommand, stdErr)

    return


def main():

    # python filterRadia.py TCGA-AB-2995 12
    # ../data/test/TCGA-AB-2995.vcf
    # ../data/test/
    # ../scripts/

    # create the usage statement
    usage = "usage: python %prog id chrom inputFile outputDir scriptDir [Opts]"
    i_cmdLineParser = OptionParser(usage=usage)

    # a,b,c,d,e,f,g,l,n,p,r,s,t,

    i_cmdLineParser.add_option(
        "-b", "--blacklistDir",
        dest="blacklistDir", metavar="BLACKLIST_DIR",
        help="the path to the blacklist directory")
    i_cmdLineParser.add_option(
        "-t", "--targetDir",
        dest="targetDir", metavar="TARGET_DIR",
        help="the path to the exon capture targets directory")
    i_cmdLineParser.add_option(
        "-d", "--dbSnpDir",
        dest="dbSnpDir", metavar="SNP_DIR",
        help="the path to the dbSNP directory")
    i_cmdLineParser.add_option(
        "-r", "--retroGenesDir",
        dest="retroGenesDir", metavar="RETRO_DIR",
        help="the path to the retrogenes directory")
    i_cmdLineParser.add_option(
        "-p", "--pseudoGenesDir",
        dest="pseudoGenesDir", metavar="PSEUDO_DIR",
        help="the path to the pseudogenes directory")
    i_cmdLineParser.add_option(
        "-c", "--cosmicDir",
        dest="cosmicDir", metavar="COSMIC_DIR",
        help="the path to the cosmic directory")
    i_cmdLineParser.add_option(
        "-a", "--radarDir",
        dest="radarDir", metavar="RADAR_DIR",
        help="the path to the radar directory")
    i_cmdLineParser.add_option(
        "-n", "--darnedDir",
        dest="darnedDir", metavar="DARNED_DIR",
        help="the path to the darned directory")
    i_cmdLineParser.add_option(
        "-s", "--snpEffDir",
        dest="snpEffDir", metavar="SNP_EFF_DIR",
        help="the path to the snpEff directory")
    i_cmdLineParser.add_option(
        "-e", "--snpEffGenome",
        dest="snpEffGenome", default="GRCh37.75", metavar="SNP_EFF_GENOME",
        help="the snpEff Genome, %default by default")
    i_cmdLineParser.add_option(
        "", "--canonical",
        action="store_true", default=False,
        dest="canonical", metavar="CANONICAL",
        help="include this argument if only the canonical transcripts from " +
             "snpEff should be used, %default by default")
    i_cmdLineParser.add_option(
        "-f", "--blatFastaFilename",
        dest="blatFastaFilename", metavar="FASTA_FILE",
        help="the fasta file that can be used during the BLAT filtering, " +
             "default is the one specified in the VCF header")
    i_cmdLineParser.add_option(
        "-o", "--outputFilename",
        dest="outputFilename", metavar="OUTPUT_FILE",
        help="the name of the output file, otherwise a file will be " +
             "automatically created in the outputDir with the following " +
             "format:  patientId + '_chr' + chrom + '.vcf')")

    i_cmdLineParser.add_option(
        "", "--rnaGeneBlckFile",
        dest="rnaGeneBlckFile", metavar="RNA_GENE_FILE",
        help="the RNA gene blacklist file")
    i_cmdLineParser.add_option(
        "", "--rnaGeneFamilyBlckFile",
        dest="rnaGeneFamilyBlckFile", metavar="RNA_GENE_FAMILY_FILE",
        help="the RNA gene family blacklist file")
    i_cmdLineParser.add_option(
        "", "--prefix",
        dest="prefix", metavar="UNIQUE_FILE_PREFIX", default=None,
        help="a prefix to be added to all temp and output files to ensure " +
             "they are unique, otherwise all files will be automatically " +
             "created in the outputDir with the following format:  " +
             "patientId + '_chr' + chrom + '.vcf'")

    # we do all filtering by default, so it's better for the user to specify
    # --no flags to disable some filters, but internally, the code is nicer
    # if we can avoid the double negatives, so store true by default and
    # drop the "no" in the flag name
    i_cmdLineParser.add_option(
        "", "--noBlacklist",
        dest="blacklist", action="store_false", default=True,
        help="include this argument if the accessibility mask info/filter " +
             "should not be applied")
    i_cmdLineParser.add_option(
        "", "--noTargets",
        dest="targets", action="store_false", default=True,
        help="include this argument if the target info/filter " +
             "should not be applied")
    i_cmdLineParser.add_option(
        "", "--noDbSnp",
        dest="dbSnp", action="store_false", default=True,
        help="include this argument if the dbSNP info/filter " +
             "should not be applied")
    i_cmdLineParser.add_option(
        "", "--noRetroGenes",
        dest="retroGenes", action="store_false", default=True,
        help="include this argument if the retrogene info/filter " +
             "should not be applied")
    i_cmdLineParser.add_option(
        "", "--noPseudoGenes",
        dest="pseudoGenes", action="store_false", default=True,
        help="include this argument if the pseudogene info/filter " +
             "should not be applied")
    i_cmdLineParser.add_option(
        "", "--noCosmic",
        dest="cosmic", action="store_false", default=True,
        help="include this argument if the cosmic annotation " +
             "should not be applied")
    i_cmdLineParser.add_option(
        "", "--noRadar",
        dest="radar", action="store_false", default=True,
        help="include this argument if the RADAR info/filter " +
             "should not be applied")
    i_cmdLineParser.add_option(
        "", "--noDarned",
        dest="darned", action="store_false", default=True,
        help="include this argument if the DARNED info/filter " +
             "should not be applied")
    i_cmdLineParser.add_option(
        "", "--noBlat",
        dest="blat", action="store_false", default=True,
        help="include this argument if the blat filter should not be applied")
    i_cmdLineParser.add_option(
        "", "--noRnaBlacklist",
        dest="rnaBlacklist", action="store_false", default=True,
        help="include this argument if the RNA exclude gene filter " +
             "should not be applied")
    i_cmdLineParser.add_option(
        "", "--noSnpEff",
        dest="snpEff", action="store_false", default=True,
        help="include this argument if the snpEff annotation should not be " +
             "applied (without the snpEff annotation, filtering of RNA "
             "excluded genes will also not be applied")
    i_cmdLineParser.add_option(
        "", "--targetsInfo",
        dest="targetsInfo", action="store_true", default=False,
        help="include this argument if the targets should be added " +
             "to the INFO instead of the FILTER")
    i_cmdLineParser.add_option(
        "", "--ignoreScriptsDir",
        dest="ignoreScriptsDir", action="store_true", default=False,
        help="include this argument if the scriptsDir should be " +
             "ignored (e.g. when RADIA is installed via a wheel")
    i_cmdLineParser.add_option(
        "", "--dnaOnly",
        dest="dnaOnly", action="store_true", default=False,
        help="include this argument if you only have DNA or " +
             "filtering should only be done on the DNA")
    i_cmdLineParser.add_option(
        "", "--rnaOnly",
        dest="rnaOnly", action="store_true", default=False,
        help="include this argument if the filtering should " +
             "only be done on the RNA")
    i_cmdLineParser.add_option(
        "", "--gzip",
        dest="gzip", action="store_true", default=False,
        help="include this argument if the VCF should be compressed with gzip")
    i_cmdLineParser.add_option(
        "", "--transcriptNameTag",
        dest="transcriptNameTag", metavar="TX_NAME_TAG",
        help="the INFO key where the original transcript name can be found")
    i_cmdLineParser.add_option(
        "", "--transcriptCoordinateTag",
        dest="transcriptCoordinateTag", metavar="TX_COORDINATE_TAG",
        help="the INFO key where the original transcript coordinate " +
             "can be found")
    i_cmdLineParser.add_option(
        "", "--transcriptStrandTag",
        dest="transcriptStrandTag", metavar="TX_STRAND_TAG",
        help="the INFO key where the original transcript strand can be found")
    i_cmdLineParser.add_option(
        "", "--rnaIncludeSecondaryAlignments",
        action="store_true", default=False,
        dest="rnaIncludeSecondaryAlignments",
        help="if you align the RNA to transcript isoforms, " +
             "then you may want to include RNA secondary " +
             "alignments in the samtools mpileups")
    i_cmdLineParser.add_option(
        "", "--readSupportMinMapQual",
        type="int", default=int(10),
        dest="readSupportMinMapQual", metavar="READ_SUPPORT_MIN_MAP_QUAL",
        help="the minimum mapping quality for reads supporting the ALT, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--rnaMpileupMinMapQual",
        type="int", default=int(15),
        dest="rnaMpileupMinMapQual", metavar="RNA_MPILEUP_MIN_MAP_QUAL",
        help="at least 1 ALT read needs this minimum mapping quality, " +
             "%default by default")
    i_cmdLineParser.add_option(
        "", "--rnaMpileupMinAvgMapQual",
        type="int", default=int(20),
        dest="rnaMpileupMinAvgMapQual",
        metavar="RNA_MPILEUP_SUPPORT_MIN_MAP_QUAL",
        help="the minimum average mapping quality for the ALT reads, " +
             "%default by default")

    i_cmdLineParser.add_option(
        "-l", "--log",
        dest="logLevel", default="WARNING", metavar="LOG",
        help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), " +
             "%default by default")
    i_cmdLineParser.add_option(
        "-g", "--logFilename",
        dest="logFilename", metavar="LOG_FILE",
        help="the name of the log file, STDOUT by default")

    # range(inclusiveFrom, exclusiveTo, by)
    i_possibleArgLengths = range(5, 69, 1)
    i_argLength = len(sys.argv)

    # check if this is one of the possible correct commands
    if (i_argLength not in i_possibleArgLengths):
        i_cmdLineParser.print_help()
        sys.exit(1)

    # try to use the exact python version that the user
    # specifies on the command line
    if (sys.executable is not None and sys.executable != ""):
        i_pythonExecutable = sys.executable
    else:
        i_pythonExecutable = "python"

    # get the required parameters
    (cmdLineOpts, cmdLineArgs) = i_cmdLineParser.parse_args()
    i_id = str(cmdLineArgs[0])
    i_chr = str(cmdLineArgs[1])
    i_inputFilename = str(cmdLineArgs[2])
    i_outputDir = str(cmdLineArgs[3])
    i_scriptsDir = str(cmdLineArgs[4])

    # get the optional params with default values
    i_blacklistFlag = cmdLineOpts.blacklist
    i_rnaBlacklistFlag = cmdLineOpts.rnaBlacklist
    i_targetsFlag = cmdLineOpts.targets
    i_dbSnpFlag = cmdLineOpts.dbSnp
    i_retroGenesFlag = cmdLineOpts.retroGenes
    i_pseudoGenesFlag = cmdLineOpts.pseudoGenes
    i_cosmicFlag = cmdLineOpts.cosmic
    i_radarFlag = cmdLineOpts.radar
    i_darnedFlag = cmdLineOpts.darned
    i_blatFlag = cmdLineOpts.blat
    i_snpEffFlag = cmdLineOpts.snpEff
    i_targetsInfo = cmdLineOpts.targetsInfo
    i_ignoreScriptsDir = cmdLineOpts.ignoreScriptsDir
    i_dnaOnlyFlag = cmdLineOpts.dnaOnly
    i_rnaOnlyFlag = cmdLineOpts.rnaOnly
    i_logLevel = cmdLineOpts.logLevel
    i_gzip = cmdLineOpts.gzip
    i_snpEffGenome = cmdLineOpts.snpEffGenome
    i_snpEffCanonical = cmdLineOpts.canonical
    i_rnaIncludeSecAligns = cmdLineOpts.rnaIncludeSecondaryAlignments
    i_readSupportMinMapQual = cmdLineOpts.readSupportMinMapQual
    i_rnaMpileupMinMapQual = cmdLineOpts.rnaMpileupMinMapQual
    i_rnaMpileupMinAvgMapQual = cmdLineOpts.rnaMpileupMinAvgMapQual

    # try to get any optional parameters with no defaults
    i_prefix = i_id
    i_outputFilename = None
    i_blacklistDir = None
    i_targetDir = None
    i_dbSnpDir = None
    i_retroGenesDir = None
    i_pseudoGenesDir = None
    i_cosmicDir = None
    i_radarDir = None
    i_darnedDir = None
    i_joblistDir = None
    i_shebang = None
    i_logFilename = None
    i_blatFastaFilename = None
    i_snpEffDir = None
    i_rnaGeneBlckFilename = None
    i_rnaGeneFamilyBlckFilename = None
    i_transcriptNameTag = None
    i_transcriptCoordinateTag = None
    i_transcriptStrandTag = None
    readFilenameList = [i_inputFilename]
    writeFilenameList = [i_outputDir]
    dirList = [i_scriptsDir]
    if (cmdLineOpts.blacklistDir is not None):
        i_blacklistDir = str(cmdLineOpts.blacklistDir)
        dirList += [i_blacklistDir]
    if (cmdLineOpts.targetDir is not None):
        i_targetDir = str(cmdLineOpts.targetDir)
        dirList += [i_targetDir]
    if (cmdLineOpts.dbSnpDir is not None):
        i_dbSnpDir = str(cmdLineOpts.dbSnpDir)
        dirList += [i_dbSnpDir]
    if (cmdLineOpts.retroGenesDir is not None):
        i_retroGenesDir = str(cmdLineOpts.retroGenesDir)
        dirList += [i_retroGenesDir]
    if (cmdLineOpts.pseudoGenesDir is not None):
        i_pseudoGenesDir = str(cmdLineOpts.pseudoGenesDir)
        dirList += [i_pseudoGenesDir]
    if (cmdLineOpts.cosmicDir is not None):
        i_cosmicDir = str(cmdLineOpts.cosmicDir)
        dirList += [i_cosmicDir]
    if (cmdLineOpts.radarDir is not None):
        i_radarDir = str(cmdLineOpts.radarDir)
        dirList += [i_radarDir]
    if (cmdLineOpts.darnedDir is not None):
        i_darnedDir = str(cmdLineOpts.darnedDir)
        dirList += [i_radarDir]
    # if (cmdLineOpts.joblistDir is not None):
    #    i_joblistDir = str(cmdLineOpts.joblistDir)
    #    dirList += [i_joblistDir]
    if (cmdLineOpts.snpEffDir is not None):
        i_snpEffDir = str(cmdLineOpts.snpEffDir)
        dirList += [i_snpEffDir]
    # if (cmdLineOpts.shebang is not None):
    #    i_shebang = str(cmdLineOpts.shebang)
    if (cmdLineOpts.outputFilename is not None):
        i_outputFilename = str(cmdLineOpts.outputFilename)
        writeFilenameList += [i_outputFilename]
    if (cmdLineOpts.logFilename is not None):
        i_logFilename = str(cmdLineOpts.logFilename)
        writeFilenameList += [i_logFilename]
    if (cmdLineOpts.blatFastaFilename is not None):
        i_blatFastaFilename = str(cmdLineOpts.blatFastaFilename)
        writeFilenameList += [i_blatFastaFilename]
    if (cmdLineOpts.rnaGeneBlckFile is not None):
        i_rnaGeneBlckFilename = str(cmdLineOpts.rnaGeneBlckFile)
        readFilenameList += [i_rnaGeneBlckFilename]
    if (cmdLineOpts.rnaGeneFamilyBlckFile is not None):
        i_rnaGeneFamilyBlckFilename = str(cmdLineOpts.rnaGeneFamilyBlckFile)
        readFilenameList += [i_rnaGeneFamilyBlckFilename]
    if (cmdLineOpts.transcriptNameTag is not None):
        i_transcriptNameTag = cmdLineOpts.transcriptNameTag
    if (cmdLineOpts.transcriptCoordinateTag is not None):
        i_transcriptCoordinateTag = cmdLineOpts.transcriptCoordinateTag
    if (cmdLineOpts.transcriptStrandTag is not None):
        i_transcriptStrandTag = cmdLineOpts.transcriptStrandTag
    if (cmdLineOpts.prefix is not None):
        i_prefix = cmdLineOpts.prefix

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
        logging.debug("pythonExecutable=%s", i_pythonExecutable)
        logging.debug("id=%s", i_id)
        logging.debug("chr=%s", i_chr)
        logging.debug("inputFilename=%s", i_inputFilename)
        logging.debug("outputFilename=%s", i_outputFilename)
        logging.debug("outputDir=%s", i_outputDir)
        logging.debug("scriptsDir=%s", i_scriptsDir)
        logging.debug("logLevel=%s", i_logLevel)
        logging.debug("gzip=%s", i_gzip)
        logging.debug("logFile=%s", i_logFilename)
        logging.debug("prefix=%s", i_prefix)
        logging.debug("blatfastaFile=%s", i_blatFastaFilename)
        logging.debug("rnaGeneExcludeFile=%s", i_rnaGeneBlckFilename)
        logging.debug("rnaGeneFamExcludeFile=%s", i_rnaGeneFamilyBlckFilename)
        logging.debug("blacklistFlag? %s", i_blacklistFlag)
        logging.debug("rnaBlacklistFlag? %s", i_rnaBlacklistFlag)
        logging.debug("targetsFlag? %s", i_targetsFlag)
        logging.debug("targetsInfo? %s", i_targetsInfo)
        logging.debug("ignoreScriptsDir? %s", i_ignoreScriptsDir)
        logging.debug("dbSnpFlag? %s", i_dbSnpFlag)
        logging.debug("snpEffFlag? %s", i_snpEffFlag)
        logging.debug("retroGenesFlag? %s", i_retroGenesFlag)
        logging.debug("pseudoGenesFlag? %s", i_pseudoGenesFlag)
        logging.debug("cosmicFlag? %s", i_cosmicFlag)
        logging.debug("radarFlag? %s", i_radarFlag)
        logging.debug("darnedFlag? %s", i_darnedFlag)
        logging.debug("blatFlag? %s", i_blatFlag)
        logging.debug("dnaOnlyFlag? %s", i_dnaOnlyFlag)
        logging.debug("rnaOnlyFlag? %s", i_rnaOnlyFlag)
        logging.debug("blacklistDir %s", i_blacklistDir)
        logging.debug("targetDir %s", i_targetDir)
        logging.debug("dbSnpDir %s", i_dbSnpDir)
        logging.debug("snpEffDir %s", i_snpEffDir)
        logging.debug("snpEffGenome %s", i_snpEffGenome)
        logging.debug("snpEffCanonical %s", i_snpEffCanonical)
        logging.debug("retroGenesDir %s", i_retroGenesDir)
        logging.debug("pseudoGenesDir %s", i_pseudoGenesDir)
        logging.debug("cosmicDir %s", i_cosmicDir)
        logging.debug("radarDir %s", i_radarDir)
        logging.debug("darnedDir %s", i_darnedDir)
        logging.debug("joblistDir %s", i_joblistDir)
        logging.debug("transcriptNameTag %s", i_transcriptNameTag)
        logging.debug("transcriptCoordinateTag %s", i_transcriptCoordinateTag)
        logging.debug("transcriptStrandTag %s", i_transcriptStrandTag)
        logging.debug("rnaInclSecAlign=%s" % i_rnaIncludeSecAligns)
        logging.debug("readSupportMinMapQual=%s" % i_readSupportMinMapQual)
        logging.debug("rnaMpileupMinMapQual=%s" % i_rnaMpileupMinMapQual)
        logging.debug("rnaMpileupMinAvgMapQual=%s" % i_rnaMpileupMinAvgMapQual)

    if (i_dnaOnlyFlag):
        i_rnaBlacklistFlag = False

    if (i_blacklistFlag):
        if (i_blacklistDir is None):
            logging.critical("No blacklist directory has been specified.")
            sys.exit(1)

    if (i_dbSnpFlag):
        if (i_dbSnpDir is None):
            logging.critical("No dbSNP directory has been specified.")
            sys.exit(1)

    if (i_retroGenesFlag):
        if (i_retroGenesDir is None):
            logging.critical("No retrogenes directory has been specified.")
            sys.exit(1)

    if (i_pseudoGenesFlag):
        if (i_pseudoGenesDir is None):
            logging.critical("No pseudogenes directory has been specified.")
            sys.exit(1)

    if (i_cosmicFlag):
        if (i_cosmicDir is None):
            logging.critical("No COSMIC directory has been specified.")
            sys.exit(1)

    if (i_radarFlag):
        if (i_radarDir is None):
            logging.critical("No RADAR directory has been specified.")
            sys.exit(1)

    if (i_darnedFlag):
        if (i_darnedDir is None):
            logging.critical("No DARNED directory has been specified.")
            sys.exit(1)

    if (i_targetsFlag):
        if (i_targetDir is None):
            logging.critical("No target directory has been specified.")
            sys.exit(1)

    if (i_snpEffFlag):
        if (i_snpEffDir is None):
            logging.critical("No snpEff directory has been specified.")
            sys.exit(1)

    if (i_rnaBlacklistFlag):
        if (i_rnaGeneBlckFilename is None):
            logging.critical("No RNA gene exclude list has been specified.")
            sys.exit(1)
        if (i_rnaGeneFamilyBlckFilename is None):
            logging.critical("No RNA gene family exclude list " +
                             "has been specified.")
            sys.exit(1)

    # check to see if the files exist
    if (not radiaUtil.check_for_argv_errors(dirList,
                                            readFilenameList,
                                            writeFilenameList)):
        sys.exit(1)

    i_joblistFileHandler = None
    if (i_joblistDir is not None):
        joblistFile = os.path.join(i_joblistDir, i_id + "_chr" + i_chr + ".sh")
        i_joblistFileHandler = radiaUtil.get_write_fileHandler(joblistFile)
        if (i_shebang is not None):
            i_joblistFileHandler.write(i_shebang + "\n")

    previousFilename = i_inputFilename
    rmTmpFilesList = list()

    if (i_dnaOnlyFlag):
        # filter by blacklist
        if (i_blacklistFlag):
            previousFilename = filter_blacklist(i_pythonExecutable,
                                                i_id,
                                                i_chr,
                                                previousFilename,
                                                i_outputDir,
                                                i_prefix,
                                                i_blacklistDir,
                                                i_scriptsDir,
                                                i_ignoreScriptsDir,
                                                i_joblistFileHandler,
                                                i_gzip,
                                                i_debug)
            rmTmpFilesList.append(previousFilename)

        # flag snp
        if (i_dbSnpFlag):
            previousFilename = flag_dbSnp(i_pythonExecutable,
                                          i_id,
                                          i_chr,
                                          previousFilename,
                                          i_outputDir,
                                          i_prefix,
                                          i_dbSnpDir,
                                          i_scriptsDir,
                                          i_ignoreScriptsDir,
                                          i_joblistFileHandler,
                                          i_gzip,
                                          i_debug)
            rmTmpFilesList.append(previousFilename)

        # flag retro genes
        if (i_retroGenesFlag):
            previousFilename = flag_retroGenes(i_pythonExecutable,
                                               i_id,
                                               i_chr,
                                               previousFilename,
                                               i_outputDir,
                                               i_prefix,
                                               i_retroGenesDir,
                                               i_scriptsDir,
                                               i_ignoreScriptsDir,
                                               i_joblistFileHandler,
                                               i_gzip,
                                               i_debug)
            rmTmpFilesList.append(previousFilename)

        # flag pseudo genes
        if (i_pseudoGenesFlag):
            previousFilename = flag_pseudoGenes(i_pythonExecutable,
                                                i_id,
                                                i_chr,
                                                previousFilename,
                                                i_outputDir,
                                                i_prefix,
                                                i_pseudoGenesDir,
                                                i_scriptsDir,
                                                i_ignoreScriptsDir,
                                                i_joblistFileHandler,
                                                i_gzip,
                                                i_debug)
            rmTmpFilesList.append(previousFilename)

        # flag cosmic
        if (i_cosmicFlag):
            previousFilename = flag_cosmic(i_pythonExecutable,
                                           i_id,
                                           i_chr,
                                           previousFilename,
                                           i_outputDir,
                                           i_prefix,
                                           i_cosmicDir,
                                           i_scriptsDir,
                                           i_ignoreScriptsDir,
                                           i_joblistFileHandler,
                                           i_gzip,
                                           i_debug)
            rmTmpFilesList.append(previousFilename)

        # flag radar
        if (i_radarFlag):
            previousFilename = flag_radar(i_pythonExecutable,
                                          i_id,
                                          i_chr,
                                          previousFilename,
                                          i_outputDir,
                                          i_prefix,
                                          i_radarDir,
                                          i_scriptsDir,
                                          i_ignoreScriptsDir,
                                          i_joblistFileHandler,
                                          i_gzip,
                                          i_debug)
            rmTmpFilesList.append(previousFilename)

        # flag darned
        if (i_darnedFlag):
            previousFilename = flag_darned(i_pythonExecutable,
                                           i_id,
                                           i_chr,
                                           previousFilename,
                                           i_outputDir,
                                           i_prefix,
                                           i_darnedDir,
                                           i_scriptsDir,
                                           i_ignoreScriptsDir,
                                           i_joblistFileHandler,
                                           i_gzip,
                                           i_debug)
            rmTmpFilesList.append(previousFilename)

        # filter targets
        if (i_targetsFlag):
            previousFilename = filter_targets(i_pythonExecutable,
                                              i_id,
                                              i_chr,
                                              previousFilename,
                                              i_outputDir,
                                              i_prefix,
                                              i_targetDir,
                                              i_scriptsDir,
                                              i_ignoreScriptsDir,
                                              i_joblistFileHandler,
                                              i_gzip,
                                              i_targetsInfo,
                                              i_debug)
            rmTmpFilesList.append(previousFilename)

        # filter mpileup
        previousFilename = filter_mpileupSupport_dna(i_pythonExecutable,
                                                     i_id,
                                                     i_chr,
                                                     previousFilename,
                                                     None,
                                                     True,
                                                     i_outputDir,
                                                     i_prefix,
                                                     i_scriptsDir,
                                                     i_ignoreScriptsDir,
                                                     i_joblistFileHandler,
                                                     i_gzip,
                                                     i_debug)
        rmTmpFilesList.append(previousFilename)

    else:
        # filter by blacklist
        if (i_blacklistFlag):
            previousFilename = filter_blacklist(i_pythonExecutable,
                                                i_id,
                                                i_chr,
                                                previousFilename,
                                                i_outputDir,
                                                i_prefix,
                                                i_blacklistDir,
                                                i_scriptsDir,
                                                i_ignoreScriptsDir,
                                                i_joblistFileHandler,
                                                i_gzip,
                                                i_debug)
            rmTmpFilesList.append(previousFilename)

        # filter by dbsnp
        if (i_dbSnpFlag):
            previousFilename = flag_dbSnp(i_pythonExecutable,
                                          i_id,
                                          i_chr,
                                          previousFilename,
                                          i_outputDir,
                                          i_prefix,
                                          i_dbSnpDir,
                                          i_scriptsDir,
                                          i_ignoreScriptsDir,
                                          i_joblistFileHandler,
                                          i_gzip,
                                          i_debug)
            rmTmpFilesList.append(previousFilename)

        # flag retro genes
        if (i_retroGenesFlag):
            previousFilename = flag_retroGenes(i_pythonExecutable,
                                               i_id,
                                               i_chr,
                                               previousFilename,
                                               i_outputDir,
                                               i_prefix,
                                               i_retroGenesDir,
                                               i_scriptsDir,
                                               i_ignoreScriptsDir,
                                               i_joblistFileHandler,
                                               i_gzip,
                                               i_debug)
            rmTmpFilesList.append(previousFilename)

        # flag pseudo genes
        if (i_pseudoGenesFlag):
            previousFilename = flag_pseudoGenes(i_pythonExecutable,
                                                i_id,
                                                i_chr,
                                                previousFilename,
                                                i_outputDir,
                                                i_prefix,
                                                i_pseudoGenesDir,
                                                i_scriptsDir,
                                                i_ignoreScriptsDir,
                                                i_joblistFileHandler,
                                                i_gzip,
                                                i_debug)
            rmTmpFilesList.append(previousFilename)

        # flag cosmic
        if (i_cosmicFlag):
            previousFilename = flag_cosmic(i_pythonExecutable,
                                           i_id,
                                           i_chr,
                                           previousFilename,
                                           i_outputDir,
                                           i_prefix,
                                           i_cosmicDir,
                                           i_scriptsDir,
                                           i_ignoreScriptsDir,
                                           i_joblistFileHandler,
                                           i_gzip,
                                           i_debug)
            rmTmpFilesList.append(previousFilename)

        # flag radar
        if (i_radarFlag):
            previousFilename = flag_radar(i_pythonExecutable,
                                          i_id,
                                          i_chr,
                                          previousFilename,
                                          i_outputDir,
                                          i_prefix,
                                          i_radarDir,
                                          i_scriptsDir,
                                          i_ignoreScriptsDir,
                                          i_joblistFileHandler,
                                          i_gzip,
                                          i_debug)
            rmTmpFilesList.append(previousFilename)

        # flag darned
        if (i_darnedFlag):
            previousFilename = flag_darned(i_pythonExecutable,
                                           i_id,
                                           i_chr,
                                           previousFilename,
                                           i_outputDir,
                                           i_prefix,
                                           i_darnedDir,
                                           i_scriptsDir,
                                           i_ignoreScriptsDir,
                                           i_joblistFileHandler,
                                           i_gzip,
                                           i_debug)
            rmTmpFilesList.append(previousFilename)

        # filter targets
        if (i_targetsFlag):
            previousFilename = filter_targets(i_pythonExecutable,
                                              i_id,
                                              i_chr,
                                              previousFilename,
                                              i_outputDir,
                                              i_prefix,
                                              i_targetDir,
                                              i_scriptsDir,
                                              i_ignoreScriptsDir,
                                              i_joblistFileHandler,
                                              i_gzip,
                                              i_targetsInfo,
                                              i_debug)
            rmTmpFilesList.append(previousFilename)

        # filter RNA mpileup
        # the output file contains all filters for all possible
        # mod types and no final mod type is chosen
        rnaFilename = filter_mpileupSupport_rna(i_pythonExecutable,
                                                i_id,
                                                i_chr,
                                                previousFilename,
                                                True,
                                                i_rnaMpileupMinMapQual,
                                                i_rnaMpileupMinAvgMapQual,
                                                i_outputDir,
                                                i_prefix,
                                                i_scriptsDir,
                                                i_ignoreScriptsDir,
                                                i_joblistFileHandler,
                                                i_gzip,
                                                i_debug)
        rmTmpFilesList.append(rnaFilename)

        # filter DNA mpileup
        # the output file contains all filters for all possible
        # mod types and no final mod type is chosen
        dnaFilename = filter_mpileupSupport_dna(i_pythonExecutable,
                                                i_id,
                                                i_chr,
                                                previousFilename,
                                                None,
                                                True,
                                                i_outputDir,
                                                i_prefix,
                                                i_scriptsDir,
                                                i_ignoreScriptsDir,
                                                i_joblistFileHandler,
                                                i_gzip,
                                                i_debug)
        rmTmpFilesList.append(dnaFilename)

        # compare the rna and dna
        # calls that pass in both the DNA and RNA will be in the overlaps file
        # calls that don't pass in the DNA but pass in the RNA are in the
        # non-overlaps file - these are the RNA Rescue and RNA Editing calls
        (overlapFilename,
         nonoverlapFilename) = radia_compare(i_pythonExecutable,
                                             i_id,
                                             i_chr,
                                             rnaFilename,
                                             dnaFilename,
                                             i_outputDir,
                                             i_prefix,
                                             i_scriptsDir,
                                             i_ignoreScriptsDir,
                                             i_joblistFileHandler,
                                             i_gzip,
                                             i_debug)
        rmTmpFilesList.append(overlapFilename)
        rmTmpFilesList.append(nonoverlapFilename)

        # filter DNA mpileup
        # filter the RNA Rescue and RNA Editing calls based on
        # the DNA to get rid of any possible germline calls
        previousFilename = filter_mpileupSupport_dna(i_pythonExecutable,
                                                     i_id,
                                                     i_chr,
                                                     nonoverlapFilename,
                                                     dnaFilename,
                                                     False,
                                                     i_outputDir,
                                                     i_prefix,
                                                     i_scriptsDir,
                                                     i_ignoreScriptsDir,
                                                     i_joblistFileHandler,
                                                     i_gzip,
                                                     i_debug)
        rmTmpFilesList.append(previousFilename)

        # filter out possible germline calls
        previousFilename = filter_rnaOnly(i_chr,
                                          previousFilename,
                                          i_outputDir,
                                          i_prefix,
                                          i_joblistFileHandler,
                                          i_gzip,
                                          i_debug)
        rmTmpFilesList.append(previousFilename)

        # if we should blat and we have something to blat
        if (i_blatFlag and
            os.path.isfile(previousFilename) and
            os.stat(previousFilename).st_size > 20):

            # create blat input
            blatInputFilename = create_blat_input(i_pythonExecutable,
                                                  i_id,
                                                  i_chr,
                                                  previousFilename,
                                                  i_transcriptNameTag,
                                                  i_transcriptCoordinateTag,
                                                  i_transcriptStrandTag,
                                                  i_rnaIncludeSecAligns,
                                                  i_outputDir,
                                                  i_prefix,
                                                  i_scriptsDir,
                                                  i_ignoreScriptsDir,
                                                  i_joblistFileHandler,
                                                  i_debug)
            rmTmpFilesList.append(blatInputFilename)

            # filter by BLAT
            (blatOutputFilename,
             previousFilename) = filter_blat(i_pythonExecutable,
                                             i_id,
                                             i_chr,
                                             previousFilename,
                                             rnaFilename,
                                             blatInputFilename,
                                             i_blatFastaFilename,
                                             i_transcriptNameTag,
                                             i_transcriptCoordinateTag,
                                             i_transcriptStrandTag,
                                             i_rnaIncludeSecAligns,
                                             i_outputDir,
                                             i_prefix,
                                             i_scriptsDir,
                                             i_ignoreScriptsDir,
                                             i_joblistFileHandler,
                                             i_gzip,
                                             i_debug)
            rmTmpFilesList.append(blatOutputFilename)
            rmTmpFilesList.append(previousFilename)

        # if RNA only, just merge the RNA Confirmation and RNA Rescue calls
        if (i_rnaOnlyFlag):
            # the dnaFilename and --dnaHeaderOnly=True means that we only
            # extract the header from the dnaFilename and ignore the rest
            previousFilename = merge_rnaAndDna(i_pythonExecutable,
                                               i_id,
                                               i_chr,
                                               dnaFilename,
                                               rnaFilename,
                                               overlapFilename,
                                               previousFilename,
                                               True,
                                               i_outputDir,
                                               i_prefix,
                                               i_scriptsDir,
                                               i_ignoreScriptsDir,
                                               i_joblistFileHandler,
                                               i_gzip,
                                               i_debug)
            rmTmpFilesList.append(previousFilename)
        else:
            # merge RNA and DNA
            # the dnaFilename and --dnaHeaderOnly=False means that we
            # merge the header and the results in the dnaFilename
            previousFilename = merge_rnaAndDna(i_pythonExecutable,
                                               i_id,
                                               i_chr,
                                               dnaFilename,
                                               rnaFilename,
                                               overlapFilename,
                                               previousFilename,
                                               False,
                                               i_outputDir,
                                               i_prefix,
                                               i_scriptsDir,
                                               i_ignoreScriptsDir,
                                               i_joblistFileHandler,
                                               i_gzip,
                                               i_debug)
            rmTmpFilesList.append(previousFilename)

    if (i_snpEffFlag):
        # keep track of the orginal file
        preSnpEffFilename = previousFilename

        # extracting passing calls
        previousFilename = extract_passing(i_chr,
                                           previousFilename,
                                           i_outputDir,
                                           i_prefix,
                                           i_joblistFileHandler,
                                           i_gzip,
                                           i_debug)
        rmTmpFilesList.append(previousFilename)

        previousFilename = run_snpEff(i_chr,
                                      previousFilename,
                                      i_snpEffDir,
                                      i_snpEffGenome,
                                      i_snpEffCanonical,
                                      i_outputDir,
                                      i_prefix,
                                      i_joblistFileHandler,
                                      i_debug)
        rmTmpFilesList.append(previousFilename)

        if (not i_dnaOnlyFlag and i_rnaBlacklistFlag):
            # filter RNA by geneNames/Families
            previousFilename = filter_rnaBlacklist(i_pythonExecutable,
                                                   i_chr,
                                                   previousFilename,
                                                   i_rnaGeneBlckFilename,
                                                   i_rnaGeneFamilyBlckFilename,
                                                   i_outputDir,
                                                   i_prefix,
                                                   i_scriptsDir,
                                                   i_ignoreScriptsDir,
                                                   i_joblistFileHandler,
                                                   i_gzip,
                                                   i_debug)
            rmTmpFilesList.append(previousFilename)

        # merge passing with snpEff back with originals
        previousFilename = merge_passingAndOriginals(i_pythonExecutable,
                                                     i_chr,
                                                     previousFilename,
                                                     preSnpEffFilename,
                                                     i_outputDir,
                                                     i_prefix,
                                                     i_scriptsDir,
                                                     i_ignoreScriptsDir,
                                                     i_joblistFileHandler,
                                                     i_gzip,
                                                     i_debug)
        rmTmpFilesList.append(previousFilename)

    # everything gets run through the read support filter
    # startTime = time.time()
    previousFilename = filter_readSupport(i_pythonExecutable,
                                          i_chr,
                                          previousFilename,
                                          i_transcriptNameTag,
                                          i_transcriptCoordinateTag,
                                          i_transcriptStrandTag,
                                          i_rnaIncludeSecAligns,
                                          i_readSupportMinMapQual,
                                          i_outputDir,
                                          i_prefix,
                                          i_outputFilename,
                                          i_scriptsDir,
                                          i_ignoreScriptsDir,
                                          i_joblistFileHandler,
                                          i_gzip,
                                          i_debug)
    # stopTime = time.time()
    '''
    logging.warning("filter_readSupport: prefix=%s: Total time=%s hrs, " +
                    "%s mins, %s secs", i_prefix,
                    ((stopTime-startTime)/(3600)),
                    ((stopTime-startTime)/60),
                    (stopTime-startTime))
    '''
    # if we aren't debugging, then remove all the tmp files
    if (not i_debug):
        # remove all the temp files
        remove_tmpFiles(rmTmpFilesList, i_joblistFileHandler, i_debug)

    if (i_joblistDir is not None):
        i_joblistFileHandler.close()

    return


main()
sys.exit(0)
