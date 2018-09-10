#!/usr/bin/env python

from collections import OrderedDict
import re


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


class VCFFormatError(Exception):

    def __init__(self, text):
        Exception.__init__(self, text)

class Filter:
    def __init__(self, anId, aDescription):
        self.id = anId
        self.description = aDescription

class Info:
    def __init__(self, anId, aNumber, aType, aDescription):
        self.id = anId
        self.number = aNumber
        self.type = aType
        self.description = aDescription
        
    def __str__(self):
        return '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">' % (self.id, self.number, self.type, self.description)


class Sample:
    def __init__(self, anId, anIndividual, aDescription, aFile, aPlatform, aSource, anAccession=None):
        self.id = anId
        self.individual = anIndividual
        self.description = aDescription
        self.file = aFile
        self.platform = aPlatform
        self.source = aSource
        self.accession = anAccession
    
    def __init__(self, aSampleDict):
        if "ID" in aSampleDict:
            self.id = aSampleDict["ID"]
        else:
            self.id = ""
        if "Individual" in aSampleDict:
            self.individual = aSampleDict["Individual"].strip('"')
        else:
            self.individual = ""
        if "Description" in aSampleDict:
            self.description = aSampleDict["Description"].strip('"')
        else:
            self.description=""
        if "File" in aSampleDict:
            self.file = aSampleDict["File"].strip('"')
        else:
            self.file = ""
        if "Platform" in aSampleDict:
            self.platform = aSampleDict["Platform"].strip('"')
        else:
            self.platform = ""
        if "Source" in aSampleDict:
            self.source = aSampleDict["Source"].strip('"')
        else:
            self.source = ""
        if "Accession" in aSampleDict:
            self.accession = aSampleDict["Accession"].strip('"')
        else:
            self.accession = ""

    def __str__(self):
        if self.accession:
            return '##SAMPLE=<ID=%s,Individual="%s",Description="%s",File="%s",Platform="%s",Source="%s",Accession=%s>'  % (self.id, self.individual, self.description, self.file, self.platform, self.source, self.accession)  
        else:
            return '##SAMPLE=<ID=%s,Individual="%s",Description="%s",File="%s",Platform="%s",Source="%s">' % (self.id, self.individual, self.description, self.file, self.platform, self.source)


def init_from_match(anId, anIndividual, aDescription, aFile, aPlatform, aSource, anAccession=None):
    initDict = {}
    initDict["ID"] = anId
    initDict["Individual"] = anIndividual
    initDict["Description"] = aDescription
    initDict["File"] = aFile
    initDict["Platform"] = aPlatform
    initDict["Source"] = aSource
    if anAccession:
        initDict["Accession"] = anAccession
    return Sample(initDict)


def parse_info(anInfoField):
    infoDict = OrderedDict()
    for item in anInfoField.split(";"):
        if "=" in item:
            try:
                (tag, value) = item.split("=")
            except:
                # some protein change annotations have an '=' in them (e.g. ProtCh=p.T394=;VC=Silent)
                pcList = item.split("=")
                tag = pcList[0]
                value = pcList[1] + "="
            infoDict[tag] = value
        else:
            infoDict[item] = True
    return infoDict


def format_info(anInfoDict):
    lineList = []
    for key, value in anInfoDict.items():
        if value is True:
            lineList.append(key)
        else:
            lineList.append("%s=%s" % (key, value))
    return ";".join(lineList)
    

class Data:
    def __init__(self, aChrom, aPos, anId, aRef, anAltField, aQual, aFilterField, anInfoField, aGenotypeData = None):
        self.chrom = aChrom
        self.pos = int(aPos)
        self.id = anId
        self.ref = aRef
        self.alt = anAltField.split(",")
        self.qual = aQual
        self.filter = aFilterField.split(";")
        self.info = parse_info(anInfoField)
        self.genotype = aGenotypeData

    def __str__(self):
        return "\t".join(map(str, [self.chrom, self.pos, self.id, self.ref, ",".join(self.alt), self.qual, ";".join(self.filter), format_info(self.info), "\t".join(self.genotype)]))


class VCF:
    def __init__(self):
        self.meta = []
        self.infos = []
        self.filters = []
        self.formats = []
        self.headers = []
        self.data = []


    def make_info(self, aValue):
        #make re match spec
        match = re.match(r"""<ID=(.*),Number=(.*),Type=(.*),Description=['"](.*)['"]""", aValue)
        if match:
            return Info(*match.groups())
        else:
            raise VCFFormatError("Improperly formatted INFO metadata: " + aValue)
            
    
    def make_sample(self, aValue):
        match = re.match(r"""<ID=(.*),Individual="(.*)",Description="(.*)",File="(.*)",Platform="(.*)",Source="(.*)">""", aValue)
        if match:
            return init_from_match(*match.groups())
        else:
            match = re.match(r"""<ID=(.*),Individual="(.*)",Description="(.*)",File="(.*)",Platform="(.*)",Source="(.*)",Accession=(.*)>""", aValue)
            if match:
                return Sample(*match.groups())
            else:
                sampleDict = {}
                for pair in aValue.strip("<>").split(","):
                    key, val = pair.split("=")
                    sampleDict[key] = val
                return Sample(sampleDict)
                
        raise VCFFormatError("Improperly formatted SAMPLE metadata: " + aValue)
        
        
    def add_info(self, aValue):
        self.infos.append(self.make_info(aValue))
    
    def add_filter(self, aValue):
        #make re match spec
        match = re.match(r"""<ID=(.*),Description=['"](.*)['"]>""", aValue)
        if match:
            self.filters.append(Filter(*match.groups()))
        else:
            raise VCFFormatError("Improperly formatted FILTER metadata: " + aValue)
    

    def set_headers(self, aHeaderList):
        self.headers = aHeaderList


    def make_data(self, aDataList):
        baseData = aDataList[:8]
        genotypeData = aDataList[8:]
        data = baseData + [genotypeData]
        return Data(*data)


    def add_data(self, data_list):
        self.data.append(self.make_data(data_list))


def make_accession(aFilename):
    accessionFile = open(aFilename, 'r')
    
    #ignore first line
    accessionFile.readline()
    
    accession = {}
    for line in accessionFile:
        data = line.strip().split(',')
        srs = data[3]
        match = re.search(r'TCGA-\w{2}-\w{4}-\w{3}-\w{3}-\w{4}-\w{2}', data[4])
        if match:
            filename = match.group()
            accession[filename] = srs
    return accession

