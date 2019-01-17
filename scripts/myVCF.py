#!/usr/bin/env python

from collections import OrderedDict
import re
import sys
from itertools import izip

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
        return ('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">' %
                (self.id, self.number, self.type, self.description))


class Sample:

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
            self.description = ""
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
            return ('##SAMPLE=<ID=%s,Individual="%s",Description="%s",' +
                    'File="%s",Platform="%s",Source="%s",Accession=%s>' %
                    (self.id, self.individual, self.description, self.file,
                     self.platform, self.source, self.accession))
        else:
            return ('##SAMPLE=<ID=%s,Individual="%s",Description="%s",' +
                    'File="%s",Platform="%s",Source="%s">' %
                    (self.id, self.individual, self.description,
                     self.file, self.platform, self.source))


def init_from_match(anId, anIndividual, aDescription, aFile,
                    aPlatform, aSource, anAccession=None):
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
            except ValueError:
                # some protein change annotations have an
                # '=' in them (e.g. ProtCh=p.T394=;VC=Silent)
                pcList = item.split("=")
                tag = pcList[0]
                value = pcList[1] + "="
            infoDict[tag] = value.split(",")
        else:
            infoDict[item] = True
    return infoDict


def format_info(anInfoDict):
    lineList = []
    for key, value in anInfoDict.items():
        if value is True:
            lineList.append(key)
        else:
            lineList.append("%s=%s" % (key, ",".join(value)))
    return ";".join(lineList)


def format_sample_data(aFormatField, aDataDict):

    # GT:DP:AD:AF:INS:DEL:START:STOP:MQ0:MMQ:MQA:BQ:SB:MMP
    # 0/1:7:2,5:0.29,0.71:0:0:0:0:1,3:1,1:1,0:82,71:1.0,1.0

    # a dict with:
    #    - the format item as key (e.g. DP)
    #      and single items as values or
    #    - the format item as key (e.g. AD)
    #      and a list of allele specific items as values
    # e.g. rnaTumorDict["DP"] = [100]
    # e.g. rnaTumorDict["AD"] = [50,50]

    dataFieldList = list()
    for formatItem in aFormatField.split(":"):
        if (formatItem == "GT"):
            sep = "/"
        else:
            sep = ","

        dataFieldList.append(sep.join(aDataDict[formatItem]))

    return ":".join(dataFieldList)


def parse_sample_data(aFormatField, aDataField, aNumAlleles):

    # GT:DP:AD:AF:INS:DEL:START:STOP:MQ0:MMQ:MQA:BQ:SB:MMP
    # 0/1:7:2,5:0.29,0.71:0:0:0:0:1,3:1,1:1,0:82,71:1.0,1.0

    # a dict with:
    #    - the format item as key (e.g. DP)
    #      and single items as values or
    #    - the format item as key (e.g. AD)
    #      and a list of allele specific items as values
    # e.g. rnaTumorDict["DP"] = [100]
    # e.g. rnaTumorDict["AD"] = [50,50]

    dataDict = dict()
    for (formatItem, dataItem) in izip(aFormatField.split(":"),
                                       aDataField.split(":")):
        if (formatItem == "GT"):
            sep = "/"
        else:
            sep = ","

        if (dataItem != "."):
            dataDict[formatItem] = dataItem.split(sep)
        else:
            dataDict[formatItem] = ["."] * aNumAlleles

    return dataDict


class Data:
    def __init__(self, aHeadersList, aChrom, aPos, anId, aRef,
                 anAltField, aQual, aFilterField, anInfoField,
                 aFormatField, aSampleData=None):

        self.chrom = aChrom
        self.pos = int(aPos)
        self.id = anId
        self.ref = aRef
        self.altList = anAltField.split(",")
        self.allelesList = [self.ref] + self.altList
        lenAltList = len(self.allelesList)
        self.qual = aQual
        self.filterList = aFilterField.split(";")
        self.infoDict = parse_info(anInfoField)
        self.formatField = aFormatField
        self.sampleData = aSampleData
        self.dnaNormalDict = None
        self.rnaNormalDict = None
        self.dnaTumorDict = None
        self.rnaTumorDict = None

        # a dict with:
        #    - the format item as key (e.g. DP)
        #      and single items as values or
        #    - the format item as key (e.g. AD)
        #      and a list of allele specific items as values
        # e.g. rnaTumorDict["DP"] = [100]
        # e.g. rnaTumorDict["AD"] = [50,50]
        try:
            # if we have a 9th column, figure out which dataset it is
            if (len(aHeadersList) > 9):
                if (aHeadersList[9] == "DNA_NORMAL"):
                    self.dnaNormalDict = parse_sample_data(aFormatField,
                                                           self.sampleData[0],
                                                           lenAltList)
                elif (aHeadersList[9] == "RNA_NORMAL"):
                    self.rnaNormalDict = parse_sample_data(aFormatField,
                                                           self.sampleData[0],
                                                           lenAltList)
                elif (aHeadersList[9] == "DNA_TUMOR"):
                    self.dnaTumorDict = parse_sample_data(aFormatField,
                                                          self.sampleData[0],
                                                          lenAltList)
                elif (aHeadersList[9] == "RNA_TUMOR"):
                    self.rnaTumorDict = parse_sample_data(aFormatField,
                                                          self.sampleData[0],
                                                          lenAltList)
            # if we have a 10th column, figure out which dataset it is
            if (len(aHeadersList) > 10):
                if (aHeadersList[10] == "RNA_NORMAL"):
                    self.rnaNormalDict = parse_sample_data(aFormatField,
                                                           self.sampleData[1],
                                                           lenAltList)
                elif (aHeadersList[10] == "DNA_TUMOR"):
                    self.dnaTumorDict = parse_sample_data(aFormatField,
                                                          self.sampleData[1],
                                                          lenAltList)
                elif (aHeadersList[10] == "RNA_TUMOR"):
                    self.rnaTumorDict = parse_sample_data(aFormatField,
                                                          self.sampleData[1],
                                                          lenAltList)
            # if we have a 11th column, figure out which dataset it is
            if (len(aHeadersList) > 11):
                if (aHeadersList[11] == "DNA_TUMOR"):
                    self.dnaTumorDict = parse_sample_data(aFormatField,
                                                          self.sampleData[2],
                                                          lenAltList)
                elif (aHeadersList[11] == "RNA_TUMOR"):
                    self.rnaTumorDict = parse_sample_data(aFormatField,
                                                          self.sampleData[2],
                                                          lenAltList)
            # if we have a 12th column, figure out which dataset it is
            if (len(aHeadersList) > 12):
                if (aHeadersList[12] == "RNA_TUMOR"):
                    self.rnaTumorDict = parse_sample_data(aFormatField,
                                                          self.sampleData[3],
                                                          lenAltList)
        except:
            print >> sys.stderr, ("Problem parsing line: ", self.chrom,
                                  self.pos, aFormatField, aSampleData)
            raise

    def __str__(self):
        dnaNormalString = None
        rnaNormalString = None
        dnaTumorString = None
        rnaTumorString = None

        sampleList = list()
        if self.dnaNormalDict is not None:
            dnaNormalString = format_sample_data(self.formatField,
                                                 self.dnaNormalDict)
            sampleList.append(dnaNormalString)
        if self.rnaNormalDict is not None:
            rnaNormalString = format_sample_data(self.formatField,
                                                 self.rnaNormalDict)
            sampleList.append(rnaNormalString)
        if self.dnaTumorDict is not None:
            dnaTumorString = format_sample_data(self.formatField,
                                                self.dnaTumorDict)
            sampleList.append(dnaTumorString)
        if self.rnaTumorDict is not None:
            rnaTumorString = format_sample_data(self.formatField,
                                                self.rnaTumorDict)
            sampleList.append(rnaTumorString)

        return "\t".join(map(str, [self.chrom, self.pos, self.id, self.ref,
                                   ",".join(self.altList), self.qual,
                                   ";".join(self.filterList),
                                   format_info(self.infoDict),
                                   self.formatField, "\t".join(sampleList)]))


class VCF:
    def __init__(self):
        self.meta = []
        self.infos = []
        self.filters = []
        self.formats = []
        self.headers = []
        self.data = []

    def make_info(self, aValue):
        # make re match spec
        match = re.match(r"""<ID=(.*),Number=(.*),Type=(.*),""" +
                         """Description=['"](.*)['"]""", aValue)
        if match:
            return Info(*match.groups())
        else:
            raise VCFFormatError("Improperly formatted INFO metadata: " +
                                 aValue)

    def make_sample(self, aValue):
        match = re.match(r"""<ID=(.*),Individual="(.*)",""" +
                         """Description="(.*)",File="(.*)",""" +
                         """Platform="(.*)",Source="(.*)">""", aValue)
        if match:
            return init_from_match(*match.groups())
        else:
            match = re.match(r"""<ID=(.*),Individual="(.*)",""" +
                             """Description="(.*)",File="(.*)",""" +
                             """Platform="(.*)",Source="(.*)",""" +
                             """Accession=(.*)>""", aValue)
            if match:
                return init_from_match(*match.groups())
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
        # make re match spec
        match = re.match(r"""<ID=(.*),Description=['"](.*)['"]>""", aValue)
        if match:
            self.filters.append(Filter(*match.groups()))
        else:
            raise VCFFormatError("Improperly formatted FILTER metadata: " +
                                 aValue)

    def set_headers(self, aHeaderList):
        self.headers = aHeaderList

    def make_data(self, aDataList):
        baseData = aDataList[:9]
        genotypeData = aDataList[9:]
        data = baseData + [genotypeData]

        return Data(self.headers, *data)

    def add_data(self, data_list):
        self.data.append(self.make_data(data_list))
