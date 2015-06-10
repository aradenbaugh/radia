#!/usr/bin/env python

from collections import OrderedDict
import re


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


class VCFFormatError(Exception):

    def __init__(self, text):
        Exception.__init__(self, text)

class Filter:
    def __init__(self, id, description):
        self.id = id
        self.description = description

class Info:
    def __init__(self, id, number, type, description):
        self.id = id
        self.number = number
        self.type = type
        self.description = description
        
    def __str__(self):
        return '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">' % (self.id, self.number, self.type, self.description)

class Sample:
    def __init__(self, id, individual, description, file, platform, source, accession=None):
        self.id = id
        self.individual = individual
        self.description = description
        self.file = file
        self.platform = platform
        self.source = source
        self.accession = accession
    
    def __init__(self, sampledict):
        if "ID" in sampledict:
            self.id = sampledict["ID"]
        else:
            self.id = ""
        if "Individual" in sampledict:
            self.individual = sampledict["Individual"].strip('"')
        else:
            self.individual = ""
        if "Description" in sampledict:
            self.description = sampledict["Description"].strip('"')
        else:
            self.description=""
        if "File" in sampledict:
            self.file = sampledict["File"].strip('"')
        else:
            self.file = ""
        if "Platform" in sampledict:
            self.platform = sampledict["Platform"].strip('"')
        else:
            self.platform = ""
        if "Source" in sampledict:
            self.source = sampledict["Source"].strip('"')
        else:
            self.source = ""
        if "Accession" in sampledict:
            self.accession = sampledict["Accession"].strip('"')
        else:
            self.accession = ""

    def __str__(self):
        if self.accession:
            return '##SAMPLE=<ID=%s,Individual="%s",Description="%s",File="%s",Platform="%s",Source="%s",Accession=%s>'  % (self.id, self.individual, self.description, self.file, self.platform, self.source, self.accession)  
        else:
            return '##SAMPLE=<ID=%s,Individual="%s",Description="%s",File="%s",Platform="%s",Source="%s">' % (self.id, self.individual, self.description, self.file, self.platform, self.source)

def init_from_match(id, individual, description, file, platform, source, accession=None):
    initdict = {}
    initdict["ID"] = id
    initdict["Individual"] = individual
    initdict["Description"] = description
    initdict["File"] = file
    initdict["Platform"] = platform
    initdict["Source"] = source
    if accession:
        initdict["Accession"] = accession
    return Sample(initdict)


def parse_info(info):
    info_dict = OrderedDict()
    for i in info.split(";"):
        if "=" in i:
            id,value = i.split("=")
            info_dict[id]=value
        else:
            info_dict[i]=True
    return info_dict

def format_info(info_dict):
    line_list = []
    for key, value in info_dict.items():
        if value is True:
            line_list.append(key)
        else:
            line_list.append("%s=%s" % (key, value))
    return ";".join(line_list)
    
class Data:
    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, genotype_data = None, extra_headers = None):
        self.chrom = chrom
        self.pos = int(pos)
        self.id = id
        self.ref = ref
        self.alt = alt.split(",")
        self.qual = qual
        self.filter = filter.split(";")
        self.info = parse_info(info)
        self.genotype = genotype_data

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

    def make_info(self, value):
        #make re match spec
        match = re.match(r"""<ID=(.*),Number=(.*),Type=(.*),Description=['"](.*)['"]""", value)
        if match:
            return Info(*match.groups())
        else:
            raise VCFFormatError("Improperly formatted INFO metadata: " + value)
            
    def make_sample(self, value):
        match = re.match(r"""<ID=(.*),Individual="(.*)",Description="(.*)",File="(.*)",Platform="(.*)",Source="(.*)">""", value)
        if match:
            return init_from_match(*match.groups())
        else:
            match = re.match(r"""<ID=(.*),Individual="(.*)",Description="(.*)",File="(.*)",Platform="(.*)",Source="(.*)",Accession=(.*)>""", value)
            if match:
                return Sample(*match.groups())
            else:
                sampledict = {}
                for pair in value.strip("<>").split(","):
                    key, value = pair.split("=")
                    sampledict[key] = value
                return Sample(sampledict)
                
        
        raise VCFFormatError("Improperly formatted SAMPLE metadata: " + value)
        
    def add_info(self, value):
        self.infos.append(self.make_info(value))
    
    def add_filter(self, value):
        #make re match spec
        match = re.match(r"""<ID=(.*),Description=['"](.*)['"]>""", value)
        if match:
            self.filters.append(Filter(*match.groups()))
        else:
            raise VCFFormatError("Improperly formatted FILTER metadata: " + value)
        #self.filters.append(Filter(values[0], values[1]))

    def set_headers(self, header_list):
        self.headers = header_list

    def make_data(self, data_list):
        baseData = data_list[:8]
        genotypeData = data_list[8:]
        data = baseData+[genotypeData]+[self.headers[8:]]
        return Data(*data)
        
    def add_data(self, data_list):
        self.data.append(self.make_data(data_list))


def make_accession(filename):
    file = open(filename, 'r')
    
    #ignore first line
    file.readline()
    
    accession = {}
    for line in file:
        data = line.strip().split(',')
        srs = data[3]
        match = re.search(r'TCGA-\w{2}-\w{4}-\w{3}-\w{3}-\w{4}-\w{2}', data[4])
        if match:
            file = match.group()
            accession[file] = srs
    return accession

