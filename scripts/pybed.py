#!/usr/bin/env python
import sys
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
    

def overlap(t, q, buf = 0):
    if len(t) == 0 or len(q) == 0:
        return []

    tuples = []
    for tst, tsp, tv in t:
        for qst, qsp, qv in q:
            if tst <= qsp + buf and tsp >= qst - buf:
                tuples.append( (tst, tsp, tv) )
                break
            
    return tuples

class pybed:
    def __init__(self, binsize=100000):
        self.paired  = False
        self.binsize = binsize
        self.chroms  = {}
        self.rchroms = {}
        for i in range(1, 23):
            self.chroms[str(i)] = i
            self.chroms['chr' + str(i)] = i
            self.rchroms[i] = 'chr' + str(i)
        
        i += 1
        self.chroms['X'] = i
        self.chroms['chrX'] = i
        self.rchroms[i] = 'chrX'

        i += 1
        self.chroms['Y'] = i
        self.chroms['chrY'] = i
        self.rchroms[i] = 'chrY'

        i += 1
        self.chroms['M'] = i
        self.chroms['chrM'] = i
        self.rchroms[i] = 'chrM'

        i += 1
        self.chroms['MT'] = i
        self.chroms['chrMT'] = i
        self.rchroms[i] = 'chrMT'

        numchroms = max(self.chroms.values())

        self.data = {}
        for i in range(1, numchroms+1):
            if i not in self.data:
                self.data[i] = {}

    def length(self):
        sd = self.data

        count = 0
        for c in sd:
            for b in sd[c]:
                count += len(sd[c][b])
        return count
    
    def output(self, f=sys.stdout, delim='\t'):
        sd = self.data
        
        for c in sd:
            chrom = self.rchroms[c]
            for b in sd[c]:
                for st, sp, v in sd[c][b]:
                    v = [chrom, str(st), str(sp), v]
                    print >> f, delim.join(v)

    def intersect(self, other, buffer=0):
        new = pybed()
        sd = self.data
        od = other.data

        for c in sd:
            chrom = self.rchroms[c]
            for b in sd[c]:
                if b not in od[c]:
                    continue
                ovtuples = overlap(sd[c][b], od[c][b], buffer)
                for ov in ovtuples:
                    st, sp, v = ov
                    new.loadtuple( (chrom, st, sp, v) )

        return new
    
    def findbin(self, pos):
        return int( round(float(pos) / float(self.binsize)) * self.binsize )

    def overlapswith(self, tuple, anIncludeCount, buf=0):
        chrom, st, sp = tuple
        #print chrom, st, sp

        try:
            c = self.chroms[chrom]
        except:
            #print "no chrom"
            return (False, "", 0)       

        b = self.findbin(st)
        #print self.data
        #print "bin", b

        if b not in self.data[c]:
            #print "bin not in data"
            return (False, "", 0)
        
        count=0
        for qst, qsp, qv in self.data[c][b]:
            #if st <= qsp + buf and sp >= qst - buf:
            #print "loop", st, qst, qsp, qv
            if st >= qst + buf and sp <= qsp - buf:
                if (anIncludeCount):
                    count += 1
                else:
                    return (True, qv, 0)
                #print st, ">=", qst, sp, "<=", qsp, "count=", count
            # if the start is past the bin stop, then stop looping
            if (qst > sp):
                #print qst, ">", sp, "break"
                break;
                
        if (count > 0):
            return (True, qv, count)

        #print "didn't overlap"
        return (False, "", 0)
        
    def loadtuple(self, tuple):
        chrom, st, sp, v = tuple
        #print "trying to load", chrom, st, sp

        try:
            c = self.chroms[chrom]
        except:
            #print "problem with chrom"
            return
        
        startBin = self.findbin(st)
        stopBin = self.findbin(sp)
        #print "loading bin" , startBin, stopBin

        for currentBin in range(startBin, (stopBin + self.binsize), self.binsize):
            if currentBin not in self.data[c]:
                #print currentBin, "not in data"
                self.data[c][currentBin] = []

            self.data[c][currentBin].append( (st, sp, v) )
            #print "appending", st, sp

    def loadfromfile(self, fname, ci=0, sti=1, spi=2, vi=3):
            
        inFile = get_read_fileHandler(fname)

        for line in inFile:
            data = line[:-1].split('\t')

            c  = data[ci]
            st = int(data[sti])
            sp = int(data[spi])

            if len(data) < 4:
                v = ''
            else:
                v = data[vi]
            
            self.loadtuple( (c, st, sp, v) )

        inFile.close()

    
