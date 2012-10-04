#!/usr/bin/env python2.7
import sys

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
            chr = self.rchroms[c]
            for b in sd[c]:
                for st, sp, v in sd[c][b]:
                    v = [chr, str(st), str(sp), v]
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

    def overlapswith(self, tuple, buf=0):
        chrom, st, sp = tuple
        #print chrom, st, sp

        try:
            c = self.chroms[chrom]
        except:
            #print "no chrom"
            return (False, "")       

        b = self.findbin(st)
        #print self.data
        #print "bin", b

        if b not in self.data[c]:
            #print "bin not in data"
            return (False, "")
        
        for qst, qsp, qv in self.data[c][b]:
            #if st <= qsp + buf and sp >= qst - buf:
            if st >= qst + buf and sp <= qsp - buf:
                return (True, qv)

        #print "didn't overlap"
        return (False, "")
        
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

        for bin in range(startBin, (stopBin + self.binsize), self.binsize):
            if bin not in self.data[c]:
                #print bin, "not in data"
                self.data[c][bin] = []

            self.data[c][bin].append( (st, sp, v) )
            #print "appending", st, sp

        '''
        pb = b
        b = self.findbin(sp)
        if b == pb:
            return

        if b not in self.data[c]:
            self.data[c][b] = []
       
        self.data[c][b].append( (st, sp, v) )
        '''

    def loadfromfile(self, fname, ci=0, sti=1, spi=2, vi=3):
            
        inFile = open(fname)

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

    
