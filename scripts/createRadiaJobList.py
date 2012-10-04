#!/usr/bin/env python2.7
import sys, os, subprocess, re

'''
'    python2.7 createRadiaJobList.py A1,A2,A7,A8,AC,AN,AO,AQ,AR,B6,BH,C7,C8,D8,E2,E9,EW,GI,GM,HN,JL,K3,LD,LL,LQ,M2 ../cancer/brca/dirLists/allDirList.txt /inside/grotto/users/aradenba/data/ brca ../cancer/brca/
'    python2.7 createRadiaJobList.py A1,A2,A7,A8,AC,AN,AO,AQ,AR,B6,BH,C7,C8,D8,E2,E9,EW,GI,GM,HN,JL,K3,LD,LL,LQ,M2 ../cancer/brca/dirLists/rnaUuidDirList.txt /inside/grotto/users/aradenba/data/ brca ../cancer/brca/
'    python2.7 createRadiaJobList.py A1,A2,A7,A8,AC,AN,AO,AQ,AR,B6,BH,C7,C8,D8,E2,E9,EW,GI,GM,HN,JL,K3,LD,LL,LQ,M2 ../cancer/brca/dirLists/finalDirList.txt /inside/grotto/users/aradenba/data/ brca ../cancer/brca/
'
'''

# ref sequences... must update as new sequences are brought into fold
faDir = '/inside/depot/fa/'
    
refs = {}
refs['GRCh37-lite']                     = ('', 'hg19', "MT", os.path.join(faDir, 'GRCh37-lite.fa'))
refs['GRCh37-lite.fa']                  = ('', 'hg19', "MT", os.path.join(faDir, 'GRCh37-lite.fa'))
refs['GRCh37-lite_w_chr_prefix']        = ('chr', 'hg19', "M", os.path.join(faDir, 'GRCh37-lite_w_chr_prefix.fa'))
refs['GRCh37-lite_w_chr_prefix.fa']     = ('chr', 'hg19', "M", os.path.join(faDir, 'GRCh37-lite_w_chr_prefix.fa'))
refs['NCBI36_WUGSC_variant']            = ('', 'hg18', "MT", os.path.join(faDir, 'NCBI36_WUGSC_variant.fasta'))
refs['NCBI36_WUGSC_variant.fasta']      = ('', 'hg18', "MT", os.path.join(faDir, 'NCBI36_WUGSC_variant.fasta'))
refs['NCBI36-HG18_Broad_variant']       = ('chr', 'hg18', "M", os.path.join(faDir, 'NCBI36-HG18_Broad_variant.fasta'))
refs['NCBI36-HG18_Broad_variant.fasta'] = ('chr', 'hg18', "M", os.path.join(faDir, 'NCBI36-HG18_Broad_variant.fasta'))
refs['Homo_sapiens_assembly18']         = ('chr', 'hg18', "M", os.path.join(faDir, 'Homo_sapiens_assembly18.fasta'))
refs['Homo_sapiens_assembly18.fasta']   = ('chr', 'hg18', "M", os.path.join(faDir, 'Homo_sapiens_assembly18.fasta'))
refs['Homo_sapiens_assembly19']         = ('', 'hg19', "MT", os.path.join(faDir, 'Homo_sapiens_assembly19.fasta'))
refs['Homo_sapiens_assembly19.fasta']   = ('', 'hg19', "MT", os.path.join(faDir, 'Homo_sapiens_assembly19.fasta'))
refs['all_sequences']                   = ('', 'hg18', "MT", os.path.join(faDir, 'all_sequences.fasta'))
refs['all_sequences.fasta']             = ('', 'hg18', "MT", os.path.join(faDir, 'all_sequences.fasta'))

# file from BCM could not be indexed by samtools faidx, using sequence from Broad, which appears equivalent.
refs['hsap36.1-hg18.fasta']             = ('chr', 'hg18', "M", os.path.join(faDir, 'hsap36.1-hg18.fasta'))
refs['hsap_36.1_hg18.fa']               = ('chr', 'hg18', "M", os.path.join(faDir, 'hsap36.1-hg18.fasta'))
refs['hsap36.fasta']                    = ('chr', 'hg18', "M", os.path.join(faDir, 'hsap36.1-hg18.fasta'))

refs['AS:hg19']                         = ('chr', 'hg19', "M", os.path.join(faDir, 'hg19.fasta'))
refs['LB:hg19']                         = ('chr', 'hg19', "M", os.path.join(faDir, 'hg19.fasta'))
refs['LB:HG19']                         = ('chr', 'hg19', "M", os.path.join(faDir, 'hg19.fasta'))

skip = ['bai', 'ncbi_enc', 'md5', 'mirna']

normal_suffices = ['10', '11']
tumor_suffices = ['01', '03']

suffices = ['_SOLiD_whole_exome_extensions', '_SOLiD', 
            '_whole', '_capture',
            '_IlluminaGA-DNASeq_whole', '_IlluminaGA-DNASeq_capture', '_IlluminaGA-DNASeq_exome',
            '_DNASeq_whole', '_DNASeq_exome',
            'trimmed.annotated.translated_to_genomic']

# this regular expression will match full TCGA sample Ids, e.g. TCGA-AG-A016-01A-01R or TCGA-37-4133-10A-01D
i_tcgaNameRegEx = re.compile("TCGA-(\\w){2}-(\\w){4}-(\\w){3}-(\\w){3}")


def get_ids(anIdFilename):
    '''
    '  Return all of the Ids found in the given file.
    '  The file is formatted with one line of comma separated Ids.
    '  Returns a set of Ids with the format "TCGA-XX-YYYY" where X 
    '  is a TCGA cancer type such as "AB" and Y is alphanumeric.  
    '  For example, "TCGA-AB-1234".
    '''
    
    i_fileHandler = open(anIdFilename, "r")
    ids = list()
    
    for line in i_fileHandler:
        line = line.rstrip("\r\n")
        ids += line.split(",")
    
    i_fileHandler.close()
    return ids


def getChroms():
    
    #chromIds = range(19, 20)
    #chromIds = map(str, chromIds)
    #chromIds.append("X")
    #chromIds.append("Y")
    #chromIds.append("M")     # if the chrom is "MT", then set the bam specific flag in radia
    
    #chromIds = ["5"]
    #chromIds.append("6")
    
    chromIds = ["6"]
    #chromIds = ["9"]
    
    return chromIds


def runCommand(cmd):
    # run a command in the shell
    p = subprocess.Popen([cmd], shell=True,
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    
    #lines = p.stderr.readlines()
    #lines.extend(p.stdout.readlines())
    lines = p.stdout.readlines()
    return lines


def getDirList(aDirFilename):
    '''
    '  Return a list of directories that are found in the directory file that is given. 
    '''
    dirList = list()
    fileHandler = open(aDirFilename, "r")
    for line in fileHandler:
        line = line.rstrip("\r\n")
        dirList.append(line)   
    fileHandler.close()
    return dirList


def getManifestInfo(aFilename):
    # key=filename, id=legacyId
    idDict = {}
    fileHandler = open(aFilename, "r")
    for line in fileHandler:
        line = line.rstrip("\r\n")
        splitLine = line.split("\t")
        idDict[splitLine[4]] = splitLine[10]   
    fileHandler.close()
    return idDict


def fileHasUniqueKey(aFile, aUniqueKeyList):
    for key in aUniqueKeyList:
        if (key in aFile):
            return True
        
    return False
        

def getBamsFromDirs(aCancerTypesList, aDirList, aWgUniqueKeyList, anExomeUniqueKeyList, anRnaUniqueKeyList):
    '''
    '  Return all of the unique Ids from the files found in the list of directories given.
    '  Returns a set of whole-genome or exome-capture Ids with the format "TCGA-XX-YYYY" where X is a 
    '  TCGA cancer type such as "AB" and Y is alphanumeric.  For example, "TCGA-AB-1234".
    '''

    wgIds = set()
    exomeIds = set()
    rnaIds = set()
    allIds = set()
    # dict[id][analyte][norm/tum][ex/wg/rna][bam]
    # dict[id][dna/rna][10/01][ex/wg/rna][bam]
    allBams = {}
    
    for directory in aDirList:
        dirList = os.listdir(directory)
    
        # for each file in the directory 
        for fileName in dirList:
            
            # if bam isn't in file, then skip it
            if fileName.find('bam') == -1:
                continue
            
            # if the file should be skipped
            toskip = False
            for s in skip:
                if fileName.find(s) != -1:
                    toskip = True
                    break
            if toskip:
                continue
            
            # get an iterator of match objects for the TCGA id
            iterator = i_tcgaNameRegEx.finditer(fileName)
            
            # for each match object in the iterator
            for match in iterator:
                # get the pattern that matched the reg ex, i.e. TCGA-AB-2973-03A-01D-0739-09_whole.bam
                tcgaId = match.group()
                
                # split the id 
                idList = tcgaId.split("-")
                tcga = idList[0]
                cancerType = idList[1]
                patientId = idList[2]
                sample = idList[3][:2]
                analyte = idList[4][2:]
                
                # if there are no cancer types, then add the id
                # if there are cancer types, then make sure the type is in the list
                if (len(aCancerTypesList) == 0 or cancerType in aCancerTypesList):
                    id = tcga + "-" + cancerType + "-" + patientId
                    
                    # add to all ids
                    allIds.add(id) 
                    
                    if (id not in allBams): 
                        allBams[id] = {}
                    
                    if (analyte == "R" or analyte == "T"):
                        analyteKey = "rna"
                    else:
                        analyteKey = "dna"
                    
                    if (analyteKey not in allBams[id]):
                        allBams[id][analyteKey] = {}
                    
                    if (sample not in allBams[id][analyteKey]):
                        allBams[id][analyteKey][sample] = {}
                        
                    if (fileHasUniqueKey(fileName, aWgUniqueKeyList)):
                        wgIds.add(id)
                                                
                        if ("wg" not in allBams[id][analyteKey][sample]):
                            allBams[id][analyteKey][sample]["wg"] = {}    
                        allBams[id][analyteKey][sample]["wg"]["bam"] = os.path.join(directory, fileName)
                        allBams[id][analyteKey][sample]["wg"]["bamFilename"] = fileName
                        
                    if (fileHasUniqueKey(fileName, anExomeUniqueKeyList)):
                        exomeIds.add(id)
                        
                        if ("exome" not in allBams[id][analyteKey][sample]):
                            allBams[id][analyteKey][sample]["exome"] = {}
                        allBams[id][analyteKey][sample]["exome"]["bam"] = os.path.join(directory, fileName)
                        allBams[id][analyteKey][sample]["exome"]["bamFilename"] = fileName
                            
                    if (fileHasUniqueKey(fileName, anRnaUniqueKeyList)):
                        rnaIds.add(id)
                        
                        if ("rna" not in allBams[id][analyteKey][sample]):
                            allBams[id][analyteKey][sample]["rna"] = {}
                        allBams[id][analyteKey][sample]["rna"]["bam"] = os.path.join(directory, fileName)
                        allBams[id][analyteKey][sample]["rna"]["bamFilename"] = fileName
                                
    return (wgIds, exomeIds, rnaIds, allIds, allBams)


def outputIds(aWgIdSet, anExomeIdSet, anRnaIdSet, anAllIdSet, anIdListOutputDir):

    print "WGIds:", len(aWgIdSet)
    print "ExomeIds:", len(anExomeIdSet)
    print "RNAIds:", len(anRnaIdSet)
    print "AllIds:", len(anAllIdSet)
    
    i_rnaWGIntersection = anRnaIdSet.intersection(aWgIdSet)
    i_rnaExomeIntersection = anRnaIdSet.intersection(anExomeIdSet)
    i_rnaDNAUnion = i_rnaWGIntersection.union(i_rnaExomeIntersection)
    #print anExomeIdSet.difference(anRnaIdSet)
    
    print >> sys.stderr, "Intersection of RNA-Seq and whole-genome:", len(i_rnaWGIntersection)
    print >> sys.stderr, "Intersection of RNA-Seq and exome-capture:", len(i_rnaExomeIntersection)
    print >> sys.stderr, "Union of RNA-Seq and DNA:", len(i_rnaDNAUnion)
    
    # print out the ids
    fileHandler = open(os.path.join(anIdListOutputDir, "wgIds.txt"), "w")
    fileHandler.write("\n".join(sorted(aWgIdSet)))   
    fileHandler.close()
    fileHandler = open(os.path.join(anIdListOutputDir, "exomeIds.txt"), "w")
    fileHandler.write("\n".join(sorted(anExomeIdSet)))   
    fileHandler.close()
    fileHandler = open(os.path.join(anIdListOutputDir, "rnaIds.txt"), "w")
    fileHandler.write("\n".join(sorted(anRnaIdSet)))   
    fileHandler.close()
    
    # if both sets aren't empty
    if (aWgIdSet and anExomeIdSet):
        i_dnaIds = aWgIdSet.union(anExomeIdSet)
    # if wg is not empty
    elif (aWgIdSet):
        i_dnaIds = aWgIdSet
    else:
        i_dnaIds = anExomeIdSet
        
    fileHandler = open(os.path.join(anIdListOutputDir, "dnaIds.txt"), "w")
    fileHandler.write("\n".join(sorted(i_dnaIds)))   
    fileHandler.close()
    
    return


def getBamHeaderInfo(anAllBamsDict, aJoblistOutputDir, anRnaUniqueKeysList):
    '''
    '    This isn't as easy as expected.  There is an entry in the header that we can check:
    '    "SO:unsorted" vs. "SO:coordinate"
    '    but this is program dependent, so it might not be accurate.
    '
    '    We will get the following error if someone tries to index a bam that is not sorted (bugfix to samtools version 0.1.12a (r862)):
    '    [bam_index_core] the alignment is not sorted: reads without coordinates prior to reads with coordinates
    '
    '    So, let's just check for "SO:coordinate" for now...
    '''
    
    sortfileHandler = open(os.path.join(aJoblistOutputDir, "sortBams.list"), "w")
        
    for patient in anAllBamsDict.iterkeys():
        for analyte in anAllBamsDict[patient].iterkeys():
            for sample in anAllBamsDict[patient][analyte].iterkeys():
                for coverage in anAllBamsDict[patient][analyte][sample].iterkeys():
                    
                    bam = anAllBamsDict[patient][analyte][sample][coverage]["bam"]
                    
                    cmd = 'samtools view -H %s' % bam

                    sorted = False
                    bamRef = None
                    for line in runCommand(cmd):
                        if line.find("SO:coordinate") != -1:
                            sorted = True
                    
                        for ref in refs:
                            if line.find(ref) != -1:
                                bamRef = ref
                                break
                            
                        if sorted and bamRef is not None:
                            break
                        
                    if (not sorted):
                        #print >> sys.stderr, 'error: bam might not be sorted %s' % bam
                        sortfileHandler.write("samtools sort " + bam + "\n")
                                    
                    # this is a hack for UNC's RNA-Seq files
                    for key in anRnaUniqueKeysList:
                        if key in bam: 
                            bamRef = "GRCh37-lite_w_chr_prefix"
                            break
                            
                    if bamRef is None:
                        print >> sys.stderr, 'error: could not find ref sequence for %s' % bam 
                    else:
                        (prefix, db, mitochon, bamRefFile) = refs[bamRef]
                        anAllBamsDict[patient][analyte][sample][coverage]["ref"] = bamRefFile
                        anAllBamsDict[patient][analyte][sample][coverage]["db"] = db
                        anAllBamsDict[patient][analyte][sample][coverage]["prefix"] = prefix
                        anAllBamsDict[patient][analyte][sample][coverage]["mitochon"] = mitochon
                            
    sortfileHandler.close()
    return anAllBamsDict


def makeIndexAndLinks(anAllBamsDict, aRadiaBasicOutputDir, aRadiaSpecificOutputDir, aJobListOutputDir):
    '''
    '
    '''
    
    linkFileHandler = open(os.path.join(aJobListOutputDir, "symbolicLinks.list"), "w")
    copyIndexFileHandler = open(os.path.join(aJobListOutputDir, "copyIndices.list"), "w")
    indexFileHandler = open(os.path.join(aJobListOutputDir, "indexBams.list"), "w")
        
    for patient in anAllBamsDict.iterkeys():
        for analyte in anAllBamsDict[patient].iterkeys():
            for sample in anAllBamsDict[patient][analyte].iterkeys():
                for coverage in anAllBamsDict[patient][analyte][sample].iterkeys():
                    
                    bam = anAllBamsDict[patient][analyte][sample][coverage]["bam"]
                    bamFilename = anAllBamsDict[patient][analyte][sample][coverage]["bamFilename"]
                    bai = bam + ".bai"
                    
                    # db must exist or something is wrong with the reference
                    if "db" in anAllBamsDict[patient][analyte][sample][coverage]:
                        db = anAllBamsDict[patient][analyte][sample][coverage]["db"]
                    else:
                        print >> sys.stderr, 'error: bam db and reference are unknown %s' % bam
                        continue
                    
                    # if the index file doesn't exist, then we have to create it
                    if not os.path.exists(bai):
                        print >> sys.stderr, 'error: no index file exists for bam %s' % bam
                        indexFileHandler.write("samtools index " + bam + " " + os.path.join(aRadiaBasicOutputDir, db, aRadiaSpecificOutputDir, coverage, bamFilename + ".bai") + "\n")
                    else:
                        # otherwise copy all of the index files to the appropriate sub-dirs
                        copyIndexFileHandler.write("cp " + bai + " " + os.path.join(aRadiaBasicOutputDir, db, aRadiaSpecificOutputDir, coverage) + "\n")
                    
                    # we need to link all the bams to the appropriate sub-dirs
                    # /inside/grotto/users/aradenba/data/hg19/brca/exome/
                    linkFileHandler.write("ln -s " + bam + " " + os.path.join(aRadiaBasicOutputDir, db, aRadiaSpecificOutputDir, coverage) + "\n")
                                                
    linkFileHandler.close()
    copyIndexFileHandler.close()
    indexFileHandler.close()
    
    return


def makeManifestLinks(aManifestDict, aRadiaBasicOutputDir, aRadiaSpecificOutputDir, aJobListOutputDir):
    '''
    '
    '''
    
    linkFileHandler = open(os.path.join(aJobListOutputDir, "symbolicLinks.list"), "w")
    copyIndexFileHandler = open(os.path.join(aJobListOutputDir, "copyIndices.list"), "w")
    indexFileHandler = open(os.path.join(aJobListOutputDir, "indexBams.list"), "w")
        
    for bamFilename in aManifestDict.iterkeys():
        # bamFilename = UNCID_1095596.27e70f84-80c5-4f6c-9d22-83c5a6e66eb5.sorted_genome_alignments.bam
        db = "hg19"
        coverage = "rna"
        tcgaBarcode = aManifestDict[bamFilename]
        tcgaBarcodeFilename = tcgaBarcode + ".sorted_genome_alignments.bam"
        tcgaBarcodeDirFilename = os.path.join(aRadiaBasicOutputDir, db, aRadiaSpecificOutputDir, coverage, tcgaBarcodeFilename)
        
        bamDirFilename = os.path.join("/inside/depot4/bams/brca/rnaseq/", bamFilename)
        bamBaiDirFilename = bamDirFilename + ".bai"
        
        #print bamBaiDirFilename
        
        # if the index file doesn't exist, then we have to create it
        if not os.path.exists(bamBaiDirFilename):
            print >> sys.stderr, 'error: no index file exists for bam %s' % bamDirFilename
            indexFileHandler.write("samtools index " + bamDirFilename + " " + tcgaBarcodeDirFilename + ".bai\n")
        else:
            # otherwise copy all of the index files to the appropriate sub-dirs
            copyIndexFileHandler.write("cp " + bamBaiDirFilename + " " + tcgaBarcodeDirFilename + ".bai\n")
        
        # we need to link all the bams to the appropriate sub-dirs
        # /inside/grotto/users/aradenba/data/hg19/brca/exome/
        linkFileHandler.write("ln -s " + bamDirFilename + " " + tcgaBarcodeDirFilename + "\n")
                                                
    linkFileHandler.close()
    copyIndexFileHandler.close()
    indexFileHandler.close()
    
    return

def makeRadiaList(anAllBamsDict, aRadiaBasicOutputDir, aRadiaSpecificOutputDir, aJobListOutputDir, anIdListOutputDir):
    '''
    '
    '''
    
    radiaFileHandler = open(os.path.join(aJobListOutputDir, "radia.list"), "w")
    
    # 01 = rna tumor tissue = 7
    # 11 = rna normal tissue = 5
    
    # 01 = dna tumorTissue = 479
    # 03 = dna tumorBlood = 0
    
    # 10 = dna normalBlood = 440
    # 11 = dna normalTissue = 70
    
    '''
    normalDNA, normalRNA, tumorDNA, tumorRNA
    normalDNA, tumorDNA, tumorRNA
    
    normalDNA, normalRNA
    tumorDNA, tumorRNA
    
    normalDNA, tumorDNA
    normalRNA, tumorRNA
    '''
    
    dnaTumorSample = "01"
    rnaTumorSample = "01"
    rnaNormalSample = "11"
    
    normalSamples = ["10"]
    coverages = ["exome"]

    includeNormalRNA = False
    includeTumorRNA = True
    includeNormalDNA = True
    includeTumorDNA = True
    
    script = "python2.7 /inside/home/aradenba/rnaEditing/scripts/radia.py"
    chromIds = getChroms()
    overlappingIds = set()
    
    #successfulIds = get_ids("/inside/home/aradenba/rnaEditing/cancer/brca/jobLists/successfulChr5All3.list")
    #successfulIds = get_ids("/inside/home/aradenba/rnaEditing/cancer/brca/jobLists/inProgressChr6.ids")
    #failedIds = get_ids("/inside/home/aradenba/rnaEditing/cancer/brca/jobLists/redoChr6.list")
       
    for patient in anAllBamsDict.iterkeys():
        '''
        # for the crashed Chr 5 and 6 jobs
        shortId = patient[:12]
        if (shortId not in failedIds):
            continue
        '''
        '''
        # for the crashed Chr 5 and 6 jobs
        shortId = patient[:12]
        if (shortId in successfulIds):
            continue
        '''
        
        #print "patient", patient
        preChrOutputList = [script, patient]
        
        normalRNA = None
        tumorRNA = None
        
        if ("rna" in anAllBamsDict[patient] and rnaNormalSample in anAllBamsDict[patient]["rna"]):
            normalRNA = anAllBamsDict[patient]["rna"][rnaNormalSample]["rna"]["bamFilename"]
            normalRNARef = anAllBamsDict[patient]["rna"][rnaNormalSample]["rna"]["ref"]
            normalRNADb = anAllBamsDict[patient]["rna"][rnaNormalSample]["rna"]["db"]
            normalRNAPrefix = anAllBamsDict[patient]["rna"][rnaNormalSample]["rna"]["prefix"]
            normalRNAMitochon = anAllBamsDict[patient]["rna"][rnaNormalSample]["rna"]["mitochon"]
            #print "normalRNA", normalRNA, anAllBamsDict[patient]
            #print "normalRNA", normalRNA
            
        if (includeNormalRNA and normalRNA is None):
            continue
            
        if ("rna" in anAllBamsDict[patient] and rnaTumorSample in anAllBamsDict[patient]["rna"]):
            tumorRNA = anAllBamsDict[patient]["rna"][rnaTumorSample]["rna"]["bamFilename"]
            tumorRNARef = anAllBamsDict[patient]["rna"][rnaTumorSample]["rna"]["ref"]
            tumorRNADb = anAllBamsDict[patient]["rna"][rnaTumorSample]["rna"]["db"]
            tumorRNAPrefix = anAllBamsDict[patient]["rna"][rnaTumorSample]["rna"]["prefix"]
            tumorRNAMitochon = anAllBamsDict[patient]["rna"][rnaTumorSample]["rna"]["mitochon"]
            #print "tumorRNA", tumorRNA, anAllBamsDict[patient]
            
        if (includeTumorRNA and tumorRNA is None):
            continue
        
        if (includeNormalDNA or includeTumorDNA):
            
            if ("dna" in anAllBamsDict[patient]):
                for coverage in coverages:
                    tumorDNA = None
                    if (dnaTumorSample in anAllBamsDict[patient]["dna"] and coverage in anAllBamsDict[patient]["dna"][dnaTumorSample]):
                        tumorDNA = anAllBamsDict[patient]["dna"][dnaTumorSample][coverage]["bamFilename"]
                        tumorDNARef = anAllBamsDict[patient]["dna"][dnaTumorSample][coverage]["ref"]
                        tumorDNADb = anAllBamsDict[patient]["dna"][dnaTumorSample][coverage]["db"]
                        tumorDNAPrefix = anAllBamsDict[patient]["dna"][dnaTumorSample][coverage]["prefix"]
                        tumorDNAMitochon = anAllBamsDict[patient]["dna"][dnaTumorSample][coverage]["mitochon"]
                        tumorDNACoverage = coverage
                        #print "tumorDNA", tumorDNA
                    
                if (includeTumorDNA and tumorDNA is None):
                        continue
                
                for normalSample in normalSamples:
                    for coverage in coverages:
                        normalDNA = None 
                        postChrOutputList = ["--batchSize", "1000000"]
                        
                        if (normalSample in anAllBamsDict[patient]["dna"] and coverage in anAllBamsDict[patient]["dna"][normalSample]):
                            normalDNA = anAllBamsDict[patient]["dna"][normalSample][coverage]["bamFilename"]
                            normalDNARef = anAllBamsDict[patient]["dna"][normalSample][coverage]["ref"]
                            normalDNADb = anAllBamsDict[patient]["dna"][normalSample][coverage]["db"]
                            normalDNAPrefix = anAllBamsDict[patient]["dna"][normalSample][coverage]["prefix"]
                            normalDNAMitochon = anAllBamsDict[patient]["dna"][normalSample][coverage]["mitochon"]
                            #print "normalDNA", normalDNA
                        
                        if (includeNormalDNA and normalDNA is None):
                            continue
                        
                        if (includeNormalDNA and normalDNA is not None):
                            postChrOutputList += ["-n", os.path.join(aRadiaBasicOutputDir, normalDNADb, aRadiaSpecificOutputDir, coverage, normalDNA)]
                            postChrOutputList += ["--dnaNormalFasta", normalDNARef]
                            if (normalDNAPrefix == "chr"):
                                postChrOutputList += ["--dnaNormalUseChr"]
                            
                        # db has to be same, ref can be different
                        if (includeNormalRNA and normalRNA is not None):
                            if (normalRNADb != normalDNADb):
                                print >> sys.stderr, "error: normal DNA db=%s while normal RNA db=%s" % (normalDNADb, normalRNADb)
                                continue
                            #elif (normalRNARef != normalDNARef):
                            #    print >> sys.stderr, "warning: normal DNA ref=%s while normal RNA ref=%s" % (normalDNARef, normalRNARef)
                            
                            postChrOutputList += ["-x", os.path.join(aRadiaBasicOutputDir, normalRNADb, aRadiaSpecificOutputDir, "rna", normalRNA)]
                            postChrOutputList += ["--rnaNormalFasta", normalRNARef]
                            if (normalRNAPrefix == "chr"):
                                postChrOutputList += ["--rnaNormalUseChr"]
                            
                        # db and ref must be the same
                        if (includeTumorDNA and tumorDNA is not None):
                            if (tumorDNADb != normalDNADb):
                                print >> sys.stderr, "error: normal DNA db=%s while tumor DNA db=%s" % (normalDNADb, tumorDNADb)
                                continue
                            elif (tumorDNARef != normalDNARef):
                                print >> sys.stderr, "error: normal DNA ref=%s while tumor DNA ref=%s" % (normalDNARef, tumorDNARef)
                                continue
                            else:
                                postChrOutputList += ["-t", os.path.join(aRadiaBasicOutputDir, tumorDNADb, aRadiaSpecificOutputDir, tumorDNACoverage, tumorDNA)]
                                postChrOutputList += ["--dnaTumorFasta", tumorDNARef]
                                if (tumorDNAPrefix == "chr"):
                                    postChrOutputList += ["--dnaTumorUseChr"]
                                
                        # db has to be same, ref can be different
                        if (includeTumorRNA and tumorRNA is not None):
                            if (tumorRNADb != normalDNADb):
                                print >> sys.stderr, "error: normal DNA db=%s while tumor RNA db=%s" % (normalDNADb, tumorRNADb)
                                continue
                            #elif (tumorRNARef != normalDNARef):
                            #    print >> sys.stderr, "warning: normal DNA ref=%s while tumor RNA ref=%s" % (normalDNARef, tumorRNARef)
                            
                            postChrOutputList += ["-r", os.path.join(aRadiaBasicOutputDir, tumorRNADb, aRadiaSpecificOutputDir, "rna", tumorRNA)]
                            postChrOutputList += ["--rnaTumorFasta", tumorRNARef]
                            if (tumorRNAPrefix == "chr"):
                                postChrOutputList += ["--rnaTumorUseChr"]
                             
                        # at least one dna is needed
                        if (normalDNA is not None or tumorDNA is not None):
                            # get ref and db from one of the dna
                            if (normalDNA is not None):
                                ref = normalDNARef
                                db = normalDNADb
                            else:
                                ref = tumorDNARef
                                db = tumorDNADb
                            
                            postChrOutputList += ["-e", db]
                            postChrOutputList += ["-u", ref]
                            postChrOutputList += ["-m", ref]
                            postChrOutputList += ["-c", os.path.join("/inside/home/aradenba/rnaEditing/data/", db, db + "_chromSizes_sorted.tab")]
                            postChrOutputList += ["-s", os.path.join(aRadiaBasicOutputDir, db, aRadiaSpecificOutputDir, "radia/stats/")]
                            
                            # keep track of the ids that overlap
                            overlappingIds.add(patient)
                            
                            for chrom in chromIds:
                                if (chrom == "M"): 
                                    if (normalDNA is not None and normalDNAMitochon == "MT"):
                                        postChrOutputList += ["--dnaNormalMitochon", "MT"]
                                    if (normalRNA is not None and normalRNAMitochon == "MT"):
                                        postChrOutputList += ["--rnaNormalMitochon", "MT"]
                                    if (tumorDNA is not None and tumorDNAMitochon == "MT"):
                                        postChrOutputList += ["--dnaTumorMitochon", "MT"]
                                    if (tumorRNA is not None and tumorRNAMitochon == "MT"):
                                        postChrOutputList += ["--rnaTumorMitochon", "MT"]    
                                        
                                outputFile = ["-o", os.path.join(aRadiaBasicOutputDir, db, aRadiaSpecificOutputDir, "radia", patient + "_chr" + chrom + ".vcf")]
                                #print " ".join(preChrOutputList) + " " + chrom + " " + " ".join(postChrOutputList) + " " + " ".join(outputFile)
                                '''
                                radiaFileHandler.write("qsub /inside/home/aradenba/rnaEditing/cancer/brca/jobLists/qsub/" + patient + "_chr" + chrom + ".sh\n")
                        
                                tempFileHandler = open("../cancer/brca/jobLists/qsub/" + patient + "_chr" + chrom + ".sh", "w")
                                tempFileHandler.write("#!/bin/bash\n")
                                tempFileHandler.write("source ~/.bashrc\n")
                                tempFileHandler.write(" ".join(preChrOutputList) + " " + chrom + " " + " ".join(postChrOutputList) + " " + " ".join(outputFile) + "\n")
                                tempFileHandler.close()
                                '''
                                radiaFileHandler.write(" ".join(preChrOutputList) + " " + chrom + " " + " ".join(postChrOutputList) + " " + " ".join(outputFile) + "\n")
                    
        elif (normalRNA is not None and tumorRNA is not None):
            
            # db and ref must be the same
            if (normalRNADb != tumorRNADb):
                print >> sys.stderr, "error: normal RNA db=%s while tumor RNA db=%s" % (normalRNADb, tumorRNADb)
                continue
            elif (normalRNARef != tumorRNARef):
                print >> sys.stderr, "error: normal RNA ref=%s while tumor RNA ref=%s" % (normalRNARef, tumorRNARef)
                continue
                        
            postChrOutputList += ["-x", os.path.join(aRadiaBasicOutputDir, normalRNADb, aRadiaSpecificOutputDir, "rna", normalRNA)]
            postChrOutputList += ["-r", os.path.join(aRadiaBasicOutputDir, tumorRNADb, aRadiaSpecificOutputDir, "rna", tumorRNA)]
            
            postChrOutputList += ["-f", normalRNARef]
            postChrOutputList += ["-e", normalRNADb]
            postChrOutputList += ["-u", normalRNARef]
            postChrOutputList += ["-m", normalRNARef]
            postChrOutputList += ["-c", os.path.join("/inside/home/aradenba/rnaEditing/data/", normalRNADb, normalRNADb + "_chromSizes_sorted.tab")]
            postChrOutputList += ["-s", os.path.join(aRadiaBasicOutputDir, normalRNADb, aRadiaSpecificOutputDir, "radia/stats/")]
            if (normalRNAPrefix == "chr"):
                postChrOutputList += ["--useChrPrefix"]  
            
            for chrom in chromIds:
                if (chrom == "M"):
                    if (normalRNAMitochon == "MT"):
                        postChrOutputList += ["--rnaNormalMitochon", "MT"]
                        postChrOutputList += ["--rnaTumorMitochon", "MT"] 
                outputFile = ["-o", os.path.join(aRadiaBasicOutputDir, normalRNADb, aRadiaSpecificOutputDir, "radia", patient + "_chr" + chrom + ".vcf")]
                print " ".join(preChrOutputList) + " " + chrom + " " + " ".join(postChrOutputList) + " " + " ".join(outputFile)
                radiaFileHandler.write(" ".join(preChrOutputList) + " " + chrom + " " + " ".join(postChrOutputList) + " " + outputFile + "\n")    
       
    overlappingIdsFileHandler = open(os.path.join(anIdListOutputDir, "overlappingIds.txt"), "w")
    overlappingIdsFileHandler.write("\n".join(sorted(overlappingIds)))   
    overlappingIdsFileHandler.close()
                                                                 
    radiaFileHandler.close()
    
    return


        
def main(cancerTypes, inputDirsFilename, radiaBasicOutputDir, radiaSpecificOutputDir, listOutputDir):
    jobListOutputDir = os.path.join(listOutputDir, "jobLists")
    idListOutputDir = os.path.join(listOutputDir, "idLists")
    if not os.path.exists(radiaBasicOutputDir):
        print >> sys.stderr, 'error: output directory: %s  does not exist.' % radiaBasicOutputDir
        sys.exit(0)
    if not os.path.exists(jobListOutputDir):
        print >> sys.stderr, 'error: joblist directory: %s  does not exist.' % jobListOutputDir
        sys.exit(0)
    if not os.path.exists(idListOutputDir):
        print >> sys.stderr, 'error: id list directory: %s  does not exist.' % idListOutputDir
        sys.exit(0)

    # get a list of all the dirs to search through
    inDirs = getDirList(inputDirsFilename)
    
    '''
    manifestDict = getManifestInfo("/inside/home/aradenba/rnaEditing/data/manifests/2012-08-21_brca_rnaseq_manifest.tsv")
    #print len(manifestDict.keys())
    
    makeManifestLinks(manifestDict, radiaBasicOutputDir, radiaSpecificOutputDir, jobListOutputDir)
    '''
    
    # hard-coded for now
    wgUniqueKeys = ["_whole"]
    exomeUniqueKeys = ["_capture", "_exome"]
    rnaUniqueKeys = ["trimmed.annotated.translated_to_genomic", "sorted_genome_alignments"]
    
    # get all the patients that belong to this study
    (i_wgIds, i_exomeIds, i_rnaIds, i_allIds, i_allBamsDict) = getBamsFromDirs(cancerTypes, inDirs, wgUniqueKeys, exomeUniqueKeys, rnaUniqueKeys)
        
    # output the ids
    outputIds(i_wgIds, i_exomeIds, i_rnaIds, i_allIds, idListOutputDir)
    print "# of unique patient ids", len(i_allBamsDict.keys())
    
    # make sure all the bams are sorted, otherwise create a joblist to sort them
    i_allBamsDict = getBamHeaderInfo(i_allBamsDict, jobListOutputDir, rnaUniqueKeys)
    #print i_allBamsDict
    
    makeIndexAndLinks(i_allBamsDict, radiaBasicOutputDir, radiaSpecificOutputDir, jobListOutputDir)
    
    # make radia list
    makeRadiaList(i_allBamsDict, radiaBasicOutputDir, radiaSpecificOutputDir, jobListOutputDir, idListOutputDir)
    
    
    
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print "createRadiaJobList.py cancerTypes inputDirsFilename radiaBasicOutputDir radiaSpecificOutputDir joblistOutputDir"
        print "  cancerTypes            = specifying the 2-alphanumeric TCGA cancer types of interest,"
        print "  inputDirsFilename      = file containing list of directories where bam files are stored,"
        print "  radiaDataOutputDir     = basic data directory where radia results will be output,"
        print "  radiaSpecificSubDir    = specific sub-directory where radia results will be output,"
        print "  listOutputDir          = where job and id lists will be output."
        sys.exit(0)

    cancerTypes             = sys.argv[1]
    inputDirsFilename       = sys.argv[2]  
    radiaDataOutputDir      = sys.argv[3]
    radiaSpecificSubDir     = sys.argv[4]
    listOutputDir           = sys.argv[5]
      
    main(cancerTypes, inputDirsFilename, radiaDataOutputDir, radiaSpecificSubDir, listOutputDir)
