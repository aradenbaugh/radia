#!/usr/bin/env python2.7

import sys
import os


'''
'    Amie Radenbaugh - 02/09/2011
'    UCSC - RNA Editing
'    Program name: "rnaEditingUtil.py"
'
'    This program provides utility functions that are used from various different scripts:
'
'    check_for_argv_errors:  Error-handling to make sure the directories and filenames provided by the user exist
'    
'''


def check_for_argv_errors(aDirList, aReadFilenameList, aWriteFilenameList):
    '''
    '  Error-handling to make sure the directories and filenames provided by the user exist
    '''
    
    if (aWriteFilenameList != None):    
        for filename in aWriteFilenameList:
            # if the file is nested in a directory
            if (filename.rfind("/") != -1):
                # get the directory and add it to the dirlist
                directory = filename[0:(filename.rfind("/")+1)]
                if (aDirList != None):
                    aDirList += [directory]
                else:
                    aDirList = [directory]    
    
    if (aDirList != None):
        try:
            for directory in aDirList:
                if (not os.path.isdir(directory)):
                    print >> sys.stderr, "Error:  Directory", directory, "does not exist."
                    return False
        except IOError:
            print >> sys.stderr, "Error with the command line arguments."
            return False
    
    if (aReadFilenameList != None):    
        for filename in aReadFilenameList:
            try:
                fileHandler = open(filename, "r")
            except IOError:
                print >> sys.stderr, "Cannot read the file, check to see if the path exists: ", filename
                return False
            else:
                fileHandler.close()
                            
    return True
