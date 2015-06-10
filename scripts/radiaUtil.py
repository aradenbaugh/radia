#!/usr/bin/env python

import sys
import os


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
            finally:
                fileHandler.close()
                            
    return True
