#! /usr/bin/env python3.4

########################################################################
#                                                                      #
#                                                                      #
#                                                                      #
# Known issues:                                                        #
#                                                                      #
#                                                                      #
#                                                                      #
# This script was written by Thomas Heavey in 2015.                    #
#        theavey@bu.edu     thomasjheavey@gmail.com                    #
#                                                                      #
# Copyright 2015 Thomas J. Heavey IV                                   #
#                                                                      #
# Licensed under the Apache License, Version 2.0 (the "License");      #
# you may not use this file except in compliance with the License.     #
# You may obtain a copy of the License at                              #
#                                                                      #
#    http://www.apache.org/licenses/LICENSE-2.0                        #
#                                                                      #
# Unless required by applicable law or agreed to in writing, software  #
# distributed under the License is distributed on an "AS IS" BASIS,    #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or      #
# implied.                                                             #
# See the License for the specific language governing permissions and  #
# limitations under the License.                                       #
#                                                                      #
########################################################################

# Import Statements
import os           # Generally good for referencing files

# todo define as a class?
def configmanager(basename, filename=".config"):
    """Function for managing a file for the current state of a
    calculation.
    .config is the file name by default
    It has a structure of keyword (whitespace) value on each line
    lines starting with # are completely ignored, but also unrecognized
    keywords are currently not interpretted as anything (I think)."""
    configexists = checkforconfig(filename)
    # todo print if one exists, and what is being done about it
    pass
    # todo Obviously need to do this



def checkforconfig(filename):
    """Return True if one exists, False if not.
    Uses os"""
    return(os.path.isfile(filename))




def makenewconfig(basename, filename, **starting_values):
    """"""
    with open(filename, 'x') as file:
        file.write('basename    {}\n'.format(basename))
        for key in starting_values:
            file.write('{}    {}'.format(key, starting_values[key]))
    print("Wrote new configuration file to {}".format(filename))



def readconfig(filename):
    """"""
    filedata = dict()
    with open(filename, 'r') as file:
        for line in file:
            # take off leading and trailing whitespace
            line = line.strip()
            # ignore commented lines
            if line.startswith("#"):
                continue
            # add first two parts of each line to the dictionary for
            # output as key: value pairs
            filedata.update(line.split()[0:1])
    print("Read configuration file {}".format(filename))
    return(filedata)




def updateconfig(filename, **added_values):
    """"""
    with open(filename, 'a') as file:
        for key in added_values:
            file.write('{}    {}'.format(key, added_values[key]))
