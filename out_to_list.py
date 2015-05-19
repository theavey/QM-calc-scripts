#! /usr/bin/env python3.4

########################################################################
#                                                                      #
#                                                 #
#                                                                      #
# Known issues:                                                        #
# None                                                                 #
#                                                                      #
# .                           #
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




# This is written to work with python 3.4 because it should be good to
# be working on the newest version of python.

import fileinput  # allows easy iteration over a file
import sys        # For importing the arguments given
import re         # RegEx package for sorting data
import os.path    # Allows for checking if file already exists
import glob       # Allows referencing file system/file names

# Really should adapt to use argparse for parsing command line
# arguments. see https://docs.python.org/2/library/argparse.html
# Probably particularly useful for setting a flag for method
# or multistate input files. I can try to detect, and just
# take it as a command line argument for ease. May not be any
# faster as I already iterate through every line of the file.


base_name     =     sys.argv[1]
in_name_list  =     glob.glob(base_name + '*.out')
num_files     = len(in_name_list)

in_name_list.sort()        # sort files alphanumerically
in_name_list.sort(key=len) # sort by length (because otherwise would
# out 1,10,11,... as opposed to 1,...,9,10,...
# if number 01,02,... They should all be the same length and the
# second sort won't do anything.

print('base file name is {}'.format(base_name))
print('number of files is {:,}'.format(num_files))
print('Energies of {}'.format(base_name))


out_name = base_name + '_energies'

# Define methods I might use that change how to process file.
methods = {'wb97xd': 'dft', 'mp2': 'mp2', 'pbe50': 'sfdft',
           'eom-ccsd': 'sf', 'hf': 'hf', 'pbe1pbe': 'gdft'}
# reading more about regex, I really should use them. I could
# define how it matches based on the method, as opposed to all
# the If Then loops I currently have.
# I could essentially have a dictionary of regular expressions and just
# reference the appropriate one for the type of output it's reading.
# Alternatively, I could just break it up into a lot of subroutines
# that could be called depending on the type of output.

i = 0
# Initialize not knowing method used: "Not Yet Determined"
method = 'nyd'

# Gaussian method regex:
# finds alphanumeric character string before a / and after a
# whitespace character. Requires that the input be in the form
# (method)/(something, normally basis set).
# Are there methods/functionals with non-alphanumeric characters?
# Maybe manual hybrid functions would be bad with this. I don't
# currently remember their syntax.
gregex = re.compile(r"\s+\w+/")
# QChem SF methods energy regex:
sfenergy = re.compile('[=:]\s*-\d*\.\d*')

with open(out_name, 'w') as out_file:
    for name in in_name_list:
        energy_list = []
        with open(name, 'r') as in_file:
            for line in in_file:
                line = line.lower().strip()
                if method is 'nyd':
                    if 'method' in line or
                      'exchange' in line:
                        # make the line into a list split by spaces
                        linelist = re.split(' +',line)
                        # could maybe shorten these next lines with
                        # a Try Except construction
                        if linelist[-1] in methods:
                            method = methods[linelist[-1]]
                        else:
                            print('Unknown method {} used'.format(
                                linelist[-1]))
                            print('Assuming output formatted as HF')
                            method = 'hf'
                    if 'entering gaussian sys' in line:
                        method = 'gaugen'
                        # Gaussian output file, method undetermined
                    continue
                if method is 'gaugen':
                    if line.startswith('#'):
                        gmethodmatch = gregex.search(line)
                        gmethod =re.search(r"\w+", gmethodmatch.group())
                        if re.search(r"\w+", gmethod.group()):
                            try:
                                method = methods[gmethod]
                            except KeyError:
                                print('unknown Gaussian method. ' +
                                      'Assuming (g)dft.')
                                method = 'gdft'
                        if re.search(r'\btd\s*[(=]', line)
                    continue
                if method is 'dft':
                    if line.startswith('total energy'):
                        # make the line into a list split by spaces
                        linelist = re.split(' +',line)
                        out_file.write(linelist[-1])
                        out_file.write('\n')
                    continue
                if method is 'mp2':
                    if 'total energy' in line:
                        # make the line into a list split by spaces
                        linelist = re.split(' +',line)
                        out_file.write(linelist[-2])
                        out_file.write('\n')
                    continue
                if method is 'sf':
                    if 'otal energy' in line:
                        if 'basis set' in line:
                            # Ignore HF energy
                            continue
                        # use RegEx to find energy in the line:
                        match = sfenergy.search(line)
                        energy_list.append(match.group()[2:])
                    continue
                if method is 'sfdft':
                    if 'otal energy' in line:
                        # use RegEx to find energy in the line:
                        match = sfenergy.search(line)
                        energy_list.append(match.group()[2:])
                    continue
        if energy_list:
            # Only true if not empty
            out_file.write(str(energy_list) + '\n')
        i += 1

print("Opened {0} files, and wrote data to {1}".format(i, out_name))


                
