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
import glob


base_name     =     sys.argv[1]
in_name_list  =     glob.glob(base_name + '*.out')
num_files     = len(in_name_list)

in_name_list.sort()
in_name_list.sort(key=len)

print('base file name is {}'.format(base_name))
print('number of files is {:,}'.format(num_files))
print('Energies of {}\n'.format(base_name))


out_name = base_name + '_energies'

i = 0

with open(out_name, 'w') as out_file:
    for name in in_name_list:
        energy_list = []
        with open(name, 'r') as in_file:
            for line in in_file:
                # First line of file, not useful
                # Note, in most places, these lines are "stripped" because
                # they come with leading spaces, which messes with
                # startswith and the split function, too I think.
                if line.strip().startswith('Total energy'):
                    # make the line into a list split by spaces
                    linelist = re.split(' {1,}',line.strip())
                    energy_list.append(linelist[-1])
        i += 1
        if energy_list:
            # Only true if not empty
            out_file.write(str(energy_list) + '\n')

print("Opened {0} files, and wrote data to {1}".format(i, out_name))


                
