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
import glob       # Iteration over files in directory


base_name     =     sys.argv[1]
spin_mult     =     sys.argv[2]
rem_file      =     sys.argv[3]
in_name_list  =     glob.glob(base_name + '*.xyz')

in_name_list.sort()
in_name_list.sort(key=len)

print('base file name is {}'.format(base_name))
print('number of files is {:,}'.format(len(in_name_list)))
print('spin multiplicity is {}'.format(spin_mult))
print('rem file name is {}'.format(rem_file))


for i, in_name in enumerate(in_name_list, start=1):
    out_name = base_name + str(i) + '.inp'
    with open(out_name, 'w') as out_file:
        with open(in_name, 'r') as in_file:
            print('opened {}'.format(in_name))
            out_file.write('$molecule\n')
            out_file.write(spin_mult + '\n')
            for line in in_file:
                # First line of file, not useful
                # Note, in most places, these lines are "stripped" because
                # they come with leading spaces, which messes with
                # startswith and the split function, too I think.
                if line.strip().isdigit():
                    # the first line is the number of atoms
                    continue
                if line.strip().startswith('Create'):
                    continue
                # else:
                out_file.write(line)
        out_file.write('\n$end\n\n')
        with open(rem_file, 'r') as rem_code:
            for line in rem_code:
                out_file.write(line)
        out_file.write('\n\n')



                
