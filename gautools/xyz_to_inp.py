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


base_name  =     sys.argv[1]
num_files  = int(sys.argv[2])
spin_mult  =     sys.argv[3]
rem_file   =     sys.argv[4]

print('base file name is {}'.format(base_name))
print('number of files is {:,}'.format(num_files))
print('spin multiplicity is {}'.format(spin_mult))
print('rem file name is {}'.format(rem_file))

in_name_list  = []
out_name_list = []

i = 1

while i < num_files + 1:
    in_name  = base_name + str(i) + ".xyz"
    out_name = base_name + str(i) + ".inp"
    in_name_list.append(in_name)
    out_name_list.append(out_name)
    i += 1

i = 0

while i < num_files:
    with open(out_name_list[i], 'w') as out_file:
        with open(in_name_list[i], 'r') as in_file:
            print('opened {}'.format(in_file))
            out_file.write('$molecule\n')
            out_file.write(spin_mult + '\n')
            for line in in_file:
                # First line of file, not useful
                # Note, in most places, these lines are "stripped" because
                # they come with leading spaces, which messes with
                # startswith and the split function, too I think.
                if line.strip().startswith('12'):
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
    i += 1


                
