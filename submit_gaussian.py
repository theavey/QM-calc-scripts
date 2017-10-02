#! /usr/bin/env python3.4

########################################################################
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

import argparse   # For parsing commandline arguments
import re         # RegEx package for sorting data
import glob       # Allows referencing file system/file names
import readline   # Allows easier file input (with tab completion?)
import subprocess # Allows for submitting commands to the shell
import fileinput  # allows easy iteration over a file

descrip = 'Create and submit a script to run a Gaussian job on SCC'

parser = argparse.ArgumentParser(description=descrip)
parser.add_argument('in_name',
                    help='Name of Gaussian input file')
parser.add_argument('-c', '--numcores', type=int, default=16,
                    help='Number of cores for job')
# I should probably check validity of this time request
# Maybe it doesn't matter so much because it just won't
# submit the job and it will give quick feedback about that?
parser.add_argument('-t', '--time',
                    help='Time required as "hh:mm:ss"',
                    default='12:00:00')
parser.add_argument('-b', '--batch', action='store_true',
                    help='create multiple scripts (batch job)')
parser.add_argument('-x', '--template', default=None,
                    help='template file for creating input from coords')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='make program more verbose')
args = parser.parse_args()

# An input function that can prefill in the text entry
def rlinput(prompt, prefill=''):
   readline.set_startup_hook(lambda: readline.insert_text(prefill))
   try:
      return raw_input(prompt)
   finally:
      readline.set_startup_hook()

def create_gau_input(coord_name, template):
    """This function takes as input a file with a set of molecular
    coordinates (the form should not matter, it will just be copied
    into the next file) and a template file that should be the header
    for the desired calculation (including charge and multiplicity),
    returns nothing, but creates a Gaussian input file ending with
    '.com' """
    out_name = coord_name.rsplit('.', 1)[0] + '.com'
    with open(out_name, 'w') as out_file:
        with open(template, 'r') as templ_file:
            if args.verbose:
                print('opened {}'.format(template))
            for line in templ_file:
                out_file.write(line)
            if '\n' not in line:
                out_file.write('\n')
        with open(coord_name, 'r') as in_file:
            if args.verbose:
                print('opened {}'.format(coord_name))
            for line in in_file:
                if line.strip().isdigit():
                    # the first line is the number of atoms
                    continue
                # XYZ files created by mathematica have a comment
                # as the second line saying something like:
                # "Created by mathematica". Obv. want to ignore that
                if line.strip().startswith('Create'):
                    continue
                # else:
                out_file.write(line)
        out_file.write('\n\n\n')
    if args.verbose:
        print('created Gaussian input file {}'.format(out_name))

in_name_list  =     glob.glob(args.in_name + '*')
in_name_list.sort()        # sort files alphanumerically
in_name_list.sort(key=len) # sort by length (because otherwise would
# put 1,10,11,... as opposed to 1,...,9,10,...
# if number 01,02,... They should all be the same length and the
# second sort won't do anything.

num_files     = len(in_name_list)

yes = ['y', 'yes', '1']

if not args.batch:
    if num_files > 1:
        print('Multiple files starting with {}'.format(args.in_name))
        if input('Did you mean to execute a batch job? ') in yes:
            args.batch = True
        else:
            print('What file name shall I use?')
            in_name_list = [ rlinput('file name: ',args.in_name)]

if args.template:
    for in_name in in_name_list:
        create_gau_input(in_name, args.template)
    in_name_list = glob.glob(args.in_name + '*.com')
    in_name_list.sort()
    in_name_list.sort(key=len)

script_list = []

for in_name in in_name_list:
    if in_name.endswith('.com'):
        short_name = in_name.rsplit('.', 1)[0]
        if not short_name + '.com' == in_name:
            raise SyntaxError('problem interpreting file name. ' +
                              'Period in file name?')
        out_name = short_name + '.out'
    elif '.' in in_name:
        short_name      = in_name.rsplit('.', 1)[0]
        input_extension = in_name.rsplit('.', 1)[-1]
        if not short_name + '.' + input_extension == in_name:
            raise SyntaxError('problem interpretting file name. ' +\
                              'Period in file name?')
        out_name = short_name + '.out'
    else:
        short_name = in_name
        in_name    = short_name + '.com'
        print('Assuming input file is {}'.format(in_name))
        out_name   = short_name + '.out'

    script_name = 'submit' + short_name + '.sh'
    script_list.append(script_name)

    with open(script_name, 'w') as script_file:
        script_file.write('#!/bin/sh \n\n')
        script_file.write('#$ -pe omp {}\n'.format(args.numcores))
        script_file.write('#$ -M theavey@bu.edu\n')
        script_file.write('#$ -m eas\n')
        script_file.write('#$ -l h_rt={}\n'.format(args.time))
        script_file.write('#$ -N {}\n'.format(short_name))
        script_file.write('#$ -j y\n\n')
        script_file.write('INPUTFILE={}\n'.format(in_name))
        script_file.write('OUTPUTFILE={}\n\n'.format(out_name))
        script_file.write('CURRENTDIR=`pwd`\n')
        script_file.write('SCRATCHDIR=/scratch/$USER\n')
        script_file.write('mkdir -p $SCRATCHDIR\n\n')
        script_file.write('cd $SCRATCHDIR\n\n')
        script_file.write('cp $CURRENTDIR/$INPUTFILE .\n\n')
        script_file.write('g09 <$INPUTFILE > $OUTPUTFILE\n\n')
        script_file.write('cp $OUTPUTFILE $CURRENTDIR/.\n\n')
        script_file.write('echo ran in /net/`hostname -s`$SCRATCHDIR\n')
        script_file.write('echo output was copied to $CURRENTDIR\n\n')

    print('script written to {}'.format(script_name))

if not len(script_list) == len(in_name_list):
    # This should never be the case as far as I know, but I would
    # like to make sure everything input gets a script and all the
    # script names are there to be submitted.
    raise IOError('num scripts dif. from num names given')

if args.batch:
    if input('submit all jobs? ') in yes:
        for script in script_list:
            bashCommand = 'qsub {}'.format(script)
            # Don't really know how this works. Copied from 
            # http://stackoverflow.com/questions/4256107/
            # running-bash-commands-in-python
            process = subprocess.Popen(bashCommand.split(),
                                       stdout=subprocess.PIPE)
            output = process.communicate()[0]
            print(output)
    else:
        print('No jobs submitted, but scripts created')
else:
    if input('submit job {}? '.format(script_list[0])) in yes:
        bashCommand = 'qsub {}'.format(script_list[0])
        # Don't really know how this works. Copied from 
        # http://stackoverflow.com/questions/4256107/
        # running-bash-commands-in-python
        process = subprocess.Popen(bashCommand.split(),
                                   stdout=subprocess.PIPE)
        output = process.communicate()[0]
        print(output)
    else:
        print('{} not submitted'.format(script_list))

