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

import glob
import argparse
import subprocess
import readline
import os
import shutil
from datetime import datetime

__version__ = 0.1

parser = argparse.ArgumentParser(description='Use PLUMED utility to sum '
                                             'HILLS and then put into '
                                             'Mathematica friendly format')
parser.add_argument('-H', '--hills', default='HILLS',
                    help='name of the HILLS file')
parser.add_argument('-s', '--stride', default=10000,
                    help='specify the stride for integrating '
                         'hills file (default 1000)')
# todo figure out how to get -pi or -3.14,-3.14 to work with --min
parser.add_argument('-i', '--min', default='',
                    help='the lower bounds for the grid')
parser.add_argument('-a', '--max', default='',
                    help='the upper bounds for the grid')
parser.add_argument('-b', '--bin', default='',
                    help='the number of bins for the grid')
parser.add_argument('-p', '--spacing', default='',
                    help='The spacing for the grid. Currently unused!')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='make script more verbose')
parser.add_argument('-f', '--folder', default='SumHills',
                    help='Folder in which this will be run. Can be '
                         'deleted automatically using -c 3.')
parser.add_argument('-t', '--temp_file', default='temp_data_file.m',
                    help='File in which to store all the data')
parser.add_argument('-n', '--var_name', default='summedHills',
                    help='Name of variable to be assigned for Mathematica')
parser.add_argument('-T', '--template', default='sumHillsTempl.m',
                    help='Output template file')
parser.add_argument('-o', '--output_name', default='summedHills.m',
                    help='Name of the file to be output')
parser.add_argument('-e', '--exists', action='store_true',
                    help='Use this argument if the fes data already exists')
parser.add_argument('-c', '--clean', type=int, default=2,
                    help='Argument for how much to clean up\n'
                         '0 does not delete or move anything\n'
                         '>0 moves output to starting folder '
                         'and deletes copy of HILLS file\n'
                         '>1 deletes temp data file\n'
                         '>2 deletes temp folder and contents\n'
                         'default is 2')
parser.add_argument('--version', action='version',
                    version='%(prog)s v{}'.format(__version__))
args = parser.parse_args()


# An input function that can prefill in the text entry
def rlinput(prompt, prefill=''):
    readline.set_startup_hook(lambda: readline.insert_text(prefill))
    try:
        # This was raw_input, but that's deprecated in py3
        return input(prompt)
    finally:
        readline.set_startup_hook()


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def run_plumed_sum_hills():
    """This function takes no arguments and returns nothing.
    Using the command line values stored in args, it will run the PLUMED
    sum_hills utility."""
    command_line = ['mpirun', '-n', '1', 'plumed', 'sum_hills']
    command_line += ['--hills', args.hills]
    if args.verbose:
        print('hills file name is {}'.format(args.hills))
    command_line += ['--stride', str(args.stride)]
    if args.verbose:
        print('data output stride is {}'.format(args.stride))
    minim = args.min
    maxim = args.max
    # todo add option for spacing instead of bin
    if args.bin:
        command_line += ['--bin', str(args.bin)]
        if args.verbose:
            print('number of bins is {}'.format(args.bin))
        if minim or maxim:
            if not (minim and maxim and args.bin):
                print('If you give a min or max, you need min, '
                      'max, and bin')
                minim = rlinput('min = ', minim)
                maxim = rlinput('max = ', maxim)
            if args.verbose:
                print('min: {}, max: {}'.format(minim, maxim))
            command_line += ['--min', str(minim), '--max', str(maxim)]
    else:
        if minim or maxim:
            print('If you give a min or max, you need min, '
                  'max, and bin')
            command_line += ['--bin', str(rlinput('bin = '))]
            if not (minim and maxim):
                minim = rlinput('min = ', minim)
                maxim = rlinput('max = ', maxim)
            if args.verbose:
                print('min: {}, max: {}'.format(minim, maxim))
            command_line += ['--min', str(minim), '--max', str(maxim)]
    command_line_str = ' '.join(command_line)
    print('command line argument is:\n{}'.format(command_line_str))
    print('Running PLUMED sum_hills utility...')
    # Run the PLUMED sum_hills utility and save/print the output as it comes
    with open('plumed_output.log', 'w') as log_file:
        with subprocess.Popen(command_line,
                              stdout=subprocess.PIPE, bufsize=1,
                              universal_newlines=True) as proc:
            for line in proc.stdout:
                log_file.write(line)
                if args.verbose:
                    print(line, end='')
    print('Done running PLUMED sum_hills utility')


def setup_folder():
    """This function takes no arguments and returns nothing.
    It will create a folder for making all this data, and
    it will copy the specified hills file there."""
    working_folder = args.folder
    # Make the working folder, don't raise error if it already exists
    os.makedirs(working_folder, exist_ok=True)
    copyout = shutil.copy(args.hills, working_folder)
    if args.verbose:
        print('HILLS file copied to {}'.format(copyout))
    os.chdir(working_folder)


def read_plumed_stuff():
    """This function takes no arguments and returns nothing.
    It will save the output data from plumed to a formatted
    temporary file that can then be read to put the data
    into the file that Mathematica can read.
    If defines the global variable num_of_cvs."""
    global num_of_cvs, formatted_data, fes_file_names
    print('Reading PLUMED output files')
    fes_file_names = glob.glob('fes*.dat')
    # Make sure the list is in the proper order:
    fes_file_names.sort()         # sort files alphanumerically
    fes_file_names.sort(key=len)  # sort files by length
    # Find number of CVs:
    #  At least in the current implementation of PLUMED, the output is
    #  the list of CV coordinates, then the height there, then the
    #  derivative with respect to each of the CVs, hence the
    #  (number of fields - 1) / 2.
    with open(fes_file_names[0], 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            num_fields = len(line.split())
            num_of_cvs = (num_fields - 1) / 2
            if num_of_cvs.is_integer():
                num_of_cvs = int(num_of_cvs)
            else:
                print('number of CVs found to be {}!'.format(num_of_cvs))
                num_of_cvs = rlinput('Real number of CVs = ',
                                     str(int(num_of_cvs)))
            break
    all_data = []
    for file in fes_file_names:
        f_data = []
        if args.verbose:
            print('Reading file {}'.format(file))
        with open(file, 'r') as crf:
            l_data = []
            for line in crf:
                if line.startswith('#'):
                    continue
                try:
                    if is_number(line.split()[0]):
                        l_data += [', '.join(line.split()[0:(num_of_cvs+1)])]
                except IndexError:
                    # Blank lines in files have length 0
                    continue
            f_data += ['},\n{'.join(l_data)]
        all_data += f_data
    formatted_data = '{{{' + '}},\n{{'.join(all_data) + '}}}\n\n'
    with open(args.temp_file, 'w') as tf:
        tf.write(formatted_data)
    print('Done reading PLUMED output data')


def data_into_mfile():
    """This function takes no arguments and returns nothing.
    It will take data saved to a temporary file from read_plumed_stuff
    and put it into the template .m file so that it can be read into
    Mathematica"""
    print('Putting data into output file...')
    about_content = []
    about_content += ['Number of CVs: {}'.format(num_of_cvs)]
    about_content += ['Number of points per time chunk: '
                      '{}'.format(args.stride)]
    about_content += ['Originally processed on {}'.format(datetime.now())]
    about_content += ['Processed with '
                      '{} v{}'.format(os.path.basename(__file__),
                                      __version__)]
    replacements = dict(varname=args.var_name, data=formatted_data,
                        numcvs=num_of_cvs, stride=args.stride)
    if args.spacing:
        # todo change this when spacing option implemented
        about_content += ['Grid spacing (unused): {}'.format(args.spacing)]
        replacements['spacing'] = args.spacing
    else:
        replacements['spacing'] = '(Print["getGridSize not currently' \
                                  'defined"]; $Failed'
    about = '{' + ', '.join(about_content) + '}'
    replacements['about'] = about
    print(replacements.keys())
    with open(args.template, 'r') as template, \
            open(args.output_name, 'w') as output:
        for line in template:
            if line.startswith('#'):
                # Remove the '#'
                line = line[1:]
                try:
                    # Unpack the dict, then do the replacements
                    output.write(line.format(**replacements))
                except KeyError as e:
                    print('Error! Key {} not found!'.format(e))
                    choice = input('(a)bort or (s)kip? ')
                    if choice in 'abort':
                        raise e
                    if choice in 'skip':
                        continue
            else:
                output.write(line)
    print('Output saved as {}'.format(args.output_name))


def clean_up():
    """This function takes no arguments and returns nothing.
    It will ask if the temp data should be deleted, and if so, will clean
    it all up. It can also move the output file back to the original
    directory.
    Default is to move output file to the starting directory, removes
    the HILLS file if copied, and removes the temp data file. Use -c
    argument to change this behavior.
    """
    print('Cleaning up...')
    if args.clean > 0:
        if args.verbose:
            print('Copying {} to {}...'.format(args.output_name, current_dir))
        shutil.copy(args.output_name, current_dir)
        if not args.exists:
            # If HILLS file was copied (only true if args.exists is false)
            # delete the copy
            if args.verbose:
                print('Removing {}...'.format(args.hills))
            os.remove(args.hills)
    if args.clean > 1:
        if args.verbose:
            print('Removing {}...'.format(args.temp_file))
        os.remove(args.temp_file)
    if args.clean > 2:
        temp_folder = current_dir + '/' + args.folder
        if args.verbose:
            print('Removing {} and contents...'.format(temp_folder))
        shutil.rmtree(temp_folder)
    print('Done cleaning up files/folders.')


current_dir = os.getcwd()
if not args.exists:
    setup_folder()
    run_plumed_sum_hills()
read_plumed_stuff()
data_into_mfile()
clean_up()
print('Done!')
