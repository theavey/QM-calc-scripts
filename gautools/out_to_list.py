#! /usr/bin/env python3.4

########################################################################
#                                                                      #
#                                                                      #
#                                                                      #
# Known issues:                                                        #
# Methods which are declared with non-alphanumeric characters          #
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




# This is written to work with python 3.4 because it should be good to
# be working on the newest version of python.

import argparse  # For parsing commandline arguments
import glob  # Allows referencing file system/file names
import re  # RegEx package for sorting data

from . import geomRegex


# Really should adapt to use argparse for parsing command line
# arguments. see https://docs.python.org/2/library/argparse.html
# Probably particularly useful for setting a flag for method
# or multistate input files. I can try to detect, and just
# take it as a command line argument for ease. May not be any
# faster as I already iterate through every line of the file.

def outtolist():
    """This function reads command line arguments (use -h for help)
    to take output from quantum chemistry calculations and put useful
    data such as energy and geometries into (a) separate file(s)."""
    descrip = ('This function takes a base file name for output files and makes'
               ' a file basename_energies that is a list of the energies from '
               'the the read output files.')

    parser = argparse.ArgumentParser(description=descrip)
    parser.add_argument('base_name', help='base name of files to read')
    # todo use argparse to check for multiple positional arguments?
    # If I pass it as list set it all as in_name_list?
    parser.add_argument('-m', '--method',
                        help=('calculation method (changes how files '
                              'are interpreted)'))
    parser.add_argument('-s', '--singlefile', action='store_true',
                        help=('use if output is single file with '
                              'multiple geometries'))
    parser.add_argument('-g', '--geometries', action='count', default=0,
                        help=('Flag for creating file of geometries in '
                              'XYZ style format.'
                              '\n-g for stationary points, -gg for all.'))
    # todo add option for saving all energies -e?
    # todo add flag for gaussian vs. qchem?
    # maybe to separate file with different name? Probably not
    parser.add_argument('-a')
    args = parser.parse_args()

    base_name     =     args.base_name
    in_name_list  =     glob.glob(base_name + '*out')
    num_files     = len(in_name_list)

    in_name_list.sort()        # sort files alphanumerically
    in_name_list.sort(key=len) # sort by length (because otherwise would
    # out 1,10,11,... as opposed to 1,...,9,10,...
    # if number 01,02,... They should all be the same length and the
    # second sort won't do anything.

    yes = ['y', 'yes', '1']
    # process in case of single file
    if args.singlefile:
        if num_files == 0:
            in_name_list = glob.glob(base_name)
            num_files    = len(in_name_list)
            print('number of files is {}'.format(num_files))
        if num_files == 1:
            in_name = in_name_list[0]
            print('(single) file name is {}'.format(in_name))
        if num_files > 1:
            print('More than one file matched input file name.')
            if input('Get energies from multiple output files?') in yes:
                args.singlefile = False
            else:
                raise SyntaxError('Okay, well files that matched your input '
                                  'are {}'.format(in_name_list))
    # Print some stuff if not a single file. Note, not "else" because
    # I might set singlefile to false above.
    if not args.singlefile:
        print('base file name is {}'.format(base_name))
        print('number of files is {:,}'.format(num_files))
    # todo add energies to comment line of geoms file
    # will be difficult/require changes, but would be good, I think
    if args.geometries == 1:
        # todo Need to figure this out obvs!!!!!!
        raise NotImplementedError('Not able to pull out single geom yet.')
    # Collect all geometries from files for -gg (or more)
    if args.geometries > 1:
        geoms_out_name = base_name + '_allgeoms'
        with open(geoms_out_name, 'w') as geom_file:
            if args.singlefile:
                geoms = geomRegex.findallgeoms(in_name)
                for geom in geoms:
                    num_atoms = len(geom)
                    geom_file.write(str(num_atoms) + '\n')
                    geom_file.write('geometry from {}\n'.format(base_name))
                    for atom in geom:
                        geom_file.write(atom + '\n')
            else:
                for file_name in in_name_list:
                    geoms = geomRegex.findallgeoms(file_name)
                    for geom in geoms:
                        num_atoms = len(geom)
                        geom_file.write(str(num_atoms) + '\n')
                        geom_file.write('geometry from {}\n'.format(base_name))
                        for atom in geom:
                            geom_file.write(atom + '\n')

    out_name = base_name + '_energies'

    # Define methods I might use that change how to process file.
    methods = {'wb97xd': 'dft', 'mp2': 'mp2', 'pbe50': 'sfdft',
               'eom-ccsd': 'sf', 'hf': 'hf', 'pbe1pbe': 'gdft',
               'cis': 'gcis'}
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
    # If given as commandline argument, set with the given method
    if args.method:
        try:
            method = methods[args.method.lower()]
        except KeyError:
            print("method {} unrecognized. Going to ".format(args.method) +
                   "try to find method based on output file.")

    # Gaussian method regex:
    # finds non-whitespace characters before a / and after a
    # whitespace character. Requires that the input be in the form
    # (method)/(something, normally basis set).
    gregex = re.compile(r"\s+\S+/")
    # QChem SF methods energy regex:
    sfenergy = re.compile(r'[=:]\s*-\d*\.\d*')
    # Gaussian DFT methods ground state energy regex:
    # (could use sfenergy, or combine into one)
    gdftenergy = re.compile(r'=\s*-\d*\.\d*')
    # Gaussian TD-DFT methods excited state energies regex:
    tdgdftenergy = re.compile(r'\s-*\d+\.\d+\s+ev')
    # Gaussian CIS ground state energy regex:
    # Goes to end of the line because given in scientific notation.
    gcisenergy = re.compile(r'eump2\s*=\s*-.+')
    # todo write geometry regex
    # todo find some way to find "stationary point" and process as needed


    # todo add geometry getting to all this multifile stuff
    with open(out_name, 'w') as out_file:
        for name in in_name_list:
            energy_list = []
            with open(name, 'r') as in_file:
                for line in in_file:
                    line = line.lower().strip()
                    if method is 'nyd':
                        if 'method' in line or 'exchange' in line:
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
                            # todo do this with only one regex using groups
                            # and non-capturing groups
                            if gmethodmatch:
                                gmethod =re.search(r"\w+", gmethodmatch.group())
                            if gmethod:
                                try:
                                    method = methods[gmethod.group()]
                                except KeyError:
                                    print('unknown Gaussian method. ' +
                                          'Assuming (g)dft.')
                                    method = 'gdft'
                            if re.search(r'\btd\s*[(=]', line):
                                method = 'td' + method
                                # Note, this will cause problems if TD is
                                # declared on a line before the
                                # functional/method.
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
                    if method.endswith('gdft'):
                        # Ground state energy search for (TD)DFT with Gauss
                        if 'scf done' in line:
                            match = gdftenergy.search(line)
                            energy_list.append(match.group()[2:])
                            continue
                        # Note on this line below: because I don't set the
                        # method name for TD methods in one step, the "is"
                        # comparison here will fail, because they point at
                        # different places, but the (slower) equality
                        # comparison will work because it will go through
                        # and actually compare each character.
                        if method == 'tdgdft':
                            if line.startswith('excited state'):
                                match = tdgdftenergy.search(line)
                                if match:
                                    energy_list.append(match.group()[:-3])
                        continue
                    if method is 'gcis':
                        if 'eump2' in line:
                            match = gcisenergy.search(line)
                            energy_list.append(match.group()[8:])
                            continue
                        if line.startswith('excited state'):
                            match = tdgdftenergy.search(line)
                            if match:
                                energy_list.append(match.group()[:-3])
                        continue
            # if energy_list:
                # Only true if not empty
            if True:
                # Using for now because some error blank files
                # should still produce a line in the output, even if blank
                out_file.write(str(energy_list) + '\n')
            i += 1
    # todo save files with desired information

    print("Opened {0} files, and wrote data to {1}".format(i, out_name))
    print('Files processed for {} method.'.format(method))
    try:
        # will only work if geoms_out_name defined above
        print('geometries written to {}'.format(geoms_out_name))
    except NameError:
        print('No geometries saved')
    # todo write statements about files that were saved

if __name__ == "__main__":
    outtolist()
                
