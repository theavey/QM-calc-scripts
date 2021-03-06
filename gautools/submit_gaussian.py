#! /usr/bin/env python3

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

# This is written to work with python 3 because it should be good to
# be working on the newest version of python.
from __future__ import print_function

import argparse   # For parsing commandline arguments
import datetime
import glob       # Allows referencing file system/file names
import os
import re
import readline   # Allows easier file input (with tab completion?)
import subprocess # Allows for submitting commands to the shell
from warnings import warn
from thtools import cd, make_obj_dir, save_obj, resolve_path

yes = ['y', 'yes', '1']


# An input function that can prefill in the text entry
# Not sure if this works in 3.5+ because raw_input is gone
def rlinput(prompt, prefill=''):
    readline.set_startup_hook(lambda: readline.insert_text(prefill))
    try:
        return input(prompt)
    finally:
        readline.set_startup_hook()


def _dir_and_file(path):
    warn('_dir_and_file is deprecated. Use os.path.split instead',
         DeprecationWarning)
    if '/' in path:
        rel_dir, f_name = path.rsplit('/', 1)
        rel_dir = rel_dir + '/'
    else:
        rel_dir = ''
        f_name = path
    return rel_dir, f_name


def create_gau_input(coord_name, template, verbose=True):
    """
    make gaussian input file by combining header and coordinates files

    This function takes as input a file with a set of molecular
    coordinates (the form should not matter, it will just be copied
    into the next file) and a template file that should be the header
    for the desired calculation (including charge and multiplicity),
    returns the name of the file, and creates a Gaussian input file ending
    with '.com'

    :param str coord_name: name of file with coordinates in a format
        Gaussian can read
    :param str template: name of file with header for Gaussian calculation
        (up to and including the charge and multiplicity)
    :param bool verbose: If True, some status messages will be printed
        (including file names)
    :return: name of the written file
    :rtype: str
    """
    if verbose:
        print('Creating Gaussian input file...')
    _out_name = coord_name.rsplit('.', 1)[0] + '.com'
    with open(_out_name, 'w') as out_file:
        with open(template, 'r') as templ_file:
            if verbose:
                print('opened {}'.format(template))
            for line in templ_file:
                out_file.write(line)
            if '\n' not in line:
                out_file.write('\n')
        with open(coord_name, 'r') as in_file:
            if verbose:
                print('opened {}'.format(coord_name))
            for i, line in enumerate(in_file):
                if i < 2:
                    # ignore first two lines
                    # number of atoms and the title/comment
                    continue
                # if line.strip().isdigit():
                #     # the first line is the number of atoms
                #     continue
                # # XYZ files created by mathematica have a comment
                # # as the second line saying something like:
                # # "Created by mathematica". Obv. want to ignore that
                # if line.strip().startswith('Create') or
                #         line.strip().startswith('generated'):
                #     continue
                # else:
                out_file.write(line)
        out_file.write('\n\n\n')
    if verbose:
        print('created Gaussian input file {}'.format(_out_name))
    return _out_name


def get_input_files(base_name, batch):
    _in_name_list = glob.glob(base_name + '*')
    _in_name_list.sort()        # sort files alphanumerically
    _in_name_list.sort(key=len) # sort by length (because otherwise would
    # put 1,10,11,... as opposed to 1,...,9,10,...
    # if number 01,02,... They should all be the same length and the
    # second sort won't do anything.
    if not batch:
        num_files = len(_in_name_list)
        if num_files > 1:
            print('Multiple files starting with {}'.format(base_name))
            if input('Did you mean to execute a batch job? ') in yes:
                batch = True
            else:
                print('What file name shall I use?')
                _in_name_list = [rlinput('file name: ', base_name)]
    return _in_name_list, batch


def use_template(template, in_names, verbose):
    made_name_list = []
    for in_name in in_names:
        out_name = create_gau_input(in_name, template, verbose=verbose)
        made_name_list.append(out_name)
        if verbose:
            print('Added {} to files to possibly submit.'.format(out_name))
    _in_name_list = made_name_list
    _in_name_list.sort()
    _in_name_list.sort(key=len)
    return _in_name_list


def write_sub_script(input_name, num_cores=16, time='12:00:00', verbose=False,
                     mem='125', executable='g09',
                     chk_file=None, copy_chk=False,
                     ln_running=None,
                     hold_jid=None, xyz=None, make_xyz=None, make_input=False,
                     ugt_dict=None):
    """
    Write submission script for (Gaussian) jobs for submission to queue

    If make_xyz is not None, the file make_xyz will be checked to exist
    first to make sure to not waste time when missing a necessary input file.

    :param str input_name: Name of the file to use as input
    :param int num_cores: Number of cores to request
    :param str time: Amount of time to request in the format 'hh:mm:ss'
    :param bool verbose: If True, print out some status messages and such
    :type mem: int or str
    :param mem: Minimum amount of memory to request
    :param str executable: Executable file to use for the job

        Example, 'g09', 'g16'

    :param str chk_file: If not None, this file will be copied back after the
        job has completed. If this is not None and make_input is True,
        this will also be passed to use_gen_template.
    :param bool copy_chk: If this is True, the script will attempt to copy
        what should be an existing checkpoint file to the scratch directory
        before running the job. `chk_file` must be not None as well.
    :param str ln_running: If not None, this will be the base name for
        linking the output file to the current directory. If chk_file is not
        None, it will also be linked with the same base name.
    :param str hold_jid: Job on which this job should depend.
        This should be the name of another job in the queuing system.
    :param str xyz: Name of an xyz file to use as input to use_gen_template
        (if make_input is True).
    :param str make_xyz: The name of a file to pass to obabel to be used to
        create an xyz file to pass to use_gen_template.
    :param bool make_input: If True, use_gen_template will be used to create
        input for the Gaussian calculation.
    :param dict ugt_dict: dict of arguments to pass to use_gen_template.

        This should not include out_file, xyz, nproc, mem, or checkpoint
        because those will all be used from other arguments to this function.
        out_file will be input_name; xyz will be xyz or a time-based name if
        make_xyz is not None; nproc will be $NSLOTS (useful if this gets
        changed after job submission); mem will be mem; and checkpoint will
        be chk_file.
    :return: The name of the script file
    :rtype: str
    """
    rel_dir, file_name = os.path.split(input_name)
    if file_name.endswith('.com'):
        short_name = os.path.splitext(file_name)[0]
        if not short_name + '.com' == file_name:
            raise SyntaxError('problem interpreting file name. ' +
                              'Period in file name?')
        out_name = short_name + '.out'
    elif '.' in file_name:
        short_name, input_extension = os.path.splitext(file_name)
        if not short_name + '.' + input_extension == file_name:
            raise SyntaxError('problem interpreting file name. ' +
                              'Period in file name?')
        out_name = short_name + '.out'
    else:
        short_name = file_name
        file_name = short_name + '.com'
        print('Assuming input file is {}'.format(file_name))
        out_name = short_name + '.out'
    job_name = re.match(r'.*?([a-zA-Z].*)', short_name).group(1)
    if len(job_name) == 0:
        job_name = 'default'
    _script_name = os.path.join(rel_dir, 'submit'+short_name+'.sh')
    temp_xyz = os.path.abspath('.temp' +
                               datetime.datetime.now().strftime('%H%M%S%f') +
                               '.xyz')
    if xyz is None or make_xyz is not None:
        n_xyz = temp_xyz
    else:
        n_xyz = resolve_path(xyz)
    temp_pkl = temp_xyz[:-4]
    if ugt_dict is not None:
        make_obj_dir()
        pkl_path = save_obj(ugt_dict, temp_pkl)
    if chk_file is not None:
        chk_line = 'checkpoint=\'{}\','.format(chk_file)
    else:
        chk_line = ''

    with open(_script_name, 'w') as script_file:
        sfw = script_file.write
        sfw('#!/bin/bash -l\n\n')
        sfw('#$ -pe omp {}\n'.format(num_cores))
        sfw('#$ -M theavey@bu.edu\n')
        sfw('#$ -m eas\n')
        sfw('#$ -l h_rt={}\n'.format(time))
        sfw('#$ -l mem_total={}G\n'.format(mem))
        sfw('#$ -N {}\n'.format(job_name))
        sfw('#$ -j y\n')
        sfw('#$ -o {}.log\n\n'.format(short_name))
        if hold_jid is not None:
            sfw('#$ -hold_jid {}\n\n'.format(hold_jid))
        if make_xyz is not None:
            sfw('if [ ! -f {} ]; then\n'.format(
                os.path.abspath(make_xyz)) +
                '    exit 17\n'
                'fi\n\n')
            sfw('module load wxwidgets/3.0.2\n')
            sfw('module load openbabel/2.4.1\n\n')
            sfw('obabel {} -O {}\n\n'.format(os.path.abspath(
                make_xyz), os.path.abspath(n_xyz)))
        if make_input:
            sfw('python -c "from gautools.tools import '
                'use_gen_template as ugt;\n'
                'from thtools import load_obj, get_node_mem;\n'
                'm = get_node_mem();\n'
                'd = load_obj(\'{}\');\n'.format(
                  os.path.abspath(pkl_path)) +
                'ugt(\'{}\',\'{}\','.format(
                    file_name, os.path.abspath(n_xyz)) +
                'nproc=$NSLOTS,mem=m,{}'.format(chk_line) +
                '**d)"\n\n')
        sfw('INPUTFILE={}\n'.format(file_name))
        sfw('OUTPUTFILE={}\n'.format(out_name))
        if chk_file is not None:
            sfw('CHECKFILE={}\n\n'.format(chk_file))
        else:
            sfw('\n')
        if ln_running is not None:
            sfw('WORKINGOUT={}.out\n'.format(ln_running))
            if chk_file is not None:
                sfw('WORKINGCHK={}.chk\n\n'.format(ln_running))
            else:
                sfw('\n')
        sfw('CURRENTDIR=`pwd`\n')
        sfw('SCRATCHDIR=/scratch/$USER\n')
        sfw('mkdir -p $SCRATCHDIR\n\n')
        sfw('cd $SCRATCHDIR\n\n')
        sfw('cp $CURRENTDIR/$INPUTFILE .\n')
        if chk_file is not None:
            sfw('# ') if not copy_chk else None
            sfw('cp $CURRENTDIR/$CHECKFILE .\n\n')
        else:
            sfw('\n')
        if ln_running is not None:
            sfw('ln -s -b /net/`hostname -s`$PWD/$OUTPUTFILE '
                '$CURRENTDIR/$WORKINGOUT\n')
            if chk_file is not None:
                sfw('ln -s -b /net/`hostname -s`$PWD/$CHECKFILE '
                    '$CURRENTDIR/$WORKINGCHK\n\n')
            else:
                sfw('\n')
        sfw('echo About to run {} in /net/`'.format(executable) +
            'hostname -s`$SCRATCHDIR\n\n')
        sfw('{} <$INPUTFILE > $OUTPUTFILE'.format(executable))
        sfw('\n\n')
        if ln_running is not None:
            sfw('rm $CURRENTDIR/$WORKINGOUT')
            if chk_file is not None:
                sfw(' $CURRENTDIR/$WORKINGCHK\n\n')
            else:
                sfw('\n\n')
        sfw('cp $OUTPUTFILE $CURRENTDIR/.\n')
        if chk_file is not None:
            sfw('cp $CHECKFILE $CURRENTDIR/.\n\n')
        else:
            sfw('\n')
        sfw('echo ran in /net/`hostname -s`$SCRATCHDIR\n')
        sfw('echo output was copied to $CURRENTDIR\n\n')

    if verbose:
        print('script written to {}'.format(_script_name))
    return _script_name


def submit_scripts(scripts, batch=False, submit=False, verbose=False):
    outputs = []
    if batch:
        if submit or input('submit all jobs? ') in yes:
            for script in scripts:
                rd, f = _dir_and_file(script)
                with cd(rd, ignore_blank=True):
                    cl = ['qsub', f]
                    # Don't really know how this works. Copied from
                    # http://stackoverflow.com/questions/4256107/
                    # running-bash-commands-in-python
                    process = subprocess.Popen(cl,
                                               stdout=subprocess.PIPE,
                                               universal_newlines=True)
                    output = process.communicate()[0]
                if verbose:
                    print(output)
                outputs.append(output)
        else:
            if verbose:
                print('No jobs submitted, but scripts created')
    else:
        if submit or input('submit job {}? '.format(scripts[0])) in yes:
            rd, f = _dir_and_file(scripts[0])
            with cd(rd, ignore_blank=True):
                cl = ['qsub', f]
                # Don't really know how this works. Copied from
                # http://stackoverflow.com/questions/4256107/
                # running-bash-commands-in-python
                process = subprocess.Popen(cl,
                                           stdout=subprocess.PIPE,
                                           universal_newlines=True)
                output = process.communicate()[0]
            if verbose:
                print(output)
            outputs.append(output)
        else:
            if verbose:
                print('{} not submitted'.format(scripts))
    _job_info = [' '.join(output.split(' ')[2:4]) for output in outputs]
    return _job_info


if __name__ == '__main__':
    description = 'Create and submit a script to run a Gaussian job on SCC'

    parser = argparse.ArgumentParser(description=description)
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
    parser.add_argument('-e', '--executable', type=str, default='g09',
                        help='name of executable to run')
    parser.add_argument('-b', '--batch', action='store_true',
                        help='create multiple scripts (batch job)')
    parser.add_argument('-x', '--template', default=None,
                        help='template file for creating input from coords')
    parser.add_argument('-s', '--submit', action='store_true',
                        help='Automatically submit jobs?')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='make program more verbose')
    parser.add_argument('-j', '--nojobinfo', action='store_false',
                        help='Do not return the submitted job information')
    parser.add_argument('-k', '--chk_file', default=None,
                        help='checkpoint file to be written and copied back')
    parser.add_argument('--copy_chk', action='store_true',
                        help='Copy check file to the scratch directory')
    parser.add_argument('-l', '--ln_running', type=str, default=None,
                        help='base name for linking output to cwd while '
                             'running')
    parser.add_argument('-d', '--hold_jid', default=None,
                        help='job on which this job should depend')
    args = parser.parse_args()

    in_name_list, args.batch = get_input_files(args.in_name, args.batch)
    if args.template:
        in_name_list = use_template(args.template, in_name_list, args.verbose)
    script_list = []
    for in_name in in_name_list:
        script_name = write_sub_script(input_name=in_name,
                                       num_cores=args.numcores,
                                       time=args.time,
                                       verbose=args.verbose,
                                       executable=args.executable,
                                       chk_file=args.chk_file,
                                       copy_chk=args.copy_chk,
                                       ln_running=args.ln_running,
                                       hold_jid=args.hold_jid)
        script_list.append(script_name)
    if not len(script_list) == len(in_name_list):
        # This should never be the case as far as I know, but I would
        # like to make sure everything input gets a script and all the
        # script names are there to be submitted.
        raise IOError('num scripts dif. from num names given')
    job_info = submit_scripts(script_list, args.batch, args.submit,
                              args.verbose)
    if job_info and args.nojobinfo:
        for job in job_info:
            print(job)
    if args.verbose:
        print('Done. Completed normally.')
