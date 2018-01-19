"""
A set of tools for working with computational chemistry files and such
"""


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

import re
from six import string_types
from paratemp import copy_no_overwrite


def fix_atom_names_xyz(xyz, dry_run=False, return_list=False):
    """
    Fix atomic symbols in xyz file

    VMD (sometimes?) writes XYZ files with the atomic symbols replaced with
    the atom names, which cause issues with other programs interpreting them.
    This should fix them by leaving only the first letter of each atom
    name/symbol.

    !! NOTE !! This is only written to work for C, H, O, N, and F currently !!

    :param str xyz: Name of the XYZ file to fix. It will be written to the
        same name but a backup will be made with '.bak.' added at the
        beginning of the name.
    :param bool dry_run: Default: False. If true, the output will not be
        written (but a file named 'test-temp.xyz' will be created/overwritten).
    :param bool return_list: Default: False. If True, the output will be
        written to 'test-temp.xyz' and the lines of the fixed file will be
        returned as a list of strings.
    :return: If return_list, a list of strings of the lines of the fixed XYZ
        file, otherwise None.
    :rtype: List[str] or None
    """
    if not dry_run and not return_list:
        xyz_c = '.bak.'+xyz
        copy_no_overwrite(xyz, xyz_c)
    else:
        xyz_c = xyz
        xyz = 'test-temp.xyz'

    def atom_name_sub(m):
        return m.group(1) + ' '*(len(m.group(0))-1)
    lines = []  # Could make this actually not write to the temp file
    with open(xyz_c, 'r') as f_in, open(xyz, 'w') as f_out:
        for i, line in enumerate(f_in):
            if i > 1:  # skip header lines
                line = re.sub(r'([CHONF])\S*',
                              atom_name_sub,
                              line)
            lines += line
            if not dry_run:
                f_out.write(line)
            else:
                print(line)
    if return_list:
        return lines


def use_gen_template(out_file, xyz, job_name='default job name',
                     checkpoint='checkpoint.cpt',
                     nproc=16, mem=125,
                     opt='opt', td=False,
                     func='wb97xd', basis='6-31g(d)',
                     charg_mult='0 1',
                     template='/projectnb/nonadmd/theavey'
                              '/qm-basics/templ-gen.txt'):
    """
    Use general template file to write Gaussian input file

    :type out_file: str or TextIOBase
    :param out_file: name of file or open file object to write output to

    :type xyz: str or list
    :param xyz: List of lines from an xyz file or string of path to an xyz
        file.

    :param str job_name: Default: 'default job name'. Name of the job to put
        into the Gaussian input.

    :param str checkpoint: Default: 'checkpoint.cpt'. File name for the
        checkpoint file.

    :type nproc: int or str
    :param nproc: Default: 16. Number of processors to tell Gaussian to use.

    :type mem: int or str
    :param mem: Default: 125. Number of gigabytes of memory to tell Gaussian
        to use.

    :param str opt: Default: 'opt'. Opt keywords to tell Gaussian.
        If True, this will be set to 'opt'.
        If this evaluates to False, it will be set to the blank string.

    :param str td: Default: False. TD keywords to tell Gaussian.
        If True, this will be set to TD.
        If this evaluates to False, it will be set to the blank string.

    :param str func: Default: 'wb97xd'. Functional for Gaussian to use.

    :param str basis: Default: '6-31g(d)'. Basis set for Gaussian to use.

    :param str charg_mult: Default: '0 1'. Charge and multiplicity line.

    :param str template: Default: '~nbth/qm-basics/templ-gen.txt'.
        The general template file to use. It should have keywords in curly
        braces with the same names as the keyword arguments to this function.
    :return: None
    """
    if opt:
        if opt is True:
            opt = 'opt'
    else:
        opt = ''
    if td:
        if td is True:
            td = 'TD'
    else:
        td = ''
    d_fill = dict(job_name=job_name,
                  checkpoint=checkpoint,
                  nproc=str(nproc), mem=str(mem),
                  opt=opt, td=td,
                  func=func, basis=basis,
                  charg_mult=charg_mult)
    if isinstance(xyz, string_types):
        xyz_lines = open(xyz, 'r').readlines()[2:]
    else:
        xyz_lines = xyz[2:]
    own_handle = False
    if isinstance(out_file, string_types):
        own_handle = True
        out_file = open(out_file, 'x')
    try:
        with open(template, 'r') as f_templ:
            line = ''  # To satisfy IDE in case of empty template
            for line in f_templ:
                line = line.format(**d_fill)
                out_file.write(line)
        if '\n' not in line:
            out_file.write('\n')
        for line in xyz_lines:
            out_file.write(line)
        out_file.write('\n\n\n')
    finally:
        if own_handle:
            out_file.close()
