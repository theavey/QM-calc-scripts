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
