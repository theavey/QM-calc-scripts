"""
A set of tools for working with computational chemistry files and such
"""


########################################################################
#                                                                      #
# This script was written by Thomas Heavey in 2019.                    #
#        theavey@bu.edu     thomasjheavey@gmail.com                    #
#                                                                      #
# Copyright 2015-2019 Thomas J. Heavey IV                              #
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
import pathlib
from six import string_types

from .oniom import OniomUniverse
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
                     checkpoint='checkpoint.chk',
                     rwf='readwrite.rwf',
                     nproc=16, mem=125,
                     opt='opt', td=False,
                     func='wb97xd', basis='6-31g(d)',
                     charg_mult='0 1',
                     footer='\n\n',
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
    :param str rwf: Default: 'readwrite.rwf'. File name for the read-write
        file.
    :type nproc: int or str
    :param nproc: Default: 16. Number of processors to tell Gaussian to use.
        Note, this now uses the newer '%cpu' syntax, and I'm not sure how
        that will work using fewer than all CPUs on the node because it says
        to use 0 to nproc-1.
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
    :param str footer: Default: '\n\n'. Footer of input file. Useful for RESP
        charge calculation jobs and such.
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
                  checkpoint=checkpoint, rwf=rwf,
                  nproc=str(int(nproc)-1), mem=str(mem),
                  opt=opt, td=td,
                  func=func, basis=basis,
                  charg_mult=charg_mult)
    xyz_lines = _get_xyz_lines(xyz)
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
        out_file.write(footer)
    finally:
        if own_handle:
            out_file.close()


def make_gaussian_input(out_file, xyz, job_name='default job name',
                        checkpoint='checkpoint.chk',
                        rwf='readwrite.rwf',
                        nproc=16, mem=125,
                        route=None,
                        opt='opt', td=False,
                        func='wb97xd', basis='6-31g(d)',
                        route_other=None,
                        charg_mult='0 1',
                        footer=None,
                        oniom: dict = None):
    """
    Write Gaussian input file

    :type out_file: str or TextIOBase
    :param out_file: name of file or open file object to write output to

    :type xyz: str or list
    :param xyz: List of lines from an xyz file or string of path to an xyz
        file.

    :param str job_name: Default: 'default job name'. Name of the job to put
        into the Gaussian input.

    :param str checkpoint: Default: 'checkpoint.cpt'. File name for the
        checkpoint file.

    :param str rwf: Default: 'readwrite.rwf'. File name for the read-write
        file.

    :type nproc: int or str
    :param nproc: Default: 16. Number of processors to tell Gaussian to use.
        Note, this now uses the newer '%cpu' syntax, and I'm not sure how
        that will work using fewer than all CPUs on the node because it says
        to use 0 to nproc-1.

    :type mem: int or str
    :param mem: Default: 125. Number of gigabytes of memory to tell Gaussian
        to use.

    :param str route: If not None, this will be the entire route section and
        the following commands will be ignored: `opt`, `td`, `func`, `basis`.

    :param str opt: Default: 'opt'. Opt keywords to tell Gaussian.
        If True, this will be set to 'opt'.
        If this evaluates to False, it will be set to the blank string.
        If something else, it will be set to the given string.

    :param str td: Default: False. TD keywords to tell Gaussian.
        If True, this will be set to TD.
        If this evaluates to False, it will be set to the blank string.
        If something else, it will be set to the given string.

    :param str func: Default: 'wb97xd'. Functional for Gaussian to use.
        If True or evaluates as false, it will be set to a blank string,
        which will likely be an error.

    :param str basis: Default: '6-31g(d)'. Basis set for Gaussian to use.
        If True or evaluates as false, it will be set to a blank string,
        which will likely be an error.

    :param str route_other: Other commands to use in the route section
        (e.g., 'SCRF=(solvent=dichloromethane) Int=Ultrafile freq')

    :param str charg_mult: Default: '0 1'. Charge and multiplicity line.

    :param str footer: Default: None. Footer of input file. Useful for RESP
        charge calculation jobs and such.

    :param dict oniom: dict to pass to :py:class:`gautools.oniom.OniomUniverse`
        constructor. The create object will then be used to make the molecule
        specification, and add the connectivity and MM parameters to the footer.
    :return: The Path to the written file
    :rtype: pathlib.Path
    """
    link0 = _make_link0(checkpoint, rwf, str(int(nproc)-1), mem)
    route_sec = _make_route(route, opt, td, func, basis, route_other)
    if oniom is not None:
        ou = OniomUniverse(**oniom)
        xyz_lines = ou.molecule_section
        bon_sec = ''.join(ou.bonds_section)
        par_sec = ''.join(ou.params_section)
        if footer is None:
            footer_list = [bon_sec, par_sec]
        else:
            footer_list = [bon_sec, footer, par_sec]
            # This should be the right order in most cases:
            # http://gaussian.com/input/
        footer = '\n'.join(footer_list)
    else:
        xyz_lines = _get_xyz_lines(xyz)
    own_handle = False
    if isinstance(out_file, string_types):
        own_handle = True
        out_file_path = pathlib.Path(out_file)
        out_file = open(out_file, 'x')
    else:
        out_file_path = pathlib.Path(out_file.name)
    try:
        out_file.write(link0)
        out_file.write(route_sec)
        out_file.write('\n')  # blank line between sections
        out_file.write(_make_newline_terminated(job_name))
        out_file.write('\n')  # blank line between sections
        out_file.write(_make_newline_terminated(charg_mult))
        line = ''  # in case empty xyz_lines
        for line in xyz_lines:
            out_file.write(line)
        if line[-1] == '\n':  # hopefully only 1 newline at most
            out_file.write('\n')  # blank line between sections
        else:
            out_file.write('\n\n')  # blank line between sections
        if footer:
            out_file.write(_make_newline_terminated(footer))
        out_file.write('\n')  # blank line before end of file
    finally:
        if own_handle:
            out_file.close()
    return out_file_path.resolve()


_link0_template_dict = {'nproc': '%cpu=0-{nproc}',
                        'mem': '%mem={mem}GB',
                        'rwf': '%rwf={rwf}\n%NoSave',
                        'checkpoint': '%chk={checkpoint}'}


def _make_link0(checkpoint, rwf, nproc, mem):
    # http://gaussian.com/link0/
    output = []
    kwarg_dict = dict()
    # want at least rwf and checkpoint to be ordered (for %NoSave),
    # so this might not be perfect in Python < 3.6
    kwarg_dict['mem'] = mem
    kwarg_dict['nproc'] = nproc
    kwarg_dict['rwf'] = rwf
    kwarg_dict['checkpoint'] = checkpoint
    for key in kwarg_dict:
        if kwarg_dict[key]:
            output.append(_link0_template_dict[key].format(**kwarg_dict))
    if output:
        return _make_newline_terminated('\n'.join(output))
    else:
        return str()


_route_template = '# {opt} {td} {func}/{basis} {route_other}'


def _make_route(route, opt, td, func, basis, route_other):
    if route:
        if not route[0] == '#':
            route = '# ' + route
        return _make_newline_terminated(route)
    kwarg_dict = dict(opt=opt, td=td, func=func,
                      basis=basis, route_other=route_other)
    defaults_dict = dict(
        opt='opt', td='TD',
        func='', basis='',  # not sure what good defaults are here
        route_other='')
    for key in kwarg_dict:
        kwarg_dict[key] = _process_keyword(kwarg_dict[key],
                                           defaults_dict[key])
    return _make_newline_terminated(_route_template.format(**kwarg_dict))


def _process_keyword(key, key_default):
    if key:
        key = key_default if key is True else key
    else:
        key = ''
    return key


def _get_xyz_lines(xyz):
    if isinstance(xyz, string_types):
        xyz_lines = open(xyz, 'r').readlines()[2:]
    else:
        xyz_lines = xyz[2:]
    return xyz_lines


def _make_newline_terminated(line):
    if line[-1] == '\n':
        return line
    else:
        return line + '\n'
