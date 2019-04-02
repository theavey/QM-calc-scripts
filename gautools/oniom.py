"""
A set of tools for setting up Gaussian ONIOM calculations
"""


########################################################################
#                                                                      #
# This script was written by Thomas Heavey in 2019.                    #
#        theavey@bu.edu     thomasjheavey@gmail.com                    #
#                                                                      #
# Copyright 2019 Thomas J. Heavey IV                                   #
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

import collections
import re
from typing import List, Tuple

import MDAnalysis
import parmed

__all__ = ['make_molecule_section', 'make_bonds_section',
           'make_params_section']

_re_element = re.compile(r'[A-Z][a-z]?')


def _get_elem(atom: MDAnalysis.core.groups.Atom) -> str:
    """
    Get element name from Atom object

    This counts on any multi-letter element being named as Ca (capital
    followed by lower case). Also, single-letter element names must be
    capitalized and not followed by a lower-case letter.

    An alternative method could be a mapping from masses to elements, or
    use the parmed Structure which knows element information.
    """
    elem_match = _re_element.match(atom.name)
    if elem_match:
        return elem_match.group(0)
    else:
        raise ValueError(f'Could not find element for atom {atom}')


def _make_atom_line(atom: MDAnalysis.core.groups.Atom,
                    level: str,
                    freeze_dict: dict = None) -> str:
    if freeze_dict is None:
        freeze_dict = {'H': 0, 'L': -1}
    elem = _get_elem(atom)
    line = (f'{elem}-'
            f'{atom.type}-'
            f'{atom.charge:3f} {freeze_dict[level]} '
            f'{atom.position[0]:4f} '
            f'{atom.position[1]:4f} '
            f'{atom.position[2]:4f} {level} \n')
    return line


def make_molecule_section(univ: MDAnalysis.Universe,
                          high_select: str,
                          low_select: str,
                          ignore_overlap: bool = False
                          ) -> Tuple[List[str], dict, int]:
    """
    Create molecule specification lines for ONIOM calculation

    :param univ: Universe including all atoms to be included (at the correct
        timestep/geometry)
    :param high_select: Selection string for the atoms to be included in the
        "high" calculation
    :param low_select: Selection string for the atoms to be included in the
        "low" calculation
    :param ignore_overlap: If True, overlap between the high and low selections
        will be ignored and anything in the overlap will be included in the
        "high" region
    :return: The lines to be written into the input, a dict mapping `Atom`s
        to atom number (line number in input), and number of atoms included
    """
    high_atoms = univ.select_atoms(high_select)
    low_atoms = univ.select_atoms(low_select)
    if high_atoms.intersection(low_atoms) and not ignore_overlap:
        raise ValueError('The selections are not mutually exclusive')
    lines = []
    line_num = 0
    atom_to_line_num = dict()
    for atom in univ.atoms:
        if atom in high_atoms:
            level = 'H'
        elif atom in low_atoms:
            level = 'L'
        else:
            continue
        line_num += 1
        lines.append(_make_atom_line(atom=atom, level=level,
                                     freeze_dict={'H': 0, 'L': -1}))
        atom_to_line_num[atom] = line_num
    sel_n_atoms = (high_atoms.n_atoms + low_atoms.n_atoms)
    if line_num != sel_n_atoms:
        raise ValueError('Number of lines and n_atoms in selections differ '
                         f'({line_num} and {sel_n_atoms})')
    return lines, atom_to_line_num, sel_n_atoms


def make_bonds_section(univ: MDAnalysis.Universe,
                       atom_to_line_num: dict,
                       n_atoms: int,  # Could likely get this from the dict...
                       ) -> List[str]:
    """
    Create bond specifications for a Gaussian job with `geom=connectivity`

    :param univ: Universe including all atoms to be included
    :param atom_to_line_num: Mapping from atoms to atom number as specified in
        the Gaussian input
    :param n_atoms: Number of atoms in the Gaussian input (not necessarily the
        total number of atoms in the Universe if some of the atoms are not
        selected)
    :return: The lines to be written in the input after the molecule
        specification
    """
    atln = atom_to_line_num
    bond_dict = collections.defaultdict(list)
    for bond in univ.bonds:
        a1, a2 = bond.atoms
        try:
            bond_dict[atln[a1]].append(f'{atln[a2]} 1.0')
        except KeyError:
            continue
    lines = []
    for i in range(n_atoms):
        bonds = ' '.join(bond_dict[i])
        lines.append(f'{i} {bonds}\n')
    return lines


def _get_atom_types(structure) -> set:
    atom_types = set()
    for atom in structure.atoms:
        atom_types.add(atom.atom_type)
    return atom_types


def _make_atomtype_line(atom_type) -> str:
    sigma = atom_type.urmin
    epsilon = atom_type.uepsilon
    sigma = sigma.value_in_unit(parmed.unit.angstrom)
    epsilon = epsilon.value_in_unit(parmed.unit.kilocalorie_per_mole)
    return f'VDW {atom_type.name: <2} {sigma:4f}  {epsilon:4f}\n'


def _get_types(instances, types) -> List:
    instance_types = list()
    for _type in types:
        for inst in instances:
            if inst.type == _type:
                instance_types.append(inst)
                break
    return instance_types


def _get_bond_types(structure: parmed.structure.Structure) -> List:
    return _get_types(structure.bonds, structure.bond_types)


def _make_bondtype_line(bond: parmed.topologyobjects.Bond) -> str:
    a1, a2 = bond.atom1.type, bond.atom2.type
    k = bond.type.uk.value_in_unit(parmed.unit.kilocalorie_per_mole /
                                   parmed.unit.angstrom ** 2)
    req = bond.type.ureq.value_in_unit(parmed.unit.angstrom)
    return f'HrmStr1 {a1:2} {a2:2}  {k: <5.1f}  {req: <5.3f}\n'


def _get_angle_types(structure: parmed.structure.Structure) -> List:
    return _get_types(structure.angles, structure.angle_types)


def _make_angletype_line(angle: parmed.topologyobjects.Angle) -> str:
    a1, a2, a3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
    k = angle.type.uk.value_in_unit(parmed.unit.kilocalorie_per_mole /
                                    parmed.unit.radian ** 2)
    thetaeq = angle.type.utheteq.value_in_unit(parmed.unit.degree)
    return f'HrmBnd1 {a1:2} {a2:2} {a3:2}  {k: >5.1f}  {thetaeq:6.2f}'


def _get_improper_types(structure: parmed.structure.Structure) -> List:
    # Somewhere along antechamber -> acpype, the impropers are stored
    # as dihedrals (of GROMACS function 1)
    return _get_types(structure.dihedrals, structure.dihedral_types)


def _make_impropertype_line(dihed: parmed.topologyobjects.Dihedral) -> str:
    a1, a2, a3, a4 = (dihed.atom1.type, dihed.atom2.type,
                      dihed.atom3.type, dihed.atom4.type)
    phi_k = dihed.type.uphi_k.value_in_unit(parmed.unit.kilocalorie_per_mole)
    phase = dihed.type.uphase.value_in_unit(parmed.unit.degree)
    per = dihed.type.per
    return (f'ImpTrs {a1:2} {a2:2} {a3:2} {a4:2}  '
            f'{phi_k: >5.1f}  {phase:5.1f} {per:3.1f}')


def _get_dihedral_types(structure: parmed.structure.Structure) -> List:
    # Somewhere along antechamber -> acpype, the impropers are stored
    # as dihedrals (of GROMACS function 1)
    # and the dihedrals get stored as Ryckaert-Bellemans dihedrals (function 3)
    return _get_types(structure.rb_torsions, structure.rb_torsion_types)


def _make_dihedraltype_line(dihed: parmed.topologyobjects.Dihedral) -> str:
    a1, a2, a3, a4 = (dihed.atom1.type, dihed.atom2.type,
                      dihed.atom3.type, dihed.atom4.type)
    dtl = parmed.DihedralTypeList.from_rbtorsion(dihed.type)
    phases = [0] * 4
    phis = [0.] * 4
    for dihed_type in dtl:
        phi_k = dihed_type.uphi_k.value_in_unit(
            parmed.unit.kilocalorie_per_mole)
        phase = dihed_type.uphase.value_in_unit(parmed.unit.degree)
        per = dihed_type.per
        phases[per], phis[per] = phase, phi_k
    output = (f'AmbTrs {a1:2} {a2:2} {a3:2} {a4:2}  ' +
              ' '.join([f'{i: >3d}' for i in phases]) + ' ' +
              ' '.join([f'{i: >6.3f}' for i in phis]) +
              ' 0.0\n')
    return output


def make_params_section(structure: parmed.structure.Structure, ) -> List[str]:
    """
    Create parameter specification for Gaussian job using Amber MM

    :param structure: Structure parameterized to include all needed bonded and
        non-bonded (VdW specifications in atom types) parameters
    :return: The lines to be included for Gaussian jobs using Amber MM and
        HardFirst, SoftFirst, or SoftOnly
    """
    # This doesn't seem perfect: the selection functions don't
    # really work because (for example) a single bond_type might
    # be used for different atom types, which is what is
    # currently assumed. Need to find a way to either find all
    # atom types for which it should be used (possibly using
    # wildcards), or just iterate over all bonds/angles/dihedrals
    # instead of iterating over *_types.
    lines = list()
    lines.append('! Van der Waals parameters\n')
    atom_types = _get_atom_types(structure)
    for at in atom_types:
        lines.append(_make_atomtype_line(at))
    lines.append('! Stretch parameters\n')
    bond_types = _get_bond_types(structure)
    for bond in bond_types:
        lines.append(_make_bondtype_line(bond))
    lines.append('! Bend parameters\n')
    angle_types = _get_angle_types(structure)
    for angle in angle_types:
        lines.append(_make_angletype_line(angle))
    lines.append('! Dihedral parameters\n')
    dihedral_types = _get_dihedral_types(structure)
    for dihed in dihedral_types:
        lines.append(_make_dihedraltype_line(dihed))
    lines.append('! Improper dihedral parameters\n')
    improper_types = _get_improper_types(structure)
    for dihed in improper_types:
        lines.append(_make_impropertype_line(dihed))
    return lines
