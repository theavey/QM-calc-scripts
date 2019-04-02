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
from typing import List

import MDAnalysis
import parmed

__all__ = ['OniomUniverse']


class OniomUniverse(object):
    """
    Object to help easily create Gaussian ONIOM input sections


    There are a few ways to instantiate this object.
    First, it can be instatiated with an existing MDAnalysis Universe
    instance and an existing parmed Structure instance:

    >>> ou = OniomUniverse(univ=univ, structure=structure)

    Alternatively, the Universe and/or Structure can be instantiated here:

    >>> ou = OniomUniverse(univ_args=['geom.pdb', 'traj.xtc'], \
                           structure_file='topology.top')

    Any combination of these methods can also be used.

    Also, `high_select` and `low_select` need to be provided during or
    after instantiation before accessing the Gaussian input sections:

    >>> ou = OniomUniverse(univ=univ, structure=structure, \
                           high_select='resid 2-3')
    >>> ou.low_select = 'protein or byres around 5 resid 2-3'

    Then, the Gaussian input sections can be created:

    >>> mol_sec = ou.molecule_section
    >>> bond_sec = ou.bonds_section
    >>> param_sec = ou.params_section

    Note, if you do not need the parameter section (e.g., only using
    already included AMBER atom types), the structure or structure files
    need not be specified.
    """

    def __init__(self, univ: MDAnalysis.Universe = None,
                 structure: parmed.Structure = None,
                 high_select: str = None,
                 low_select: str = None,
                 overlap_okay: bool = False,
                 univ_args=None, univ_kwargs=None,
                 structure_file=None,
                 structure_args=None, structure_kwargs=None):
        """
        Instantiate OniomUniverse to create Gaussian ONIOM input sections

        :param univ: Universe with the geometry and bonding information for the
            system of interest. Note, the geometry must include bonding
            information or MDAnalysis will have to be told to guess them:
            https://www.mdanalysis.org/docs/documentation_pages/topology/guessers.html#MDAnalysis.topology.guessers.guess_bonds
        :param structure: Structure with atom types, bonds, angles,
            etc. Note, this is only currently written to work with AMBER,
            and it is unclear how it will work for other force fields.
        :param high_select: Selection string for the atoms to be included in the
            "high" calculation
        :param low_select: Selection string for the atoms to be included in the
            "low" calculation
        :param overlap_okay: If True, overlap between the high and low
            selections will be ignored and anything in the overlap will be
            included in the "high" region
        :param univ_args: arguments to be provided to instantiate the Universe
        :param univ_kwargs: keyword arguments to be provided to instantiate
            the Universe
        :param structure_file: filename (first argument) to be provided to
            instantiate the Structure
        :param structure_args: arguments to be provided to instantiate the
            Structure
        :param structure_kwargs: keyword arguments to be provided to instantiate
            the Structure
        """
        if univ is None:
            self.universe = MDAnalysis.Universe(*univ_args, **univ_kwargs)
        else:
            self.universe = univ
        if structure is None:
            if structure_file is None:
                self.structure = None
            else:
                self.structure = parmed.load_file(structure_file,
                                                  *structure_args,
                                                  **structure_kwargs)
        else:
            self.structure = structure
        self.high_select = high_select
        self.low_select = low_select
        self.overlap_okay = overlap_okay
        self.atom_to_line_num = dict()
        self.n_atoms_in_input = 0

    @property
    def molecule_section(self,) -> List[str]:
        """
        Create molecule specification lines for ONIOM calculation

        This defines a dict mapping `Atom`s to atom number (line number in
        input) as `self.atom_to_line_number`, and number of atoms included in
        the input as `self.n_atoms_in_input`.

        :return: The lines to be written into the input
        """
        if self.high_select is None or self.low_select is None:
            raise self.SelectionError('Both `high_select` and `low_select` '
                                      'must be specified')
        high_atoms = self.universe.select_atoms(self.high_select)
        low_atoms = self.universe.select_atoms(self.low_select)
        if high_atoms.intersection(low_atoms) and not self.overlap_okay:
            raise ValueError('The selections are not mutually exclusive')
        lines = []
        line_num = 0
        for atom in self.universe.atoms:
            if atom in high_atoms:
                level = 'H'
            elif atom in low_atoms:
                level = 'L'
            else:
                continue
            line_num += 1
            lines.append(self._make_atom_line(atom=atom, level=level,
                                              freeze_dict={'H': 0, 'L': -1}))
            self.atom_to_line_num[atom] = line_num
        sel_n_atoms = (high_atoms.n_atoms + low_atoms.n_atoms)
        if line_num != sel_n_atoms:
            raise ValueError('Number of lines and n_atoms in selections differ '
                             f'({line_num} and {sel_n_atoms})')
        self.n_atoms_in_input = sel_n_atoms
        return lines

    @property
    def bonds_section(self, ) -> List[str]:
        """
        Create bond specifications for a Gaussian job with `geom=connectivity`

        :return: The lines to be written in the input after the molecule
            specification
        """
        if self.n_atoms_in_input == 0:
            raise ValueError('No atoms have been put into the molecule '
                             'specification yet so the bonds cannnot yet be '
                             'defined. Either run `molecule_section` '
                             'first or check your selections.')
        atln = self.atom_to_line_num
        bond_dict = collections.defaultdict(list)
        for bond in self.universe.bonds:
            a1, a2 = bond.atoms
            try:
                bond_dict[atln[a1]].append(f'{atln[a2]} 1.0')
            except KeyError:
                continue
        lines = []
        for i in range(self.n_atoms_in_input):
            bonds = ' '.join(bond_dict[i])
            lines.append(f'{i} {bonds}\n')
        return lines

    @property
    def params_section(self, ) -> List[str]:
        """
        Create parameter specification for Gaussian job using Amber MM

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
        atom_types = self._get_atom_types()
        for at in atom_types:
            lines.append(self._make_atomtype_line(at))
        lines.append('! Stretch parameters\n')
        bond_types = self._get_bond_types()
        for bond in bond_types:
            lines.append(self._make_bondtype_line(bond))
        lines.append('! Bend parameters\n')
        angle_types = self._get_angle_types()
        for angle in angle_types:
            lines.append(self._make_angletype_line(angle))
        lines.append('! Dihedral parameters\n')
        dihedral_types = self._get_dihedral_types()
        for dihed in dihedral_types:
            lines.append(self._make_dihedraltype_line(dihed))
        lines.append('! Improper dihedral parameters\n')
        improper_types = self._get_improper_types()
        for dihed in improper_types:
            lines.append(self._make_impropertype_line(dihed))
        return lines

    _re_element = re.compile(r'[A-Z][a-z]?')

    def _get_elem(self, atom: MDAnalysis.core.groups.Atom) -> str:
        """
        Get element name from Atom object

        This counts on any multi-letter element being named as Ca (capital
        followed by lower case). Also, single-letter element names must be
        capitalized and not followed by a lower-case letter.

        An alternative method could be a mapping from masses to elements, or
        use the parmed Structure which knows element information.
        """
        elem_match = self._re_element.match(atom.name)
        if elem_match:
            return elem_match.group(0)
        else:
            raise ValueError(f'Could not find element for atom {atom}')

    def _make_atom_line(self, atom: MDAnalysis.core.groups.Atom,
                        level: str,
                        freeze_dict: dict = None) -> str:
        if freeze_dict is None:
            freeze_dict = {'H': 0, 'L': -1}
        elem = self._get_elem(atom)
        line = (f'{elem}-'
                f'{atom.type}-'
                f'{atom.charge:3f} {freeze_dict[level]} '
                f'{atom.position[0]:4f} '
                f'{atom.position[1]:4f} '
                f'{atom.position[2]:4f} {level} \n')
        return line

    def _get_atom_types(self,) -> set:
        atom_types = set()
        for atom in self.structure.atoms:
            atom_types.add(atom.atom_type)
        return atom_types

    @staticmethod
    def _make_atomtype_line(atom_type) -> str:
        sigma = atom_type.urmin
        epsilon = atom_type.uepsilon
        sigma = sigma.value_in_unit(parmed.unit.angstrom)
        epsilon = epsilon.value_in_unit(parmed.unit.kilocalorie_per_mole)
        return f'VDW {atom_type.name: <2} {sigma:4f}  {epsilon:4f}\n'

    @staticmethod
    def _get_types(instances, types) -> List:
        instance_types = list()
        for _type in types:
            for inst in instances:
                if inst.type == _type:
                    instance_types.append(inst)
                    break
        return instance_types

    def _get_bond_types(self,) -> List:
        return self._get_types(self.structure.bonds, self.structure.bond_types)

    @staticmethod
    def _make_bondtype_line(bond: parmed.topologyobjects.Bond) -> str:
        a1, a2 = bond.atom1.type, bond.atom2.type
        k = bond.type.uk.value_in_unit(parmed.unit.kilocalorie_per_mole /
                                       parmed.unit.angstrom ** 2)
        req = bond.type.ureq.value_in_unit(parmed.unit.angstrom)
        return f'HrmStr1 {a1:2} {a2:2}  {k: <5.1f}  {req: <5.3f}\n'

    def _get_angle_types(self,) -> List:
        return self._get_types(self.structure.angles,
                               self.structure.angle_types)

    @staticmethod
    def _make_angletype_line(angle: parmed.topologyobjects.Angle) -> str:
        a1, a2, a3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
        k = angle.type.uk.value_in_unit(parmed.unit.kilocalorie_per_mole /
                                        parmed.unit.radian ** 2)
        thetaeq = angle.type.utheteq.value_in_unit(parmed.unit.degree)
        return f'HrmBnd1 {a1:2} {a2:2} {a3:2}  {k: >5.1f}  {thetaeq:6.2f}'

    def _get_improper_types(self,) -> List:
        # Somewhere along antechamber -> acpype, the impropers are stored
        # as dihedrals (of GROMACS function 1)
        return self._get_types(self.structure.dihedrals,
                               self.structure.dihedral_types)

    @staticmethod
    def _make_impropertype_line(dihed: parmed.topologyobjects.Dihedral
                                ) -> str:
        a1, a2, a3, a4 = (dihed.atom1.type, dihed.atom2.type,
                          dihed.atom3.type, dihed.atom4.type)
        phi_k = dihed.type.uphi_k.value_in_unit(
            parmed.unit.kilocalorie_per_mole)
        phase = dihed.type.uphase.value_in_unit(parmed.unit.degree)
        per = dihed.type.per
        return (f'ImpTrs {a1:2} {a2:2} {a3:2} {a4:2}  '
                f'{phi_k: >5.1f}  {phase:5.1f} {per:3.1f}')

    def _get_dihedral_types(self,) -> List:
        # Somewhere along antechamber -> acpype, the impropers are stored
        # as dihedrals (of GROMACS function 1)
        # and the dihedrals get stored as Ryckaert-Bellemans
        # dihedrals (function 3)
        return self._get_types(self.structure.rb_torsions,
                               self.structure.rb_torsion_types)

    @staticmethod
    def _make_dihedraltype_line(dihed: parmed.topologyobjects.Dihedral
                                ) -> str:
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

    class SelectionError(ValueError):
        pass
