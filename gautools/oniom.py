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
from typing import List, Union, Tuple, Dict

import MDAnalysis
import numpy as np
import parmed

__all__ = ['OniomUniverse']


class NoStructureException(Exception):
    pass


class OniomStructure(object):

    def __init__(self,
                 structure: parmed.structure.Structure = None,
                 structure_file: str = None,
                 structure_args: Union[Tuple, List] = None,
                 structure_kwargs: dict = None,
                 only_unique_types: bool = False,
                 only_used_terms: bool = True):
        """
        Initialize OniomStructure to create Gaussian Amber MM input section

        :param structure_file: filename (first argument) to be provided to
            instantiate the Structure
        :param structure_args: arguments to be provided to instantiate the
            Structure
        :param structure_kwargs: keyword arguments to be provided to instantiate
            the Structure
        :param only_unique_types: If False (default), all bonds, angles,
            dihedrals, and impropers will be included.
            If True, only the unique elements for each of those will be
            included, which may not define the terms for all possible
            interactions because one type may be used for several atom types.
            For example, there might be an angle_type that should be used for
            "*-C3-C3", but it might only get defined for "H1-C3-C3".
        :param only_used_terms: If True (default), the params returned will
            only include those with all atoms contained in the atoms actually
            used. This can make the section returned shorter if not all atoms
            have been selected, especially if `only_unique_types` is False.
            This will also require `atoms_used_indices` to be defined, which
            it will be if `molecule_section` is accessed first (in an associated
            OniomUniverse).
            If False, all parameters will be given.
        """
        if structure is None:
            if (structure_file is None and
                    structure_args is None and
                    structure_kwargs is None):
                raise NoStructureException
            else:
                self.structure = parmed.load_file(structure_file,
                                                  *structure_args,
                                                  **structure_kwargs)
        else:
            self.structure = structure
        self.only_unique_types = only_unique_types
        self._unique_types = {
            'bonds': self._get_bond_types_uniq(),
            'angles': self._get_angle_types_uniq(),
            'dihedrals': self._get_dihedral_types_uniq(),
            'impropers': self._get_improper_types_uniq()
        }
        self._non_unique_types = {
            'bonds': self._get_bond_types_nu(),
            'angles': self._get_angle_types_nu(),
            'dihedrals': self._get_dihedral_types_nu(),
            'impropers': self._get_improper_types_nu()}
        self._types_dict = {True: self._unique_types,
                            False: self._non_unique_types}
        self.atoms_used_indices = None
        self.only_used_terms = only_used_terms

    @property
    def params_section(self) -> List[str]:
        """
        Parameter specification for Gaussian job using Amber MM

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
        types = self._types_dict[self.only_unique_types]
        if self.only_used_terms and self.atoms_used_indices is not None:
            types = self._remove_unused_terms(types)
        lines.append('! Van der Waals parameters\n!\n')
        atom_types = self._get_atom_types()
        param_lines = set()
        for at in atom_types:
            param_lines.add(self._make_atomtype_line(at))
        lines += list(param_lines)
        lines.append('! Stretch parameters\n!\n')
        bond_types = types['bonds']
        param_lines = set()
        for bond in bond_types:
            param_lines.add(self._make_bondtype_line(bond))
        lines += list(param_lines)
        lines.append('! Bend parameters\n!\n')
        angle_types = types['angles']
        param_lines = set()
        for angle in angle_types:
            param_lines.add(self._make_angletype_line(angle))
        lines += list(param_lines)
        lines.append('! Dihedral parameters\n!\n')
        dihedral_types = types['dihedrals']
        param_lines = set()
        for dihed in dihedral_types:
            param_lines.add(self._make_dihedraltype_line(dihed))
        lines += list(param_lines)
        lines.append('! Improper dihedral parameters\n!\n')
        improper_types = types['impropers']
        param_lines = set()
        for dihed in improper_types:
            param_lines.add(self._make_impropertype_line(dihed))
        lines += list(param_lines)
        return lines

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

    def _get_bond_types_uniq(self, ) -> List:
        return self._get_types(self.structure.bonds, self.structure.bond_types)

    def _get_bond_types_nu(self) -> List:
        return self.structure.bonds

    @staticmethod
    def _make_bondtype_line(bond: parmed.topologyobjects.Bond) -> str:
        a1, a2 = bond.atom1.type, bond.atom2.type
        k = bond.type.uk.value_in_unit(parmed.unit.kilocalorie_per_mole /
                                       parmed.unit.angstrom ** 2)
        req = bond.type.ureq.value_in_unit(parmed.unit.angstrom)
        return f'HrmStr1 {a1:2} {a2:2}  {k: <5.1f}  {req: <5.3f}\n'

    def _get_angle_types_uniq(self, ) -> List:
        return self._get_types(self.structure.angles,
                               self.structure.angle_types)

    def _get_angle_types_nu(self) -> List:
        return self.structure.angles

    @staticmethod
    def _make_angletype_line(angle: parmed.topologyobjects.Angle) -> str:
        a1, a2, a3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
        k = angle.type.uk.value_in_unit(parmed.unit.kilocalorie_per_mole /
                                        parmed.unit.radian ** 2)
        thetaeq = angle.type.utheteq.value_in_unit(parmed.unit.degree)
        return f'HrmBnd1 {a1:2} {a2:2} {a3:2}  {k: >5.1f}  {thetaeq:6.2f}\n'

    def _get_improper_types_uniq(self, ) -> List:
        # Somewhere along antechamber -> acpype, the impropers are stored
        # as dihedrals (of GROMACS function 1)
        return self._get_types(self.structure.dihedrals,
                               self.structure.dihedral_types)

    def _get_improper_types_nu(self) -> List:
        # Somewhere along antechamber -> acpype, the impropers are stored
        # as dihedrals (of GROMACS function 1)
        return self.structure.dihedrals

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
                f'{phi_k: >5.1f}  {phase:5.1f} {per:3.1f}\n')

    def _get_dihedral_types_uniq(self, ) -> List:
        # Somewhere along antechamber -> acpype, the impropers are stored
        # as dihedrals (of GROMACS function 1)
        # and the dihedrals get stored as Ryckaert-Bellemans
        # dihedrals (function 3)
        return self._get_types(self.structure.rb_torsions,
                               self.structure.rb_torsion_types)

    def _get_dihedral_types_nu(self) -> List:
        # Somewhere along antechamber -> acpype, the impropers are stored
        # as dihedrals (of GROMACS function 1)
        # and the dihedrals get stored as Ryckaert-Bellemans
        # dihedrals (function 3)
        return self.structure.rb_torsions

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
                  ' -1.0\n')
        return output

    def _remove_unused_terms(self, types) -> Dict[str, list]:
        atoms_used = np.array(self.structure.atoms)[self.atoms_used_indices]
        n_atoms_by_type = {'bonds': 2, 'angles': 3, 'dihedrals': 4,
                           'impropers': 4}
        for key in types:
            n_atoms = n_atoms_by_type[key]
            return_params = []
            input_params = types[key]
            for param in input_params:
                for i in range(n_atoms):
                    i += 1
                    if not getattr(param, f'atom{i}') in atoms_used:
                        break
                else:
                    return_params.append(param)
            types[key] = return_params
        return types


class OniomUniverse(object):
    """
    Object to help easily create Gaussian ONIOM input sections


    There are a few ways to instantiate this object.
    First, it can be instantiated with an existing MDAnalysis Universe
    instance and an existing parmed Structure instance:

    >>> univ = MDAnalysis.Universe('geom.pdb', 'traj.xtc')
    >>> structure = parmed.load_file('topology.top')
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

    The interfaces between high and low are not treated specially, so link
    atoms will need to be manually treated. That can be done after writing
    to an input file, or using something like::

        interface_atom = univ.select_atoms('bynum 88')
        interface_atom_index = ou.atom_to_line_num[interface_atom] - 1
        # (because the dict gives a (Gaussian) 1-based index)
        interface_atom_line = mol_sec[interface_atom_index][:-2]+' H-H1-0.1\\n'
        # (remove newline, add link atom definition and newline)
        mol_sec[interface_atom_index] = interface_atom_line

    """

    def __init__(self, univ: MDAnalysis.Universe = None,
                 structure: parmed.Structure = None,
                 high_select: str = None,
                 low_select: str = None,
                 overlap_okay: bool = False,
                 univ_args=None, univ_kwargs=None,
                 structure_file=None,
                 structure_args=None, structure_kwargs=None,
                 freeze_dict: dict = None):
        """
        Initialize OniomUniverse to create Gaussian ONIOM input sections

        :param univ: Universe with the geometry and bonding information for the
            system of interest. Note, the geometry must include bonding
            information or MDAnalysis will have to be told to guess them:
            https://www.mdanalysis.org/docs/documentation_pages/topology/guessers.html#MDAnalysis.topology.guessers.guess_bonds
        :param structure: Structure with atom types, bonds, angles,
            etc. Note, this is only currently written to work with AMBER (or
            really GAFF as made by Antechamber/AcPype), and it is unclear
            how it will work for other force fields or implementations.
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
        :param freeze_dict: mapping from levels ('H' and 'L') to freeze
            commands (0 for unfrozen, -1 for frozen).
            Default is `{'H': 0, 'L': -1}`
        """
        if univ is None:
            self.universe = MDAnalysis.Universe(*univ_args, **univ_kwargs)
        else:
            self.universe = univ
        try:
            self.oniom_structure = OniomStructure(
                structure=structure,
                structure_file=structure_file,
                structure_args=structure_args,
                structure_kwargs=structure_kwargs)
        except NoStructureException:
            self.oniom_structure = None
        self.high_select = high_select
        self.low_select = low_select
        self.overlap_okay = overlap_okay
        self.atom_to_line_num = dict()
        self.n_atoms_in_input = 0
        self.freeze_dict = ({'H': 0, 'L': -1} if freeze_dict is None
                            else freeze_dict)

    @property
    def molecule_section(self,) -> List[str]:
        """
        Molecule specification lines for ONIOM calculation

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
        n_atoms_in_both = high_atoms.intersection(low_atoms).n_atoms
        if n_atoms_in_both and not self.overlap_okay:
            raise ValueError('The selections are not mutually exclusive. '
                             'Make mutually exclusive or set overlap_okay=True')
        atoms_used_indices = []
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
            lines.append(self._make_atom_line(atom=atom, level=level,))
            atoms_used_indices.append(atom.index)
            self.atom_to_line_num[atom] = line_num
        sel_n_atoms = (high_atoms.n_atoms + low_atoms.n_atoms - n_atoms_in_both)
        if line_num != sel_n_atoms:
            raise ValueError('Number of lines and n_atoms in selections differ '
                             f'({line_num} and {sel_n_atoms})')
        self.n_atoms_in_input = sel_n_atoms
        if self.oniom_structure is not None:
            self.oniom_structure.atoms_used_indices = atoms_used_indices
        return lines

    @property
    def bonds_section(self, ) -> List[str]:
        """
        Bond specifications for a Gaussian job with `geom=connectivity`

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
            i += 1  # use 1-based indexing
            bonds = ' '.join(bond_dict[i])
            lines.append(f'{i} {bonds}\n')
        return lines

    @property
    def params_section(self):
        if self.oniom_structure is None:
            raise NoStructureException('No Structure given for this '
                                       'OniomUniverse')
        else:
            return self.oniom_structure.params_section

    params_section.__doc__ = OniomStructure.params_section.__doc__

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
                        level: str,) -> str:
        elem = self._get_elem(atom)
        line = (f'{elem}-'
                f'{atom.type}-'
                f'{atom.charge:3f} {self.freeze_dict[level]} '
                f'{atom.position[0]:4f} '
                f'{atom.position[1]:4f} '
                f'{atom.position[2]:4f} {level}\n')
        return line

    class SelectionError(ValueError):
        pass
