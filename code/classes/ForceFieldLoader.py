from typing import Tuple

"""
Force Field Loader Module for Molecular Dynamics Simulations

Description:
    This module provides classes and functions to load, manage, and process
    force field parameters from specified file paths within 'forcefields/' directory.
    It uses the FF_atom, FF_pair_coef, FF_bond, FF_angle, FF_torsion, FF_improper, and FF_equivalence classes 
    - loads cvff data from clayff and cvff forcefields in .frc files
    - also works with only one file
    - loads atom types, equivalences, pair_coef, bond_coef, angle_coef, torsion_coef, improper_coef parameters
    - pair_coef parameters are converted to give sigma and epsilon values
    - 2k values are divided by 2 to get k values
    - cvff equivalences are used to ensure all atom types are represented in the force field parameters
    - checks for duplicate atom types and appends '[dup_#]' to the type if a duplicate is found, 
        with the second file iteration. FF_atom.type are thus unique identifiers
    - for torsion_coefs, '*' is used to represent a wildcard atom type
"""

from code.classes.SystemParts import FF_atom, FF_pair_coef, FF_bond_coef, FF_angle_coef, FF_torsion_coef, FF_improper_coef, FF_equivalence
from code.classes.ChosenSettings import ChosenSettings
from typing import List, Set, Dict, Optional

class ForceFieldLoader:
    def __init__(self, settings: ChosenSettings) -> None:
        """
        Initialize the ForceFieldLoader with the given settings.

        Args:
            settings (ChosenSettings): The settings object containing various parameters.
        """
        self.s = settings

        self.ff_atoms: List[FF_atom] = []
        self.ff_equivalences: List[FF_equivalence] = []
        self.pair_coefs: List[FF_pair_coef] = []
        self.bond_coefs: List[FF_bond_coef] = []
        self.angle_coefs: List[FF_angle_coef] = []
        self.torsion_coefs: List[FF_torsion_coef] = []
        self.improper_coefs: List[FF_improper_coef] = []
        self.duplicates: Set[str] = set()

        self.load_all_forcefield_params()

    def load_all_forcefield_params(self) -> None:
        """
        Load all force field parameters from the specified files.
        """
        i = 0
        files = [f"forcefields/{name}.frc" for name in self.s.ff_files]
        for file_path in files:
            self.read_ff_atom_types(file_path, i)
            self.read_ff_equivalences(file_path, i)

            self.read_ff_pair_coef(file_path, i)
            self.add_pair_coef_equivalences()
            self.add_pair_coefs_for_ff_types()

            self.read_ff_bond(file_path, i)
            self.add_bond_equivalences()

            self.read_ff_angle(file_path, i)
            self.add_angle_equivalences()

            self.read_ff_torsion(file_path, i)
            self.add_torsion_equivalences()

            self.read_ff_improper(file_path, i)
            self.add_improper_equivalences()

            i += 1

    def get_info_by_type(self, identifier: str) -> Dict[str, Optional[object]]:
        """
        Get information related to a specific atom type.

        Args:
            identifier (str): The atom type identifier.

        Returns:
            Dict[str, Optional[object]]: A dictionary containing related force field parameters.
        """
        results = {
            "ff_atom": None,
            "pair_coef": None,
            "bond_coef": [],
            "angle_coef": [],
            "torsion_coef": [],
            "improper_coef": []
        }

        for atom in self.ff_atoms:
            if atom.type == identifier:
                results["ff_atom"] = atom
        
        for pair_coef in self.pair_coefs:
            if pair_coef.ff_atom.type == identifier:
                results["pair_coef"] = pair_coef
        
        for bond in self.bond_coefs:
            if any(atom.type == identifier for atom in bond.ff_atoms):
                results["bond_coef"].append(bond)
        
        for angle in self.angle_coefs:
            if any(atom.type == identifier for atom in angle.ff_atoms):
                results["angle_coef"].append(angle)
        
        for torsion in self.torsion_coefs:
            if any(atom.type == identifier for atom in torsion.ff_atoms):
                results["torsion_coef"].append(torsion)
        
        for improper in self.improper_coefs:
            if any(atom.type == identifier for atom in improper.ff_atoms):
                results["improper_coef"].append(improper)
        
        return results

    def find_ff_atom(self, atom_type: str, iteration: int) -> Optional[FF_atom]:
        """
        Find an FF_atom by its type, considering duplicates.

        Args:
            atom_type (str): The atom type.
            iteration (int): The current file iteration index.

        Returns:
            Optional[FF_atom]: The found FF_atom object or None.
        """
        if iteration > 0 and atom_type in self.duplicates:
            atom_type += f'[dup{iteration}]'

        if atom_type == '*':
            ff_atom = FF_atom()
            ff_atom.type = '*'
            return ff_atom

        for ff_atom in self.ff_atoms:
            if ff_atom.type == atom_type:
                return ff_atom

        return None

    def read_ff_atom_types(self, file_path: str, iteration: int) -> None:
        """
        Read force field atom types from a file.

        Args:
            file_path (str): Path to the force field file.
            iteration (int): The current file iteration index.
        """
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#') and start_processing:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>'):
                    continue
                if line.strip().startswith('#atom_types'):
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()
                
                    ff_atom = FF_atom()
                    ff_atom.type = parts[2]
                    ff_atom.mass = float(parts[3])
                    ff_atom.element = parts[4]

                    try:
                        ff_atom.connections = int(parts[5]) if parts[5] else None
                    except (IndexError, ValueError):
                        ff_atom.connections = None

                    try:
                        charge_converted = float(parts[-1])
                        ff_atom.charge = charge_converted
                        ff_atom.description = ' '.join(parts[6:-1])
                    except ValueError:
                        ff_atom.description = ' '.join(parts[6:])

                    existing_types = {ff_atom.type for ff_atom in self.ff_atoms}
                    if ff_atom.type in existing_types:
                        self.duplicates.add(ff_atom.type)
                        ff_atom.type += f'[dup{iteration}]'
                    
                    self.ff_atoms.append(ff_atom)

    def read_ff_equivalences(self, file_path: str, iteration: int) -> None:
        """
        Read force field equivalences from a file.

        Args:
            file_path (str): Path to the force field file.
            iteration (int): The current file iteration index.
        """
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#') and start_processing:
                    return
                if line.startswith('!') or line.strip() == '' or line.startswith('>'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#equivalence' and line_parts[1] == 'cvff':
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()

                    ff_equivalence = FF_equivalence()
                    ff_equivalence.ff_atom = self.find_ff_atom(parts[2], iteration)
                    ff_equivalence.pair_coef = self.find_ff_atom(parts[3], iteration)
                    ff_equivalence.bond_coef = self.find_ff_atom(parts[4], iteration)
                    ff_equivalence.angle_coef = self.find_ff_atom(parts[5], iteration)
                    ff_equivalence.torsion_coef = self.find_ff_atom(parts[6], iteration)
                    ff_equivalence.improper_coef = self.find_ff_atom(parts[7], iteration)
                    
                    if ff_equivalence.ff_atom is not None:
                        self.ff_equivalences.append(ff_equivalence)

    def read_ff_pair_coef(self, file_path: str, iteration: int) -> None:
        """
        Read force field pair coefficients from a file.

        Args:
            file_path (str): Path to the force field file.
            iteration (int): The current file iteration index.
        """

        def calculate_sigma_epsilon(A: float, B: float) -> Tuple[float, float]:
            if A == 0 or B == 0:
                return 0.0, 0.0
            sigma = (A / B) ** (1/6)
            epsilon = B / (4 * sigma**6)
            return sigma, epsilon

        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#') and start_processing:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>') or line.startswith('@'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#nonbond(12-6)' and line_parts[1] == 'cvff':
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()
                    ff_pair_coef = FF_pair_coef()
                    ff_pair_coef.ff_atom = self.find_ff_atom(parts[2], iteration)
                    A, B = float(parts[3]), float(parts[4])
                    ff_pair_coef.sigma, ff_pair_coef.epsilon = calculate_sigma_epsilon(A, B)
                    
                    self.pair_coefs.append(ff_pair_coef)

    def read_ff_bond(self, file_path: str, iteration: int) -> None:
        """
        Read force field bond coefficients from a file.

        Args:
            file_path (str): Path to the force field file.
            iteration (int): The current file iteration index.
        """
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#') and start_processing:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>') or line.startswith('@'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#quadratic_bond' and line_parts[1] == 'cvff': 
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()

                    ff_bond = FF_bond_coef()
                    ff_bond.ff_atoms = [self.find_ff_atom(parts[2], iteration), self.find_ff_atom(parts[3], iteration)]
                    ff_bond.r0 = float(parts[4])
                    ff_bond.k = float(parts[5])
                    
                    self.bond_coefs.append(ff_bond)

    def read_ff_angle(self, file_path: str, iteration: int) -> None:
        """
        Read force field angle coefficients from a file.

        Args:
            file_path (str): Path to the force field file.
            iteration (int): The current file iteration index.
        """
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#') and start_processing:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>') or line.startswith('@'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#quadratic_angle' and line_parts[1] == 'cvff':    
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()

                    ff_angle = FF_angle_coef()
                    ff_angle.ff_atoms = [self.find_ff_atom(parts[2], iteration), self.find_ff_atom(parts[3], iteration), self.find_ff_atom(parts[4], iteration)]
                    ff_angle.theta0 = float(parts[5])
                    ff_angle.k = float(parts[6])
                    
                    self.angle_coefs.append(ff_angle)

    def read_ff_torsion(self, file_path: str, iteration: int) -> None:
        """
        Read force field torsion coefficients from a file.

        Args:
            file_path (str): Path to the force field file.
            iteration (int): The current file iteration index.
        """
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#') and start_processing:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>') or line.startswith('@'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#torsion_1' and line_parts[1] == 'cvff':
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()

                    ff_torsion = FF_torsion_coef()
                    ff_torsion.ff_atoms = [self.find_ff_atom(parts[2], iteration), self.find_ff_atom(parts[3], iteration), self.find_ff_atom(parts[4], iteration), self.find_ff_atom(parts[5], iteration)]
                    ff_torsion.kphi = float(parts[6])
                    ff_torsion.phi0 = float(parts[8])
                    ff_torsion.n = int(parts[7])
                    
                    self.torsion_coefs.append(ff_torsion)

    def read_ff_improper(self, file_path: str, iteration: int) -> None:
        """
        Read force field improper coefficients from a file.

        Args:
            file_path (str): Path to the force field file.
            iteration (int): The current file iteration index.
        """
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#') and start_processing:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>') or line.startswith('@'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#out_of_plane' and line_parts[1] == 'cvff':    
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()

                    ff_improper = FF_improper_coef()
                    ff_improper.ff_atoms = [self.find_ff_atom(parts[2], iteration), self.find_ff_atom(parts[3], iteration), self.find_ff_atom(parts[4], iteration), self.find_ff_atom(parts[5], iteration)]
                    ff_improper.kchi = float(parts[6])
                    ff_improper.n = int(parts[7])
                    ff_improper.chi0 = float(parts[8])
                    
                    self.improper_coefs.append(ff_improper)

    def add_pair_coef_equivalences(self) -> None:
        """
        Add pair coefficients based on equivalences.
        """
        equivalences = [equivalence for equivalence in self.ff_equivalences]
        non_bonds = [pair_coef.ff_atom for pair_coef in self.pair_coefs]

        for equivalence in equivalences:
            if equivalence.ff_atom not in non_bonds:
                new_ff_pair_coef = FF_pair_coef()
                new_ff_pair_coef.ff_atom = equivalence.ff_atom
                
                for pair_coef in self.pair_coefs:
                    if equivalence.pair_coef == pair_coef.ff_atom:
                        new_ff_pair_coef.epsilon = pair_coef.epsilon
                        new_ff_pair_coef.sigma = pair_coef.sigma
                        break

                if new_ff_pair_coef.sigma is not None and new_ff_pair_coef.epsilon is not None:  
                    self.pair_coefs.append(new_ff_pair_coef)
            
    def add_bond_equivalences(self) -> None:
        """
        Add bond coefficients based on equivalences.
        """
        equivalences = [equivalence for equivalence in self.ff_equivalences]
        bond_atoms = [atom for bond_coef in self.bond_coefs for atom in bond_coef.ff_atoms]

        for equivalence in equivalences:
            if equivalence.ff_atom not in bond_atoms:
                new_ff_bond = FF_bond_coef()

                for bond_coef in self.bond_coefs:
                    if bond_coef.ff_atoms is not None:
                        if equivalence.bond_coef == bond_coef.ff_atoms[0] and equivalence.bond_coef == bond_coef.ff_atoms[1]:
                            new_ff_bond.ff_atoms = [equivalence.ff_atom, equivalence.ff_atom]
                            new_ff_bond.k = bond_coef.k
                            new_ff_bond.r0 = bond_coef.r0
                            break
                        elif equivalence.bond_coef is bond_coef.ff_atoms[0]:
                            new_ff_bond.ff_atoms = [equivalence.ff_atom, bond_coef.ff_atoms[1]]
                            new_ff_bond.k = bond_coef.k
                            new_ff_bond.r0 = bond_coef.r0
                            break
                        elif equivalence.bond_coef is bond_coef.ff_atoms[1]:
                            new_ff_bond.ff_atoms = [bond_coef.ff_atoms[0], equivalence.ff_atom]
                            new_ff_bond.k = bond_coef.k
                            new_ff_bond.r0 = bond_coef.r0
                            break

                if new_ff_bond.ff_atoms is not None:           
                    self.bond_coefs.append(new_ff_bond)

    def add_angle_equivalences(self) -> None:
        """
        Add angle coefficients based on equivalences.
        """
        angle_atoms = {atom for angle_coef in self.angle_coefs for atom in angle_coef.ff_atoms}

        remaining_equivalences = [
            equivalence for equivalence in self.ff_equivalences
            if equivalence.ff_atom not in angle_atoms
        ]

        for equivalence in remaining_equivalences:
            new_ff_angle = FF_angle_coef()

            for angle_coef in self.angle_coefs:
                if angle_coef.ff_atoms:
                    match_indices = [i for i in range(3) if equivalence.angle_coef == angle_coef.ff_atoms[i]]

                    if match_indices:
                        new_atoms = angle_coef.ff_atoms[:]
                    
                        for index in match_indices:
                            new_atoms[index] = equivalence.ff_atom

                        new_ff_angle.ff_atoms = new_atoms
                        new_ff_angle.k = angle_coef.k
                        new_ff_angle.theta0 = angle_coef.theta0
                        break
            if new_ff_angle.ff_atoms:
                self.angle_coefs.append(new_ff_angle)

    def add_torsion_equivalences(self) -> None:
        """
        Add torsion coefficients based on equivalences.
        """
        torsion_atoms = {atom for torsion_coef in self.torsion_coefs for atom in torsion_coef.ff_atoms}

        remaining_equivalences = [
            equivalence for equivalence in self.ff_equivalences
            if equivalence.ff_atom not in torsion_atoms
        ]

        for equivalence in remaining_equivalences:
            new_ff_torsion = FF_torsion_coef()

            for torsion_coef in self.torsion_coefs:
                if torsion_coef.ff_atoms:
                    match_indices = [i for i in range(4) if equivalence.torsion_coef == torsion_coef.ff_atoms[i]]

                    if match_indices:
                        new_atoms = torsion_coef.ff_atoms[:]
                        
                        for index in match_indices:
                            new_atoms[index] = equivalence.ff_atom

                        new_ff_torsion.ff_atoms = new_atoms
                        new_ff_torsion.kphi = torsion_coef.kphi
                        new_ff_torsion.n = torsion_coef.n
                        new_ff_torsion.phi0 = torsion_coef.phi0
                        break

            if new_ff_torsion.ff_atoms:
                self.torsion_coefs.append(new_ff_torsion)

    def add_improper_equivalences(self) -> None:
        """
        Add improper coefficients based on equivalences.
        """
        improper_atoms = {atom for improper_coef in self.improper_coefs for atom in improper_coef.ff_atoms}
        remaining_equivalences = [
            equivalence for equivalence in self.ff_equivalences
            if equivalence.ff_atom not in improper_atoms
        ]

        for equivalence in remaining_equivalences:
            new_ff_improper = FF_improper_coef()

            for improper_coef in self.improper_coefs:
                if improper_coef.ff_atoms:
                    match_indices = [i for i in range(4) if equivalence.improper_coef == improper_coef.ff_atoms[i]]

                    if match_indices:
                        new_atoms = improper_coef.ff_atoms[:]
                        
                        for index in match_indices:
                            new_atoms[index] = equivalence.ff_atom

                        new_ff_improper.ff_atoms = new_atoms
                        new_ff_improper.kchi = improper_coef.kchi
                        new_ff_improper.n = improper_coef.n
                        new_ff_improper.chi0 = improper_coef.chi0
                        break

            if new_ff_improper.ff_atoms:
                self.improper_coefs.append(new_ff_improper)

    def add_pair_coefs_for_ff_types(self) -> None:
        """
        Assign pair coefficients to each FF_atom.
        """
        for ff_atom in self.ff_atoms:
            for pair_coef in self.pair_coefs:
                if pair_coef.ff_atom == ff_atom:
                    ff_atom.ff_pair_coef = pair_coef
                    break
