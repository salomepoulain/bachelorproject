from code.classes.SystemParts import Molecule
from code.classes.ClayBuilder import ClayBuilder
from code.classes.ForceFieldLoader import ForceFieldLoader
from typing import Dict, List, Set, Tuple
from itertools import combinations, permutations, product

"""
This script defines the SystemAllocator class, which is responsible for allocating force field
parameters to atoms in the system based on input settings. The class allocates force field
parameters for clay atoms, water molecules, and ions, and adds bonds, angles, torsions, and
impropers to the system based on the allocated parameters.
"""

class SystemAllocator(ClayBuilder):
    def __init__(self, settings) -> None:
        """
        Initialize the SystemAllocator with the given settings and allocate force field attributes.

        Args:
            settings (ChosenSettings): The settings object containing various parameters.
        """
        super().__init__(settings)
        self.s = settings
        self.forcefieldloader = ForceFieldLoader(self.s)
        self.ff_attributes = self.add_ff_attributes()
        self.preprocess_all()
        self.process_types = {'clay'}
        self.allocate_ff_atoms(self.process_types)
        self.allocate_bonds(self.process_types)

    def add_ff_attributes(self) -> Dict[str, Dict]:
        """
        Add force field attributes for all specified atom types.

        Returns:
            Dict[str, Dict]: A dictionary containing force field parameters for each atom type.
        """
        super_dictionary = {}
        for ff_type in self.s.ff_atom_types:
            dictionary = self.forcefieldloader.get_info_by_type(ff_type)
            super_dictionary[ff_type] = dictionary
        return super_dictionary

    def allocate_ff_atoms(self, process_types: Set[str]) -> None:
        """
        Allocate force field atoms based on process types.

        Args:
            process_types (Set[str]): The set of process types to allocate force field atoms for.
        """
        if 'clay' in process_types:
            self.add_clay_ff()
            self.add_clay_ff_oxygens()
        if 'water' in process_types:
            self.add_solvent()
        if 'ion' in process_types:
            self.add_ions()

    def allocate_bonds(self, process_types: Set[str]) -> None:
        """
        Allocate bonds, angles, torsions, and impropers based on process types.

        Args:
            process_types (Set[str]): The set of process types to allocate bonds for.
        """
        for process_type in process_types:
            self.add_bonds(process_type)
            self.add_angles(process_type)
            self.add_torsion(process_type)
            self.add_improper(process_type)

    def preprocess_all(self) -> None:
        """
        Preprocess bond, angle, and torsion coefficients for efficient lookup.
        """
        self.preprocess_bond_coefs()
        self.preprocess_angle_coefs()
        self.preprocess_torsion_coefs()

    def preprocess_bond_coefs(self) -> None:
        """
        Preprocess bond coefficients and store them in a dictionary.
        """
        for key, attributes in self.ff_attributes.items():
            if "bond_coef" in attributes:
                bond_dict = {}
                for bond in attributes["bond_coef"]:
                    pair = tuple(sorted((bond.ff_atoms[0].type, bond.ff_atoms[1].type)))
                    bond_dict[pair] = bond
                self.ff_attributes[key]["bond_coef"] = bond_dict

    def preprocess_angle_coefs(self) -> None:
        """
        Preprocess angle coefficients and store them in a dictionary.
        """
        for key, attributes in self.ff_attributes.items():
            if "angle_coef" in attributes:
                angle_dict = {}
                for angle in attributes["angle_coef"]:
                    atoms = angle.ff_atoms
                    if len(atoms) == 3:
                        sorted_ends = sorted([atoms[0].type, atoms[2].type])
                        angle_key = (sorted_ends[0], atoms[1].type, sorted_ends[1])
                    else:
                        continue  
                    angle_dict[angle_key] = angle
                self.ff_attributes[key]["angle_coef"] = angle_dict

    def preprocess_torsion_coefs(self) -> None:
        """
        Preprocess torsion coefficients and store them in a dictionary, including wildcard combinations.
        """
        for key, attributes in self.ff_attributes.items():
            if "torsion_coef" in attributes:
                torsion_dict = {}
                for torsion in attributes["torsion_coef"]:
                    atoms = torsion.ff_atoms

                    if len(atoms) == 4:
                        atom_types = []
                        for atom in atoms:
                            if atom.type == '*':
                                atom_types.append(list(self.ff_attributes.keys()))
                            else:
                                atom_types.append([atom.type])

                        combinations = product(*atom_types)

                        for comb in combinations:
                            sorted_ends = sorted([comb[0], comb[3]])
                            torsion_key = (sorted_ends[0], comb[1], comb[2], sorted_ends[1])
                            torsion_dict[torsion_key] = torsion

                            backward_key = tuple(reversed(torsion_key))
                            torsion_dict[backward_key] = torsion

                self.ff_attributes[key]["torsion_coef"] = torsion_dict

    def add_clay_ff(self) -> None:
        """
        Assign force field parameters to clay atoms based on their element type.
        """
        for molecule in self.molecules:
            if molecule.type == "clay":
                for atom in molecule.atoms:
                    if atom.element == "Al":
                        atom.ff_atom = self.ff_attributes["ao"]["ff_atom"]
                    if atom.element == "Mg":
                        atom.ff_atom = self.ff_attributes["mgo"]["ff_atom"]
                    if atom.element == "Si":
                        atom.ff_atom = self.ff_attributes["st"]["ff_atom"]
                    if atom.element == "H":
                        atom.ff_atom = self.ff_attributes["ho"]["ff_atom"]

    def add_clay_ff_oxygens(self) -> None:
        """
        Assign force field parameters to oxygen atoms in clay based on their position and nearby atoms.
        """
        for molecule in self.molecules:
            if molecule.type == "clay":
                min_z, max_z = self.find_min_max_z_for_molecule(molecule)
                
                height_dictionary = {i: [] for i in range(4)}  

                for atom in molecule.atoms:
                    if atom.element == "O":
                        z_position = atom.position[2]
                        row_number = int((z_position - min_z) / (max_z - min_z) * 4)
                        row_number = min(row_number, 3)
                        height_dictionary[row_number].append(atom)

                for row in height_dictionary:
                    for oxygen_atom in height_dictionary[row]:
                        close_to_magnesium = any(
                            atom.element == "Mg" and self.distance(oxygen_atom, atom) <= self.s.mg_cutoff
                            for atom in molecule.atoms
                        )
                        close_to_hydrogen = any(
                            atom.element == "H" and self.distance(oxygen_atom, atom) <= self.s.h_cutoff
                            for atom in molecule.atoms
                        )
                        ff_type = 'ob' if row in [0, 3] else (
                            'ohs' if close_to_magnesium and close_to_hydrogen else 
                            'obos' if close_to_magnesium else 
                            'oh' if close_to_hydrogen else 
                            'ob'
                        )
                        
                        oxygen_atom.ff_atom = self.ff_attributes[ff_type]["ff_atom"]

    def find_min_max_z_for_molecule(self, molecule: Molecule) -> Tuple[float, float]:
        """
        Find the minimum and maximum z-coordinates for oxygen atoms in a molecule.

        Args:
            molecule (Molecule): The molecule to analyze.

        Returns:
            Tuple[float, float]: The minimum and maximum z-coordinates.
        """
        min_z = float('inf')
        max_z = float('-inf')
        for atom in molecule.atoms:
            if atom.element == "O":
                z_position = atom.position[2]
                if z_position < min_z:
                    min_z = z_position
                if z_position > max_z:
                    max_z = z_position
        return min_z, max_z

    def add_solvent(self) -> None:
        """
        Assign force field parameters to water molecules.
        """
        for molecule in self.molecules:
            if molecule.type == "water":
                for atom in molecule.atoms:
                    if atom.element == "O":
                        atom.ff_atom = self.ff_attributes["o*"]["ff_atom"]
                    if atom.element == "H":
                        atom.ff_atom = self.ff_attributes["h*"]["ff_atom"]
                            
    def add_ions(self) -> None:
        """
        Assign force field parameters to ions.
        """
        for molecule in self.molecules:
            if molecule.type == "ion":
                if molecule.atoms[0].element == "Ca":
                    molecule.atoms[0].ff_atom = self.ff_attributes["Ca"]["ff_atom"]

    def add_bonds(self, process_type: str) -> None:
        """
        Add bonds to molecules based on force field parameters.

        Args:
            process_type (str): The type of molecules to process.
        """
        for molecule in self.molecules:
            if molecule.type != process_type:
                continue

            all_atoms = molecule.atoms
            distance_cache = {}
            bond_cache = {}

            for atom1, atom2 in combinations(all_atoms, 2):
                pair_key = (id(atom1), id(atom2))

                if pair_key not in distance_cache:
                    distance = self.distance(atom1, atom2)
                    distance_cache[pair_key] = distance
                else:
                    distance = distance_cache[pair_key]

                if distance > self.s.bond_cutoff:
                    continue

                if molecule.is_bonded(atom1, atom2):
                    continue

                bond_types = tuple(sorted((atom1.ff_atom.type, atom2.ff_atom.type)))

                if bond_types not in bond_cache:
                    bond_dict = self.ff_attributes.get(atom1.ff_atom.type, {}).get("bond_coef", {})
                    bond_coef = bond_dict.get(bond_types)
                    bond_cache[bond_types] = bond_coef
                else:
                    bond_coef = bond_cache[bond_types]

                if bond_coef:
                    molecule.add_bond(atom1, atom2, bond_coef)
        
    def add_angles(self, process_type: str) -> None:
        """
        Add angles to molecules based on force field parameters.

        Args:
            process_type (str): The type of molecules to process.
        """
        for molecule in self.molecules:
            if molecule.type != process_type:
                continue

            if len(molecule.atoms) < 3:
                continue

            atom_types_present = {atom.ff_atom.type for atom in molecule.atoms if atom.ff_atom}
            if not any(self.ff_attributes.get(atom_type, {}).get("angle_coef") for atom_type in atom_types_present):
                continue

            processed_angles = set()
            for first_atom, central_atom, last_atom in permutations(molecule.atoms, 3):
                current_permutation = (id(first_atom), id(central_atom), id(last_atom))
                if current_permutation in processed_angles:
                    continue

                processed_angles.add(current_permutation)
                if not (molecule.is_bonded(central_atom, first_atom) and molecule.is_bonded(central_atom, last_atom)):
                    continue

                if molecule.has_angle(first_atom, central_atom, last_atom):
                    continue

                angle_key = (first_atom.ff_atom.type, central_atom.ff_atom.type, last_atom.ff_atom.type)
                angle_dict = self.ff_attributes.get(central_atom.ff_atom.type, {}).get("angle_coef", {})
                angle_coef = angle_dict.get(angle_key)
                
                if angle_coef:
                    molecule.add_angle(first_atom, central_atom, last_atom, angle_coef)

    def add_torsion(self, process_type: str) -> None:
        """
        Add torsions to molecules based on force field parameters.

        Args:
            process_type (str): The type of molecules to process.
        """
        for molecule in self.molecules:
            if molecule.type != process_type:
                continue

            if len(molecule.atoms) < 4:
                continue

            atom_types_present = {atom.ff_atom.type for atom in molecule.atoms if atom.ff_atom}
            if not any(self.ff_attributes.get(atom_type, {}).get("torsion_coef") for atom_type in atom_types_present):
                continue

            processed_torsions = set()
            for torsion_atoms in permutations(molecule.atoms, 4):
                current_permutation = tuple(id(atom) for atom in torsion_atoms)
                if current_permutation in processed_torsions:
                    continue

                processed_torsions.add(current_permutation)

                # Check if A-B, B-C, and C-D are connected
                if not (molecule.is_bonded(torsion_atoms[0], torsion_atoms[1]) and
                        molecule.is_bonded(torsion_atoms[1], torsion_atoms[2]) and
                        molecule.is_bonded(torsion_atoms[2], torsion_atoms[3])):
                    continue

                if molecule.has_torsion(*torsion_atoms):
                    continue

                torsion_key = tuple(atom.ff_atom.type for atom in torsion_atoms)
                torsion_dict = self.ff_attributes.get(torsion_atoms[1].ff_atom.type, {}).get("torsion_coef", {})
                torsion_coef = torsion_dict.get(torsion_key)

                if torsion_coef:
                    molecule.add_torsion(*torsion_atoms, torsion_coef)

    def add_improper(self, process_type: str) -> None:
        """
        Add improper torsions to molecules based on force field parameters.

        Args:
            process_type (str): The type of molecules to process.
        """
        for molecule in self.molecules:
            if molecule.type != "clay":
                continue

            # Placeholder for actual implementation
            # Here you can implement the logic for adding improper torsions based on the given process_type
