from code.classes.SystemParts import Molecule
from code.classes.SystemBuilder import SystemBuilder
from code.classes.ForceFieldLoader import ForceFieldLoader
from typing import List
import math
from itertools import combinations, permutations

class SystemAllocator:
    def __init__(self,
                 input_file,
                replication,
                height,
                al_mg_ratio,
                ca_si_ratio,
                water_file,
                vdw_radii_file,
                vdw_def_radius,
                vdw_def_scale,
                ff_files,
                mg_cutoff,
                h_cutoff,
                bond_cutoff,
                ff_params):

        self.system = SystemBuilder(input_file, replication, height, al_mg_ratio, ca_si_ratio, water_file, vdw_radii_file, vdw_def_radius, vdw_def_scale)
        self.dimensions = self.system.dimensions
        self.molecules: List[Molecule] = self.system.molecules

        self.forcefieldloader = ForceFieldLoader(ff_files)
        self.params = ff_params

        self.mg_cutoff = mg_cutoff
        self.h_cutoff = h_cutoff
        self.bond_cutoff = bond_cutoff

        self.ff_attributes = {}
        for ff_type in self.params:
            dictionary = self.forcefieldloader.get_info_by_type(ff_type)
            self.ff_attributes[ff_type] = dictionary

        self.preprocess_all()

        self.pair_coeffs = []

        self.allocate_ff_atoms()
        self.allocate_bonds()


    def allocate_ff_atoms(self):
        """
        THESE FUNCTIONS ARE MANUALLY ADDED AND HARDCODED.
        SYSTEM DEPENDENT
        """
        self.add_clay_ff()
        self.add_clay_ff_oxygens()
        self.add_solvent()
        self.add_ions()

    def preprocess_all(self):
        self.preprocess_bond_coefs()
        self.preprocess_angle_coefs()

    def allocate_bonds(self):
        self.add_pair_coefficients()
        self.add_bonds()
        self.add_angles()
        self.add_torsion()
        self.add_improper()

    def preprocess_bond_coefs(self):
        for key, attributes in self.ff_attributes.items():
            if "bond_coef" in attributes:
                bond_dict = {}
                for bond in attributes["bond_coef"]:
                    pair = tuple(sorted((bond.ff_atoms[0].type, bond.ff_atoms[1].type)))
                    bond_dict[pair] = bond
                self.ff_attributes[key]["bond_coef"] = bond_dict

    def preprocess_angle_coefs(self):
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

    def find_min_max_z(self):
        min_z = float('inf')  
        max_z = float('-inf') 
        for molecule in self.molecules:
            if molecule.type == "clay":
                for atom in molecule.atoms:
                    if atom.element == "O":
                        z_position = atom.position[2]
                        if z_position < min_z:
                            min_z = z_position
                        if z_position > max_z:
                            max_z = z_position
        return min_z, max_z
    
    def distance(self, atom1, atom2):
        normal_distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(atom1.position, atom2.position)))

        pbc_distance = []
        for (dim_lo, dim_hi), (a, b) in zip(self.dimensions, zip(atom1.position, atom2.position)):
            length = dim_hi - dim_lo  
            delta = b - a
            delta -= round(delta / length) * length
            pbc_distance.append(delta ** 2)
        
        pbc_distance = math.sqrt(sum(pbc_distance))
        min_distance = min(normal_distance, pbc_distance)
        return min_distance

    def add_clay_ff(self):
        for molecule in self.molecules:
            if molecule.type == "clay":
                for atom in molecule.atoms:
                    # Al aluminiums are octahedral aluminums
                    if atom.element == "Al":
                        atom.ff_atom = self.ff_attributes["ao"]["ff_atom"]

                    # Mg magnesiums are octahedral magnesiums
                    if atom.element == "Mg":
                        atom.ff_atom = self.ff_attributes["mgo"]["ff_atom"]

                    # Si silicons are tetrahedral silicons
                    if atom.element == "Si":
                        atom.ff_atom = self.ff_attributes["st"]["ff_atom"]

                    # Hydrogens are hydroxyl hydrogens
                    if atom.element == "H":
                        atom.ff_atom = self.ff_attributes["ho"]["ff_atom"]
    
    def add_clay_ff_oxygens(self):
        min_z, max_z = self.find_min_max_z()

        for molecule in self.molecules:
            if molecule.type == "clay":
                height_dictionary = {i: [] for i in range(4)} 

                # Classify oxygen atoms into rows based on their z position
                for atom in molecule.atoms:
                    if atom.element == "O":
                        z_position = atom.position[2]
                        row_number = int((z_position - min_z) / (max_z - min_z) * 4)
                        row_number = min(row_number, 3) 
                        height_dictionary.setdefault(row_number, []).append(atom)

                for row in range(4):
                    if height_dictionary[row]:  
                        for oxygen_atom in height_dictionary[row]:
                            close_to_magnesium = any(atom.element == "Mg" and self.distance(oxygen_atom, atom) <= self.mg_cutoff for atom in molecule.atoms)
                            close_to_hydrogen = any(atom.element == "H" and self.distance(oxygen_atom, atom) <= self.h_cutoff for atom in molecule.atoms)
                            ff_type = 'ob' if row in [0, 3] else ('ohs' if close_to_magnesium and close_to_hydrogen else 'obos' if close_to_magnesium else 'oh' if close_to_hydrogen else 'ob')
                            
                            oxygen_atom.ff_atom = self.ff_attributes[ff_type]["ff_atom"]
                            
    def add_ions(self):
        for molecule in self.molecules:
            if molecule.type == "ion":
                if molecule.atoms[0].element == "Ca":
                    molecule.atoms[0].ff_atom = self.ff_attributes["Ca"]["ff_atom"]

    def add_solvent(self):
        for molecule in self.molecules:
            if molecule.type == "water":
                for atom in molecule.atoms:
                    if atom.element == "O":
                        atom.ff_atom = self.ff_attributes["o*"]["ff_atom"]
                    if atom.element == "H":
                        atom.ff_atom = self.ff_attributes["h*"]["ff_atom"]
    
    def add_pair_coefficients(self):
        used_ff = set([atom.ff_atom.type for molecule in self.molecules for atom in molecule.atoms])
        for ff_type in used_ff:
            self.pair_coeffs.append(self.ff_attributes[ff_type]["nonbond_coef"])

    def add_bonds(self):
        for molecule in self.molecules:
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

                if distance > self.bond_cutoff:
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
        
    def add_angles(self):
        for molecule in self.molecules:
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






    

                

            



            








    def add_torsion(self):
        pass

    def add_improper(self):
        pass



   

    
