from code.classes.ForcefieldPartials import FF_atom, FF_nonbond_coef, FF_bond_coef, FF_angle_coef, FF_torsion_coef, FF_improper_coef
from code.classes.ForceFieldLoader import ForceFieldLoader
from code.classes.UnitcellPartials import Molecule, Atom
from code.classes.UnitcellBuilder import UnitcellBuilder
from code.classes.SystemPartials import Bond, Angle, Torsion, Improper
from typing import List, Set
import math


class SystemAllocator:
    def __init__(self, file_name, replication, height, al_mg_ratio) -> None:

        # Before allocation:

        self.unitcell = UnitcellBuilder(file_name, replication, height, al_mg_ratio)
        self.dimensions = self.unitcell.dimensions
        self.atoms: List[Atom] = self.unitcell.atoms
        self.molecules: List[Molecule] = self.unitcell.molecules

        self.forcefieldloader = ForceFieldLoader()
        self.ff_atoms: List[FF_atom] = self.forcefieldloader.ff_atoms
        self.distance_threshold_mg = 2.0
        self.distance_threshold_h = 1.0
        self.bond_cutoff = 1.0

        self.nonbond_coefs: List[FF_nonbond_coef] = self.forcefieldloader.nonbond_coefs
        self.bond_coefs: List[FF_bond_coef] = self.forcefieldloader.bond_coefs
        self.angle_coefs: List[FF_angle_coef] = self.forcefieldloader.angle_coefs
        self.torsion_coefs: List[FF_torsion_coef] = self.forcefieldloader.torsion_coefs
        self.improper_coefs: List[FF_improper_coef] = self.forcefieldloader.improper_coefs

        # After allocation:

        self.used_ff_atoms: Set[FF_atom] = set()
        self.used_nonbond_coefs: Set[FF_nonbond_coef] = set()
        self.used_bond_coefs: Set[FF_bond_coef] = set()
        self.used_angle_coefs: Set[FF_angle_coef] = set()
        self.used_torsion_coefs: Set[FF_torsion_coef] = set()
        self.used_improper_coefs: Set[FF_improper_coef] = set()

        self.bonds: List[Bond] = []
        self.angles: List[Angle] = []
        self.torsions: List[Torsion] = []
        self.impropers: List[Improper] = []

        self.allocate_ff_atoms()
        self.allocate_bonds()

    '''
    Helper function to find the minimum and maximum z position of the clay atoms
    '''
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
    
    '''
    Helper function to calculate the distance between two atoms
    '''
    def distance(self, atom1, atom2):
        normal_distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(atom1.position, atom2.position)))
        pbc_distance = []
        for dim, (a, b) in zip(self.dimensions, zip(atom1.position, atom2.position)):
            delta = b - a
            delta -= round(delta / dim) * dim  
            pbc_distance.append(delta ** 2)
        pbc_distance = math.sqrt(sum(pbc_distance))

        return min(normal_distance, pbc_distance)
    
    def allocate_ff_atoms(self):
        self.add_clay_ff()
        self.add_clay_ff_oxygens()

    def allocate_bonds(self):
        self.add_nonbonds()
        self.add_bonds()

    def add_clay_ff(self):
        for molecule in self.molecules:
            if molecule.type == "clay":
                for atom in molecule.atoms:
                    # Al aluminiums are octahedral aluminums
                    if atom.element == "Al":
                        for ff_atom in self.ff_atoms:
                            if ff_atom.type == "ao":
                                atom.ff_atom = ff_atom

                                self.used_ff_atoms.add(ff_atom)

                    # Mg magnesiums are octahedral magnesiums
                    if atom.element == "Mg":
                        for ff_atom in self.ff_atoms:
                            if ff_atom.type == "mgo":
                                atom.ff_atom = ff_atom

                                self.used_ff_atoms.add(ff_atom)

                    # Si silicons are tetrahedral silicons
                    if atom.element == "Si":
                        for ff_atom in self.ff_atoms:
                            if ff_atom.type == "st":
                                atom.ff_atom = ff_atom

                                self.used_ff_atoms.add(ff_atom)

                    # Hydrogens are hydroxyl hydrogens
                    if atom.element == "H":
                        for ff_atom in self.ff_atoms:
                            if ff_atom.type == "ho":
                                atom.ff_atom = ff_atom

                                self.used_ff_atoms.add(ff_atom)
                                
    
    def add_clay_ff_oxygens(self):
        min_z, max_z = self.find_min_max_z()

        for molecule in self.molecules:
            if molecule.type == "clay":
                height_dictionary = {i: [] for i in range(4)}  

                for atom in molecule.atoms:
                    if atom.element == "O":
                        z_position = atom.position[2]
                        # Normalize z_position to range 0 to less than 1, then multiply by 4 and take the floor integer
                        row_number = int((z_position - min_z) / (max_z - min_z) * 4)
                        # Ensure row_number is capped at 3
                        row_number = min(row_number, 3)
                        height_dictionary.setdefault(row_number, []).append(atom)


                for row in range(4):
                    if height_dictionary[row]:
                        for oxygen_atom in height_dictionary[row]:
                            close_to_magnesium = any(atom.element == "Mg" and self.distance(oxygen_atom, atom) <= self.distance_threshold_mg for atom in molecule.atoms)
                            close_to_hydrogen = any(atom.element == "H" and self.distance(oxygen_atom, atom) <= self.distance_threshold_h for atom in molecule.atoms)
                            ff_type = 'ob' if row in [0, 3] else ('ohs' if close_to_magnesium and close_to_hydrogen else 'obos' if close_to_magnesium else 'oh' if close_to_hydrogen else 'ob')
                            for ff_atom in self.ff_atoms:
                                if ff_atom.type == ff_type:
                                    oxygen_atom.ff_atom = ff_atom
                                    self.used_ff_atoms.add(ff_atom)
                                    break

    def add_nonbonds(self):
        for nonbond_coef in self.nonbond_coefs:
            for ff_atom in self.used_ff_atoms:
                if ff_atom == nonbond_coef.ff_atoms:
                    self.used_nonbond_coefs.add(nonbond_coef)
                    break

    def add_bonds(self):    
        for molecule in self.molecules:
            # Create a dictionary to map ff_atom types to atom objects for quick lookup
            atom_dict = {}
            for atom in molecule.atoms:
                if atom.ff_atom not in atom_dict:
                    atom_dict[atom.ff_atom] = []
                atom_dict[atom.ff_atom].append(atom)

            # Loop through each bond coefficient
            for bond_coef in self.bond_coefs:
                ff_atom1, ff_atom2 = bond_coef.ff_atoms

                # Check if both atom types are present in the dictionary
                if ff_atom1 in atom_dict and ff_atom2 in atom_dict:
                    # Get all atoms of type ff_atom1
                    for atom1 in atom_dict[ff_atom1]:
                        # Get all atoms of type ff_atom2
                        for atom2 in atom_dict[ff_atom2]:
                            if atom1 != atom2 and self.distance(atom1, atom2) <= self.bond_cutoff:
                                bond = Bond()
                                bond.ff_bond_coef = bond_coef
                                bond.atoms = [atom1, atom2]
                                self.used_bond_coefs.add(bond_coef)
                                self.bonds.append(bond)     
    
    def add_angle(self):
        pass

    def add_torsion(self):
        pass

    def add_improper(self):
        pass



   

    
