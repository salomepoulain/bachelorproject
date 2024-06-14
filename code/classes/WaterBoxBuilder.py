from code.classes.SystemParts import Atom, Molecule
from code.classes.ChosenSettings import ChosenSettings
import random
import math
from typing import List, Tuple

"""
This module defines the WaterBoxBuilder class, which is responsible for building and
manipulating water molecules based on input settings. The class loads water molecules
from an XYZ file, calculates the dimensions of the water molecules, determines the
duplication factor for water molecules, replicates water molecules to fill the system
dimensions, removes water molecules outside the system dimensions, and translates water
molecules up by the height of the clay plus the water distance.
"""

class WaterBoxBuilder:
    def __init__(self, ion_count: int, settings: ChosenSettings, seed: int = 42) -> None:
        self.s = settings
        random.seed(self.s.random_seed) 

        self.water_file = 'water/' + self.s.water_file + '.xyz'
        self.molecules: List[Molecule] = self.load_water_molecules()
        self.dimensions: Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]] = self.get_water_dimensions()
        self.ion_count: int = ion_count

    def determine_water_duplication(self) -> int:
        """
        Determine the duplication factor for water molecules based on the water per ion setting.
        """
        duplication_factor = 1
        if self.s.water_per_ion == 0 or self.ion_count == 0:
            return 0

        real_water_per_ion = len(self.molecules) / self.ion_count

        if real_water_per_ion < self.s.water_per_ion:
            duplication_factor = math.ceil(self.s.water_per_ion / real_water_per_ion)

        return duplication_factor

    def add_molecules(self, atoms: List[Tuple[str, Tuple[float, float, float]]]) -> List[Molecule]:
        """
        Create water molecules from the given atoms list.
        """
        molecules = []
        counter = 0
        water = Molecule()
        water.type = "water"
        for atom in atoms:
            element, position = atom
            if element == "O":
                water.add_atom(element, position)
                counter += 1
            elif element == "H":
                water.add_atom(element, position)
                counter += 1

            if counter == 3:
                molecules.append(water)
                counter = 0
                water = Molecule()
                water.type = "water"
        
        return molecules

    def load_water_molecules(self) -> List[Molecule]:
        """
        Load water molecules from the water file.
        """
        atoms = []
        with open(self.water_file, 'r') as file:
            lines = file.readlines()
            for line in lines[2:]:  # Skip header lines
                parts = line.strip().split()
                if len(parts) != 4:
                    continue
                element = parts[0]
                position = (float(parts[1]), float(parts[2]), float(parts[3]))
                atoms.append((element, position))

        return self.add_molecules(atoms)

    def get_water_dimensions(self) -> Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
        """
        Calculate the dimensions of the water molecules.
        """
        atoms = [atom for molecule in self.molecules for atom in molecule.atoms]

        min_x = min_y = min_z = float('inf')
        max_x = max_y = max_z = float('-inf')
        
        for atom in atoms:
            x, y, z = atom.position
            min_x = min(min_x, x)
            max_x = max(max_x, x)
            min_y = min(min_y, y)
            max_y = max(max_y, y)
            min_z = min(min_z, z)
            max_z = max(max_z, z)
        
        return ((min_x, max_x), (min_y, max_y), (min_z, max_z))

    def translate_new_origin(self, new_origin_x: float, new_origin_y: float, new_origin_z: float) -> None:
        """
        Translate water molecules to a new origin.
        """
        atoms = [atom for molecule in self.molecules for atom in molecule.atoms]

        current_min_x, _ = self.dimensions[0]
        current_min_y, _ = self.dimensions[1]
        current_min_z, _ = self.dimensions[2]

        translation_x = new_origin_x - current_min_x
        translation_y = new_origin_y - current_min_y
        translation_z = new_origin_z - current_min_z

        for atom in atoms:
            x_pos, y_pos, z_pos = atom.position
            atom.position = (x_pos + translation_x, y_pos + translation_y, z_pos + translation_z)

        self.get_water_dimensions()

    def calculate_replication(self, system_dimensions: Tuple[float, float, float]) -> Tuple[int, int, int]:
        """
        Calculate the replication factors for water molecules to fill the system dimensions.
        """
        unit_x = self.dimensions[0][1] - self.dimensions[0][0]
        unit_y = self.dimensions[1][1] - self.dimensions[1][0]
        unit_z = self.determine_water_duplication()

        n_x = math.ceil(system_dimensions[0] / unit_x)
        n_y = math.ceil(system_dimensions[1] / unit_y)

        return n_x, n_y, unit_z

    def replicate_water(self, system_dimensions: Tuple[float, float, float]) -> None:
        """
        Replicate water molecules to fill the system dimensions.
        """
        n_x, n_y, n_z = self.calculate_replication(system_dimensions)
        original_molecules = list(self.molecules)
        unit_x = self.dimensions[0][1] - self.dimensions[0][0]
        unit_y = self.dimensions[1][1] - self.dimensions[1][0]
        unit_z = self.dimensions[2][1] - self.dimensions[2][0]

        self.molecules.clear()

        for i in range(n_x):
            for j in range(n_y):
                for k in range(n_z):
                    for molecule in original_molecules:
                        new_molecule = Molecule()
                        new_molecule.type = molecule.type
                        for atom in molecule.atoms:
                            x, y, z = atom.position
                            new_x = x + i * unit_x
                            new_y = y + j * unit_y
                            new_z = z + k * unit_z

                            new_atom = Atom()
                            new_atom.element = atom.element
                            new_atom.position = (new_x, new_y, new_z)
                            new_molecule.add_atom(new_atom.element, new_atom.position)
                        self.molecules.append(new_molecule)

        self.get_water_dimensions()

    def remove_outside_box(self, system_dimensions: Tuple[float, float, float]) -> None:
        """
        Remove water molecules that are outside the given system dimensions.
        """
        new_molecules = []
        for molecule in self.molecules:
            all_atoms_within_bounds = True
            for atom in molecule.atoms:
                x, y, _ = atom.position 
                if not (0 <= x <= system_dimensions[0] and 0 <= y <= system_dimensions[1]):
                    all_atoms_within_bounds = False
                    break

            if all_atoms_within_bounds:
                new_molecules.append(molecule)

        self.molecules = new_molecules
        self.get_water_dimensions()

    def translate_up(self, clay_height: float) -> None:
        """
        Translate water molecules up by the height of the clay plus the water distance.
        """
        for molecule in self.molecules:
            for atom in molecule.atoms:
                x, y, z = atom.position
                atom.position = (x, y, z + clay_height + self.s.water_distance)

        self.get_water_dimensions()

    def change_water_per_ions(self) -> None:
        """
        Adjust the number of water molecules based on the water per ion setting.
        """
        if self.s.water_per_ion is None or self.s.water_per_ion < 0.0:
            return            
        
        desired_water_count = int(self.s.water_per_ion * self.ion_count)

        waters = [molecule for molecule in self.molecules if molecule.type == 'water']
        current_water_count = len(waters)

        if current_water_count > desired_water_count:
            remove_water_count = current_water_count - desired_water_count
            if remove_water_count > 0:
                water_to_remove = random.sample(waters, remove_water_count)
                for water in water_to_remove:
                    self.molecules.remove(water)
