from code.classes.SystemParts import Molecule
from code.classes.ChosenSettings import ChosenSettings
import random
import math

class ClayBuilder():
    def __init__(self, settings):
        super().__init__()
        self.s = settings
        
        self.molecules = []
        self.dimensions = ((None, None), (None, None), (None, None))

        self.read_xyz_clay()
        self.translate_to_zero()
        self.expand_and_duplicate()
        self.mg_substitute()
        
    def read_xyz_clay(self):
        file_path = 'unitcell/' + self.s.input_file + '.xyz'
        
        clay = Molecule()
        clay.type = "clay"
        with open(file_path, 'r') as file:
            next(file)  
            next(file)
            for line in file:
                parts = line.strip().split()
                if len(parts) == 4:
                    element, x, y, z = parts
                    position = (float(x), float(y), float(z))
                    clay.add_atom(element, position)
        
        self.molecules.append(clay)
        self.get_unit_cell_dimensions()

    def get_unit_cell_dimensions(self):
        min_x = min_y = min_z = float('inf')
        max_x = max_y = max_z = float('-inf')

        for molecule in self.molecules:
            for atom in molecule.atoms:
                x, y, z = atom.position
                min_x = min(min_x, x)
                max_x = max(max_x, x)
                min_y = min(min_y, y)
                max_y = max(max_y, y)
                min_z = min(min_z, z)
                max_z = max(max_z, z)
    
        self.dimensions = ((min_x, max_x), (min_y, max_y), (min_z, max_z))

    def translate_to_zero(self):
        min_x, min_y, min_z = self.molecules[0].atoms[0].position

        for molecule in self.molecules:
            for atom in molecule.atoms:
                x, y, z = atom.position
                min_x = min(min_x, x)
                min_y = min(min_y, y)
                min_z = min(min_z, z)

        for molecule in self.molecules:
            for atom in molecule.atoms:
                x, y, z = atom.position
                atom.position = (x - min_x, y - min_y, z - min_z)

        self.get_unit_cell_dimensions()


    def expand_and_duplicate(self):
        original_molecules = [molecule for molecule in self.molecules if molecule.type == "clay"]
        unit_x = self.dimensions[0][1] - self.dimensions[0][0]
        unit_y = self.dimensions[1][1] - self.dimensions[1][0]
        n_x = self.s.replication[0] 
        n_y = self.s.replication[1]
        self.molecules.clear()

        new_molecule = Molecule()
        new_molecule.type = original_molecules[0].type
        
        for i in range(n_x):
            for j in range(n_y):

                for atom in original_molecules[0].atoms: 
                    x, y, z = atom.position
                    new_x = x + i * unit_x  
                    new_y = y + j * unit_y 

                    new_molecule.add_atom(atom.element, (new_x, new_y, z))
                
        self.molecules.append(new_molecule)
        self.get_unit_cell_dimensions()

    def write_xyz(self, output_file):
        output_path = 'output/' + output_file + '.xyz'
        with open(output_path, 'w') as file:
            file.write(f"{len([atom for molecule in self.molecules for atom in molecule.atoms])}\n")
            file.write("Generated by Unitcell.write_xyz from " + self.input_file + ".xyz\n")
            for molecule in self.molecules:
                for atom in molecule.atoms:
                    element = atom.element
                    x, y, z = atom.position
                    file.write(f"{element} {x} {y} {z}\n")

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
    
    def custom_round(self, number):
        if (number - math.floor(number)) < 0.5:
            return int(math.floor(number))
        else:
            return int(math.ceil(number))

    def mg_substitute(self):
        random.seed(self.s.random_seed)

        all_atoms = [atom for molecule in self.molecules if molecule.type == 'clay' for atom in molecule.atoms]
        al_atoms = [atom for atom in all_atoms if atom.element == 'Al']

        substitution_ratio = self.s.al_mg_ratio  
        total_substitutions = math.ceil((len(al_atoms) / substitution_ratio))

        random.shuffle(al_atoms)
        
        substituted_atoms = []

        for al_atom in al_atoms:
            if len(substituted_atoms) >= total_substitutions:
                break
            if all(self.distance(al_atom, mg_atom) > self.s.clay_sub_cutoff for mg_atom in substituted_atoms):
                al_atom.element = 'Mg'
                substituted_atoms.append(al_atom)
        
        return



    
