from code.classes.SystemParts import Molecule
from code.classes.WaterBuilder import WaterBuilder
from typing import Tuple
import math

class SystemBuilder:
    def __init__(self, input_file, replication, height, al_mg_ratio, ca_si_ratio, water_file, vdw_radii_file, vdw_def_radius, vdw_def_scale):
        self.input_file = input_file
        self.height = height
        self.al_mg_ratio = al_mg_ratio
        self.ca_si_ratio = ca_si_ratio
        self.replication: Tuple[int, int] = replication
        
        self.dimensions = ((None, None), (None, None), (None, None))
        self.molecules = []

        self.waterbox = WaterBuilder(water_file, vdw_radii_file, vdw_def_radius, vdw_def_scale)

        self.read_xyz_clay(input_file)
        self.translate_to_zero()
        # self.rescale_coordinates()  # Not sure if needed ?????
        self.expand_and_duplicate()
        self.mg_substitute()
        self.expand_simulation_height()
        self.translate_simulation_height()
        self.add_solvent()
        self.add_ions_uniformly()
        
    def read_xyz_clay(self, input_file):
        file_path = 'unitcell/' + input_file + '.xyz'
        
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
        n_x, n_y = self.replication
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

    def mg_substitute(self):
        all_atoms = [atom for molecule in self.molecules if molecule.type == 'clay' for atom in molecule.atoms]
        substitution_ratio = self.al_mg_ratio  
        fraction = 1 / substitution_ratio
        
        cumulative_fraction = 0.0  
        substitution_count = 0

        for atom in all_atoms:
            if atom.element == 'Al':
                cumulative_fraction += fraction

                if cumulative_fraction >= 1:
                    atom.element = 'Mg'
                    substitution_count += 1
                    cumulative_fraction -= 1 

        self.write_xyz(self.input_file + "_mg")
        return

    def rescale_coordinates(self):
        target_heigth = 9.6 
        current_heigth = self.dimensions[2]
        scale_factor = target_heigth / current_heigth

        all_atoms = [atom for molecule in self.molecules for atom in molecule.atoms]
        for atom in all_atoms:
            x, y, z = atom.position
            atom.position = (x * scale_factor, y * scale_factor, z * scale_factor)

        self.get_unit_cell_dimensions()

    def expand_simulation_height(self):
        self.dimensions = (self.dimensions[0], self.dimensions[1], (self.dimensions[2][0], self.dimensions[2][0] + self.height))

    def get_clay_height(self):
        zhi = None
        zlo = None
        
        for molecule in self.molecules:
            if molecule.type == 'clay':
                for atom in molecule.atoms:
                    _, _, z = atom.position
                    
                    if zhi is None or zlo is None:
                        zhi = zlo = z
                    else:
                        zhi = max(zhi, z)
                        zlo = min(zlo, z)
            return zhi, zlo

    def translate_simulation_height(self):
        zhi, zlo = self.get_clay_height()
        distance = zhi - zlo
        half_distance = distance / 2   

        self.dimensions = (self.dimensions[0], self.dimensions[1], (self.dimensions[2][0] + half_distance, self.dimensions[2][1] + half_distance))

    def add_solvent(self):
        self.waterbox.translate_new_origin(self.dimensions[0][0], self.dimensions[1][0], self.dimensions[2][0])
        self.waterbox.replicate_water((self.dimensions[0][1] - self.dimensions[0][0], self.dimensions[1][1] - self.dimensions[1][0], self.dimensions[2][1] - self.dimensions[2][0]))
        self.waterbox.remove_outside_box((self.dimensions[0][1] - self.dimensions[0][0], self.dimensions[1][1] - self.dimensions[1][0], self.dimensions[2][1] - self.dimensions[2][0]))
        self.waterbox.remove_overlaying_water(self.molecules, self.dimensions)

        for molecule in self.waterbox.molecules:
            self.molecules.append(molecule)
        self.write_xyz(self.input_file + "_solvent")

    def custom_round(self, number):
        if (number - math.floor(number)) < 0.5:
            return int(math.floor(number))
        else:
            return int(math.ceil(number))

    def add_ions_uniformly(self):
        ion_ratio = self.ca_si_ratio

        if any(molecule.type == 'clay' for molecule in self.molecules):
            total_si = sum(1 for molecule in self.molecules if molecule.type == 'clay' for atom in molecule.atoms if atom.element == 'Si')
            ion_count = self.custom_round(total_si * ion_ratio)

            if ion_count == 0:
                return

            ca_height = (self.dimensions[2][0] + self.dimensions[2][1]) / 2

            n_x = int(ion_count ** 0.5)
            n_y = n_x if n_x * n_x == ion_count else n_x + 1

            x_spacing = (self.dimensions[0][1] - self.dimensions[0][0]) / (n_x if n_x > 1 else 1)
            y_spacing = (self.dimensions[1][1] - self.dimensions[1][0]) / (n_y if n_y > 1 else 1)

            ca_ions_added = 0

            for j in range(n_y):
                for i in range(n_x):
                    if ca_ions_added >= ion_count:
                        break
                    x = self.dimensions[0][0] + (i * x_spacing) % (self.dimensions[0][1] - self.dimensions[0][0])
                    y = self.dimensions[1][0] + (j * y_spacing) % (self.dimensions[1][1] - self.dimensions[1][0])
                    
                    ion = Molecule()
                    ion.type = 'ion'
                    ion.add_atom('Ca', (x, y, ca_height))
                    self.molecules.append(ion)

                    ca_ions_added += 1

            self.write_xyz(self.input_file + "_mg_ion")


