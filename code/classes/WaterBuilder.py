from code.classes.SystemParts import Atom, Molecule
import math

class WaterBuilder:
    def __init__(self, water_file, vdw_radii_file, vdw_def_radius, vdw_def_scale):
        self.water_file = 'water/' + water_file + '.gro'
        self.vdw_radii_file = 'water/' + vdw_radii_file + '.dat'
        self.vdw_def_radius = vdw_def_radius
        self.vdw_def_scale = vdw_def_scale
        self.vdw_radii = {}

        self.dimensions = ((0, 0), (0, 0), (0, 0))
        self.molecules = []

        self.load_water_molecules()
        self.load_atom_radii()
        self.get_water_dimensions()

    def load_water_molecules(self):
        atoms = []
        with open(self.water_file, 'r') as file:
            lines = file.readlines()
            for line in lines[2:]:
                if len(line.strip().split()) != 6:
                    continue
            
                element = line[10:15].strip()
                if element == "OW":
                    element = "O"
                elif element == "HW1" or element == "HW2":
                    element = "H"
                x = float(line[20:28].strip()) * 10
                y = float(line[28:36].strip()) * 10
                z = float(line[36:44].strip()) * 10
                position = (x, y, z)
                atoms.append((element, position))

        self.add_molecules(atoms)
        
    def add_molecules(self, atoms):
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
                self.molecules.append(water)
                counter = 0
                water = Molecule()
                water.type = "water"

    def load_atom_radii(self):
        with open(self.vdw_radii_file, 'r') as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip().split()
                if line[0] == ";":
                    continue
                element, radius = line[1], float(line[2])
                self.vdw_radii[element] = radius
    
    def get_water_dimensions(self):
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
        
        self.dimensions = ((min_x, max_x), (min_y, max_y), (min_z, max_z))

    def translate_new_origin(self, new_origin_x, new_origin_y, new_origin_z):
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

    def calculate_replication(self, system_dimensions: tuple):
        unit_x = self.dimensions[0][1] - self.dimensions[0][0]
        unit_y = self.dimensions[1][1] - self.dimensions[1][0]
        unit_z = self.dimensions[2][1] - self.dimensions[2][0]

        n_x = math.ceil(system_dimensions[0] / unit_x)
        n_y = math.ceil(system_dimensions[1] / unit_y)
        n_z = math.ceil(system_dimensions[2] / unit_z)

        return n_x, n_y, n_z

    def replicate_water(self, system_dimensions: tuple):
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
    
    def remove_outside_box(self, system_dimensions: tuple):
        new_molecules = []
        for molecule in self.molecules:
            all_atoms_within_bounds = True
            for atom in molecule.atoms:
                x, y, z = atom.position
                if not (0 <= x <= system_dimensions[0] and 0 <= y <= system_dimensions[1] and 0 <= z <= system_dimensions[2]):
                    all_atoms_within_bounds = False
                    break

            if all_atoms_within_bounds:
                new_molecules.append(molecule)

        self.molecules = new_molecules
        self.get_water_dimensions()

    def distance(self, position1, position2, system_dimensions: tuple):
        normal_distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(position1, position2)))

        pbc_distance = []
        for (dim_lo, dim_hi), (a, b) in zip(system_dimensions, zip(position1, position2)):
            length = dim_hi - dim_lo
            delta = b - a
            delta -= round(delta / length) * length
            pbc_distance.append(delta ** 2)

        pbc_distance = math.sqrt(sum(pbc_distance))
        return min(normal_distance, pbc_distance)
    
    def remove_overlaying_water(self, other_molecules: list, system_dimensions: tuple):
        other_atoms_info = [
            (atom.position, self.vdw_radii.get(atom.element, self.vdw_def_radius) * self.vdw_def_scale)
            for molecule in other_molecules for atom in molecule.atoms
        ]

        filtered_molecules = []
        for molecule in self.molecules:
            keep_molecule = True
            for atom in molecule.atoms:
                atom_radius = self.vdw_radii.get(atom.element, self.vdw_def_radius) * self.vdw_def_scale
                if any(self.distance(atom.position, other_position, system_dimensions) < (atom_radius + other_radius)
                    for other_position, other_radius in other_atoms_info):
                    keep_molecule = False
                    break

            if keep_molecule:
                filtered_molecules.append(molecule)

        self.molecules = filtered_molecules
        self.get_water_dimensions()





    

    


    
