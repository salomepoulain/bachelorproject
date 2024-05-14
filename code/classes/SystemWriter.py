from code.classes.SystemParts import Bond, Angle, Torsion, Improper
from code.classes.VerticalDuplicator import VerticalDuplicator
from typing import List

class SystemWriter(VerticalDuplicator):
    def __init__(self, settings):
        super().__init__(settings)
        self.s = settings

        # Load system in from previous classes
        self.atoms = [atom for molecule in self.molecules for atom in molecule.atoms]

        self.ff_atoms = list(set([atom.ff_atom for molecule in self.molecules for atom in molecule.atoms]))

        self.bond_coeffs = list(set([bond.ff_bond_coef for molecule in self.molecules for bond in molecule.bonds]))
        self.angle_coeffs = list(set([angle.ff_angle_coef for molecule in self.molecules for angle in molecule.angles]))
        self.torsion_coeffs = list(set([torsion.ff_torsion_coef for molecule in self.molecules for torsion in molecule.torsions]))
        self.improper_coeffs = list(set([improper.ff_improper_coef for molecule in self.molecules for improper in molecule.impropers]))

        self.bonds: List[Bond] = [bond for molecule in self.molecules for bond in molecule.bonds]
        self.angles: List[Angle] = [angle for molecule in self.molecules for angle in molecule.angles]
        self.torsions: List[Torsion] = [torsion for molecule in self.molecules for torsion in molecule.torsions]
        self.impropers: List[Improper] = [improper for molecule in self.molecules for improper in molecule.impropers]

        # Store all in correct order
        self.stored_description = []
        self.stored_masses = []
        self.stored_pair_coeffs = []
        self.stored_bond_coeffs = []
        self.stored_angle_coeffs = []
        self.stored_torsion_coeffs = []
        self.stored_improper_coeffs = []

        self.stored_atoms = []
        self.stored_bonds = []
        self.stored_angles = []
        self.stored_torsions = []
        self.stored_impropers = []
        
        self.stored_footer = []

        self.data_file = None

        # Add and write all data to file
        self.give_all_ids()
        self.store_all()
        self.write_data_file()

        # Write out groups to clay_top.dat
        # self.write_out_groups()

    def give_all_ids(self):
        self.give_atom_id()
        self.give_bond_coeff_id()
        self.give_angle_coeff_id()
        self.give_torsion_coeff_id()
        self.give_improper_coeff_id()
        self.give_molecule_id()
        self.give_ff_atom_type_id()
        self.give_bond_id()
        self.give_angle_id()
        self.give_torsion_id()
        self.give_improper_id()


    def give_atom_id(self):
        atom_id_counter = 1
        sorted_molecules = sorted(self.molecules, key=lambda molecule: molecule.type)

        for molecule in sorted_molecules:
            sorted_atoms = sorted(molecule.atoms, key=lambda atom: (atom.element, atom.position))

            for atom in sorted_atoms:
                atom.id = atom_id_counter
                atom_id_counter += 1


    def give_bond_coeff_id(self):
        if len(self.bond_coeffs) == 0:
            return
        self.bond_coeffs.sort(key=lambda x: (x.ff_atoms[0].type, x.ff_atoms[1].type))
        for index, bond_coeff in enumerate(self.bond_coeffs):
            bond_coeff.id = index + 1

    def give_angle_coeff_id(self):
        if len(self.angle_coeffs) == 0:
            return
        self.angle_coeffs.sort(key=lambda x: (x.ff_atoms[0].type, x.ff_atoms[1].type, x.ff_atoms[2].type))
        for index, angle_coeff in enumerate(self.angle_coeffs):
            angle_coeff.id = index + 1

    def give_torsion_coeff_id(self):
        if len(self.torsion_coeffs) == 0:
            return
        self.torsion_coeffs.sort(key=lambda x: (x.ff_atoms[0].type, x.ff_atoms[1].type, x.ff_atoms[2].type, x.ff_atoms[3].type))
        for index, torsion_coeff in enumerate(self.torsion_coeffs):
            torsion_coeff.id = index + 1

    def give_improper_coeff_id(self):
        if len(self.improper_coeffs) == 0:
            return
        self.improper_coeffs.sort(key=lambda x: (x.ff_atoms[0].type, x.ff_atoms[1].type, x.ff_atoms[2].type, x.ff_atoms[3].type))
        for index, improper_coeff in enumerate(self.improper_coeffs):
            improper_coeff.id = index + 1

    def give_molecule_id(self):
        if len(self.molecules) == 0:
            return
        self.molecules.sort(key=lambda molecule: min(atom.id for atom in molecule.atoms))
        for index, molecule in enumerate(self.molecules):
            molecule.id = index + 1

    def give_ff_atom_type_id(self):
        sorted_ff_atoms = sorted(self.ff_atoms, key=lambda x: (x.mass, x.description))
        for index, ff_atom in enumerate(sorted_ff_atoms):
            ff_atom.id = index + 1

    def give_bond_id(self):
        if len(self.bonds) == 0:
            return
        self.bonds.sort(key=lambda x: (x.ff_bond_coef.id, (min(x.atoms[0].id, x.atoms[1].id), max(x.atoms[0].id, x.atoms[1].id))))
        for index, bond in enumerate(self.bonds):
            bond.id = index + 1
    
    def give_angle_id(self):
        if len(self.angles) == 0:
            return
        self.angles.sort(key=lambda x: (x.ff_angle_coef.id, (x.atoms[0].id, x.atoms[1].id, x.atoms[2].id)))
        for index, angle in enumerate(self.angles):
            angle.id = index + 1
        
    def give_torsion_id(self):
        if len(self.torsions) == 0:
            return
        self.torsions.sort(key=lambda x: (x.ff_torsion_coef.id, (x.atoms[0].id, x.atoms[1].id, x.atoms[2].id, x.atoms[3].id)))
        for index, torsion in enumerate(self.torsions):
            torsion.id = index + 1

    def give_improper_id(self):
        if len(self.impropers) == 0:
            return
        self.impropers.sort(key=lambda x: (x.ff_improper_coef.id, (x.atoms[0].id, x.atoms[1].id, x.atoms[2].id, x.atoms[3].id)))
        for index, improper in enumerate(self.impropers):
            improper.id = index + 1

    def store_all(self):
        self.store_description()
        self.store_masses()
        self.store_pair_coeffs()
        self.store_bond_coeffs()
        self.store_angle_coeffs()
        self.store_dihedral_coeffs()
        self.store_improper_coeffs()
        self.store_atoms()
        self.store_bonds()
        self.store_angles()
        self.store_dihedrals()
        self.store_impropers()
        self.store_footer()
    
    def generate_footer(self):
        extra_spaces = 4 
        labels = [
            "replication", 
            "al_mg_ratio",
            "net_charge",
            "water_per_ion",
            "ff_atom_types",
            "water_distance",
            "mg_cutoff",
            "h_cutoff",
            "bond_cutoff",
            "input_file",
            "ff_files",
            "water_file"
        ] 
            # "Input file", 
            # "Replication", 
            # "Al/Mg Ratio",
            # "Ca/Si Ratio", 
            # "Water file used", 
            # "Water distance used",
            # "Water duplication", 
            # "Water per ion",
            # "Force fields used", 
            # "Mg cutoff distance", 
            # "H cutoff distance",
            # "Bond cutoff distance", 
            # "Force field parameters"

        max_label_length = max(len(label) for label in labels) + extra_spaces
        footer = [
            f"# {'replication:'.ljust(max_label_length)} {self.s.replication[0]}x{self.s.replication[1]}x{self.s.replication[2]}",
            f"# {'al_mg_ratio:'.ljust(max_label_length)} {self.s.al_mg_ratio}",
            f"# {'net_charge:'.ljust(max_label_length)} {self.s.net_charge}",
            f"# {'water_per_ion:'.ljust(max_label_length)} {self.write_water_per_ion()}",
            f"# {'water_distance:'.ljust(max_label_length)} {self.s.water_distance} Angstroms",
            f"# {'mg_cutoff:'.ljust(max_label_length)} {self.s.mg_cutoff} Angstroms",
            f"# {'h_cutoff:'.ljust(max_label_length)} {self.s.h_cutoff} Angstroms",
            f"# {'bond_cutoff:'.ljust(max_label_length)} {self.s.bond_cutoff} Angstroms",
            f"# {'input_file:'.ljust(max_label_length)} {self.s.input_file}",
            f"# {'ff_files:'.ljust(max_label_length)} {', '.join(self.s.ff_files)}",
            f"# {'water_file:'.ljust(max_label_length)} {self.s.water_file}",
            "# Force field parameters:"
        ]

        params_start = max_label_length + 2

        for param in self.s.ff_atom_types:
            footer.append(f"#{' ':<{params_start}}- {param}")

        footer_length = 0
        with open('misc/banner.txt', 'r') as file:
            lines = file.readlines()
            for line in lines:
                footer.append("# " + line.rstrip('\n'))
                footer_length = len(line) if len(line) > footer_length else footer_length

        dashed_line = "# " + "-" * (footer_length - 1)
        footer.insert(0, dashed_line)
        footer.append(dashed_line)

        return "\n".join(footer)

    def write_water_per_ion(self):
        if self.s.water_per_ion is None or self.s.water_per_ion > 23.125:
            return "Default by szscerba: 23.125"

        return self.s.water_per_ion

    def store_footer(self):
        self.stored_footer.append(self.generate_footer())

    def title_line_writer(self, item):
        attribute = getattr(self, item, [])
        amount = len(attribute)
        return f"{amount} {item}"

    def store_description(self):
        self.stored_description.append("Generated DATA file using 'LAMMPS MT DATA FILE GENERATOR'")

        items = ["atoms", "bonds", "angles", "dihedrals", "impropers"]
        for item in items:
            line = self.title_line_writer(item)
            self.stored_description.append(line)

        self.stored_description.append(f"\n{len(self.ff_atoms)} atom types")
        self.stored_description.append(f"{len(self.bond_coeffs)} bond types")
        self.stored_description.append(f"{len(self.angle_coeffs)} angle types")
        self.stored_description.append(f"{len(self.torsion_coeffs)} dihedral types")
        self.stored_description.append(f"{len(self.improper_coeffs)} improper types\n")

        width = 20

        self.stored_description.append(f"{self.dimensions[0][0]:<{width}} {self.dimensions[0][1]:<{width}} xlo xhi")
        self.stored_description.append(f"{self.dimensions[1][0]:<{width}} {self.dimensions[1][1]:<{width}} ylo yhi")
        self.stored_description.append(f"{self.dimensions[2][0]:<{width}} {self.dimensions[2][1]:<{width}} zlo zhi\n")

    def store_masses(self):
        self.stored_masses.append("Masses\n")
        sorted_atoms = sorted(self.ff_atoms, key=lambda atom: atom.id)
        for ff_atom in sorted_atoms:
            description = "# " + ff_atom.description
            self.stored_masses.append((ff_atom.id, ff_atom.mass, description))

    def store_pair_coeffs(self):
        self.stored_pair_coeffs.append("Pair Coeffs\n")
        sorted_atoms = sorted(self.ff_atoms, key=lambda atom: atom.id)
        for ff_atom in sorted_atoms:
            self.stored_pair_coeffs.append((ff_atom.id, ff_atom.ff_pair_coef.epsilon, ff_atom.ff_pair_coef.sigma))

    def store_bond_coeffs(self):
        if len(self.bond_coeffs) == 0:
            return
        self.stored_bond_coeffs.append("Bond Coeffs\n")
        sorted_bond_coeffs = sorted(self.bond_coeffs, key=lambda x: x.id)
        for bond_coef in sorted_bond_coeffs:
            self.stored_bond_coeffs.append((bond_coef.id, bond_coef.k, bond_coef.r0))

    def store_angle_coeffs(self):
        if len(self.angle_coeffs) == 0:
            return
        self.stored_angle_coeffs.append("Angle Coeffs\n")
        sorted_angle_coeffs = sorted(self.angle_coeffs, key=lambda x: x.id)
        for angle_coef in sorted_angle_coeffs:
            self.stored_angle_coeffs.append((angle_coef.id, angle_coef.k, angle_coef.theta0))

    def store_dihedral_coeffs(self):
        if len(self.torsion_coeffs) == 0:
            return
        self.stored_torsion_coeffs.append("Dihedral Coeffs\n")
        sorted_torsion_coeffs = sorted(self.torsion_coeffs, key=lambda x: x.id)
        for torsion_coef in sorted_torsion_coeffs:
            self.stored_torsion_coeffs.append((torsion_coef.id, torsion_coef.kphi, torsion_coef.n, torsion_coef.phi0))

    def store_improper_coeffs(self):
        if len(self.improper_coeffs) == 0:
            return
        self.stored_improper_coeffs.append("Improper Coeffs\n")
        sorted_improper_coeffs = sorted(self.improper_coeffs, key=lambda x: x.id)
        for improper_coef in sorted_improper_coeffs:
            self.stored_improper_coeffs.append((improper_coef.id, improper_coef.kchi, improper_coef.n, improper_coef.chi0))

    def store_atoms(self):
        if len(self.atoms) == 0:
            return
        self.stored_atoms.append("Atoms\n")
        for molecule in self.molecules:
            sorted_atoms = sorted(molecule.atoms, key=lambda atom: atom.id)
            for atom in sorted_atoms:
                self.stored_atoms.append((
                    atom.id,
                    molecule.id,
                    atom.ff_atom.id,
                    atom.ff_atom.charge,
                    atom.position[0],
                    atom.position[1],
                    atom.position[2]
                ))

    def store_bonds(self):
        if len(self.bonds) == 0:
            return
        self.stored_bonds.append("Bonds\n")
        sorted_bonds = sorted(self.bonds, key=lambda x: x.id)
        for bond in sorted_bonds:
            self.stored_bonds.append((bond.id, bond.ff_bond_coef.id, bond.atoms[0].id, bond.atoms[1].id))
    
    def store_angles(self):
        if len(self.angles) == 0:
            return
        self.stored_angles.append("Angles\n")
        sorted_angles = sorted(self.angles, key=lambda x: x.id)
        for angle in sorted_angles:
            self.stored_angles.append((angle.id, angle.ff_angle_coef.id, angle.atoms[0].id, angle.atoms[1].id, angle.atoms[2].id))

    def store_dihedrals(self):
        if len(self.torsions) == 0:
            return
        self.stored_torsions.append("Dihedrals\n")
        sorted_dihedrals = sorted(self.torsions, key=lambda x: x.id)
        for torsion in sorted_dihedrals:
            self.stored_torsions.append((torsion.id, torsion.ff_torsion_coef.id, torsion.atoms[0].id, torsion.atoms[1].id, torsion.atoms[2].id, torsion.atoms[3].id))

    def store_impropers(self):
        if len(self.impropers) == 0:
            return
        self.stored_impropers.append("Impropers\n")
        sorted_impropers = sorted(self.impropers, key=lambda x: x.id)
        for improper in sorted_impropers:
            self.stored_impropers.append((improper.id, improper.ff_improper_coef.id, improper.atoms[0].id, improper.atoms[1].id, improper.atoms[2].id, improper.atoms[3].id))

    def write_section(self, section):
        if len(section) == 0:
            return

        def calculate_column_widths(section):
            column_widths = [0] * len(section[0])
            for line in section:
                for i, item in enumerate(line):
                    if isinstance(item, int) or isinstance(item, float):
                        if i >= len(column_widths):
                            column_widths.append(len(str(item)))
                        else:
                            column_widths[i] = max(column_widths[i], len(str(item)))
            return column_widths

        column_widths = calculate_column_widths(section)

        for line in section:
            if all(isinstance(item, str) for item in line):
                self.data_file.write(line + "\n")
            else:
                formatted_line = ""
                for i, item in enumerate(line):
                    if isinstance(item, str):
                        formatted_line += item
                    else:
                        formatted_line += "{:<{width}}".format(item, width=column_widths[i])
                    if i < len(line) - 1:
                        formatted_line += "\t"
                self.data_file.write(formatted_line + "\n")

        self.data_file.write("\n")

    def write_data_file(self):
        output_name = self.s.output_file if self.s.output_file is not None else self.s.input_file + str(self.s.replication)

        self.data_file = open("output/" + output_name + ".data", "w")
        self.write_section(self.stored_description)
        
        self.write_section(self.stored_masses)
        
        self.write_section(self.stored_pair_coeffs)
        self.write_section(self.stored_bond_coeffs)
        self.write_section(self.stored_angle_coeffs)
        self.write_section(self.stored_torsion_coeffs)
        self.write_section(self.stored_improper_coeffs)

        self.write_section(self.stored_atoms)
        self.write_section(self.stored_bonds)
        self.write_section(self.stored_angles)
        self.write_section(self.stored_torsions)
        self.write_section(self.stored_impropers)

        self.write_section(self.stored_footer)
        self.data_file.close()

    def write_out_groups(self):
        with open('clay_top.dat', 'a') as file:
            current_line_length = 0  
            for molecule in self.molecules:
                if molecule.type == 'clay':
                    for atom in molecule.atoms:
                        if atom.element == 'O':
                            if 6.069 < atom.position[2] < 10:
                                output = str(atom.id) + ' '
                                output_length = len(output)
                                
                                if current_line_length + output_length > 60:
                                    file.write('&'+ '\n')
                                    current_line_length = 0  

                                file.write(output)
                                current_line_length += output_length

        with open('clay_bottom.dat', 'a') as file:
            current_line_length = 0  
            for molecule in self.molecules:
                if molecule.type == 'clay':
                    for atom in molecule.atoms:
                        if atom.element == 'O':
                            if -2 < atom.position[2] < 0.66:
                                output = str(atom.id) + ' '
                                output_length = len(output)
                                
                                if current_line_length + output_length > 60:
                                    file.write('&'+ '\n')
                                    current_line_length = 0  

                                file.write(output)
                                current_line_length += output_length

        with open('clay_left.dat', 'a') as file:
            current_line_length = 0  
            for molecule in self.molecules:
                if molecule.type == 'clay':
                    for atom in molecule.atoms:
                        if atom.element == 'O':
                            if atom.position[2] < 10:
                                if 8.156 < atom.position[1] < 8.853:
                                    output = str(atom.id) + ' '
                                    output_length = len(output)
                                    
                                    if current_line_length + output_length > 60:
                                        file.write('&' + '\n')
                                        current_line_length = 0  

                                    file.write(output)
                                    current_line_length += output_length

        with open('clay_right.dat', 'a') as file:
            current_line_length = 0  
            for molecule in self.molecules:
                if molecule.type == 'clay':
                    for atom in molecule.atoms:
                        if atom.element == 'O':
                            if atom.position[2] < 10:
                                if atom.position[1] > 52.320:
                                    output = str(atom.id) + ' '
                                    output_length = len(output)
                                    
                                    if current_line_length + output_length > 60:
                                        file.write('&'+ '\n')
                                        current_line_length = 0  

                                    file.write(output)
                                    current_line_length += output_length