from code.classes.SystemAllocator import SystemAllocator
from code.classes.ForcefieldPartials import FF_atom, FF_nonbond_coef, FF_bond_coef, FF_angle_coef, FF_torsion_coef, FF_improper_coef
from code.classes.SystemPartials import Bond, Angle, Torsion, Improper
from code.classes.UnitcellPartials import Atom, Molecule
from code.banner_writer import print_banner
from typing import List

class SystemWriter:
    def __init__(self, file_name, replication, height, al_mg_ratio, ca_si_ratio, ff_params):

        print_banner('misc/banner.txt')

        # Load system in from previous classes
        self.file_name = file_name
        self.replication = replication
        self.system = SystemAllocator(file_name, replication, height, al_mg_ratio, ca_si_ratio, ff_params)

        self.dimensions = self.system.dimensions 
        self.molecules = self.system.molecules
        self.atoms = [atom for molecule in self.molecules for atom in molecule.atoms]

        self.pair_coeffs = self.system.pair_coeffs
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

        self.data_file = None

        # Add everything

        self.give_all_ids()
        self.store_all()
        self.write_data_file()

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
        for atom in self.atoms:
            atom.id = self.atoms.index(atom) + 1

    def give_bond_coeff_id(self):
        if len(self.bond_coeffs) == 0:
            return
        for bond_coeff in self.bond_coeffs:
            bond_coeff.id = self.bond_coeffs.index(bond_coeff) + 1

    def give_angle_coeff_id(self):
        if len(self.angle_coeffs) == 0:
            return
        for angle_coeff in self.angle_coeffs:
            angle_coeff.id = self.angle_coeffs.index(angle_coeff) + 1

    def give_torsion_coeff_id(self):
        if len(self.torsion_coeffs) == 0:
            return
        for torsion_coeff in self.torsion_coeffs:
            torsion_coeff.id = self.torsion_coeffs.index(torsion_coeff) + 1

    def give_improper_coeff_id(self):
        if len(self.improper_coeffs) == 0:
            return
        for improper_coeff in self.improper_coeffs:
            improper_coeff.id = self.improper_coeffs.index(improper_coeff) + 1

    def give_molecule_id(self):
        if len(self.molecules) == 0:
            return
        for molecule in self.molecules:
            molecule.id = self.molecules.index(molecule) + 1

    def give_ff_atom_type_id(self):
        sorted_atoms = sorted(self.ff_atoms, key=lambda atom: atom.mass)
        for index, ff_atom in enumerate(sorted_atoms, start=1):
            ff_atom.id = index

    def give_bond_id(self):
        if len(self.bonds) == 0:
            return
        for bond in self.bonds:
            bond.id = self.bonds.index(bond) + 1
    
    def give_angle_id(self):
        if len(self.angles) == 0:
            return
        for angle in self.angles:
            angle.id = self.angles.index(angle) + 1
        
    def give_torsion_id(self):
        if len(self.torsions) == 0:
            return
        for torsion in self.torsions:
            torsion.id = self.torsions.index(torsion) + 1

    def give_improper_id(self):
        if len(self.impropers) == 0:
            return
        for improper in self.impropers:
            improper.id = self.impropers.index(improper) + 1

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

    def title_line_writer(self, item):
        attribute = getattr(self, item, [])
        amount = len(attribute)
        return f"{amount} {item}"

    def store_description(self):
        dimensions = " " + str(self.replication[0]) + " by " + str(self.replication[1]) + " replicated"
        self.stored_description.append("Generated DATA file using OOP: " + self.file_name + ".xyz" + dimensions + "\n")

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
        if len(self.pair_coeffs) == 0:
            return
        
        sorted_pair_coeffs = sorted(self.pair_coeffs, key=lambda x: x.ff_atoms.id)

        self.stored_pair_coeffs.append("Pair Coeffs\n")
        for nonbond_coef in sorted_pair_coeffs:
            self.stored_pair_coeffs.append((nonbond_coef.ff_atoms.id, nonbond_coef.epsilon, nonbond_coef.sigma))

    def store_bond_coeffs(self):
        if len(self.bond_coeffs) == 0:
            return
        self.stored_bond_coeffs.append("Bond Coeffs\n")
        for bond_coef in self.bond_coeffs:
            self.stored_bond_coeffs.append((bond_coef.id, bond_coef.k, bond_coef.r0))

    def store_angle_coeffs(self):
        if len(self.angle_coeffs) == 0:
            return
        self.stored_angle_coeffs.append("Angle Coeffs\n")
        for angle_coef in self.angle_coeffs:
            self.stored_angle_coeffs.append((angle_coef.id, angle_coef.k, angle_coef.theta0))

    def store_dihedral_coeffs(self):
        if len(self.torsion_coeffs) == 0:
            return
        self.stored_torsion_coeffs.append("Dihedral Coeffs\n")
        for torsion_coef in self.torsion_coeffs:
            self.stored_torsion_coeffs.append((torsion_coef.id, torsion_coef.kphi, torsion_coef.n, torsion_coef.phi0))

    def store_improper_coeffs(self):
        if len(self.improper_coeffs) == 0:
            return
        self.stored_improper_coeffs.append("Improper Coeffs\n")
        for improper_coef in self.improper_coeffs:
            self.stored_improper_coeffs.append((improper_coef.id, improper_coef.kchi, improper_coef.n, improper_coef.chi0))

    def store_atoms(self):
        if len(self.atoms) == 0:
            return
        self.stored_atoms.append("Atoms\n")
        for molecule in self.molecules:
            for atom in molecule.atoms:
                self.stored_atoms.append((atom.id, molecule.id, atom.ff_atom.id, atom.ff_atom.charge, atom.position[0], atom.position[1], atom.position[2]))

    def store_bonds(self):
        if len(self.bonds) == 0:
            return
        self.stored_bonds.append("Bonds\n")
        for bond in self.bonds:
            self.stored_bonds.append((bond.id, bond.ff_bond_coef.id, bond.atoms[0].id, bond.atoms[1].id))
    
    def store_angles(self):
        if len(self.angles) == 0:
            return
        self.stored_angles.append("Angles\n")
        for angle in self.angles:
            self.stored_angles.append((angle.id, angle.ff_angle_coef.id, angle.atoms[0].id, angle.atoms[1].id, angle.atoms[2].id))

    def store_dihedrals(self):
        if len(self.torsions) == 0:
            return
        self.stored_torsions.append("Dihedrals\n")
        for torsion in self.torsions:
            self.stored_torsions.append((torsion.id, torsion.ff_torsion_coef.id, torsion.atoms[0].id, torsion.atoms[1].id, torsion.atoms[2].id, torsion.atoms[3].id))

    def store_impropers(self):
        if len(self.impropers) == 0:
            return
        self.stored_impropers.append("Impropers\n")
        for improper in self.impropers:
            self.impropers.append((improper.id, improper.ff_improper_coef.id, improper.atoms[0].id, improper.atoms[1].id, improper.atoms[2].id, improper.atoms[3].id))

    def write_section(self, section):
        if len(section) == 0:
            return

        # Function to calculate width of each column
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

        # Get column widths
        column_widths = calculate_column_widths(section)

        # Write section to file with formatted numbers
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
        self.data_file = open("output/" + self.file_name + str(self.replication) + ".data", "w")
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
        self.data_file.close()

    





    
    






    