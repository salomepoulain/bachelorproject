from code.classes.SolventIonAdder import SolventIonAdder
from code.classes.SystemParts import Molecule, Atom

class VerticalDuplicator(SolventIonAdder):
    def __init__(self, settings):
        super().__init__(settings)
        self.s = settings

        self.replicate_z_axis()
        
    def replicate_z_axis(self):
        unit_z = self.dimensions[2][1] - self.dimensions[2][0]
        n_z = self.s.replication[2]
        original_molecules = list(self.molecules)

        for k in range(1, n_z):
            for molecule in original_molecules:
                new_molecule = Molecule()
                new_molecule.type = molecule.type
                
                atom_map = {}
                for atom in molecule.atoms:
                    x, y, z = atom.position
                    new_z = z + k * unit_z
                    new_atom = Atom()
                    new_atom.element = atom.element
                    new_atom.position = (x, y, new_z)
                    new_atom.ff_atom = atom.ff_atom
                    new_molecule.atoms.append(new_atom)
                    atom_map[atom] = new_atom

                for bond in molecule.bonds:
                    new_atom1 = atom_map[bond.atoms[0]]
                    new_atom2 = atom_map[bond.atoms[1]]
                    new_molecule.add_bond(new_atom1, new_atom2, bond.ff_bond_coef)

                for angle in molecule.angles:
                    new_atom1 = atom_map[angle.atoms[0]]
                    new_atom2 = atom_map[angle.atoms[1]]
                    new_atom3 = atom_map[angle.atoms[2]]
                    new_molecule.add_angle(new_atom1, new_atom2, new_atom3, angle.ff_angle_coef)

                for torsion in molecule.torsions:
                    new_atom1 = atom_map[torsion.atoms[0]]
                    new_atom2 = atom_map[torsion.atoms[1]]
                    new_atom3 = atom_map[torsion.atoms[2]]
                    new_atom4 = atom_map[torsion.atoms[3]]
                    new_molecule.add_torsion(new_atom1, new_atom2, new_atom3, new_atom4, torsion.ff_torsion_coef)

                for improper in molecule.impropers:
                    new_atom1 = atom_map[improper.atoms[0]]
                    new_atom2 = atom_map[improper.atoms[1]]
                    new_atom3 = atom_map[improper.atoms[2]]
                    new_atom4 = atom_map[improper.atoms[3]]
                    new_molecule.add_improper(new_atom1, new_atom2, new_atom3, new_atom4, improper.ff_improper_coef)

                self.molecules.append(new_molecule)

        self.get_unit_cell_dimensions()
        self.translate_simulation_height()
