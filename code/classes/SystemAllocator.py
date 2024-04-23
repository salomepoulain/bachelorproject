from ForcefieldPartials import FF_atom, FF_nonbond, FF_bond, FF_angle, FF_torsion, FF_improper
import ForceFieldLoader
import Atom
from typing import List
import ForceFieldLoader
import Unitcell

class SystemAllocator:
    def __init__(self, Unitcell, ForceFieldLoader):
        self.dimensions = Unitcell.dimensions
        self.atoms: List[Atom] = Unitcell.atoms

        self.ff_atoms: List[FF_atom] = ForceFieldLoader.ff_atoms

        self.pair_coeffs: = []
        self.bond_coeffs: = []
        self.angle_coeffs: = []
        self.torsion_coeffs: = []
        self.improper_coeffs: = []

        self.bonds: = []
        self.angles: = []
        self.torsions: = []
        self.impropers: = []


    # Manual function I think
    def give_atom_ff_type(self):
        pass
        # Here we will allocate the system based on the unitcell and forcefield partials
        # This will involve creating a system, and adding atoms, bonds, angles, torsions, and impropers
        # We will also need to calculate the nonbonded interactions

        for atom in self.atoms:
            for ff_atom in self.ff_atoms:
                if atom.element == ff_atom.element:
                    atom.ff_atom = ff_atom
                    break

    # Manual function I think
    def give_atom_bonds(self):
        for atom in self.atoms:
            for bond in self.bonds:
                if atom.ff_atom == bond.ff_atoms[0]:
                    pass
    
    # Manual function I think               
    def give_atom_angles(self):
        pass
        # Here we will allocate the angles to the atoms

    # Manual function I think
    def give_atom_torsions(self):
        pass
        # Here we will allocate the torsions to the atoms
    
    # Manual function I think
    def give_atom_impropers(self):
        pass
        # Here we will allocate the impropers to the atoms


    def duplicate_unitcell(self):
        pass
        # Watch out for id's!!!

        # Here we will duplicate the unitcell and store it in a new object
        # This will involve copying the atom positions and elements

   

    
