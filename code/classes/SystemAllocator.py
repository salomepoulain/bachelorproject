from ForcefieldPartials import FF_atom, FF_nonbond_coef, FF_bond_coef, FF_angle_coef, FF_torsion_coef, FF_improper_coef
from ForceFieldLoader import ForceFieldLoader
from Atom import Atom
from typing import List
from Unitcell import Unitcell
from SystemPartials import Bond, Angle, Torsion, Improper



class SystemAllocator:
    def __init__(self, Unitcell, ForceFieldLoader):
        self.dimensions = Unitcell.dimensions
        self.atoms: List[Atom] = Unitcell.atoms

        self.ff_atoms: List[FF_atom] = ForceFieldLoader.ff_atoms

        self.bonds: List[Bond] = []
        self.angles: List[Angle] = []
        self.torsions: List[Torsion] = []
        self.impropers: List[Improper] = []


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


    def add_bond(self, Bond):
        pass
    
    def add_angle(self, Angle):
        pass

    def add_torsion(self, Torsion):
        pass

    def add_improper(self, Improper):
        pass

    # Not sure what order this should be done in, depends if above functions are partial manual or automatic
    def rescale_unitcell(self):
        pass
    

   

    
