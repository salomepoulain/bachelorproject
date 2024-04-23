from SystemAllocator import SystemAllocator
from ForcefieldPartials import FF_atom, FF_nonbond_coef, FF_bond_coef, FF_angle_coef, FF_torsion_coef, FF_improper_coef
from SystemPartials import Bond, Angle, Torsion, Improper
from ForceFieldLoader import ForceFieldLoader
from Atom import Atom
from Unitcell import Unitcell
from typing import List

class SystemWriter:
    def __init__(self, SystemAllocator):
        self.dimensions = SystemAllocator.dimensions 
        self.atoms  = SystemAllocator.atoms
        self.ff_atoms: List[FF_atom] = SystemAllocator.ff_atoms
       
        # atom_type, mass
        self.masses: List[FF_atom.id, FF_atom.mass] = []
        # atom_type1, atom_type2, epsilon, sigma
        self.pair_coeffs: List[FF_nonbond_coef.ff_atoms[0].id, FF_nonbond_coef.ff_atoms[1].id, FF_nonbond_coef.epsilon, FF_nonbond_coef.sigma] = []
        # bond_type, k, r0
        self.bond_coeffs: List[FF_bond_coef.id, FF_bond_coef.k, FF_bond_coef.r0] = []
        # angle_type, k, theta0
        self.angle_coeffs: List[FF_angle_coef.id, FF_angle_coef.k, FF_angle_coef.theta0] = []
        # torsion_type, kphi, n, phi0
        self.torsion_coeffs: List[FF_torsion_coef.id, FF_torsion_coef.kphi, FF_torsion_coef.n, FF_torsion_coef.phi0] = []
        # improper_type, kchi, n, chi0
        self.improper_coeffs: List[FF_improper_coef.id, FF_improper_coef.kchi, FF_improper_coef.n, FF_improper_coef.chi0] = []

        # atom_id, atom_type, x, y, z
        self.atoms: List[Atom.id, Atom.ff_atom.id, Atom.position[0], Atom.position[1], Atom.position[2]] = []
        # bond_id, bond_type, atom1, atom2
        self.bonds: List[Bond.id, Bond.ff_bond_coef.id, Bond.atoms[0].id, Bond.atoms[1].id] = []
        # angle_id, angle_type, atom1, atom2, atom3
        self.angles: List[Angle.id, Angle.ff_angle_coef.id, Angle.atoms[0].id, Angle.atoms[1].id, Angle.atoms[2].id] = []
        # torsion_id, torsion_type, atom1, atom2, atom3, atom4
        self.torsions: List[Torsion.id, Torsion.ff_torsion_coef.id, Torsion.atoms[0].id, Torsion.atoms[1].id, Torsion.atoms[2].id, Torsion.atoms[3].id] = []
        # improper_id, improper_type, atom1, atom2, atom3, atom4
        self.impropers: List[Improper.id, Improper.ff_improper_coef.id, Improper.atoms[0].id, Improper.atoms[1].id, Improper.atoms[2].id, Improper.atoms[3].id] = []


    def give_atom_id(self):
        for atom in self.atoms:
            atom.id = self.atoms.index(atom)

    def give_ff_atom_type_id(self):
        used_ff_atom_types = []
        for atom in self.atoms:
            if atom.ff_atom.type not in used_ff_atom_types:
                used_ff_atom_types.append(atom.ff_atom.type)
                atom.ff_atom.id = used_ff_atom_types.index(atom.ff_atom.type)
    
    def give_bond_id(self):
        for bond in self.bonds:
            bond.id = self.bonds.index(bond)
    
    def give_angle_id(self):
        for angle in self.angles:
            angle.id = self.angles.index(angle)
        
    def give_torsion_id(self):
        for torsion in self.torsions:
            torsion.id = self.torsions.index(torsion)

    def give_improper_id(self):
        for improper in self.impropers:
            improper.id = self.impropers.index(improper)

    def write_data_file(self):
        pass

    
    






    