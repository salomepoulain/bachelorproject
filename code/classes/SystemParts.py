from typing import List, Tuple, Union, Dict, Optional

"""
This module defines classes to model the fundamental components of a molecular dynamics simulation. 
These include molecules, atoms, and various force field parameters such as bonds, angles, torsions, and improper torsions. 
Each class encapsulates the properties necessary for defining molecular structures and interactions based on force fields.
"""

class Molecule:
    def __init__(self):   
        self.id: Optional[int] = None
        self.type: Optional[str] = None
        self.atoms: List['Atom'] = []
        
        self.bonds: List['Bond'] = []
        self.angles: List['Angle'] = []
        self.torsions: List['Torsion'] = []
        self.impropers: List['Improper'] = []

    def add_atom(self, element: str, position: Tuple[float, float, float]):
        atom = Atom()
        atom.element = element
        atom.position = position
        self.atoms.append(atom)

    def add_bond(self, atom1: 'Atom', atom2: 'Atom', bond_coef: 'FF_bond_coef'):
        bond = Bond()
        bond.atoms = [atom1, atom2]
        bond.ff_bond_coef = bond_coef
        self.bonds.append(bond)
    
    def add_angle(self, atom1: 'Atom', atom2: 'Atom', atom3: 'Atom', angle_coef: 'FF_angle_coef'):
        angle = Angle()
        angle.atoms = [atom1, atom2, atom3]
        angle.ff_angle_coef = angle_coef
        self.angles.append(angle)

    def add_torsion(self, atom1: 'Atom', atom2: 'Atom', atom3: 'Atom', atom4: 'Atom', torsion_coef: 'FF_torsion_coef'):
        torsion = Torsion()
        torsion.atoms = [atom1, atom2, atom3, atom4]
        torsion.ff_torsion_coef = torsion_coef
        self.torsions.append(torsion)

    def add_improper(self, atom1: 'Atom', atom2: 'Atom', atom3: 'Atom', atom4: 'Atom', improper_coef: 'FF_improper_coef'):
        improper = Improper()
        improper.atoms = [atom1, atom2, atom3, atom4]
        improper.ff_improper_coef = improper_coef
        self.impropers.append(improper)

    def return_atom_dict(self) -> Dict[str, List['Atom']]:
        atom_dict = {}
        for atom in self.atoms:
            atom_type_key = atom.ff_atom.type if atom.ff_atom else None
            if atom_type_key:
                if atom_type_key not in atom_dict:
                    atom_dict[atom_type_key] = []
                atom_dict[atom_type_key].append(atom)
    
        return atom_dict
    
    def is_bonded(self, atom1: 'Atom', atom2: 'Atom') -> bool:
        return any(set([atom1, atom2]) == set(bond.atoms) for bond in self.bonds)
    
    def has_angle(self, atom1: 'Atom', atom2: 'Atom', atom3: 'Atom') -> bool:
        normal_orientation = [atom1, atom2, atom3]
        reversed_orientation = [atom3, atom2, atom1]

        return any(
            (angle.atoms == normal_orientation or angle.atoms == reversed_orientation)
            for angle in self.angles
        )

class Atom:
    def __init__(self):
        self.id: Optional[int] = None
        self.element: Optional[str] = None
        self.position: Tuple[Optional[float], Optional[float], Optional[float]] = (None, None, None)
        self.ff_atom: Optional['FF_atom'] = None

class Bond:
    def __init__(self) -> None:
        self.id: Optional[int] = None
        self.ff_bond_coef: Optional['FF_bond_coef'] = None
        self.atoms: Optional[List['Atom']] = None

    def other_atom(self, current_atom: 'Atom') -> Optional['Atom']:
        if current_atom == self.atoms[0]:
            return self.atoms[1]
        elif current_atom == self.atoms[1]:
            return self.atoms[0]
        else:
            return None

class Angle:
    def __init__(self) -> None:
        self.id: Optional[int] = None
        self.ff_angle_coef: Optional['FF_angle_coef'] = None
        self.atoms: Optional[List['Atom']] = None

class Torsion:
    def __init__(self) -> None:
        self.id: Optional[int] = None
        self.ff_torsion_coef: Optional['FF_torsion_coef'] = None
        self.atoms: Optional[List['Atom']] = None

class Improper:
    def __init__(self) -> None:
        self.id: Optional[int] = None
        self.ff_improper_coef: Optional['FF_improper_coef'] = None
        self.atoms: Optional[List['Atom']] = None

class FF_atom:
    def __init__(self) -> None:
        self.id: Optional[int] = None
        self.type: Optional[str] = None # unique identifier
        self.element: Optional[str] = None
        self.description: Optional[str] = None
        self.connections: Optional[int] = None
        self.mass: Optional[float] = None
        self.charge: Optional[float] = None

        self.ff_pair_coef: Optional['FF_pair_coef'] = None

class FF_pair_coef:
    def __init__(self) -> None:
        self.ff_atom: Optional['FF_atom'] = None
        self.epsilon: Optional[float] = None
        self.sigma: Optional[float] = None

class FF_bond_coef:
    def __init__(self) -> None:
        self.id: Optional[int] = None
        self.ff_atoms: Optional[List['FF_atom']] = None
        self.k: Optional[float] = None
        self.r0: Optional[float] = None

class FF_angle_coef:
    def __init__(self) -> None:
        self.id: Optional[int] = None
        self.ff_atoms: Optional[List['FF_atom']] = None
        self.k: Optional[float] = None
        self.theta0: Optional[float] = None

class FF_torsion_coef:
    def __init__(self) -> None:
        self.id: Optional[int] = None
        self.ff_atoms: Optional[List['FF_atom']] = None
        self.kphi: Optional[float] = None
        self.n: Optional[int] = None
        self.phi0: Optional[float] = None
    
class FF_improper_coef:
    def __init__(self) -> None:
        self.id: Optional[int] = None
        self.ff_atoms: Optional[List['FF_atom']] = None
        self.kchi: Optional[float] = None
        self.n: Optional[int] = None     
        self.chi0: Optional[float] = None

class FF_equivalence:
    def __init__(self) -> None:
        self.ff_atom: Optional['FF_atom'] = None
        self.pair_coef: Optional['FF_atom'] = None
        self.bond_coef: Optional['FF_atom'] = None
        self.angle_coef: Optional['FF_atom'] = None
        self.torsion_coef: Optional['FF_atom'] = None
        self.improper_coef: Optional['FF_atom'] = None
