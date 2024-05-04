class Molecule:
    def __init__(self):   
        self.id = None
        self.type = None
        self.atoms = []
        
        self.bonds = []
        self.angles = []
        self.torsions = []
        self.impropers = []

    def add_atom(self, element: str, position: tuple):
        atom = Atom()
        atom.element = element
        atom.position = position
        self.atoms.append(atom)

    def add_bond(self, atom1, atom2, bond_coef):
        bond = Bond()
        bond.atoms = [atom1, atom2]
        bond.ff_bond_coef = bond_coef
        self.bonds.append(bond)
        # print(f"Bond added between atom {atom1.element} and {atom2.element}.")
    
    def add_angle(self, atom1, atom2, atom3, angle_coef):
        angle = Angle()
        angle.atoms = [atom1, atom2, atom3]
        angle.ff_angle_coef = angle_coef
        self.angles.append(angle)
    #     print(f"Angle added between atoms with IDs {atom1.id}, {atom2.id}, and {atom3.id}.")

    # def add_torsion(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom, torsion_coef: FF_torsion_coef = None):
    #     torsion = Torsion()
    #     torsion.atoms = [atom1, atom2, atom3, atom4]
    #     torsion.ff_torsion_coef = torsion_coef
    #     self.torsions.append(torsion)
    #     print(f"Torsion added between atoms with IDs {atom1.id}, {atom2.id}, {atom3.id}, and {atom4.id}.")

    # def add_improper(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom, improper_coef: FF_improper_coef = None):
    #     improper = Improper()
    #     improper.atoms = [atom1, atom2, atom3, atom4]
    #     improper.ff_improper_coef = improper_coef
    #     self.impropers.append(improper)
    #     print(f"Improper added between atoms with IDs {atom1.id}, {atom2.id}, {atom3.id}, and {atom4.id}.")

    def return_atom_dict(self, bond_cutoff):
        atom_dict = {}
        for atom in self.atoms:
            atom_type_key = atom.ff_atom.type if atom.ff_atom else None
            if atom_type_key:
                if atom_type_key not in atom_dict:
                    atom_dict[atom_type_key] = []
                atom_dict[atom_type_key].append(atom)
    
        return atom_dict
    
    def is_bonded(self, atom1, atom2):
        """ Helper method to determine if two atoms are bonded in a molecule """
        return any(set([atom1, atom2]) == set(bond.atoms) for bond in self.bonds)
    
    def has_angle(self, atom1, atom2, atom3):
        normal_orientation = [atom1, atom2, atom3]
        reversed_orientation = [atom3, atom2, atom1]

        return any(
            (angle.atoms == normal_orientation or angle.atoms == reversed_orientation)
            for angle in self.angles
        )


class Atom:
    def __init__(self):
        self.id = None
        self.element = None
        self.position = (None, None, None)
        self.ff_atom = None


class Bond:
    def __init__(self) -> None:
        self.id: int = None
        self.ff_bond_coef: FF_bond_coef = None
        self.atoms: list[Atom, Atom] = None

    def other_atom(self, current_atom):
        if current_atom == self.atoms[0]:
            return self.atoms[1]
        elif current_atom == self.atoms[1]:
            return self.atoms[0]
        else:
            return None

class Angle:
    def __init__(self) -> None:
        self.id: int = None
        self.ff_angle_coef: FF_angle_coef = None
        self.atoms: list[Atom, Atom, Atom] = None

class Torsion:
    def __init__(self) -> None:
        self.id: int = None
        self.ff_torsion_coef: FF_torsion_coef = None
        self.atoms: list[Atom, Atom, Atom, Atom] = None

class Improper:
    def __init__(self) -> None:
        self.id: int = None
        self.ff_improper_coef: FF_improper_coef = None
        self.atoms: list[Atom, Atom, Atom, Atom] = None



class FF_atom:
    def __init__(self) -> None:
        self.id: int = None
        self.type: str = None #uniqe identifier
        self.element: str = None
        self.description: str = None
        self.connections: int = None
        self.mass: float = None
        self.charge: float = None

    def __str__(self):
        return f"Type: {self.type}, Mass: {self.mass}, Element: {self.element}, " \
               f"Connections: {self.connections}, Description: {self.description}, Charge: {self.charge}"

class FF_nonbond_coef:
    def __init__(self) -> None:
        self.ff_atoms: FF_atom = None
        self.epsilon = None
        self.sigma = None

class FF_bond_coef:
    def __init__(self) -> None:
        self.id: int = None
        self.ff_atoms: list[FF_atom, FF_atom] = None
        self.k = None
        self.r0 = None

class FF_angle_coef:
    def __init__(self) -> None:
        self.id: int = None
        self.ff_atoms: list[FF_atom, FF_atom, FF_atom] = None
        self.k = None
        self.theta0 = None

class FF_torsion_coef:
    def __init__(self) -> None:
        self.id: int = None
        self.ff_atoms: list[FF_atom, FF_atom, FF_atom, FF_atom] = None
        self.kphi = None
        self.n = None
        self.phi0 = None
    
class FF_improper_coef:
    def __init__(self) -> None:
        self.id: int = None
        self.ff_atoms: list[FF_atom, FF_atom, FF_atom, FF_atom] = None
        self.kchi = None
        self.n = None     
        self.chi0 = None

class FF_equivalence:
    def __init__(self) -> None:
        self.ff_atom: FF_atom = None
        self.nonbond_coef: FF_atom = None
        self.bond_coef: FF_atom = None
        self.angle_coef: FF_atom = None
        self.torsion_coef: FF_atom  = None
        self.improper_coef: FF_atom = None