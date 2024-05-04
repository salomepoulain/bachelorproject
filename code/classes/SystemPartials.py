from code.classes.ForcefieldPartials import FF_bond_coef, FF_angle_coef, FF_torsion_coef, FF_improper_coef
from code.classes.UnitcellPartials import Atom

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