from ForcefieldPartials import FF_bond_coef, FF_angle_coef, FF_torsion_coef, FF_improper_coef
import Atom

class Bond:
    self.id: int = None
    self.bond_type: FF_bond_coef = None
    self.atoms: list[Atom, Atom] = None

class Angle:
    self.id: int = None
    self.angle_type: FF_angle_coef = None
    self.atoms: list[Atom, Atom, Atom] = None

class Torsion:
    self.id: int = None
    self.torsion_type: FF_torsion_coef = None
    self.atoms: list[Atom, Atom, Atom, Atom] = None

class Improper:
    self.id: int = None
    self.improper_type: FF_improper_coef = None
    self.atoms: list[Atom, Atom, Atom, Atom] = None