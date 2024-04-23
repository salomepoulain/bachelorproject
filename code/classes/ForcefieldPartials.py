
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

class FF_equivalence:
    def __init__(self) -> None:
        self.ff_atom: FF_atom = None
        self.nonbond_coef: FF_atom = None
        self.bond_coef: FF_atom = None
        self.angle_coef: FF_atom = None
        self.torsion_coef: FF_atom  = None
        self.improper_coef: FF_atom = None

class FF_nonbond_coef:
    def __init__(self) -> None:
        self.id: int = None
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