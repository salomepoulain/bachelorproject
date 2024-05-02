class Atom:
    def __init__(self):
        self.id = None
        self.element = None
        self.position = (None, None, None)
        self.ff_atom = None

class Molecule:
    def __init__(self):   
        self.id = None
        self.type = None
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.torsions = []
        self.impropers = []