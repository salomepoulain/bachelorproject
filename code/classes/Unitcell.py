from code.classes.atom import Atom

class Unitcell:
    def __init__(self, file_name) -> None:
        self.file_name = file_name
        self.dimensions = (None, None, None)
        self.atoms = []
        self.read_xyz(file_name)


    def read_xyz(self, file_name):
        file_path = 'input/' + file_name + '.xyz'
        
        with open(file_path, 'r') as file:
            next(file)  # Skip the first line (atom count)
            next(file)  # Skip the second line (comment)
            for line in file:
                parts = line.strip().split()
                if len(parts) == 4:
                    element, x, y, z = parts
                    atom = Atom()
                    atom.element = element
                    atom.position = (float(x), float(y), float(z))
                    self.atoms.append(atom)
    