from ForcefieldPartials import FF_atom, FF_nonbond_coef, FF_bond_coef, FF_angle_coef, FF_torsion_coef, FF_improper_coef, FF_equivalence

"""
Force Field Loader Module for Molecular Dynamics Simulations

Description:
    This module provides classes and functions to load, manage, and process
    force field parameters from specified file paths within 'forcefields/' directory.
    It uses the FF_atom, FF_nonbond_coef, FF_bond, FF_angle, FF_torsion, FF_improper, and FF_equivalence classes 
    
    - loads cvff data from clayff and cvff forcefields in .frc files
    - also works with only one file
    - loads atom types, equivalences, nonbond_coef, bond_coef, angle_coef, torsion_coef, improper_coef parameters
    - nonbond_coef parameters are converted to give sigma and epsilon values
    - 2k values are devided by 2 to get k values
    - cvff equivalences are used to ensure all atom types are represented in the force field parameters
    - checks for duplicate atom types and appends '[dup]' to the type if a duplicate is found, 
        with the second file iteration. FF_atom.type are thus unique identifiers
    - for torsion_coefs, '*' is used to represent a wildcard atom type
"""

class ForceFieldLoader:
    def __init__(self) -> None:
        self.file_paths = ['forcefields/clayff.frc', 'forcefields/cvff.frc']
        self.ff_atoms: List[FF_atom] = []
        self.ff_equivalences: List[FF_equivalence] = []
        self.nonbond_coefs: List[FF_nonbond_coef] = []
        self.bond_coefs: List[FF_bond] = []
        self.angle_coefs: List[FF_angle] = []
        self.torsion_coefs: List[FF_torsion] = []
        self.improper_coefs: List[FF_improper] = []
        self.duplicates: Set[str] = set()

        self.load_all_forcefield_params()

    def load_all_forcefield_params(self):
        i = 0
        for file_path in self.file_paths:
            self.read_ff_atom_types(file_path, i)
            self.read_ff_equivalences(file_path, i)

            self.read_ff_nonbond_coef(file_path, i)
            self.add_nonbond_coef_equivalences()

            self.read_ff_bond(file_path, i)
            self.add_bond_equivalences()

            self.read_ff_angle(file_path, i)
            self.add_angle_equivalences()

            self.read_ff_torsion(file_path, i)
            self.add_torsion_equivalences()

            self.read_ff_improper(file_path, i)
            self.add_improper_equivalences()

            i += 1

    def find_ff_atom(self, str, i):
        if i > 0 and str in self.duplicates:
            str += '[dup]'

        if str == '*':
            ff_atom = FF_atom()
            ff_atom.type = '*'
            return ff_atom

        for ff_atom in self.ff_atoms:
            if ff_atom.type == str:
                return ff_atom
            
    def read_ff_atom_types(self, file_path, i):
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:

                if line.startswith('#') and start_processing == True:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>'):
                    continue
                if line.strip().startswith('#atom_types'):
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()
                
                    ff_atom = FF_atom()
                    ff_atom.type = parts[2]
                    ff_atom.mass = float(parts[3])
                    ff_atom.element = parts[4]

                    try:
                        ff_atom.connections = int(parts[5]) if parts[5] else None
                    except IndexError:  
                        ff_atom.connections = None
                    except ValueError:
                        ff_atom.connections = None

                    try:
                        charge_converted = float(parts[-1])
                        ff_atom.charge = charge_converted
                        ff_atom.description = ' '.join(parts[6:-1])
                    except ValueError:
                        ff_atom.description = ' '.join(parts[6:])
                        pass
                    
                    existing_types = {ff_atom.type for ff_atom in self.ff_atoms}
                    if ff_atom.type in existing_types:
                        self.duplicates.add(ff_atom.type)
                        ff_atom.type += '[dup]'
                    
                    self.ff_atoms.append(ff_atom)

    def read_ff_equivalences(self, file_path, i):
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#') and start_processing:
                    return
                if line.startswith('!') or line.strip() == '' or line.startswith('>'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#equivalence' and line_parts[1] == 'cvff':
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()

                    ff_equivalence = FF_equivalence()
                    ff_equivalence.ff_atom = self.find_ff_atom(parts[2], i)
                    ff_equivalence.nonbond_coef = self.find_ff_atom(parts[3], i)
                    ff_equivalence.bond_coef = self.find_ff_atom(parts[4], i)
                    ff_equivalence.angle_coef = self.find_ff_atom(parts[5], i)
                    ff_equivalence.torsion_coef = self.find_ff_atom(parts[6], i)
                    ff_equivalence.improper_coef = self.find_ff_atom(parts[7], i)
                    
                    if ff_equivalence.ff_atom is not None:
                        self.ff_equivalences.append(ff_equivalence)

    def read_ff_nonbond_coef(self, file_path, i):

        def calculate_sigma_epsilon(A, B):
            if A == 0 or B == 0:
                return 0, 0  
            sigma = (A / B) ** (1/6)
            epsilon = B / (4 * sigma**6)
            return sigma, epsilon

        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:

                if line.startswith('#') and start_processing == True:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>') or line.startswith('@'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#nonbond_coef(12-6)' and line_parts[1] == 'cvff':
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()

                    ff_nonbond_coef = FF_nonbond_coef()
                    ff_nonbond_coef.ff_atoms = self.find_ff_atom(parts[2], i)
                    A, B = float(parts[3]), float(parts[4])
                    ff_nonbond_coef.sigma, ff_nonbond_coef.epsilon = calculate_sigma_epsilon(A, B)
                    
                    self.nonbond_coefs.append(ff_nonbond_coef)

    def read_ff_bond(self, file_path, i):
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:

                if line.startswith('#') and start_processing == True:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>') or line.startswith('@'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#quadratic_bond' and line_parts[1] == 'cvff': 
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()

                    ff_bond = FF_bond_coef()
                    ff_bond.ff_atoms = [self.find_ff_atom(parts[2], i), self.find_ff_atom(parts[3], i)]
                    ff_bond.r0 = float(parts[4])
                    ff_bond.k = float(parts[5]) / 2
                    
                    self.bond_coefs.append(ff_bond)

    def read_ff_angle(self, file_path, i):
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:

                if line.startswith('#') and start_processing == True:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>') or line.startswith('@'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#quadratic_angle' and line_parts[1] == 'cvff':    
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()

                    ff_angle = FF_angle_coef()
                    ff_angle.ff_atoms = [self.find_ff_atom(parts[2],i), self.find_ff_atom(parts[3],i), self.find_ff_atom(parts[4],i)]
                    ff_angle.theta0 = float(parts[5])
                    ff_angle.k = float(parts[6]) / 2
                    
                    self.angle_coefs.append(ff_angle)

    def read_ff_torsion(self, file_path, i):
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:

                if line.startswith('#') and start_processing == True:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>') or line.startswith('@'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#torsion_1' and line_parts[1] == 'cvff':
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()

                    ff_torsion = FF_torsion_coef()
                    ff_torsion.ff_atoms = [self.find_ff_atom(parts[2],i), self.find_ff_atom(parts[3],i), self.find_ff_atom(parts[4],i), self.find_ff_atom(parts[5],i)]
                    ff_torsion.kphi = float(parts[6])
                    ff_torsion.phi0 = float(parts[8])
                    ff_torsion.n = int(parts[7])
                    
                    self.torsion_coefs.append(ff_torsion)

    def read_ff_improper(self, file_path, i):
        start_processing = False
        with open(file_path, 'r') as file:
            for line in file:

                if line.startswith('#') and start_processing == True:
                    return

                if line.startswith('!') or line.strip() == '' or line.startswith('>') or line.startswith('@'):
                    continue
                line_parts = line.strip().split()
                if len(line_parts) > 1 and line_parts[0] == '#out_of_plane' and line_parts[1] == 'cvff':    
                    start_processing = True
                    continue
                if start_processing:
                    parts = line.strip().split()

                    ff_improper = FF_improper_coef()

                    ff_improper.ff_atoms = [self.find_ff_atom(parts[2],i), self.find_ff_atom(parts[3],i), self.find_ff_atom(parts[4],i), self.find_ff_atom(parts[5],i)]
                    ff_improper.kchi = float(parts[6])
                    ff_improper.n = int(parts[7])
                    ff_improper.chi0 = float(parts[8])
                    
                    self.improper_coefs.append(ff_improper)
    

    def add_nonbond_coef_equivalences(self):
        equivalences = [equivalence for equivalence in self.ff_equivalences]
        non_bonds = [nonbond_coef.ff_atoms for nonbond_coef in self.nonbond_coefs]

        for equivalence in equivalences:
            if equivalence.ff_atom not in non_bonds:
                new_ff_nonbond_coef = FF_nonbond_coef()
                new_ff_nonbond_coef.ff_atoms = equivalence.ff_atom
                
                for nonbond_coef in self.nonbond_coefs:
                    if equivalence.nonbond_coef == nonbond_coef.ff_atoms:
                        new_ff_nonbond_coef.epsilon = nonbond_coef.epsilon
                        new_ff_nonbond_coef.sigma = nonbond_coef.sigma
                        break

                if new_ff_nonbond_coef.sigma is not None and new_ff_nonbond_coef.epsilon is not None:  
                    self.nonbond_coefs.append(new_ff_nonbond_coef)
            
    def add_bond_equivalences(self):
        equivalences = [equivalence for equivalence in self.ff_equivalences]
        bond_atoms = [atom for bond_coef in self.bond_coefs for atom in bond_coef.ff_atoms]

        for equivalence in equivalences:
            if equivalence.ff_atom not in bond_atoms:
                new_ff_bond = FF_bond_coef()

                for bond_coef in self.bond_coefs:
                    if bond_coef.ff_atoms is not None:
                        if equivalence.bond_coef == bond_coef.ff_atoms[0] and equivalence.bond_coef == bond_coef.ff_atoms[1]:
                            new_ff_bond.ff_atoms = [equivalence.ff_atom, equivalence.ff_atom]
                            new_ff_bond.k = bond_coef.k
                            new_ff_bond.r0 = bond_coef.r0
                            break
                        elif equivalence.bond_coef is bond_coef.ff_atoms[0]:
                            new_ff_bond.ff_atoms = [equivalence.ff_atom, bond_coef.ff_atoms[1]]
                            new_ff_bond.k = bond_coef.k
                            new_ff_bond.r0 = bond_coef.r0
                            break
                        elif equivalence.bond_coef is bond_coef.ff_atoms[1]:
                            new_ff_bond.ff_atoms = [bond_coef.ff_atoms[0], equivalence.ff_atom]
                            new_ff_bond.k = bond_coef.k
                            new_ff_bond.r0 = bond_coef.r0
                            break

                if new_ff_bond.ff_atoms is not None:           
                    self.bond_coefs.append(new_ff_bond)

    def add_angle_equivalences(self):
        angle_atoms = {atom for angle_coef in self.angle_coefs for atom in angle_coef.ff_atoms}

        remaining_equivalences = [
            equivalence for equivalence in self.ff_equivalences
            if equivalence.ff_atom not in angle_atoms
        ]

        for equivalence in remaining_equivalences:
            new_ff_angle = FF_angle_coef()

            for angle_coef in self.angle_coefs:
                if angle_coef.ff_atoms:
                    match_indices = [i for i in range(3) if equivalence.angle_coef == angle_coef.ff_atoms[i]]

                    if match_indices:
                        new_atoms = angle_coef.ff_atoms[:]
                    
                        for index in match_indices:
                            new_atoms[index] = equivalence.ff_atom

                        new_ff_angle.ff_atoms = new_atoms
                        new_ff_angle.k = angle_coef.k
                        new_ff_angle.theta0 = angle_coef.theta0
                        break
            if new_ff_angle.ff_atoms:
                self.angle_coefs.append(new_ff_angle)

    def add_torsion_equivalences(self):
        torsion_atoms = {atom for torsion_coef in self.torsion_coefs for atom in torsion_coef.ff_atoms}

        remaining_equivalences = [
            equivalence for equivalence in self.ff_equivalences
            if equivalence.ff_atom not in torsion_atoms
        ]

        for equivalence in remaining_equivalences:
            new_ff_torsion = FF_torsion_coef()

            for torsion_coef in self.torsion_coefs:
                if torsion_coef.ff_atoms:
                    match_indices = [i for i in range(4) if equivalence.torsion_coef == torsion_coef.ff_atoms[i]]

                    if match_indices:
                        new_atoms = torsion_coef.ff_atoms[:]
                        
                        for index in match_indices:
                            new_atoms[index] = equivalence.ff_atom

                        new_ff_torsion.ff_atoms = new_atoms
                        new_ff_torsion.kphi = torsion_coef.kphi
                        new_ff_torsion.n = torsion_coef.n
                        new_ff_torsion.phi0 = torsion_coef.phi0
                        break

            if new_ff_torsion.ff_atoms:
                self.torsion_coefs.append(new_ff_torsion)

    def add_improper_equivalences(self):
        improper_atoms = {atom for improper_coef in self.improper_coefs for atom in improper_coef.ff_atoms}
        remaining_equivalences = [
            equivalence for equivalence in self.ff_equivalences
            if equivalence.ff_atom not in improper_atoms
        ]

        for equivalence in remaining_equivalences:
            new_ff_improper = FF_improper_coef()

            for improper_coef in self.improper_coefs:
                if improper_coef.ff_atoms:
                    match_indices = [i for i in range(4) if equivalence.improper_coef == improper_coef.ff_atoms[i]]

                    if match_indices:
                        new_atoms = improper_coef.ff_atoms[:]
                        
                        for index in match_indices:
                            new_atoms[index] = equivalence.ff_atom

                        new_ff_improper.ff_atoms = new_atoms
                        new_ff_improper.kchi = improper_coef.kchi
                        new_ff_improper.n = improper_coef.n
                        new_ff_improper.chi0 = improper_coef.chi0
                        break

            if new_ff_improper.ff_atoms:
                self.improper_coefs.append(new_ff_improper)


        seen_angles = set()
        duplicates = False

        for angle_coef in self.angle_coefs:
            if len(angle_coef.ff_atoms) == 3:
                ordered_atoms = (angle_coef.ff_atoms[0].type, angle_coef.ff_atoms[1].type, angle_coef.ff_atoms[2].type)

                if ordered_atoms in seen_angles:
                    print(f"Duplicate angle_coef found between atoms {ordered_atoms[0]}, {ordered_atoms[1]}, and {ordered_atoms[2]}")
                    duplicates = True
                else:
                    seen_angles.add(ordered_atoms)

        if not duplicates:
            print("No duplicate angle_coefs found.")



def test_forcefield_params():
    # Create an instance of the ForceFieldLoader
    ff_params = ForceFieldLoader()

    # Check if any atom types have been loaded
    if len(ff_params.ff_atoms) > 0:
        print(f"###################Successfully loaded {len(ff_params.ff_atoms)} atom types:")
        for atom in ff_params.ff_atoms:
            print(atom)
    else:
        print("No atom types were loaded. Check file paths and file content.")

    # Check if any nonbond_coef parameters have been loaded
    if len(ff_params.nonbond_coefs) > 0:
        print(f"###################Successfully loaded {len(ff_params.nonbond_coefs)} nonbond_coef parameters:")
        for nonbond_coef in ff_params.nonbond_coefs:
            print(f"Types: {nonbond_coef.ff_atoms.type}, Epsilon: {nonbond_coef.epsilon}, Sigma: {nonbond_coef.sigma}")
    else:
        print("No nonbond_coef parameters were loaded.")

    # Check if any bond_coef parameters have been loaded
    if len(ff_params.bond_coefs) > 0:
        print(f"###################Successfully loaded {len(ff_params.bond_coefs)} bond_coef parameters:")
        for bond_coef in ff_params.bond_coefs:
            if bond_coef.ff_atoms:
                # Using a list comprehension to extract 'type' from each atom in the list
                atom_types = ', '.join(atom.type for atom in bond_coef.ff_atoms if atom is not None)
                print(f"Types: {atom_types}, K: {bond_coef.k}, R0: {bond_coef.r0}")
            else:
                print("Warning: No atoms found in this bond_coef.")
    else:
        print("No bond_coef parameters were loaded.")

    # Check if any angle_coef parameters have been loaded
    if len(ff_params.angle_coefs) > 0:
        print(f"###################Successfully loaded {len(ff_params.angle_coefs)} angle_coef parameters:")
        for angle_coef in ff_params.angle_coefs:
            if angle_coef.ff_atoms:
                atom_types = ', '.join(atom.type for atom in angle_coef.ff_atoms if atom is not None)
                print(f"Types: {atom_types}, K: {angle_coef.k}, Theta0: {angle_coef.theta0}")
            else:
                print("Warning: No atoms found in this angle_coef.")
    else:
        print("No angle_coef parameters were loaded.")

    # Check if any torsion_coef parameters have been loaded
    if len(ff_params.torsion_coefs) > 0:
        print(f"###################Successfully loaded {len(ff_params.torsion_coefs)} torsion_coef parameters:")
        for torsion_coef in ff_params.torsion_coefs:
            if torsion_coef.ff_atoms:
                atom_types = ', '.join(atom.type for atom in torsion_coef.ff_atoms if atom is not None)
                print(f"Types: {atom_types}, Kphi: {torsion_coef.kphi}, N: {torsion_coef.n}, Phi0: {torsion_coef.phi0}")
            else:
                print("Warning: No atoms found in this torsion_coef.")
    else:
        print("No torsion_coef parameters were loaded.")

    # Check if any improper_coef parameters have been loaded
    if len(ff_params.improper_coefs) > 0:
        print(f"###################Successfully loaded {len(ff_params.improper_coefs)} improper_coef parameters:")
        for improper_coef in ff_params.improper_coefs:
            if improper_coef.ff_atoms:
                atom_types = ', '.join(atom.type for atom in improper_coef.ff_atoms if atom is not None)
                print(f"Types: {atom_types}, Kchi: {improper_coef.kchi}, N: {improper_coef.n}, Chi0: {improper_coef.chi0}")
            else:
                print("Warning: No atoms found in this improper_coef.")
    else:
        print("No improper_coef parameters were loaded.")

    # Print the equivalence parameters to the console for verification
    if ff_params.ff_equivalences:
        print(f"################### Successfully loaded {len(ff_params.ff_equivalences)} equivalence parameters:")
        for equivalence in ff_params.ff_equivalences:
            print(f"FF Atom: {equivalence.ff_atom.type}, "
                f"nonbond_coef: {equivalence.nonbond_coef.type}, "
                f"bond_coef: {equivalence.bond_coef.type}, "
                f"angle_coef: {equivalence.angle_coef.type}, "
                f"torsion_coef: {equivalence.torsion_coef.type}, "
                f"improper_coef: {equivalence.improper_coef.type}")
    else:
        print("No equivalence parameters were loaded.")


    if len(ff_params.duplicates) > 0:
        print(f"{len(ff_params.duplicates)} Detected duplicate atom types:")
        for duplicate_type in ff_params.duplicates:
            # Print all atoms that have the duplicated type
            for atom in ff_params.ff_atoms:
                if atom.type == duplicate_type:
                    print(f"Duplicate Type: {atom.type}, Description: {atom.description}")
    else:
        print("No duplicates detected.")

 
# Call the test function
test_forcefield_params()

