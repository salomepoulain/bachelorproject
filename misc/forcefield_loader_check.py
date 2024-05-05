
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
# test_forcefield_params()