from collections import OrderedDict
import argparse
import sys

"""
With this script, you can combine two LAMMPS data files.
"""

###############################################################################

'''
Usage: 
    Modify the code manually or use:
    - python combine.py -b <basefile> -a <added file> -o <output file>

    example:
    - python combine.py -b data.data -a lysine.data -o combined.data

    The added file is added to the base file. 
    The header and system box of the base file is used.

    - The script is dumb and just adds numbers and
    does not take duplicates into account
    - Removes comments after and between headers

'''

delta_molecules = 10000 # not really necessary

base_file = 'base.data'
added_file = 'added.data'

output_file = 'combined.data'

###############################################################################


def modify_mass(file_path):
    start_processing = False
    skip_empty = True
    masses = []

    with open(file_path, 'r') as file:
        for line in file:
            if start_processing and line.strip() == '' and skip_empty == False:
                return masses

            if line.strip().startswith('Masses'):
                start_processing = True
                continue

            if start_processing:
                if skip_empty:
                    skip_empty = False
                    continue

                parts = line.strip().split()
                
                new_line = []

                first = int(parts.pop(0)) + delta_atom_type
                new_line.append(first)

                for part in parts:
                    new_line.append(part)
        
                masses.append(new_line)

    return masses
    
def modify_pair_coef(file_path):
    start_processing = False
    skip_empty = True
    pair_coeffs = []

    with open(file_path, 'r') as file:
        for line in file:
            if start_processing and line.strip() == '' and skip_empty == False:
                return pair_coeffs

            if line.strip().startswith('Pair Coeffs'):
                start_processing = True
                continue

            if start_processing:
                if skip_empty:
                    skip_empty = False
                    continue

                parts = line.strip().split()
                
                new_line = []

                first = int(parts.pop(0)) + delta_atom_type
                new_line.append(first)

                for part in parts:
                    new_line.append(part)
        
                pair_coeffs.append(new_line)

    return pair_coeffs

def modify_bond_coef(file_path):
    start_processing = False
    skip_empty = True
    bond_coeffs = []

    with open(file_path, 'r') as file:
        for line in file:
            if start_processing and line.strip() == '' and skip_empty == False:
                return bond_coeffs

            if line.strip().startswith('Bond Coeffs'):
                start_processing = True
                continue

            if start_processing:
                if skip_empty:
                    skip_empty = False
                    continue

                parts = line.strip().split()
                
                parts = line.strip().split()
                
                new_line = []

                first = int(parts.pop(0)) + delta_bond_coeff
                new_line.append(first)

                for part in parts:
                    new_line.append(part)
        
                bond_coeffs.append(new_line)

    return bond_coeffs

def modify_angle_coef(file_path):
    start_processing = False
    skip_empty = True
    angle_coeffs = []

    with open(file_path, 'r') as file:
        for line in file:
            if start_processing and line.strip() == '' and skip_empty == False:
                return angle_coeffs

            if line.strip().startswith('Angle Coeffs'):
                start_processing = True
                continue

            if start_processing:
                if skip_empty:
                    skip_empty = False
                    continue

                parts = line.strip().split()
                
                new_line = []

                first = int(parts.pop(0)) + delta_angle_coeff
                new_line.append(first)

                for part in parts:
                    new_line.append(part)
        
                angle_coeffs.append(new_line)

    return angle_coeffs

def modify_dihedral_coef(file_path):
    start_processing = False
    skip_empty = True
    dihedral_coeffs = []

    with open(file_path, 'r') as file:
        for line in file:
            if start_processing and line.strip() == '' and skip_empty == False:
                return dihedral_coeffs

            if line.strip().startswith('Dihedral Coeffs'):
                start_processing = True
                continue

            if start_processing:
                if skip_empty:
                    skip_empty = False
                    continue

                parts = line.strip().split()
                
                new_line = []

                first = int(parts.pop(0)) + delta_dihedral_coeff
                new_line.append(first)

                for part in parts:
                    new_line.append(part)
        
                dihedral_coeffs.append(new_line)

    return dihedral_coeffs

def modify_improper_coef(file_path):
    start_processing = False
    skip_empty = True
    improper_coeffs = []

    with open(file_path, 'r') as file:
        for line in file:
            if start_processing and line.strip() == '' and skip_empty == False:
                return improper_coeffs

            if line.strip().startswith('Improper Coeffs'):
                start_processing = True
                continue

            if start_processing:
                if skip_empty:
                    skip_empty = False
                    continue

                parts = line.strip().split()
                
                new_line = []

                first = int(parts.pop(0)) + delta_improper_coeff
                new_line.append(first)

                for part in parts:
                    new_line.append(part)
        
                improper_coeffs.append(new_line)

    return improper_coeffs

def modify_atoms(file_path):
    start_processing = False
    skip_empty = True
    atoms = []

    with open(file_path, 'r') as file:
        for line in file:
            if start_processing and line.strip() == '' and skip_empty == False:
                return atoms

            if line.strip().startswith('Atoms'):
                start_processing = True
                continue

            if start_processing:
                if skip_empty:
                    skip_empty = False
                    continue

                parts = line.strip().split()
                
                new_line = []

                first = int(parts.pop(0)) + delta_atoms
                new_line.append(first)
                
                second = int(parts.pop(0)) + delta_molecules
                new_line.append(second)

                third = int(parts.pop(0)) + delta_atom_type
                new_line.append(third)

                for part in parts:
                    new_line.append(part)
        
                atoms.append(new_line)

    return atoms

def modify_bonds(file_path):
    start_processing = False
    skip_empty = True
    bonds = []

    with open(file_path, 'r') as file:
        for line in file:
            if start_processing and line.strip() == '' and skip_empty == False:
                return bonds

            if line.strip().startswith('Bonds'):
                start_processing = True
                continue

            if start_processing:
                if skip_empty:
                    skip_empty = False
                    continue

                parts = line.strip().split()
                
                new_line = []

                first = int(parts.pop(0)) + delta_bonds
                new_line.append(first)
                
                second = int(parts.pop(0)) + delta_bond_coeff
                new_line.append(second)

                third = int(parts.pop(0)) + delta_atoms
                new_line.append(third)

                fourth = int(parts.pop(0)) + delta_atoms
                new_line.append(fourth)

                for part in parts:
                    new_line.append(part)
        
                bonds.append(new_line)

    return bonds

def modify_angles(file_path):
    start_processing = False
    skip_empty = True
    angles = []

    with open(file_path, 'r') as file:
        for line in file:
            if start_processing and line.strip() == '' and skip_empty == False:
                return angles

            if line.strip().startswith('Angles'):
                start_processing = True
                continue

            if start_processing:
                if skip_empty:
                    skip_empty = False
                    continue

                parts = line.strip().split()
                
                new_line = []

                first = int(parts.pop(0)) + delta_angles
                new_line.append(first)
                
                second = int(parts.pop(0)) + delta_angle_coeff
                new_line.append(second)

                third = int(parts.pop(0)) + delta_atoms
                new_line.append(third)

                fourth = int(parts.pop(0)) + delta_atoms
                new_line.append(fourth)

                fifth = int(parts.pop(0)) + delta_atoms
                new_line.append(fifth)

                for part in parts:
                    new_line.append(part)
        
                angles.append(new_line)

    return angles

def modify_dihedrals(file_path):
    start_processing = False
    skip_empty = True
    dihedrals = []

    with open(file_path, 'r') as file:
        for line in file:
            if start_processing and line.strip() == '' and skip_empty == False:
                return dihedrals

            if line.strip().startswith('Dihedrals'):
                start_processing = True
                continue

            if start_processing:
                if skip_empty:
                    skip_empty = False
                    continue

                parts = line.strip().split()
                
                new_line = []

                first = int(parts.pop(0)) + delta_dihedrals
                new_line.append(first)
                
                second = int(parts.pop(0)) + delta_dihedral_coeff
                new_line.append(second)

                third = int(parts.pop(0)) + delta_atoms
                new_line.append(third)

                fourth = int(parts.pop(0)) + delta_atoms
                new_line.append(fourth)

                fifth = int(parts.pop(0)) + delta_atoms
                new_line.append(fifth)

                sixth = int(parts.pop(0)) + delta_atoms
                new_line.append(sixth)

                for part in parts:
                    new_line.append(part)
        
                dihedrals.append(new_line)
    
    return dihedrals

def modify_impropers(file_path):
    start_processing = False
    skip_empty = True
    impropers = []

    with open(file_path, 'r') as file:
        for line in file:

            if start_processing == True and skip_empty == False and line.strip() == '':
                return impropers

            if line.strip().startswith('Impropers'):
                start_processing = True
                continue

            if start_processing:
                if skip_empty:
                    skip_empty = False
                    continue

                parts = line.strip().split()
                
                new_line = []

                first = int(parts.pop(0)) + delta_impropers
                new_line.append(first)
                
                second = int(parts.pop(0)) + delta_improper_coeff
                new_line.append(second)

                third = int(parts.pop(0)) + delta_atoms
                new_line.append(third)

                fourth = int(parts.pop(0)) + delta_atoms
                new_line.append(fourth)

                fifth = int(parts.pop(0)) + delta_atoms
                new_line.append(fifth)

                sixth = int(parts.pop(0)) + delta_atoms
                new_line.append(sixth)

                for part in parts:
                    new_line.append(part)
        
                impropers.append(new_line)

    return impropers

def extract_info(file_path, string):
    with open(file_path, 'r') as file:
        for line in file:
            if string in line:
                return int(line.split()[0])

def write_info(file_path, output_file, items_dict):
    # Ensure items_dict is a dictionary with the correct format
    with open(file_path, 'r') as file, open(output_file, 'w') as out_file:
        for line in file:
            modified = False
            for string, item in items_dict.items():
                if string in line:
                    number = int(line.split()[0]) + len(item)
                    modified_line = f"{number} {' '.join(line.split()[1:])}"
                    out_file.write(modified_line + '\n')
                    modified = True
                    break
            if not modified:
                out_file.write(line)

def ensure_proper_comment_format(line):
    if '#' in line:
        parts = line.split('#', 1)
        parts[1] = '\t# ' + parts[1].lstrip()
        return parts[0] + parts[1]
    return line

def format_line(inner_array):
    if '#' in inner_array:
        comment_index = inner_array.index('#')
        before_comment = inner_array[:comment_index]
        comment = inner_array[comment_index:]
        new_line = '\t'.join(map(str, before_comment)) + '\t' + ' '.join(map(str, comment)) + '\n'
    else:
        new_line = '\t'.join(map(str, inner_array)) + '\n'
    return new_line

def write_lines(file_path, headers_dict, added_file, base_file):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    start_processing = False
    header_lines = []
    sections = {}
    current_section = None
    section_lines = []

    # Read lines and process the header section
    for line in lines:
        if not start_processing:
            header_lines.append(line)
            if 'zhi' in line.strip():
                start_processing = True
            continue

        stripped_line = line.split('#')[0].strip()

        if stripped_line.startswith('#'):
            continue
        elif stripped_line == '':
            continue
        elif stripped_line in headers_dict.keys():
            if current_section:
                sections[current_section] = section_lines
            current_section = stripped_line
            section_lines = [line]
        else:
            section_lines.append(line)

    if current_section:
        sections[current_section] = section_lines

    # Write the header section back to the file
    with open(file_path, 'w') as file:
        header_lines[0] = f"Combined datafile of {added_file} added to {base_file}\n\n"
        file.writelines(header_lines)
        file.write('\n\n')

        # Write the sections with modifications
        for header, data in headers_dict.items():
            if header in sections:
                out_lines = sections[header]
                file.write(f'\n{header}\n\n')
                # Ensure proper formatting of existing lines
                formatted_lines = [ensure_proper_comment_format(line) for line in out_lines[1:]]
                file.writelines(formatted_lines)  # Skip the original header line from out_lines
                # Append the new data at the end of the existing section
                if data:
                    for inner_array in data:
                        new_line = format_line(inner_array)
                        file.write(new_line)
            elif data:
                # If the section is missing, add it with its data
                file.write(f'\n{header}\n\n')
                for inner_array in data:
                    new_line = format_line(inner_array)
                    file.write(new_line)
                file.write('\n')

def cleanup_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Replace double newlines with a single newline
    cleaned_content = content.replace('\n\n\n', '\n\n')
    
    with open(file_path, 'w') as file:
        file.write(cleaned_content)

def write_all(base_file, output_file, added_file, masses, pair_coeffs, bond_coeffs, angle_coeffs, dihedral_coeffs, improper_coeffs, atoms, bonds, angles, dihedrals, impropers):
    items_dict = {
    'atom types': masses,   
    'bond types': bond_coeffs,
    'angle types': angle_coeffs,
    'dihedral types': dihedral_coeffs,
    'improper types': improper_coeffs,
    'atoms': atoms,
    'bonds': bonds,
    'angles': angles,
    'dihedrals': dihedrals,
    'impropers': impropers
    }

    headers_dict = OrderedDict([
    ('Masses', masses),
    ('Pair Coeffs', pair_coeffs),
    ('Bond Coeffs', bond_coeffs),
    ('Angle Coeffs', angle_coeffs),
    ('Dihedral Coeffs', dihedral_coeffs),
    ('Improper Coeffs', improper_coeffs),
    ('Atoms', atoms),
    ('Bonds', bonds),
    ('Angles', angles),
    ('Dihedrals', dihedrals),
    ('Impropers', impropers)
    ])

    # Step 1: Modify and write the input file to the output file
    write_info(base_file, output_file, items_dict)

    # Step 2: Add the new data to the output file
    write_lines(output_file, headers_dict, added_file, base_file)
    # for header, data in headers_dict.items():
    #     write_lines(output_file, data, header)

    # Step 3: Clean up the file
    cleanup_file(output_file)

    print(f'\nFile {output_file} has been created with {added_file} added to {base_file}.')

if len(sys.argv) > 2:
    parser = argparse.ArgumentParser(description="LAMMPS data file combiner")
    parser.add_argument('-b', '--basefile', required=True, help="Path to the base file.")
    parser.add_argument('-a', '--addedfile', required=True, help="Path to the added file.")
    parser.add_argument('-o', '--outputfile', required=True, help="Path to the output file.")

    args = parser.parse_args()

    if args.basefile and args.addedfile and args.outputfile:
        base_file = args.basefile
        added_file = args.addedfile
        output_file = args.outputfile

delta_atoms = extract_info(base_file, 'atoms')
delta_bonds = extract_info(base_file, 'bonds')
delta_angles = extract_info(base_file, 'angles')
delta_dihedrals = extract_info(base_file, 'dihedrals')
delta_impropers = extract_info(base_file, 'impropers')

delta_atom_type = extract_info(base_file, 'atom types')
delta_bond_coeff = extract_info(base_file, 'bond types')
delta_angle_coeff = extract_info(base_file, 'angle types')
delta_dihedral_coeff = extract_info(base_file, 'dihedral types')
delta_improper_coeff = extract_info(base_file, 'improper types')

masses = modify_mass(added_file)
pair_coeffs = modify_pair_coef(added_file)
bond_coeffs = modify_bond_coef(added_file)
angle_coeffs = modify_angle_coef(added_file)
dihedral_coeffs = modify_dihedral_coef(added_file)
improper_coeffs = modify_improper_coef(added_file)

atoms = modify_atoms(added_file)
bonds = modify_bonds(added_file)
angles = modify_angles(added_file)
dihedrals = modify_dihedrals(added_file)
impropers = modify_impropers(added_file)

write_all(base_file, 
            output_file,
            added_file,
            masses,
            pair_coeffs,
            bond_coeffs,
            angle_coeffs,
            dihedral_coeffs,
            improper_coeffs,
            atoms,
            bonds,
            angles,
            dihedrals,
            impropers)



