import MDAnalysis as mda

# Path to your LAMMPS dump file
dump_file = 'dump.data'

# Create a Universe object from the dump file
u = mda.Universe(dump_file, format='LAMMPSDUMP')

# Now you can access various attributes and methods of the Universe object
# For example, to access all the atoms you can use:
atoms = u.atoms
print(atoms)

# To access positions of atoms:
positions = atoms.positions
print(positions)

# You can also perform analysis, like calculating the center of mass:
center_of_mass = atoms.center_of_mass()
print("Center of Mass:", center_of_mass)
