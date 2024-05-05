'''
Function to trim the structure, used for creating the first unit cell
'''
def trim_structure(self):
    x_trim_start = -2.876
    x_trim_end = 2.316
    y_trim_start = -7.664
    y_trim_end = 1.351
    
    self.atoms = [atom for atom in self.atoms if x_trim_start <= atom.position[0] < x_trim_end and y_trim_start <= atom.position[1] < y_trim_end]
    
    unit_x = (x_trim_start, x_trim_end)
    unit_y = (y_trim_start, y_trim_end)
    unit_z = self.dimensions[2]  # Assumes z-dimension is unaffected by trimming
    
    self.dimensions = (unit_x, unit_y, unit_z)

    print(f"New dimensions: {self.dimensions}")
    print(f"Number of atoms after trimming: {len(self.atoms)}")
    
    output_file_name = self.file_name + "_trim.xyz"
    self.write_xyz(output_file_name)