
"""
This class is used to store the chosen settings for the simulation.
"""

class ChosenSettings:
    def __init__(self, 
                 replication, 
                 al_mg_ratio, 
                 net_charge, 
                 water_per_ion, 
                 output_file,
                 ff_atom_types, 
                 water_distance, 
                 clay_sub_cutoff,
                 random_seed,
                 mg_cutoff, 
                 h_cutoff, 
                 bond_cutoff, 
                 input_file, 
                 ff_files, 
                 water_file):
        
        self.replication = replication
        self.al_mg_ratio = al_mg_ratio
        self.net_charge = net_charge
        self.water_per_ion = water_per_ion
        self.output_file = output_file

        self.ff_atom_types = ff_atom_types
        
        self.water_distance = water_distance
        self.clay_sub_cutoff = clay_sub_cutoff
        self.random_seed = random_seed
        self.mg_cutoff = mg_cutoff
        self.h_cutoff = h_cutoff
        self.bond_cutoff = bond_cutoff

        self.input_file = input_file
        self.ff_files = ff_files
        self.water_file = water_file

        if water_per_ion < 0:  # Define the threshold as needed
            raise ValueError("Invalid value for water_per_ion. It must be non-negative.")

        if water_distance < 0: # Define the threshold as needed
            raise ValueError("Invalid value for water_distance. It must be non-negative.")
        

