"""
USAGE:
    python main.py <replication number>
    write the replication number as an integer
    if no replication number is provided, the default value is in this file

The main function is called with the following parameters:
    - file_name: the name of the input unit cell file
    - replication: a tuple with the replication factor in x and y
    - height: the height of the system in Angstrom
    - al_mg_ratio: the Al/Mg ratio in the system
    - ca_si_ratio: the Ca/Si ratio in the system
    - water_file: the name of the water model file
    - vdw_radii_file: the name of the Van der Waals radii file
    - vdw_def_radius: the default Van der Waals radius
    - vdw_def_scale: the default Van der Waals scale factor
    - ff_file_paths: a list with the force field file paths
    - mg_cutoff: the cutoff distance for Mg atoms for oxygen allocation in MT
    - h_cutoff: the cutoff distance for H atoms for oxygen allocation in MT
    - bond_cutoff: the cutoff distance for bonds
    - ff_params: a list with the force field parameters called ff_types

Note: 
    Adding a ff_param, requires to manually add this allocatoin in SystemAllocator.py
    Adding ff parameters requires knowledge about the system and manually finding the correct type in the force field file(s)
    Duplicate ff_types will be called "{ff_type}[dup]" for the second instance
    A manual check for the .data file is recommended to ensure that the correct parameters are allocated to the atoms
"""

from code.main_runner import main


if __name__ == "__main__":
    main(input_file =       'STx_prot',
         replication =      (4,4), 
         height =           30,
         al_mg_ratio =      7.1702509,      # Based on STx1b data (Castellini 2017)
         ca_si_ratio =      0.055125,       # Based on STx1b data (Castellini 2017)
         
         water_file =       'spc216',
         vdw_radii_file =   'vdwradaii',
         vdw_def_radius =    1.05,          # default if not in file
         vdw_def_scale =     5.70,          # helping to achieve a density close to 1000 g/l for proteins
                                            # ^ this is to be changed.

         ff_files =         ['clayff', 'cvff'],
         mg_cutoff     =    2.0,        
         h_cutoff      =    1.0,       
         bond_cutoff   =    1.5,      

         ff_params =        ['ao', 
                             'mgo', 
                             'st', 
                             'ho', 
                             'ob', 
                             'ohs', 
                             'obos', 
                             'oh', 
                             'Ca', 
                             'o*',          # SPC water model
                             'h*'])         # SPC water model

"""
https://link.springer.com/article/10.1346/CCMN.2017.064065 (Casteleinni 2017)
"""