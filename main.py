"""                              
 _        _    __  __ __  __ ____  ____    __  __ _____  
| |      / \  |  \/  |  \/  |  _ \/ ___|  |  \/  |_   _| 
| |     / _ \ | |\/| | |\/| | |_) \___ \  | |\/| | | |   
| |___ / ___ \| |  | | |  | |  __/ ___) | | |  | | | |   
|_____/_/_  \_\_|  |_|_|  |_|_| __|____/  |_|__|_| |_|   
|  _ \  / \|_   _|/ \    |  ___|_ _| |   | ____|         
| | | |/ _ \ | | / _ \   | |_   | || |   |  _|           
| |_| / ___ \| |/ ___ \  |  _|  | || |___| |___          
|____/_/___\_\_/_/ __\_\ |_|_  |___|_____|_____| ____    
 / ___| ____| \ | | ____|  _ \    / \|_   _/ _ \|  _ \   
| |  _|  _| |  \| |  _| | |_) |  / _ \ | || | | | |_) |  
| |_| | |___| |\  | |___|  _ <  / ___ \| || |_| |  _ <   
 \____|_____|_| \_|_____|_| \_\/_/   \_\_| \___/|_| \_\  
 author:  Salomé Poulain (2024)

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
    - water_distance: the distance between the clay and the water in Angstrom
    - water_duplication: duplicates the szerba water box
    - water_per_ion: the number of water molecules per ion
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

README.md contains all information about
"""

from code.main_runner import main


if __name__ == "__main__":
    main(replication =      (6,6,1), 
         al_mg_ratio =      6.15232975,      # Based on STx1b data (Castellini 2017) [2]
         net_charge =       0,
         water_per_ion =    12,  

         ff_atom_types =    ['ao', 
                             'mgo', 
                             'st', 
                             'ho', 
                             'ob', 
                             'ohs', 
                             'obos', 
                             'oh', 
                             'Ca', 
                             'o*',           # SPC water model
                             'h*',           # SPC water model,
                            ],

         water_distance =   4,
         mg_cutoff =        2.0,
         h_cutoff =         1.0,                   
         bond_cutoff =      1.5,

         input_file =       'STx_prot',     # Contains protonated unit cell from www.charmm-gui.org [1]
         ff_files =         ['clayff',
                             'cvff'],
         water_file =       'szscerba'
                            )           
"""
[1] https://www.charmm-gui.org/?doc=input
[2] https://link.springer.com/article/10.1346/CCMN.2017.064065 (Casteleinni 2017)
"""