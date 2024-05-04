from code.classes.SystemWriter import SystemWriter

'''
HARD CODED VALUES

reading from clayff and cvff force field in the forcefield folder
manually making functions to allocate the values of this forcefield that need to be used (i.e. 'oh' and 'al')
if a new ff_param is added, or removed, this should be manually taken care of in the SystemAllocator class
solvates system based on gromacs method using spc216 water model box
ions are added in a horizontal plane, spaced out with pbc in mind

by enlargening the system, the accuracy of the system will be increased

The following thresholds for bond defenitions:
- distance_threshold_mg = 2.0
- distance_threshold_h = 1.0
- bond_cutoff = 1.5

for the solvent system, the following thresholds are used:
- default_radius = 1.05
- self.default_scale = 5.7 (10x gromacs default)

'''


if __name__ == '__main__':

    system = SystemWriter(file_name     =     'STx_prot', 
                          replication   =     (4,4), 
                          height        =     30,
                          al_mg_ratio   =     7.1702509,    # Based on STx1b data (Castellini 2017)
                          ca_si_ratio   =     0.055125,     # Based on STx1b data (Castellini 2017)
                          ff_params     =     ['ao', 'mgo', 'st', 'ho', 'ob', 'ohs', 'obos', 'oh', 'Ca', 'o*', 'h*']
                          )        

    

'''
https://link.springer.com/article/10.1346/CCMN.2017.064065 (Casteleinni 2017)
'''