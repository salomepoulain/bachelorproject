from code.classes.SystemWriter import SystemWriter

'''
HARD CODED VALUES

reading from clayff and cvff force field in the forcefield folder
manually making functions to allocate the values of this forcefield that need to be used (i.e. 'oh' and 'al')
The following thresholds for bond defenitions:
- distance_threshold_mg = 2.0
- distance_threshold_h = 1.0
- bond_cutoff = 1.0
'''


if __name__ == '__main__':

    system = SystemWriter(file_name     =     'STx_prot', 
                          replication   =     (6,6), 
                          height        =     30,
                          al_mg_ratio   =     7.1702509,    # Based on STx1b data (Castellini 2017)
                          ca_si_ratio   =     0.055125       # Based on STx1b data (Castellini 2017)
                          )        

    

'''
https://link.springer.com/article/10.1346/CCMN.2017.064065 (Casteleinni 2017)
'''