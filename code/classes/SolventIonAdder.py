from code.classes.SystemAllocator import SystemAllocator
from code.classes.WaterBoxBuilder import WaterBoxBuilder
from code.classes.SystemParts import Molecule
import math
import numpy as np
from typing import List, Tuple

"""
This script defines the SolventIonAdder class, which is responsible for adding solvent
and ions to the system based on input settings. The class calculates the system charge,
the number of ions needed to neutralize the system, and adds ions uniformly throughout
the simulation box. The class also allocates force field parameters and bonds for ions
and water molecules.
"""

class SolventIonAdder(SystemAllocator):
    def __init__(self, settings) -> None:
        """
        Initialize the SolventIonAdder with the given settings.
        
        Args:
            settings (ChosenSettings): The settings object containing various parameters.
        """
        super().__init__(settings)
        self.s = settings

        self.system_charge = self.calculate_system_charge()
        self.ion_count = self.calculate_ion_count()

        self.add_solvent_box()
        self.add_ions_uniformly()

        self.allocate_ions_and_water()

        self.final_charge = self.calculate_system_charge()

    def calculate_system_charge(self) -> float:
        """
        Calculate the total charge of the system based on the charges of all atoms.
        
        Returns:
            float: The total system charge.
        """
        charge = 0.0
        for molecule in self.molecules:
            for atom in molecule.atoms:
                charge += atom.ff_atom.charge
        
        return charge

    def calculate_ion_count(self) -> int:
        """
        Calculate the number of ions needed to neutralize the system.
        
        Returns:
            int: The calculated number of ions.
        """
        ion_charge = self.ff_attributes['Ca']['ff_atom'].charge
        ion_count = self.custom_round(abs(self.s.net_charge - self.system_charge) / ion_charge)
        return ion_count

    def get_clay_height(self) -> Tuple[float, float]:
        """
        Get the height (z-dimension) of the clay layer.
        
        Returns:
            Tuple[float, float]: The minimum and maximum z-coordinates of the clay atoms.
        """
        zhi = None
        zlo = None
        
        for molecule in self.molecules:
            if molecule.type == 'clay':
                for atom in molecule.atoms:
                    _, _, z = atom.position
                    
                    if zhi is None or zlo is None:
                        zhi = zlo = z
                    else:
                        zhi = max(zhi, z)
                        zlo = min(zlo, z)
        return zhi, zlo

    def translate_simulation_height(self) -> None:
        """
        Translate the height of the simulation box to accommodate the clay and solvent.
        """
        clay_height = (self.get_clay_height()[0] - self.get_clay_height()[1]) / 2.0

        if self.s.water_per_ion == 0:
            # Add additional height for systems without water per ion
            self.dimensions = (
                self.dimensions[0], self.dimensions[1], 
                (self.dimensions[2][0], self.dimensions[2][1] + 20)  # Example extra distance
            )
            self.dimensions = (
                self.dimensions[0], self.dimensions[1], 
                (self.dimensions[2][0] + clay_height, self.dimensions[2][1] + clay_height)
            )
            return
        
        self.dimensions = (
            self.dimensions[0], self.dimensions[1], 
            (self.dimensions[2][0], self.dimensions[2][1] + self.s.water_distance)
        )
        self.dimensions = (
            self.dimensions[0], self.dimensions[1], 
            (self.dimensions[2][0] + clay_height, self.dimensions[2][1] + clay_height)
        )

    def add_solvent_box(self) -> None:
        """
        Add a solvent box to the system based on the water per ion ratio.
        """
        if self.s.water_per_ion == 0:
            return

        clay_height = self.get_clay_height()[0] - self.get_clay_height()[1]

        waterbox = WaterBoxBuilder(self.ion_count, self.s)

        waterbox.translate_new_origin(self.dimensions[0][0], self.dimensions[1][0], self.dimensions[2][0])
        waterbox.replicate_water((
            self.dimensions[0][1] - self.dimensions[0][0], 
            self.dimensions[1][1] - self.dimensions[1][0], 
            self.dimensions[2][1] - self.dimensions[2][0]
        ))
        waterbox.remove_outside_box((
            self.dimensions[0][1] - self.dimensions[0][0], 
            self.dimensions[1][1] - self.dimensions[1][0]
        ))
        waterbox.translate_up(clay_height)
        waterbox.change_water_per_ions()

        for molecule in waterbox.molecules:
            self.molecules.append(molecule)

    def add_ions_uniformly(self) -> None:
        """
        Add ions uniformly throughout the simulation box.
        """
        self.get_unit_cell_dimensions()
        self.translate_simulation_height()

        if any(molecule.type == 'clay' for molecule in self.molecules):
            if self.ion_count == 0:
                return

            ca_height = (self.dimensions[2][0] + self.dimensions[2][1]) / 2
            ion_count = self.ion_count

            def is_prime(n: int) -> bool:
                if n <= 1:
                    return False
                for i in range(2, int(n**0.5) + 1):
                    if n % i == 0:
                        return False
                return True

            if is_prime(ion_count):
                ion_count -= 1
                place_last_ion_in_center = True
            else:
                place_last_ion_in_center = False

            sqrt_ion_count = int(ion_count ** 0.5)

            if sqrt_ion_count ** 2 == ion_count:
                n_x = n_y = sqrt_ion_count
            else:
                factors = [(i, ion_count // i) for i in range(1, sqrt_ion_count + 1) if ion_count % i == 0]
                n_x, n_y = min(factors, key=lambda x: abs(x[0] - x[1]))

            x_spacing = (self.dimensions[0][1] - self.dimensions[0][0]) / n_x
            y_spacing = (self.dimensions[1][1] - self.dimensions[1][0]) / n_y

            for j in range(n_y):
                for i in range(n_x):
                    x = self.dimensions[0][0] + i * x_spacing
                    y = self.dimensions[1][0] + j * y_spacing

                    ion = Molecule()
                    ion.type = 'ion'
                    ion.add_atom('Ca', (x, y, ca_height))
                    self.molecules.append(ion)

            if place_last_ion_in_center:
                mid_x = (self.dimensions[0][0] + self.dimensions[0][1]) / 2
                mid_y = (self.dimensions[1][0] + self.dimensions[1][1]) / 2
                ion = Molecule()
                ion.type = 'ion'
                ion.add_atom('Ca', (mid_x, mid_y, ca_height))
                self.molecules.append(ion)

    def allocate_ions_and_water(self) -> None:
        """
        Allocate force field parameters and bonds for ions and water molecules.
        """
        self.process_types = {'water', 'ion'}
        self.allocate_ff_atoms(self.process_types)
        self.allocate_bonds(self.process_types)