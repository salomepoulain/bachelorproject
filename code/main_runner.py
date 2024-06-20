import threading
import time
import sys
from code.classes.SystemWriter import SystemWriter
from code.classes.ChosenSettings import ChosenSettings
import argparse
from typing import Tuple, List, Optional

"""
This script processes command-line arguments to configure and run a system writer 
for molecular simulations. Also adds a loading animation to the process.
"""

def main(replication: Tuple[int, int, int], 
         al_mg_ratio: float, 
         net_charge: float, 
         water_per_ion: float, 
         output_file: str,
         ff_atom_types: List[str], 
         water_distance: float, 
         clay_sub_cutoff: float,
         random_seed: Optional[int],
         mg_cutoff: float, 
         h_cutoff: float, 
         bond_cutoff: float, 
         input_file: str, 
         ff_files: List[str], 
         water_file: str) -> None:

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument('-r', '--replication_factor', type=int, nargs='+', help='Replication factor for dimensions (up to 3 values)')
    parser.add_argument('-w', '--water_per_ion', type=float, help='Water per ion')
    parser.add_argument('-o', '--output_file', type=str, help='Output file name')
    args = parser.parse_args()

    # Update replication factor if provided
    if args.replication_factor:
        if len(args.replication_factor) > 3:
            raise ValueError("Replication factor can have a maximum of 3 values.")
        elif len(args.replication_factor) == 1:
            replication = (args.replication_factor[0], args.replication_factor[0], 1)
        elif len(args.replication_factor) == 2:
            replication = (args.replication_factor[0], args.replication_factor[1], 1)
        elif len(args.replication_factor) == 3:
            replication = tuple(args.replication_factor)

    # Update water_per_ion if provided
    if args.water_per_ion is not None:
        water_per_ion = args.water_per_ion
    
    # Update output_file if provided
    output_file = args.output_file if args.output_file else output_file

    def loading_animation(min_display_time: int = 3) -> None:
        """
        Display a loading animation for a minimum amount of time.

        Args:
            min_display_time (int): Minimum display time in seconds.
        """
        start_time = time.time()
        chars = "/—\\|"
        while time.time() - start_time < min_display_time or run_animation.is_set():
            for char in chars:
                if not run_animation.is_set() and time.time() - start_time >= min_display_time:
                    break
                sys.stdout.write('\r' + 'Loading... ' + char)
                sys.stdout.flush()
                time.sleep(0.1)
        sys.stdout.write('\r' + final_text + '\n\n')

    def print_banner(file_path: str) -> None:
        """
        Print the contents of a banner file.

        Args:
            file_path (str): Path to the banner file.
        """
        try:
            with open(file_path, 'r') as file:
                content = file.read()
                print(content)
        except FileNotFoundError:
            print("The file was not found.")
        except Exception as e:
            print(f"An error occurred: {e}")

    # Print the banner
    print_banner('misc/banner.txt')

    run_animation = threading.Event()
    if replication[0] >= 4 and replication[1] >= 4:
        run_animation.set()
        thread = threading.Thread(target=loading_animation)
        thread.start()

    try:
        # Initialize settings and system writer
        settings = ChosenSettings(replication, 
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
                                  water_file)
        
        runner = SystemWriter(settings)

        # Prepare the final text to display
        final_text = (
            f"Finished writing .data file for {replication[0]}x{replication[1]}x{replication[2]} system.\n\n"
            f"○ Output file:            {runner.output_name}.data\n"
            f"○ Water per ion:          {water_per_ion}\n"
            f"○ Sheet charge:           {runner.system_charge:.3f}\n"
            f"○ Ions per sheet added:   {runner.ion_count}\n"
            f"○ Total water added:      {len([molecule for molecule in runner.molecules if molecule.type == 'water'])}\n"
            f"○ Total final charge:     {runner.final_charge:.3f}"
        )

    finally:
        # Stop the loading animation and print the final text
        if run_animation.is_set():
            run_animation.clear()
            thread.join()
        else:
            sys.stdout.write('\r' + final_text + '\n')
