import threading
import time
import sys
from code.classes.SystemWriter import SystemWriter
from code.classes.ChosenSettings import ChosenSettings
import argparse

from code.classes.ClayBuilder import ClayBuilder

def main(replication, 
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
    
    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument('replication_factor', type=int, nargs='?', default=None, help='Replication factor for x and y dimensions')
    parser.add_argument('-o', '--output_file', type=str, help='Output file name')
    args = parser.parse_args()
    
    if args.replication_factor is not None:
        replication_factor = args.replication_factor
        replication = (replication_factor, replication_factor, 1)
    
    output_file = args.output_file if args.output_file else output_file

    def loading_animation(min_display_time=3):
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
    

    def print_banner(file_path):
        try:
            with open(file_path, 'r') as file:
                content = file.read()
                print(content)
        except FileNotFoundError:
            print("The file was not found.")
        except Exception as e:
            print(f"An error occurred: {e}")

    print_banner('misc/banner.txt')

    run_animation = threading.Event()
    if replication[0] >= 4 and replication[1] >= 4:
        run_animation.set()
        thread = threading.Thread(target=loading_animation)
        thread.start()

    try:
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

        final_text = (
        f"Finished writing .data file for {replication[0]}x{replication[1]}x{replication[2]} system.\n\n"
        f"○ Output file:    {output_file}.data\n"
        f"○ Water per ion:  {water_per_ion}\n"
        f"○ System charge:  {runner.system_charge:.3f}\n"
        f"○ Ions added:     {runner.ion_count}\n"
        f"○ Final charge:   {runner.final_charge:.3f}"
        )


        if replication[0] < 4 and replication[1] < 4:
            sys.stdout.write('\r' + final_text)
    finally:
        if run_animation.is_set():
            run_animation.clear()
            thread.join()