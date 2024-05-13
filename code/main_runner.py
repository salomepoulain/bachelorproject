import threading
import time
import sys
from code.classes.SystemWriter import SystemWriter
from code.classes.ChosenSettings import ChosenSettings
from code.classes.SystemAllocator import SystemAllocator
from code.classes.SolventIonAdder import SolventIonAdder
from code.classes.VerticalDuplicator import VerticalDuplicator

from code.classes.ClayBuilder import ClayBuilder

def main(replication, 
        al_mg_ratio, 
        net_charge, 
        water_per_ion, 
        ff_atom_types, 
        water_distance, 
        mg_cutoff, 
        h_cutoff, 
        bond_cutoff, 
        input_file, 
        ff_files, 
        water_file):
    
    if len(sys.argv) > 1:
        replication_factor = int(sys.argv[1])
        replication = (replication_factor, replication_factor, 1)

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
    final_text = f"Finished writing .data file for {replication[0]}x{replication[1]}x{replication[2]} system."

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
                                ff_atom_types, 
                                water_distance, 
                                mg_cutoff, 
                                h_cutoff, 
                                bond_cutoff, 
                                input_file, 
                                ff_files, 
                                water_file)
        
        # vertical = VerticalDuplicator(settings)

        #allocator = SystemAllocator(settings)
        
        # clay = ClayBuilder(settings)
        
        runner = SystemWriter(settings)

        if replication[0] < 4 and replication[1] < 4:
            sys.stdout.write('\r' + final_text + '\n\n')
    finally:
        if run_animation.is_set():
            run_animation.clear()
            thread.join()