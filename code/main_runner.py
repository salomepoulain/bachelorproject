import threading
import time
import sys
from code.classes.SystemWriter import SystemWriter

def main(input_file,
         replication,
         height,
         al_mg_ratio,
         ca_si_ratio,
         water_file,
         vdw_radii_file,
         vdw_def_radius,
         vdw_def_scale,
         ff_files,
         mg_cutoff,
         h_cutoff,
         bond_cutoff,
         ff_params):

    if len(sys.argv) > 1:
        replication_factor = int(sys.argv[1])
        replication = (replication_factor, replication_factor)
    
    final_text = f"Finished writing .data file for {replication[0]}x{replication[1]} system."
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
        system = SystemWriter(input_file,
                                replication,
                                height,
                                al_mg_ratio,
                                ca_si_ratio,
                                water_file,
                                vdw_radii_file,
                                vdw_def_radius,
                                vdw_def_scale,
                                ff_files,
                                mg_cutoff,
                                h_cutoff,
                                bond_cutoff,
                                ff_params)

        if replication[0] < 4 and replication[1] < 4:
            sys.stdout.write('\r' + final_text + '\n\n')
    finally:
        # Stop the loading animation if it was started
        if run_animation.is_set():
            run_animation.clear()
            thread.join()