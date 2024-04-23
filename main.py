from code.classes.unitcell import Unitcell
from code.classes.forcefieldparams import ForceFieldParams

if __name__ == '__main__':
    try:
        # Ensure that 'chim222' is passed as a string
        unitcell = Unitcell('chim222')

        # Check if any atoms were read
        if not unitcell.atoms:
            print("No atoms were read from the file.")

        # Print out the elements and positions
        for atom in unitcell.atoms:
            print(f'Element: {atom.element}, Position: {atom.position}')

    except Exception as e:
        print(f"An error occurred: {e}")

