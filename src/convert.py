import gemmi
import sys
import os

filepath = sys.argv[1]
cif_file = os.path.join(filepath, 'model.cif')
structure = gemmi.read_structure(cif_file)
structure.write_pdb(os.path.join(filepath, 'model.pdb'))