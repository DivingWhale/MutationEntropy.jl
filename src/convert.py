import gemmi
import sys
import os

if len(sys.argv) < 2:
    print("Usage: python convert.py <cif_file_path> [output_pdb_path]")
    sys.exit(1)

cif_file = sys.argv[1]
if len(sys.argv) >= 3:
    output_pdb = sys.argv[2]
else:
    # Default: same directory and name as CIF but with .pdb extension
    output_pdb = os.path.splitext(cif_file)[0] + '.pdb'
    # If output has prefix, use model.pdb instead
    if '_model.cif' in cif_file:
        output_pdb = os.path.join(os.path.dirname(cif_file), 'model.pdb')
    elif cif_file.endswith('model.cif'):
        output_pdb = os.path.join(os.path.dirname(cif_file), 'model.pdb')

structure = gemmi.read_structure(cif_file)
structure.write_pdb(output_pdb)