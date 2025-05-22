import gemmi
import sys

filepath = sys.argv[1]
cif_file = f'{filepath}/model.cif'
structure = gemmi.read_structure(cif_file)
structure.write_pdb(f'{filepath}/model.pdb')