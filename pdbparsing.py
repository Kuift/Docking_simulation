#get the atoms coordinates
#calculate rmsd values
#get the poses of the atoms from vina
#retrieves atoms from the poses
#calculate rmsd values from it
#reproduce vina's rmsd values
from Bio.PDB import PDBParser
import deepchem as dc
import numpy as np

import struct
p = PDBParser()
s = p.get_structure("5tgz", "5tgz.pdb")                    

# for chains in s:
#     for chain in chains:
#         for residue in chain:                             
#             for atom in residue:
#                 print(f"{atom.get_coord()}{atom.get_name()}")

vpg = dc.dock.VinaPoseGenerator(pocket_finder=None)

poses, scores = vpg.generate_poses(
    (str("5tgz.pdb"), str("THCA6.sdf")),
    centroid=None,
    box_dims=np.array([20.0, 20.0, 20.0]),
    exhaustiveness=15,
    num_modes=25,
    out_dir=str("temp"),
    generate_scores=True,
)

print(f"{type(poses[0][0])}")
print(f"{dir(poses[0][0])}")
print(f"{type(list(poses[0][0].GetAtoms())[0])}")
print(f"{dir(list(poses[0][0].GetAtoms())[0])}") 
print(f"{list(poses[0][0].GetAtoms())[3].GetAtomicNum()}")
#if you don't get same rdms value, might be because you don't get only the m_num_movable_atoms
for tuple in poses:
    for mol in tuple:
        atoms_list = list(mol.GetAtoms())
        for conformer in mol.GetConformers():
            for i, position in enumerate(conformer.GetPositions()):
                print(f"{position}{atoms_list[i].GetAtomicNum()}")
        # for atom in mol.GetAtoms():
        #     #print(f"{atom.get_coord()}{atom.get_name()}")
        # break
    break