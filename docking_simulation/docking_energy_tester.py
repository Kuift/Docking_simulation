#go to wsl
#conda activate docking

#python3 dock.py
#pymol cbdas_substrate.sdf
# import matplotlib.pyplot as plt
# import seaborn as sns
# import plotly.express as px

# import numpy as np
# import pandas as pd
# import torch

# from sklearn.linear_model import RidgeCV
# from sklearn.manifold import TSNE

# from protein_design.constants import AA, device
# from protein_design.data import to_tensor
# from protein_design.splitter import random_split
# from protein_design.sequence import integer_to_seqs, seqs_to_onehot, read_fasta
# from protein_design.discriminative import MLP, SequenceGP
# from protein_design.trainer import train
# from protein_design.evaluator import regression_metrics

# fname = "/content/2Q8A_binding.txt"

# df = pd.read_csv(fname, sep='\t', skiprows=1)
# print(df)
# df_filtered = df.loc[df["Best?"] == True]
# df_filtered = df_filtered.drop_duplicates("Slide")
# SEQ_NB = 24
# seqs = df_filtered["Slide"].to_list()
# #X = seqs_to_onehot(seqs, flatten=True)

# y = energies = df_filtered["Energy"].to_numpy()



# plt.rcParams["figure.dpi"] = 100

# sns.kdeplot(x=energies, fill=True, label="2Q8A")
# plt.xlabel("Energy")
# plt.legend()
# plt.show()



# violin plot
# rmsd plot
from pathlib import Path
from typing import Optional
import deepchem as dc
import numpy as np
import typer
import matplotlib.pyplot as plt
import os
import time

vpg = dc.dock.VinaPoseGenerator(pocket_finder=None)

def dock(out: Path = "results", centroid: Optional[str] = None, box: Optional[str] = None, exhaustiveness: int = 40, num_modes: int = 40):
    """Dock a ligand to a protein using the deepchem interface to AutoDock Vina

    Parameters
    ----------
    protein
        File path of target protein structure in .pdb format (receptor)
    ligand
        File path of ligand to dock in .sdf format
    out
        Output folder for docked structures and intermediate files
    centroid, optional
        Center of search space for docking, by default None
    box, optional
        Dimensions of a cubic box defining the search space, by default None
    exhaustiveness, optional
        Intensity of the search for docked poses, by default 1
    num_modes, optional
        Number of docked structures, by default 20
    """
    centroidArray = [] 
    centroidArray.append(np.array([  35.763,  24.716, 292.937]))
    centroidArray.append(np.array([  43.898,  27.552, 318.184]))
    centroidArray.append(np.array([  34.676,  -6.526, 267.112]))
    centroidArray.append(np.array([  30.849,  48.148, 291.007]))
    centroidArray.append(np.array([  38.553,  34.443, 292.315]))
    centroidArray.append(np.array([  26.670,  43.479, 298.391]))
    centroidArray.append(np.array([  35.913,  45.738, 300.592]))
    centroidArray.append(np.array([  19.494,  20.135, 277.747]))

    if box is not None:
        x, y, z = tuple(box.split(","))
        box = np.array([float(x), float(y), float(z)])
    else:
        box = np.array([20.0, 20.0, 20.0])
    structures = []
    substrates = []
    for file in os.listdir():
        if file.endswith(".pdb"):
            structures.append(str(file))
    for file in os.listdir():
        if file.endswith(".sdf"):
            substrates.append(str(file))
    substrates_scores = []
    one_substrate_scores = [] 
    to_be_averaged_scores = []

    for structure in structures:
        for substrate in substrates:
            for i in range(110):
                print("now working on substrate named : " + substrate)
                poses, scores = vpg.generate_poses(
                    (str(structure), str(substrate)),
                    centroid=centroidArray[1],
                    box_dims=box,
                    exhaustiveness=exhaustiveness,
                    num_modes=num_modes,
                    out_dir=str(out),
                    generate_scores=True,
                )
                print(f"SCORES : {scores}\n CENTROID:{1}")
                print(f"iterations: {i}")
                one_substrate_scores.append(scores)
                to_be_averaged_scores.append(scores[0])
            f = open(f"{str(substrate)}-INDEX{1}-{str(structure)}-{exhaustiveness}-{num_modes}-{time.time()}.txt", "a")
            f.write(f"{sum(to_be_averaged_scores)/len(to_be_averaged_scores)}\n")
            for element in one_substrate_scores:
                f.write(f"{element}\n")
            f.close() 
            substrates_scores.append(sum(to_be_averaged_scores)/len(to_be_averaged_scores))
            one_substrate_scores = []
            to_be_averaged_scores = []
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,1])
        ax.bar(substrates,substrates_scores)
        plt.xticks(rotation='vertical')
        plt.savefig(f"{str(structure)}-index-{1}.jpg", bbox_inches='tight')
        substrates_scores.clear()

dock()