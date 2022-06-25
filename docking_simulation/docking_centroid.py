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

from pathlib import Path
from typing import Optional
import deepchem as dc
import numpy as np
import typer
import matplotlib.pyplot as plt
import os
import time

vpg = dc.dock.VinaPoseGenerator(pocket_finder=None)

def dock(out: Path = "results", centroid: Optional[str] = None, box: Optional[str] = None, exhaustiveness: int = 15, num_modes: int = 25):
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
    CsAP2L1_centroid = [] 
    CsAP2L1_centroid.append(np.array([  10.720,  15.910, -13.490]))
    CsAP2L1_centroid.append(np.array([  -14.993, -34.156, -25.767]))
    CsAP2L1_centroid.append(np.array([  -8.482,  -7.515,   1.325]))
    CsAP2L1_centroid.append(np.array([  -0.961,  -4.373,  -5.470]))
    CsAP2L1_centroid.append(np.array([  -1.852,  -2.086,   1.031]))
    CsAP2L1_centroid.append(np.array([  -4.599,   5.153,  -5.851]))
    CsAP2L1_centroid.append(np.array([  -3.413,   2.465,  -2.817]))
    CsAP2L1_centroid.append(np.array([  -11.478,   7.421,  -8.557]))
    CsAP2L1_centroid.append(np.array([  -9.342,  12.249,  -4.366]))
    CsAP2L1_centroid.append(np.array([  -0.622,   6.913,   2.897]))
    CsAP2L1_centroid.append(np.array([  1.863,   6.028,  22.018]))
    CsAP2L1_centroid.append(np.array([  18.586, -10.067,  -1.520]))

    CsWRKY1_centroid = []
    #sWRKY1_centroid.append(np.array([  -33.193, -22.610,  31.789]))
    #CsWRKY1_centroid.append(np.array([  -1.918,  36.886,   8.255]))
    #CsWRKY1_centroid.append(np.array([  -3.999,  38.820, -41.432]))
    #CsWRKY1_centroid.append(np.array([  -20.246,  16.499, -55.063]))
    #CsWRKY1_centroid.append(np.array([  35.626,  -6.674, -54.139]))
    #CsWRKY1_centroid.append(np.array([  8.323, -40.449, -20.743]))
    #CsWRKY1_centroid.append(np.array([  -6.586, -29.827, -46.453]))
    #CsWRKY1_centroid.append(np.array([  -14.427,   3.870, -23.045]))
    #CsWRKY1_centroid.append(np.array([  -16.173,  14.154,  -7.521]))
    #CsWRKY1_centroid.append(np.array([  -14.470,  13.710,  -6.267]))
    #CsWRKY1_centroid.append(np.array([  -11.585,   8.348,  -2.961]))
    #CsWRKY1_centroid.append(np.array([  -1.390,   4.242,  -4.586]))
    #CsWRKY1_centroid.append(np.array([  2.562,  -5.871,  -2.014]))
    #CsWRKY1_centroid.append(np.array([  -0.166,  -5.049,  -1.808]))
    #CsWRKY1_centroid.append(np.array([  -0.886,   1.843,  -2.969]))
    #CsWRKY1_centroid.append(np.array([  13.147,  -2.379,  11.962]))
    #CsWRKY1_centroid.append(np.array([  32.755,   4.842,  18.688]))
    #CsWRKY1_centroid.append(np.array([  6.606,  17.026,  28.291]))
    #CsWRKY1_centroid.append(np.array([  15.443,   4.722,  68.874]))
    #CsWRKY1_centroid.append(np.array([  36.185, -28.286,  15.014]))
    #CsWRKY1_centroid.append(np.array([  27.414,  12.215,  -8.369]))

    CsMYB1_centroid = []


    centroids = {}
    centroids['CsAP2L1'] = CsAP2L1_centroid
    centroids['CsWRKY1'] = CsWRKY1_centroid
    centroids['CsMYB1'] = CsMYB1_centroid 


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
        structure_name = os.path.splitext(structure)[0]
        for cent in centroids[structure_name]:
            for substrate in substrates:
                for i in range(110):
                    print("now working on substrate named : " + substrate)
                    poses, scores = vpg.generate_poses(
                        (str(structure), str(substrate)),
                        centroid=cent,
                        box_dims=box,
                        exhaustiveness=exhaustiveness,
                        num_modes=num_modes,
                        out_dir=str(out),
                        generate_scores=True,
                    )
                    print(f"SCORES : {scores}\n CENTROID:{j}")
                    print(f"iterations: {i}")
                    one_substrate_scores.append(scores)
                    to_be_averaged_scores.append(scores[0])
                f = open(f"{str(structure)}-{str(substrate)}-{str(exhaustiveness)}-{str(num_modes)}-INDEX{j}-{time.time()}.txt", "a")
                f.write(f"{sum(to_be_averaged_scores)/len(to_be_averaged_scores)}\n")
                for element in one_substrate_scores:
                    f.write(f"{element}\n")
                f.close() 
                substrates_scores.append(sum(to_be_averaged_scores)/len(to_be_averaged_scores))
                one_substrate_scores = []
                to_be_averaged_scores = []
            fig = plt.figure()
            ax = fig.add_axes([0,0,1,1])
            ax.bar(structures,substrates_scores)
            plt.xticks(rotation='vertical')
            plt.savefig(f"{str(structure)}-index-{j}.jpg", bbox_inches='tight')
            substrates_scores.clear()

dock()