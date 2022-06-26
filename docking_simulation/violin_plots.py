import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd


substrate_binding_raw_data_files = []
substrates_name = []
substrate_data = []
for file in os.listdir():
    if file.endswith(".txt"):
        substrate_binding_raw_data_files.append(str(file))
permlist= []
for file in substrate_binding_raw_data_files:
    templist = []
    if not (file in permlist):
        templist.append(file)
    else:
        continue
    removelist = []
    permlist.append(file) 
    for file2 in substrate_binding_raw_data_files:
        if str(file).split(".")[0] == str(file2).split(".")[0] and str(file) != str(file2) and not (file2 in permlist):
            templist.append(file2)
            permlist.append(file2)

    temp_array = []
    substrates_name.append(file.split(".")[0])
    for file in templist:
        f = open(file, "r")
        i = 0
        for line in f:
            if i == 0:
                i+=1
                continue
            # for number in line.strip()[1:-1].split(",")[:]:
            #     temp_array.append(number)
            temp_array.append(float(line[1:].split(",")[0]))
    substrate_data.append(temp_array)
 
# for file in substrate_binding_raw_data_files:
    # f = open(file, "r")
    # substrates_name.append(file.split(".")[0])
    # i = 0
    # temp_array = []
    # for line in f:
        # if i == 0:
            # i+=1
            # continue
        # result=float(line[1:].split(",")[0])
        # temp_array.append(float(line[1:].split(",")[0]))
    # substrate_data.append(temp_array)


sns.set(rc={'figure.figsize':(20,10)})
ax = sns.violinplot(data=substrate_data)
ax.set_xticklabels(substrates_name)
plt.title("5u09.pdb - CB1 receptor binding affinity of multiple substrates")
plt.xticks(rotation='vertical')
plt.xlabel("Substrates")
plt.ylabel("Energy")
plt.savefig("_violin_plot.jpg", bbox_inches='tight')