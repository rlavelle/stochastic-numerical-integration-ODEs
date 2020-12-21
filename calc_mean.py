import numpy as np
import matplotlib.pyplot as plt
import csv
import math
import sys
import time

def read_proteins(file):
    proteins = []
    f = open(file,'r')
    with f:
        reader = csv.reader(f,quotechar='"')
        for row in reader:
            l = []
            for r in row: l.append(float(r))
            proteins.append(l)
    f.close()
    return proteins

def calc_mean(proteins):
    mean_array = []
    for j in range(0, proteincount):
        p_avg = []
        for i in range(0, trialcount):
            p_avg.append(proteins[i][j])
        p_avg = np.asarray(p_avg)
        p_avg = np.transpose(p_avg)
        p_avg = p_avg.mean(axis=1)
        mean_array.append(p_avg)
    
    return mean_array

"""
run program with:
    python3 calc_mean.py trials protein-folder-path time output-file
"""
if __name__ == "__main__":
    if(len(sys.argv) != 5):
        raise(Exception("Error: expected trials protein-folder-path time output-file"))

    trialcount = int(sys.argv[1])
    protiens_file_path = sys.argv[2]
    t = int(float(int(sys.argv[3])/0.001))
    output_file = sys.argv[4]
    proteincount = 19

    proteins = [read_proteins(f'{protiens_file_path}/proteins-trial-{str(i)}.csv') for i in range(0,trialcount)]

    # calculate mean over trials
    mean = calc_mean(proteins)
    
    means = []
    for j in range(0, proteincount):
        means.append(mean[j][t])
    # write mean to a file
    with open(output_file, 'w') as file:
        file.writelines([f'{m}\n' for m in means])


