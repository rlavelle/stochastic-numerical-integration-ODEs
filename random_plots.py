import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import interp1d
from scipy.integrate import odeint
import math
import sys
import time
import random

# significant proteins
PROTEINS_SIG = [1,0,6,4,11,15]

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
    for j in range(0, sig_proteincount):
        p_avg = []
        for i in range(0, trialcount):
            p_avg.append(proteins[i][j])
        p_avg = np.asarray(p_avg)
        p_avg = np.transpose(p_avg)
        p_avg = p_avg.mean(axis=1)
        mean_array.append(p_avg)
    
    return mean_array

def calc_error(proteins,mean):
    proteins = np.asarray(proteins)
    mean = np.asarray(mean)

    return np.sqrt((((proteins-mean)**2).sum(axis=0))/(trialcount-1))

if __name__ == "__main__":
    if(len(sys.argv) != 4):
        raise(Exception("Error: expected a number of trials, proteins folder path and trial number"))

    trialcount = int(sys.argv[1])
    protiens_file_path = sys.argv[2]
    trial_num = sys.argv[3]
    runtime = 660 
    proteincount = 19
    sig_proteincount = 7

    proteins = [read_proteins(f'{protiens_file_path}/proteins-trial-{str(i)}.csv') for i in range(0,trialcount)]

    t_data = np.arange(0, runtime, 0.001)    

    # create significant proteins list
    sig_proteins = []
    for i in range(0,trialcount):
        proteins_sig = []
        for j in PROTEINS_SIG:
            proteins_sig.append(proteins[i][j])
        special = []
        for j in range(0,len(t_data)):
            special.append(proteins[i][18][j]-proteins[i][15][j])
        proteins_sig.append(special)
        sig_proteins.append(proteins_sig)

    # randomly sample 8 trials
    #sample = random.sample(sig_proteins, 8)
    sig_mean = calc_mean(sig_proteins)
    sig_error = calc_error(sig_proteins, sig_mean)

    # graph and save means of protiens
    # for j in range(0, sig_proteincount):
    #     for i in range(0, len(sample)):
    #         plt.plot(t_data, sample[i][j])
    # plt.title(label=f'{trial_num} 8 Stochastic samples')
    # plt.savefig(f'/N/u/rowlavel/Carbonate/stochastic-numerical-integration-ODEs/analysis_images/sample-{trial_num}.png')

    for j in range(0, sig_proteincount):
        plt.errorbar(t_data[::50], sig_mean[j][::50], sig_error[j][::50])
        plt.plot(t_data, sig_mean[j])
        plt.savefig(f'/N/u/rowlavel/Carbonate/stochastic-numerical-integration-ODEs/analysis_images/singe-protien-{j}-{trial_num}.png')
    