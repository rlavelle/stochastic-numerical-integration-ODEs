import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import interp1d
from scipy.integrate import odeint
import math
import sys
import time
import random


# protein names
names = ['clb1','clb3','cdc20T','cdc20','clb4','sp','cdc5t','cdc5a','ndd1t','ndd1a','hcm1','ndt80','sum1iIme2','sum1iCdk1','sum1iRC','ama1p','rc','dsb','ama1t']
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

"""
run program with:
    python3 random_plots.py trials start_time end_time protein-folder-path trial-label
"""
if __name__ == "__main__":
    if(len(sys.argv) != 6):
        raise(Exception("Error: expected a number of trials, proteins folder path and trial number"))

    analysis_image_path = '/N/u/rowlavel/Carbonate/stochastic-numerical-integration-ODEs/analysis_images'
    trialcount = int(sys.argv[1])
    start_time = int(sys.argv[2])
    end_time = int(sys.argv[3])
    protiens_file_path = sys.argv[4]
    trial_label = sys.argv[5]
    proteincount = 19
    sig_proteincount = 7

    proteins = [read_proteins(f'{protiens_file_path}/proteins-trial-{str(i)}.csv') for i in range(0,trialcount)]

    t_data = np.arange(start_time, end_time, 0.001)    

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

    for j in range(0, sig_proteincount):
        plt.figure()
        plt.errorbar(t_data[::50], sig_mean[j][::50], sig_error[j][::50])
        plt.plot(t_data, sig_mean[j])
        try:
            plt.title(label=f'protien {names[PROTEINS_SIG[j]]} mean and error bars {trial_label}')
            plt.savefig(f'{analysis_image_path}/single-protien-{names[PROTEINS_SIG[j]]}-{trial_label}.png')
        except:
            plt.title(label=f'protien ama1t-ama1p mean and error bars {trial_label}')
            plt.savefig(f'{analysis_image_path}/single-protien-ama1t-ama1p-{trial_label}.png')