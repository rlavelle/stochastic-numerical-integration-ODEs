import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import interp1d
from scipy.integrate import odeint
import math
import sys
import time

# significant proteins
PROTEINS_SIG = [1,0,6,4,11,15]

"""
Helper functions for analysis
"""
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
    python3 analysis.py trials protein-folder-path trial-label
"""
if __name__ == "__main__":
    if(len(sys.argv) != 4):
        raise(Exception("Error: expected a number of trials, proteins folder path and trial number"))

    analysis_image_path = '/home/cvanoeve/stochastic-numerical-integration-ODEs/analysis_images'
    trialcount = int(sys.argv[1])
    protiens_file_path = sys.argv[2]
    trial_label = sys.argv[3]
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

    # calculate mean over trials
    mean = calc_mean(sig_proteins)

    # graph and save means of protiens
    for j in range(0, sig_proteincount):
        plt.plot(t_data, mean[j])
    plt.title(label=f'{trial_label} Mean Through {str(trialcount)} Trials')
    plt.savefig(f'{analysis_image_path}/mean-{trial_label}.png')

    # calculate error over trials
    y_err = calc_error(sig_proteins,mean)

    # graph and save plot with error bars
    for j in range(0, sig_proteincount):
        plt.errorbar(t_data[::50], mean[j][::50], y_err[j][::50])
        plt.plot(t_data, mean[j])
    plt.title(label=f'{trial_label} Error Bars Over Mean Through {str(trialcount)} Trials')
    plt.savefig(f'{analysis_image_path}/error_bars-{trial_label}.png')
