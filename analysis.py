import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import interp1d
from scipy.integrate import odeint
import math
import sys

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
    for j in range(0, proteincount):
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

    y_err = np.sqrt( ( ((proteins-mean)**2).sum(axis=0) )/(trialcount-1) )

    return y_err

if __name__ == "__main__":
    if(len(sys.argv) != 2):
        raise(Exception("Error: expected a number of trials"))

    trialcount = int(sys.argv[1])
    runtime = 660 
    proteincount = 19

    proteins = [read_proteins('protein-trials/proteins-trial-'+str(i)+'.csv') for i in range(0,trialcount)]

    t_data = np.arange(0, runtime, 0.001)

    # calculate mean over trials
    mean = calc_mean(proteins)

    # graph and save mean
    for j in range(0, proteincount):
        plt.plot(t_data, mean[j])
    plt.title(label='Mean Through '+str(trialcount)+' Trials')
    plt.savefig('analysis_images/mean.png')

    # calculate error over trials
    y_err = calc_error(proteins,mean)

    # graph and save plot with error bars
    for j in range(0, proteincount):
        plt.errorbar(t_data, mean[j], yerr = y_err[j])
        plt.plot(t_data, mean[j])
    plt.savefig('analysis_images/error_bars.png')