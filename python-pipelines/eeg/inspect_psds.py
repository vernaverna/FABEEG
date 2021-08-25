#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 12:33:10 2021

@author: heikkiv
"""

import numpy as np
import argparse
import h5py 
from mne import find_layout, pick_info, pick_types
from pathlib import Path
from config_eeg import get_all_fnames, fname, bads, age_df


# Define the frequency bands (heuristic approach as in BRRR)
f_bands = [(0,1), (1,3), (3,5.2), (5.2,7.6), (7.6,10.2), (10.2, 13), (13,16),
           (16,19.2), (19.2,22.6), (22.6,26.2), (26.2,30), (30,34), (34,38.2), (38.2,42.6)]

# Helper function to calculate bandpower for each sensor
def bandpower(psd, f, fmin, fmax): 
    min_index = np.argmax(f > fmin) - 1
    max_index = np.argmax(f > fmax) -1
    
    return np.trapz(psd[:,min_index: max_index], f[min_index: max_index], axis=1)


# Relative or absolute band power?
relative = False


# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', help='The subject to process')
args = parser.parse_args()


# Try to open psds; if not successful, write the missing psd into a txt file
no_psds = []

try:
    subj_psds = fname.psds(subject=args.subject)

    f = h5py.File(subj_psds, "r") 
    list(f.keys()) #name of the dataset 
    group = f['h5io']
    freqs = np.array(group['key_freqs'])
    data_n1 = np.array(group['key_sleep N1'])
    data_n2 = np.array(group['key_sleep N2'])
    #info = np.array(group['key_info']
    f.close()
    
    # Calculate the average relative band power (RBP):
    # RBP = 100*ABP/sum(ABP)
    # i.e. sum of spectral power over all bands is 100 for each channel and subject
    # TODO: make sure that this is the desired approach!
    n1_bandpower = []
    n2_bandpower = []
    
    for band in f_bands:
        fmin, fmax = band[0], band[1]
        bandpwr1 = bandpower(psd=data_n1, f=freqs, fmin=fmin, fmax=fmax)
        bandpwr2 = bandpower(psd=data_n2, f=freqs, fmin=fmin, fmax=fmax)
    
        n1_bandpower.append(bandpwr1)
        n2_bandpower.append(bandpwr2)
    
    
    #sums over the absolute power of all freq. bands; channel-wise info
    abs_power_n1 = np.sum(n1_bandpower, axis=0) 
    abs_power_n2 = np.sum(n2_bandpower, axis=0)
    
    
    #TODO: handle possible divisions by zero! BETTER
    #normalizes data
    if relative:
        if abs_power_n1 > 0e-5:
            n1_bandpower = np.array(n1_bandpower)/abs_power_n1 #did not scale with 100 yet
            n2_bandpower = np.array(n2_bandpower)/abs_power_n2 
    else:
        n1_bandpower = np.array(n1_bandpower)
        n2_bandpower = np.array(n2_bandpower)        
    
    #Create a directory to save the .csv?? files
    parent_dir = "/projects/FABEEG/Data2R/absolute_spectra/"
    subj_dir = parent_dir + args.subject
    Path(subj_dir).mkdir(parents=True, exist_ok=True)
    
    np.savetxt(subj_dir+'/n1.csv', n1_bandpower, delimiter=',') #save N1 & N2 sleep in separate files
    np.savetxt(subj_dir+'/n2.csv', n2_bandpower, delimiter=',')

except:
    print("Did not find psds -file!")
    no_psds.append(args.subject)

with open('corrupted_PSDs.txt', 'w') as f:
    for bad_fname in no_psds:
        f.write(bad_fname)
        f.write('\n')
