#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate absolute / relative bandpower from PSDs.
@author: heikkiv
"""

import numpy as np
import argparse
import h5py 
import mne
from mne import find_layout, pick_info, pick_types
from pathlib import Path
from config_eeg import get_all_fnames, fname, bads, age_df


# Define the frequency bands (heuristic approach as in BRRR)
f_bands = [(1,3), (3,5.2), (5.2,7.6), (7.6,10.2), (10.2, 13), (13,16),
           (16,19.2), (19.2,22.6), (22.6,26.2), (26.2,30), (30,34), (34,38.2), (38.2,42.6)]

#TODO: change these!

# Helper function to calculate bandpower for each sensor
def bandpower(psd, f, fmin, fmax): 
    min_index = max(0, np.argmax(f > fmin) - 1)
    max_index = np.argmax(f > fmax) -1
    
    return np.trapz(psd[:,min_index: max_index], f[min_index: max_index], axis=1)


# Relative or absolute band power?
relative = True

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', help='The subject to process')
args = parser.parse_args()


#Create a directory to save the .csv files
if not relative:
    parent_dir = "/projects/FABEEG/Data2R/absolute_spectra_1min/"
else:
    parent_dir = "/projects/FABEEG/Data2R/relative_spectra_1min/"
subj_dir = parent_dir + args.subject
Path(subj_dir).mkdir(parents=True, exist_ok=True)

# Try to open psds; 
# if not successful, write the missing psd into a txt file

try:
    subj_psds = fname.psds(subject=args.subject)

    f = h5py.File(subj_psds, "r") 
    group = f['h5io']
    #info = np.array(group['key_info']
    freqs = np.array(group['key_freqs'])    
    keys = [key for key in group.keys()]
        
    #hold the data 
    dataset = {}
        
    for data_key in keys[:-2]: #loop through all psds
        data = np.array(group[data_key])
        seq_identifier = data_key.split(' ')[1] #to get N1a, N2b etc
        
        dataset[seq_identifier] = data
     
        
    f.close() #close file
 
    # Calculate the average relative band power (RBP):
    # RBP = 100*ABP/sum(ABP)
    # i.e. sum of spectral power over all bands is 100 for each channel and subject
    # TODO: make sure that this is the desired approach?

    for data_obj in list(dataset.keys()): #calculate bandpower for all PSDs 
        data_bandpower = []
        
        for band in f_bands:
             fmin, fmax = band[0], band[1]
             bandpwr = bandpower(psd=dataset[data_obj], f=freqs, fmin=fmin, fmax=fmax)
             data_bandpower.append(bandpwr)
             
        #sums over the absolute power of all freq. bands; channel-wise info
        abs_power = np.sum(data_bandpower, axis=0)
        #normalizes data
        if relative:
            if all (abs_power < -1e-5):
                data_bandpower = np.array(data_bandpower)/abs_power #did not scale with 100 yet
            else:
                print('Division by zero, did not calculate bandpower!')
        else:
            data_bandpower = np.array(data_bandpower)
        
        np.savetxt(subj_dir+'/'+data_obj+'.csv', data_bandpower, delimiter=',') #save spectra in separate files

    
        
    
except:
    print(f'Did not find psds -file from {args.subject}!')
    bad_fname = args.subject

    with open('corrupted_PSDs.txt', 'a') as f:
        f.write(bad_fname)
        f.write('\n')
