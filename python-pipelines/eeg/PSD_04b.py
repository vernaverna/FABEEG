#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 11:39:08 2023

@author: heikkiv
"""

import numpy as np
import argparse
import h5py 
import mne
from mne import find_layout, pick_info, pick_types
from pathlib import Path
from config_eeg import get_all_fnames, fname, bads, age_df


# Relative or absolute psd?
relative = True

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', help='The subject to process')
args = parser.parse_args()


#Create a directory to save the .csv files
if not relative:
    parent_dir = "/projects/FABEEG/Data2R/absolute_PSD_1min/"
else:
    parent_dir = "/projects/FABEEG/Data2R/relative_PSD_1min/"
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
        data_power = dataset[data_obj]
        
        #sums over the absolute power of all freq. bands; channel-wise info
        abs_power = np.sum(data_power, axis=0)
        #normalizes data
        if relative:
            if all (abs_power < -1e-5):
                data_power = np.array(data_power)/abs_power #did not scale with 100 yet
            else:
                print('Division by zero, did not calculate bandpower!')
        else:
            data_power = np.array(data_power)
        
        np.savetxt(subj_dir+'/'+data_obj+'.csv', data_power, delimiter=',') #save spectra in separate files

    
        
    
except:
    print(f'Did not find psds -file from {args.subject}!')
    bad_fname = args.subject

    with open('corrupted_PSDs.txt', 'a') as f:
        f.write(bad_fname)
        f.write('\n')
