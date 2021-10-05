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
#no_psds = []

try:
    subj_psds = fname.psds(subject=args.subject)

    f = h5py.File(subj_psds, "r") 
    list(f.keys()) #name of the dataset 
    group = f['h5io']
    freqs = np.array(group['key_freqs'])
    data_n1 = np.array(group['key_sleep N1'])
    data_n2 = np.array(group['key_sleep N2'])
    data_n2_2 = np.array(group['key_sleep N2 (2)'])
    data_n2_3 = np.array(group['key_sleep N2 (3)'])
    data_n2_4 = np.array(group['key_sleep N2 (4)'])
    data_n2_5 = np.array(group['key_sleep N2 (5)'])    
    
    #info = np.array(group['key_info']
    f.close()
    
    #Create a directory to save the .csv?? files
    parent_dir = "/projects/FABEEG/Data2R/absolute_spectra/"
    subj_dir = parent_dir + args.subject
    Path(subj_dir).mkdir(parents=True, exist_ok=True)
    
    # Calculate the average relative band power (RBP):
    # RBP = 100*ABP/sum(ABP)
    # i.e. sum of spectral power over all bands is 100 for each channel and subject
    # TODO: make sure that this is the desired approach!
    
    dataset = {'n1':data_n1, 'n2':data_n2, 'n2_2':data_n2_2, 'n2_3':data_n2_3,
               'n2_4':data_n2_4, 'n2_5':data_n2_5}
    
    for data_obj in list(dataset.keys()):
        data_bandpower = []
        
        for band in f_bands:
             fmin, fmax = band[0], band[1]
             bandpwr = bandpower(psd=dataset[data_obj], f=freqs, fmin=fmin, fmax=fmax)
             data_bandpower.append(bandpwr)
             
        #sums over the absolute power of all freq. bands; channel-wise info
        abs_power = np.sum(data_bandpower, axis=0)
        #TODO: handle possible divisions by zero! BETTER
        #normalizes data
        if relative:
            if abs_power > 0e-5:
                data_bandpower = np.array(data_bandpower)/abs_power #did not scale with 100 yet
        else:
            data_bandpower = np.array(data_bandpower)
        
        
    
        np.savetxt(subj_dir+'/'+data_obj+'.csv', data_bandpower, delimiter=',') #save spectra in separate files
   
    
except:
    print("Did not find psds -file!")
    #no_psds.append(args.subject)

#with open('corrupted_PSDs.txt', 'w') as f:
#    for bad_fname in no_psds:
#        f.write(bad_fname)
#        f.write('\n')
