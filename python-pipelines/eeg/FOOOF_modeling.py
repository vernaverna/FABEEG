#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model PSDs using periodic and aperiodic components as given by FOOOF.
Good for group-lvl analysis.

Created on Wed Feb  1 10:34:17 2023

@author: heikkiv
"""

import numpy as np
import h5py 
import mne
from mne import find_layout, pick_info, pick_types
from pathlib import Path
from fooof import FOOOF
from config_eeg import get_all_fnames, fname, bads, age_df


#First, loop through every subject and get their data N1/N2 + freqs in an array


subjects = np.genfromtxt('EEG_subjects.txt', dtype=None,encoding="utf8")
spectra_n1a = []
spectra_n2a = []

individuals = []

for subj in subjects:
    
    try:
        subj_psds = fname.psds(subject=subj)
    
        f = h5py.File(subj_psds, "r") 
        list(f.keys()) #name of the dataset 
        group = f['h5io']
        freqs = np.array(group['key_freqs'])
        data_n1 = np.array(group['key_sleep N1 (1)'])
        data_n1_2 = np.array(group['key_sleep N1 (2)'])
        data_n2 = np.array(group['key_sleep N2 (1)'])
        data_n2_2 = np.array(group['key_sleep N2 (2)'])
        data_n2_3 = np.array(group['key_sleep N2 (3)'])
        data_n2_4 = np.array(group['key_sleep N2 (4)'])
        data_n2_5 = np.array(group['key_sleep N2 (5)'])    
        
        #info = np.array(group['key_info']
        f.close()
        
        spectra_n1a.append(data_n1)
        spectra_n2a.append(data_n2)
        individuals.append(subj)
        
    except:
        print(f'Did not find psds -file of {subj}!')


#TODO: what channels are we checking? Or are we interested just in GA?

#Calculate also global averages per each subject
spectra_n1a_GA = [np.mean(subj_dat, axis=0) for subj_dat in spectra_n1a]
spectra_n2a_GA = [np.mean(subj_dat, axis=0) for subj_dat in spectra_n2a]


# Trying FOOOF on one subject, different inputs
fm = FOOOF(max_n_peaks=8, aperiodic_mode='knee')
fm.print_settings(description=True)

freq_range = [2, 40]
#spectrum = spectra_n1a_GA[57]
#fm.report(freqs, spectrum, freq_range)

for i in range(0,19):
    spectrum = spectra_n2a[5][i,:]

    fm.report(freqs, spectrum, freq_range)


#Then, fit the FOOF model on a group of spectra?










