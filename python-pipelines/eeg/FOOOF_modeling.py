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

for subj in subjects:
    
    try:
        subj_psds = fname.psds(subject=args.subject)
    
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






#Then, fit the FOOF model on a group level?
# i prolly want to do some ind lvl too, just to check.


















