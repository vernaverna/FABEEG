#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 11:39:08 2023

@author: heikkiv
"""

import numpy as np
import pandas as pd
import mne
from pathlib import Path
from config_eeg import get_all_fnames, fname, bads, age_df, all_subjects


EMG_rms_dict = {}

for subj in all_subjects:
    subj = subj.replace('sub-', '')
    try:
        emg_evokeds = mne.read_evokeds(fname.emg_evoked(subject=subj))
        
        subj_info = {}
        for obj in emg_evokeds:
            
            description = obj.comment
            description = description.replace('PSD', 'RMS')
            
            EMG_data = obj.get_data().mean(axis=0)
            RMS_value = EMG_data.mean(axis=-1)*1e6 #data is in volts => convert to microvolts
            subj_info[description] = RMS_value
        
        
        EMG_rms_dict[subj]=subj_info
    
    except:
        print(f'No emg data for subject {subj}')
    
emg_dataframe = pd.DataFrame.from_dict(EMG_rms_dict, orient='index')

emg_dataframe.to_csv('emg_RMS_values.csv')