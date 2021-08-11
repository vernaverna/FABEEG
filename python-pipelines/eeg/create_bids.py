#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:45:44 2019

@author: mhusberg
"""

import os
import os.path as op
import mne
import pandas as pd
from mne_bids import (write_raw_bids, BIDSPath, print_dir_tree)

#from config_eeg import tasks

# Set paths
data_path = '/net/theta/fishpool/projects/FABEEG/childEEG_data/'
raw_path = '/net/theta/fishpool/projects/FABEEG/FABEEG/'
output_path = op.join(data_path, 'bids')

# Subject name
subj_file = op.join(data_path, 'subj_table.txt')
subjects_in = pd.read_csv(subj_file, header=None)
subjects_in = list(subjects_in.iloc[:,0])


subjects_in = [x.strip('.edf') for x in subjects_in] #remove file ending

corrupted = [] #list for corrupted .edf files
channels = [] #list for channels per subject
bids_subjs = [] #list for subjects that have BIDS created

# Define events, CHECK TRIGGER CODES
rest_events = {}

def get_filenames(subject_in): #Write new
    
    #TODO: add functionality for selecting a subset of tasks
    
    # For this dataset
    file_endswith = '.edf'
    
    # Works for this dataset, make a function for each dataset
    filenames = ['{}{}'.format(subject_in, file_endswith)]

    
    return filenames

def get_subject_out(subject_in):
    
    #For this dataset subject_out is NOT the same as subject_in
    subject_out=subject_in.replace('-', '')
    
    return subject_out

def create_BIDS(filename, subject_in, subject_out, event_id):

    # Read in data
    raw_fname = op.join(raw_path, filename)
    
    try:
        raw = mne.io.read_raw_edf(raw_fname, preload=False)
        raw.info['line_freq'] = 50  # specify power line frequency as required by BIDS
        #raw.rename_channels(lambda x: x.strip('.'))  # remove dots from channel names, FIXME

        #events_data = mne.find_events(raw, min_duration=2/1000.)

        bids_path = BIDSPath(subject=subject_out,  
                         root=output_path)
    
        write_raw_bids(raw, 
                       bids_path, 
                       #events_data=events_data,
                       #event_id=event_id, 
                       overwrite=True)

        print_dir_tree(output_path)
        bids_subjs.append(filename)
    except:
        print("Could not read {}".format(filename))
        corrupted.append(filename)
    
    #change file permissions     
    #os.system('chmod -R ug+rwx ' + output_path + '/sub-' + subject_out )


# #######################################################################

for subject_in in subjects_in:    

    filenames = get_filenames(subject_in)
    subject_out = get_subject_out(subject_in)
    print("BIDSifying file")
    print(subject_in)
    
    for filename in filenames:

        #print('Processing run: ' + run)
        create_BIDS(filename, subject_in, subject_out, rest_events)

with open('eeg_subjects.txt', 'w') as f:
    for good_fname in bids_subjs:
        f.write(get_subject_out(good_fname.strip('.edf')))
        f.write('\n')

with open('corrupted_EEGs.txt', 'w') as f:
    for bad_fname in corrupted:
        f.write(bad_fname)
        f.write('\n')