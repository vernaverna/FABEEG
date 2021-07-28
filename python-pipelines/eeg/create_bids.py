#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:45:44 2019

@author: mhusberg
"""

import os
import os.path as op
import mne
from mne_bids import (write_raw_bids, BIDSPath, print_dir_tree)

from config_eeg import tasks

# Set paths
data_path = '/m/nbe/scratch/rubberboot/test/'
raw_path = op.join(data_path, 'data', 'MNE-eegbci-data/files/eegmmidb/1.0.0/')
output_path = op.join(data_path, 'bids')

# Subject name
subjects_in = ['S001','S002','S003','S004','S005']

# Tasks should match filenames
#tasks = ['restEO', 'restEC']

# Only one run
run = '01'

# Session 1
session = '01'

#Maxfiltered or not
proc = 'raw' # proc = 'tsss'

# Define events, CHECK TRIGGER CODES
rest_events = {}

def get_filenames(subject_in):
    
    #TODO: add functionality for selecting a subset of tasks
    
    # For this dataset
    file_endswith = ['R01.edf',  # eyes closed
                     'R02.edf']  # eyes open
    
    # Works for this dataset, make a function for each dataset
    filenames = ['{}{}'.format(subject_in, file_endswith[0]),
                 '{}{}'.format(subject_in, file_endswith[1])]
    
    return filenames

def get_subject_out(subject_in):
    
    #For this dataset subject_out is the same as subject_in
    subject_out=subject_in
    
    return subject_out

def create_BIDS(filename, subject_in, subject_out, session, run, proc, task, event_id):

    # Read in data
    raw_fname = op.join(raw_path, subject_in, filename)
    raw = mne.io.read_raw_edf(raw_fname, preload=False)
    
    # Specific for test eeg dataset
    raw.info['line_freq'] = 50  # specify power line frequency as required by BIDS
    raw.rename_channels(lambda x: x.strip('.'))  # remove dots from channel names, FIXME

    #events_data = mne.find_events(raw, min_duration=2/1000.)

    bids_path = BIDSPath(subject=subject_out, 
                         session=session, 
                         task=task, 
                         run=run, 
                         processing=proc, 
                         root=output_path)
    
    write_raw_bids(raw, 
                   bids_path, 
                   #events_data=events_data,
                   #event_id=event_id, 
                   overwrite=True)

    print_dir_tree(output_path)
    
    #change file permissions     
    #os.system('chmod -R ug+rwx ' + output_path + '/sub-' + subject_out )


# #######################################################################

for subject_in in subjects_in:    

    filenames = get_filenames(subject_in)
    subject_out = get_subject_out(subject_in)
    
    for filename, task in zip(filenames,tasks):

        print('Processing run: ' + run)
        create_BIDS(filename, subject_in, subject_out, session, run, proc, task, rest_events)
