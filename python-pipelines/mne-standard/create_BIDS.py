#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:45:44 2019

Simple example for how to use mne_bids to save your data in BIDS format

Assumes that the raw data is saved as /m/myproject/meg/subject_name/YYMMDD/rawfile_for_mytask_raw.fif
Change according to your needs

@author: mhusberg
"""

import os.path as op
import mne
from mne_bids import write_raw_bids, read_raw_bids, make_bids_basename
from mne_bids.utils import print_dir_tree


data_path = '/m/myproject/'

subject_in = 'subject_name'

# Ouput file
output_path = op.join(data_path, 'BIDS')
subject_out = '01' # subject-id: '01' -> 'sub-01' in BIDS structure


# Define event IDs according to your project
events = {'Category1':1, 'Category2': 2, 'Other event': 3}
rest_events = {}

def create_BIDS(filename, subject_in, subject_out, date, session, run, proc, task, event_id):
# Read in data
    raw_fname = op.join(data_path, 'meg', subject_in, date, filename) # Folder structure for recorded raw data
    raw = mne.io.read_raw_fif(raw_fname)
    events_data = mne.find_events(raw, min_duration=2/1000.)


    bids_basename = make_bids_basename(subject=subject_out, session=session,
                                   task=task, processing=proc,run=run)
    write_raw_bids(raw, bids_basename, output_path, events_data=events_data,
               event_id=event_id, overwrite=True)

    print_dir_tree(output_path)


#######################################################################
# Session 1
date = 'YYMMDD' # date when measured, if in folder structure
session = '01'

#Maxfiltered or not?
proc = 'raw' #'tsssmctrans'


###########
# Task: mytask
task = 'mytask'

# Filename, run 1
filename='rawfile_for_mytask_run1_raw.fif' # Fill in name of rawdata-file 1
run = '01'

create_BIDS(filename, subject_in, subject_out, date, session, run, proc, task, events)

# Filename, run 2
filename='rawfile_for_mytask_run2_raw.fif' # Fill in name of rawdata-file 2
run = '02'

create_BIDS(filename, subject_in, subject_out, date, session, run, proc, task, events)


############
# Task: rest
task='restEO'

# Filename, run 1
filename='rawfile_for_mytask_run1_raw.fif' # Fill in name of rest rawdata-file
run='01'

create_BIDS(filename, subject_in, subject_out, date, session, run, proc, task, rest_events)
