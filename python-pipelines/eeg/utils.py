#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 16:54:30 2021

@author: mhusberg
"""
import pandas as pd

def download_eegbci_data(subject=1, runs=[1,2]):
    
    from mne.datasets import eegbci
    from mne.io import concatenate_raws, read_raw_edf
    
    # download data 
    # Subjects can be in the range of 1-109 (inclusive).
    #runs = [1,2]  # use only eynes open and eyes closed data

    fnames = eegbci.load_data(subject, runs, path = '/m/nbe/scratch/rubberboot/test/data/')
    raws = [read_raw_edf(f, preload=True) for f in fnames]
    raw = concatenate_raws(raws)
    print(raw.info)

def write_bad_channels(channel_fname, info, bads):
    
    # Set up dataframe with sub-columns from bids specification
    df = pd.DataFrame(columns=['name','status','status_description'])

    # Default values
    names = info['ch_names']
    df['name'] = names
    df['status'] = 'good' #TODO: check if bad channels were marked during recording
    df['status_description'] = 'n/a'
    
    # Mark bad channels
    for bad in bads:
        df.loc[df['name'] == bad,'status'] = 'bad'
        df.loc[df['name'] == bad,'status_description'] = 'manually marked bad'
    
    try:
        df.to_csv(channel_fname,index=False)
        print('Saving bad channels')
    except:
        print('Saving bad channels failed')
        
def read_channels(channel_fname):
    
    print('To be implemented')
    # Read csv-file into dataframe
    # Find bad channels
    # Transform into list
    
    return
    
    