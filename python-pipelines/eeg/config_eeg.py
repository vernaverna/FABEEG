"""
These are all the relevant parameters that are unique to the EEG analysis
pipeline.
"""

# Some relevant files are in the parent folder.
import sys
sys.path.append('../../python-pipelines')

from fnames import FileNames
import os
import mne
import pandas as pd
import numpy as np

from config_common import (raw_data_dir, processed_data_dir, figures_dir,
                           reports_dir)


# Parameters

# Highpass filter above 1Hz. This is needed for the ICA to perform well
# later on. Lowpass filter below 100Hz to get rid of the signal produced by
# the cHPI coils. Notch filters at 50Hz and 100Hz to get rid of powerline.
fmin = 1
fmax = 45
fnotch = [50, 100]

# Computation of the PSDs
n_fft = 1024  # Higher number means more resolution at the lower frequencies

# Subjects
all_subjects = os.listdir(raw_data_dir)
for s in all_subjects:
    s_path = os.path.join(raw_data_dir, s)
    if not os.path.isdir(s_path):
        all_subjects.remove(s)

#Ages of subjects
age_df = pd.read_csv('ages.csv')
#age_df['Cap'] = np.nan
#age_df = age_df[['File', 'Sex', 'Age']]
#age_files = [x.replace('-', '') for x in list(age_df['File'])]
#age_df['File'].replace(dict( zip(list(age_df['File']), age_files) ), inplace=True)
#age_df.to_csv('ages.csv', index=False)

# Subjects removed from the EEG analysis because of some problem
bad_subjects = ['derivatives']

# Analysis is performed on these subjects
subjects = [subject for subject in all_subjects if subject not in bad_subjects]

# Define the frequency bands (heuristic approach as in BRRR)
f_bands = [(1,3), (3,5.2), (5.2,7.6), (7.6,10.2), (10.2, 13), (13,16),
           (16,19.2), (19.2,22.6), (22.6,26.2), (26.2,30), (30,34), (34,38.2), (38.2,42.6)]

#TODO: change these!

# Helper function to calculate bandpower for each sensor
def bandpower(psd, f, fmin, fmax): 
    min_index = max(0, np.argmax(f > fmin) - 1)
    max_index = np.argmax(f > fmax) -1
    
    return np.trapz(psd[:,min_index: max_index], f[min_index: max_index], axis=1)

# Bad EEG channels for each subject.
# Made-up list: this is work in progress
# TODO: save bad channels in files, BIDS compliant or MNE python compliant?
bads = {}

eog_channel = 'PIETSO' #Change?

###############################################################################
# Templates for filenames
#
# This part of the config file uses the FileNames class. It provides a small
# wrapper around string.format() to keep track of a list of filenames.
# See fnames.py for details on how this class works.
fname = FileNames()

# Some directories
fname.add('raw_data_dir', raw_data_dir)
fname.add('processed_data_dir', processed_data_dir)


### TO DO: Change these (tasks and runs out, bids)

# Continuous data
fname.add('raw', '{raw_data_dir}/sub-{subject}/eeg/sub-{subject}_task-sleep_eeg.edf')
fname.add('filt', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_filt_eeg.fif')
fname.add('clean', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_clean_eeg.fif')

# Files used during EOG and ECG artifact suppression
fname.add('bads', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_channels.csv')
fname.add('annotations', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_annotations.csv')
fname.add('ica', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_ica.h5')

# PSD files
fname.add('psds', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_psds.h5')
fname.add('psds_fif', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_psds.fif')
fname.add('evoked', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_evoked.fif')


# Filenames for MNE reports
fname.add('reports_dir', f'{reports_dir}/eeg')
fname.add('report', '{reports_dir}/sub-{subject}-report.h5')
fname.add('report_html', '{reports_dir}/sub-{subject}-report.html')

# Filenames for figures
fname.add('figures_dir', f'{figures_dir}/eeg')
fname.add('figure_psds', '{figures_dir}/psds.pdf')


def get_all_fnames(subject, kind, tasks=None, runs=None, exclude=None): #CHANGED
    """Get all filenames for a given subject of a given kind.

    Not all subjects have exactly the same files. For example, subject 1 does
    not have an emptyroom recording, while subject 4 has 2 runs of emptyroom.
    Use this function to get a list of the the files that are present for a
    given subject. It will check which raw files there are and based on that
    will generate a list with corresponding filenames of the given kind.

    You can exclude the recordings for one or more tasks with the ``exclude``
    parameter. For example, to skip the emptyroom recordings, set
    ``exclude='emptyroom'``.

    Parameters
    ----------
    subject : str
        The subject to get the names of the raw files for.
    kind : 'raw' | 'tsss' | 'filt' | 'eog_ecg_events' | 'ica'
        The kind of files to return the filenames for.
    exclude : None | str | list of str
        The tasks to exclude from the list.
        Defaults to not excluding anything.

    Returns
    -------
    all_fnames : list of str
        The names of the files of the given kind.
    """
    import os.path as op

    if exclude is None:
        exclude = []
    elif type(exclude) == str:
        exclude = [exclude]
    elif type(exclude) != list:
        raise TypeError('The `exclude` parameter should be None, str or list')

    all_fnames = list()
    print('Looking for: ' + str(fname.raw(subject=subject))) #ADDED
    if op.exists(fname.raw(subject=subject)) and kind == None:
        all_fnames.append(fname.raw(subject=subject)) #raw file
        
    elif op.exists(fname.raw(subject=subject)):
        all_fnames.append(fname.files()[f'{kind}'](subject=subject))
        
   # for task in tasks:
   #     if task in exclude:
   #         continue
   #     for run in runs:
   #         print('Looking for: ' + str(fname.raw(subject=subject))) #CHANGED
   #         if op.exists(fname.raw(subject=subject)):
   #             all_fnames.append(fname.files()[f'{kind}'](subject=subject))
    return all_fnames


def task_from_fname(fname):
    """Extract task name from a BIDS filename."""
    import re
    match = re.search(r'task-([^_]+)_run-(\d\d)', str(fname))
    task = match.group(1)
    run = int(match.group(2))
    if task == 'pasat':
        return f'{task}_run{run}'
    else:
        return task

#Function to change the metadata & pick common channels
def change_metadata(raw):
    """   
    Change the metadata, set correct channel types & pick common channels from subjects' raw data
    
    Parameters
    ----------
    raw : mne.io.Raw
        Raw eeg data (FIF)

    Returns
    -------
    raw : mne.io.Raw
        Raw eeg data (FIF) with modified metadata
    cap_status : str
        'FT' or '-', infdicating if larger cap was used

    """
    
    #For testdata only, remove . in channel-names
    raw.rename_channels(lambda x: x.strip('.'))  # remove dots from channel names, FIXME
    raw.rename_channels(lambda x: x.upper())  # capitalize the ch names
    
    # Get reference channel; either Pz or Cz
    ch_info_df = raw.describe(data_frame=True)
    cz_deviation = abs(ch_info_df.iloc[17, 6]-ch_info_df.iloc[17, 4]) #Q3-Q1
    pz_deviation = abs(ch_info_df.iloc[18, 6]-ch_info_df.iloc[18, 4])
    
    if cz_deviation < pz_deviation: #which channel has lower deviation ~> zero signal
        reference = 'CZ'
    else:
        reference = 'PZ'
    
    
    #Drop channels that are not in the shared_channels.txt .file
    allowed_chs = []
    with open('shared_channels.txt') as f:
         allowed_chs = f.read().splitlines()
        
    if ('FT9' or 'FT10') in raw.info['ch_names']:
        cap_status = 'FT'
    else:
        cap_status = '-'
    
    for ch_name in raw.info['ch_names']:
        if ch_name not in allowed_chs:
            raw.drop_channels(ch_name)

    # Drop photic
    raw.drop_channels("PHOTIC")
    #raw.drop_channels("PIETSO") # or keep pietso just in case?
    #raw.drop_channels("EKG")
    
    # Reset sensor types 
    mapping = dict( [(ch_name, 'ecg') if ch_name=='EKG' else (ch_name, 'eog') if ch_name=='PIETSO' 
                     else (ch_name, 'eeg') for ch_name in raw.info['ch_names']] )
    raw.set_channel_types(mapping)
    
    # Sensor locations not provided with data - use standard layout
    ten_twenty_montage = mne.channels.make_standard_montage('standard_1020')
    ten_twenty_montage.ch_names = [CH_NAME.upper() for CH_NAME in ten_twenty_montage.ch_names]

    raw = raw.copy().set_montage(ten_twenty_montage, on_missing='ignore') #for non-eeg chs?
    
    #Set the reference as well
    #raw.set_eeg_reference(ref_channels=[reference])
    raw.set_eeg_reference('average')
    
    return raw, cap_status



















