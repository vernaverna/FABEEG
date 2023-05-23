"""
These are all the relevant parameters that are unique to the EEG analysis
pipeline.
"""

# Some relevant files are in the parent folder.
import sys
sys.path.append('../../python-pipelines')

from fnames import FileNames
import os
import pandas as pd
import numpy as np

from config_common import (raw_data_dir, processed_data_dir, figures_dir,
                           reports_dir)


# Parameters

# Highpass filter above 1Hz. This is needed for the ICA to perform well
# later on. Lowpass filter below 100Hz to get rid of the signal produced by
# the cHPI coils. Notch filters at 50Hz and 100Hz to get rid of powerline.
fmin = 1
fmax = 40
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
bad_subjects = []

# Analysis is performed on these subjects
subjects = [subject for subject in all_subjects if subject not in bad_subjects]

# Tasks 
#tasks = ['restEO', 'restEC']

#runs = [1, 2]

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
fname.add('raw', '{raw_data_dir}/sub-{subject}/eeg/sub-{subject}_eeg.edf')
fname.add('filt', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_filt_eeg.fif')
fname.add('clean', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_clean_eeg.fif')

# Files used during EOG and ECG artifact suppression
fname.add('bads', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_channels.csv')
fname.add('annotations', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_annotations.csv')
fname.add('ica', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_ica.h5')

# PSD files
fname.add('psds', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_psds.h5')
fname.add('psds', '{processed_data_dir}/sub-{subject}/eeg/sub-{subject}_psds.fif')
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