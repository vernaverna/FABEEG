"""
These are all the relevant parameters that are unique to the EEG analysis
pipeline.
"""

# Some relevant files are in the parent folder.
import sys
sys.path.append('../')

from fnames import FileNames

from config_common import (raw_data_dir, processed_data_dir, figures_dir,
                           reports_dir)


# Parameters

# Highpass filter above 1Hz. This is needed for the ICA to perform well
# later on. Lowpass filter below 100Hz to get rid of the signal produced by
# the cHPI coils. Notch filters at 50Hz and 100Hz to get rid of powerline.
fmin = 1
fmax = 70
fnotch = [50]

# Computation of the PSDs
n_fft = 1024  # Higher number means more resolution at the lower frequencies

# Subjects

all_subjects = ['S001','S002','S003','S004','S005']



# Subjects removed from the MEG analysis because of some problem
bad_subjects = [
    24, 29,  # No EEG data present
]

# Analysis is performed on these subjects
subjects = [subject for subject in all_subjects if subject not in bad_subjects]

# Tasks 
tasks = ['restEO', 'restEC']

runs = [1, 2]

# Bad EEG channels for each subject.
# Made-up list: this is work in progress
# TODO: save bad channels in files, BIDS compliant or MNE python compliant?
bads = {
    'S001': ['Fc5.','Fc1.'],
    'S002': ['Fcz.', 'Fc1.'],
    'S003': [],
    'S004': [],
    'S005': [],
}

eog_channel = 'Fp1.' #Change!

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

# Continuous data
fname.add('raw', '{raw_data_dir}/sub-{subject}/ses-01/eeg/sub-{subject}_ses-01_task-{task}_run-{run:02d}_proc-raw_eeg.edf')
fname.add('filt', '{processed_data_dir}/sub-{subject}/ses-01/eeg/sub-{subject}_ses-01_task-{task}_run-{run:02d}_proc-filt_eeg.fif')
fname.add('clean', '{processed_data_dir}/sub-{subject}/ses-01/eeg/sub-{subject}_ses-01_task-{task}_run-{run:02d}_proc-clean_eeg.fif')

# Files used during EOG and ECG artifact suppression
fname.add('bads', '{processed_data_dir}/sub-{subject}/ses-01/eeg/sub-{subject}_ses-01_task-{task}_run-{run:02d}_channels.csv')
fname.add('annotations', '{processed_data_dir}/sub-{subject}/ses-01/eeg/sub-{subject}_ses-01_task-{task}_run-{run:02d}_annotations.csv')
fname.add('ica', '{processed_data_dir}/sub-{subject}/ses-01/eeg/sub-{subject}_ses-01_task-{task}_run-{run:02d}_ica.h5')

# PSD files
fname.add('psds', '{processed_data_dir}/sub-{subject}/ses-01/eeg/sub-{subject}_psds.h5')

# Filenames for MNE reports
fname.add('reports_dir', f'{reports_dir}/eeg')
fname.add('report', '{reports_dir}/sub-{subject}-report.h5')
fname.add('report_html', '{reports_dir}/sub-{subject}-report.html')

# Filenames for figures
fname.add('figures_dir', f'{figures_dir}/eeg')
fname.add('figure_psds', '{figures_dir}/psds.pdf')


def get_all_fnames(subject, kind, tasks, runs, exclude=None):
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
    for task in tasks:
        if task in exclude:
            continue
        for run in runs:
            print('Looking for: ' + str(fname.raw(subject=subject, task=task, run=run)))
            if op.exists(fname.raw(subject=subject, task=task, run=run)):
                all_fnames.append(fname.files()[f'{kind}'](subject=subject, task=task, run=run))
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