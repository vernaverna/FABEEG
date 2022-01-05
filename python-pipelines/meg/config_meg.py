"""
These are all the relevant parameters that are unique to the MEG analysis
pipeline.
"""

# Some relevant files are in the parent folder.
import sys
sys.path.append('../')

from fnames import FileNames

from config_common import (raw_data_dir, processed_data_dir, figures_dir,
                           reports_dir, all_subjects, tasks)


###############################################################################
# Parameters that should be mentioned in the paper

# Highpass filter above 1Hz. This is needed for the ICA to perform well
# later on. Lowpass filter below 100Hz to get rid of the signal produced by
# the cHPI coils. Notch filters at 50Hz and 100Hz to get rid of powerline.
fmin = 1
fmax = 100
fnotch = [50, 100]

# Computation of the PSDs
n_fft = 1024  # Higher number means more resolution at the lower frequencies


###############################################################################
# Parameters pertaining to the subjects

# Subjects removed from the MEG analysis because of some problem
bad_subjects = [
    12, 14, 19, 35, 36,  # These subjects have some kind of interference.
    34, 39,  # These subjects have problems with the HPI coils
]

# Analysis is performed on these subjects
subjects = [subject for subject in all_subjects if subject not in bad_subjects]

# Bad MEG channels for each subject.
# Manually marked by Marijn van Vliet.
bads = {
    1: ['MEG0533', 'MEG1112'],
    2: ['MEG0533', 'MEG1112', 'MEG0111', 'MEG1812'],
    3: ['MEG0533', 'MEG1112', 'MEG1033'],
    4: ['MEG0533', 'MEG1112'],
    5: ['MEG0533', 'MEG1112', 'MEG2141'],
    6: ['MEG0533', 'MEG1112', 'MEG0642'],
    7: ['MEG0533', 'MEG1112', 'MEG1731', 'MEG1732', 'MEG0121'],
    8: ['MEG0533', 'MEG1112', 'MEG1033', 'MEG2133', 'MEG2543', 'MEG2533', 'MEG2532', 'MEG2542', 'MEG2531', 'MEG2541', 'MEG2131', 'MEG2141'],
    9: ['MEG0533', 'MEG1112', 'MEG0313', 'MEG0311', 'MEG2121', 'MEG2512'],
    10: ['MEG0533', 'MEG1112', 'MEG1931'],
    11: ['MEG0533', 'MEG1112'],
    12: ['MEG0533', 'MEG1112'],
    13: ['MEG0533', 'MEG1112', 'MEG2232'],
    14: ['MEG0112', 'MEG0113', 'MEG0141', 'MEG0142', 'MEG0143', 'MEG1423', 'MEG0111', 'MEG1112'],
    15: ['MEG0533', 'MEG1112'],
    16: ['MEG0533', 'MEG1112'],
    17: ['MEG0533', 'MEG1112'],
    18: ['MEG0533', 'MEG1112', 'MEG1722'],
    19: ['MEG0533', 'MEG1112'],
    20: ['MEG0533', 'MEG1112', 'MEG2331'],
    21: ['MEG0533', 'MEG1112', 'MEG2542'],
    22: ['MEG0221', 'MEG0321', 'MEG0533', 'MEG1112', 'MEG0122', 'MEG0123', 'MEG0222', 'MEG0223', 'MEG0323', 'MEG2132'],
    23: ['MEG0533', 'MEG1112'],
    24: ['MEG0533', 'MEG1112', 'MEG2611', 'MEG2612', 'MEG2613'],
    25: ['MEG0533', 'MEG1112'],
    26: ['MEG0533', 'MEG1112', 'MEG1721', 'MEG2542'],
    27: ['MEG0533', 'MEG1112', 'MEG2542', 'MEG0211'],
    28: ['MEG0533', 'MEG1112', 'MEG2542', 'MEG0211'],
    29: ['MEG0533', 'MEG1112', 'MEG2331', 'MEG2131'],
    30: ['MEG0533', 'MEG0741', 'MEG0742', 'MEG1821', 'MEG1822', 'MEG1823', 'MEG1112', 'MEG2331'],
    31: ['MEG0533', 'MEG1112'],
    32: ['MEG0533', 'MEG1112', 'MEG0142', 'MEG2542', 'MEG2331'],
    33: ['MEG0533', 'MEG2523', 'MEG1112', 'MEG2542', 'MEG2331', 'MEG2131'],
    34: ['MEG0533', 'MEG1112'],
    35: ['MEG0533', 'MEG1112'],
    36: ['MEG0533', 'MEG1112'],
    37: ['MEG0533', 'MEG2142', 'MEG1112', 'MEG2543', 'MEG2533'],
    38: ['MEG0431', 'MEG0533', 'MEG1112'],
    39: ['MEG0533', 'MEG1112'],
    40: ['MEG0533', 'MEG1112', 'MEG2331'],
    41: ['MEG0533', 'MEG2111', 'MEG2441', 'MEG1112'],
    42: ['MEG0533', 'MEG1112'],
    43: ['MEG0533', 'MEG1112'],
    44: ['MEG1112', 'MEG0533', 'MEG2331'],
}

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
fname.add('raw', '{raw_data_dir}/sub-{subject:02d}/ses-01/meg/sub-{subject:02d}_ses-01_task-{task}_run-{run:02d}_meg.fif')
fname.add('tsss', '{processed_data_dir}/sub-{subject:02d}/ses-01/meg/sub-{subject:02d}_ses-01_task-{task}_run-{run:02d}_tsss.fif')
fname.add('filt', '{processed_data_dir}/sub-{subject:02d}/ses-01/meg/sub-{subject:02d}_ses-01_task-{task}_run-{run:02d}_filt.fif')
fname.add('clean', '{processed_data_dir}/sub-{subject:02d}/ses-01/meg/sub-{subject:02d}_ses-01_task-{task}_run-{run:02d}_clean.fif')

# Files used during EOG and ECG artifact suppression
fname.add('ica', '{processed_data_dir}/sub-{subject:02d}/ses-01/meg/sub-{subject:02d}_ses-01_task-{task}_run-{run:02d}_ica.h5')

# Calibration files used during tSSS filtering
fname.add('cal', '{raw_data_dir}/sub-{subject:02d}/ses-01/meg/sub-{subject:02d}_ses-01_acq-calibration_meg.dat')
fname.add('ctc', '{raw_data_dir}/sub-{subject:02d}/ses-01/meg/sub-{subject:02d}_ses-01_acq-crosstalk_meg.fif')

# PSD files
fname.add('psds', '{processed_data_dir}/sub-{subject:02d}/ses-01/meg/sub-{subject:02d}_psds.h5')

# Filenames for MNE reports
fname.add('reports_dir', f'{reports_dir}/meg')
fname.add('report', '{reports_dir}/sub-{subject:02d}-report.h5')
fname.add('report_html', '{reports_dir}/sub-{subject:02d}-report.html')

# Filenames for figures
fname.add('figures_dir', f'{figures_dir}/meg')
fname.add('figure_psds', '{figures_dir}/psds.pdf')


def get_all_fnames(subject, kind, exclude=None):
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
    subject : int
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
        for run in [1, 2]:
            if op.exists(fname.raw(subject=subject, task=task, run=run)):
                all_fnames.append(fname.files()[f'{kind}'](subject=subject, task=task, run=run))
    return all_fnames


def task_from_fname(fname):
    """Extract task name from a BIDS filename."""
    import re
    return re.search(r'task-([^_]+)', str(fname)).group(1)
