"""
Do-it script to execute the entire pipeline using the doit tool:
http://pydoit.org

All the filenames are defined in config.py
"""
import sys
sys.path.append('..')

from config_meg import subjects, get_all_fnames, fname

# Configuration for the "doit" tool.
DOIT_CONFIG = dict(
    # While running scripts, output everything the script is printing to the
    # screen.
    verbosity=2,

    # When the user executes "doit list", list the tasks in the order they are
    # defined in this file, instead of alphabetically.
    sort='definition',
)


def task_tsss():
    """Step 01: Perform TSSS maxfilter"""
    for subject in subjects:
        yield dict(
            name=f'sub-{subject:02d}',
            file_dep=get_all_fnames(subject, 'raw') + ['01_tsss.py'],
            targets=get_all_fnames(subject, 'tsss'),
            actions=[f'python 01_tsss.py {subject}'],
        )

def task_filt():
    """Step 02: Perform frequency filtering"""
    for subject in subjects:
        yield dict(
            name=f'sub-{subject:02d}',
            file_dep=get_all_fnames(subject, 'tsss') + ['02_freqfilt.py'],
            targets=get_all_fnames(subject, 'filt'),
            actions=[f'python 02_freqfilt.py {subject}'],
        )

def task_ica():
    """Step 03: Remove blink (EOG) and heart-beat (ECG) artifacts using ICA"""
    for subject in subjects:
        yield dict(
            name=f'sub-{subject:02d}',
            file_dep=(get_all_fnames(subject, 'raw', exclude='emptyroom') + 
                      get_all_fnames(subject, 'filt', exclude='emptyroom') +
                      ['03_ica.py']),
            targets=(get_all_fnames(subject, 'ica', exclude='emptyroom') +
                     get_all_fnames(subject, 'clean', exclude='emptyroom')),
            actions=[f'python 03_ica.py {subject}'],
        )

def task_psds():
    """Step 04: Compute the Power Spectral Density (PSD) for each recording."""
    for subject in subjects:
        yield dict(
            name=f'sub-{subject:02d}',
            file_dep=[
                fname.clean(subject=subject, task='eyesopen', run=1),
                fname.clean(subject=subject, task='eyesclosed', run=1),
                fname.clean(subject=subject, task='pasat', run=1),
                fname.clean(subject=subject, task='pasat', run=2),
                '04_psds.py',
            ],
            targets=[fname.psds(subject=subject)],
            actions=[f'python 04_psds.py {subject}'],
        )
