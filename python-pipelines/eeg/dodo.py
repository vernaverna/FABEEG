"""
Do-it script to execute the entire pipeline using the doit tool:
http://pydoit.org

All the filenames are defined in config.py
"""
import sys
sys.path.append('..')

from config_eeg import fname, subjects, get_all_fnames

# Configuration for the "doit" tool.
DOIT_CONFIG = dict(
    # While running scripts, output everything the script is printing to the
    # screen.
    verbosity=2,

    # When the user executes "doit list", list the tasks in the order they are
    # defined in this file, instead of alphabetically.
    sort='definition',
)

def task_filt():
    """Step 01: Perform frequency filtering"""
    for subject in subjects:
        yield dict(
            name=f'sub-{subject:02d}',
            file_dep=get_all_fnames(subject, 'raw') + ['freqfilt_01.py'],
            targets=get_all_fnames(subject, 'filt'),
            actions=[f'python 01_freqfilt.py {subject}'],
        )

#def task_ica():
#    """Step 02: Remove blink (EOG) artifacts using ICA"""
#    for subject in subjects:
#        yield dict(
#            name=f'sub-{subject:02d}',
#            file_dep=(get_all_fnames(subject, 'raw', exclude=['emptyroom', 'eyesclosed']) +
#                      get_all_fnames(subject, 'filt', exclude=['emptyroom', 'eyesclosed']) +
#                      ['02_ica.py']),
#            targets=(get_all_fnames(subject, 'ica', exclude=['emptyroom', 'eyesclosed']) +
#                     get_all_fnames(subject, 'clean', exclude=['emptyroom', 'eyesclosed'])),
#            actions=[f'python 02_ica.py {subject}'],
#        )

def task_psds():
    """Step 03: Compute the Power Spectral Density (PSD) for each recording."""
    for subject in subjects:
        yield dict(
            name=f'sub-{subject:02d}',
            file_dep=[
                fname.filt(subject=subject),
                'psds_03.py',
            ],
            targets=[fname.psds(subject=subject)],
            actions=[f'python psds_03.py {subject}'],
        )
        
def task_bandpower():
    """ Step 04: Compute the bandpower and save results as .csv"""
    for subject in subjects:
    	yield dict(
    	    name=f'sub-{subject:02d}',
    	    file_dep=[fname.psds(subject=subject), 'bandpower_04.py',]
    	    actions=[f'python bandpower_04.py {subject}']
    	)
       
       
