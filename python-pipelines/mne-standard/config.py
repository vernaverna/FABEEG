"""
===========
Config file
===========

Configuration parameters for the study.
"""

import os
from fnames import FileNames


user = os.environ['USER']

study_path = '/m/nbe/project/aging/recordings/'

N_JOBS = 4

# Preprocessing with maxfilter
proc = 'tsssmctrans' 

# Rejection limits
reject = dict(mag=4e-12, grad=3000e-13)#, eog=150e-6) 

# Default baseline window
default_baseline = (-0.2,0)

# Source spacing 
spacing = 'ico4'

sessions=['01','02']

tasks = ['naming','N400','envspeech']

session1 = {'naming'   :['01','02'], 
         'N400'     :['01','02'], 
         'envspeech':['01'],
         'restEO'   :['01'],
         'restEC'   :['01'],
         'emptyroom':['01']
         }

session2 = {'naming'   :['01','02'], 
         'N400'     :['01','02'], 
         'envspeech':['01'],
         'restEO'   :['01','02'],
         'CKC'   :['01'],
         'emptyroom':['01']
         }
###############################################################################
# Folders (TODO: decide on folder structure)

data_path = os.path.join(study_path, 'BIDS/')
subjects_dir = os.path.join(data_path, 'derivatives/', 'freesurfer/')
processed_dir = os.path.join(study_path,'processed/')
#meg_dir = os.path.join(study_path,'derivatives/','meg/')
os.environ["SUBJECTS_DIR"] = subjects_dir


###############################################################################
# Subjects

subjects=['sub-s01y1', 'sub-s02y2','sub-s03y3','sub-s04y4','sub-s05y1','sub-s06y2','sub-s07y3','sub-s08y4','sub-s09o1']

# Subject-codes for freesurfing
subject_info = {
        'sub-s01y1': 'S01Y1/S8790/00002',
        'sub-s02y2': 'S02Y2/S8848/00002',
        'sub-s03y3': 'S03Y3/S8804/00002',
        'sub-s04y4': 'S04Y04/S8822/00002'
        }

# or: read all subjects from folder-names
# or: place subjects in task-specific settings

###############################################################################
# Templates for filenames
fname = FileNames()

# Some directories
fname.add('data_path', data_path)
fname.add('megbids_dir','{data_path}/{subject}/ses-{ses}/meg/')
fname.add('subjects_dir', subjects_dir)
fname.add('processed_dir', processed_dir)
fname.add('megprocessed_dir','{processed_dir}{subject}/ses-{ses}/meg/')

# Add these so we can use them in the filenames
fname.add('spacing', spacing)

# Maxfilter
fname.add('pos', '{megbids_dir}/{subject}_ses-{ses}_task-{task}_run-{run}_movecomp.pos')
fname.add('tsss_log', '{megbids_dir}/{subject}_ses-{ses}_task-{task}_run-{run}_tsss_log.log')

# Sensor-level files
fname.add('raw', '{megbids_dir}/{subject}_ses-{ses}_task-{task}_run-{run}_proc-{proc}_meg.fif')
fname.add('epochs', '{megprocessed_dir}/epochs/{subject}_ses-{ses}_task-{task}-epo.fif')
fname.add('epochs_proc', '{megprocessed_dir}/epochs/{subject}_ses-{ses}_task-{task}_proc-{proc}-epo.fif')
fname.add('evoked', '{megprocessed_dir}/evoked/{subject}_ses-{ses}_task-{task}-ave.fif')

# Artefact removal
fname.add('ica', '{megprocessed_dir}/ica/{subject}_ses-{ses}_task-{task}_{method}-ica.fif')
          
# Coregistration files (WIP)          
fname.add('fid', '{subjects_dir}/{subject}/bem/{subject}-fiducials.fif')
fname.add('trans_file', '{megprocessed_dir}/forward/{subject}_ses-{ses}_coreg-transfile.fif')
fname.add('trans_file_old', '{subjects_dir}/{subject}/mri/T1-neuromag/sets/COR-transfile.fif')
fname.add('coreg_log', '{megprocessed_dir}/forward/coreg_log_{subject}_ses-{ses}_hs.csv')

# Source space and forward model          
fname.add('bem_sol', '{megprocessed_dir}/forward/{subject}-{ntri}-bem-sol.fif')
fname.add('fwd', '{megprocessed_dir}/forward/{subject}_ses-{ses}_task-{task}_run-{run}-{spacing}-fwd.fif')
fname.add('src', '{megprocessed_dir}/forward/{subject}-{spacing}-src.fif')

# Inverse solution
fname.add('noise_cov', '{megprocessed_dir}/noise_cov/{subject}_ses-{ses}_task-{task}_cov.fif')
fname.add('inv', '{megprocessed_dir}/inverse/{subject}_ses-{ses}_task-{task}-{spacing}-inv.fif')

# STC-files
fname.add('stc_hemi', '{megprocessed_dir}/stcs/{subject}_ses-{ses}_task-{task}-{condition}-{hemi}.stc')
fname.add('stc', '{megprocessed_dir}/stcs/{subject}_ses-{ses}_task-{task}-{condition}')

