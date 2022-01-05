# -*- coding: utf-8 -*-

import argparse
import mne
import os
import numpy as np

from config import (fname, proc,  reject)


# Handle command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', help='The subject to process')
parser.add_argument('-ses', help='The session to process. ')
parser.add_argument('-task', help='The task to process. ')
args = parser.parse_args()
subject = args.subject
ses = args.ses
task = args.task

print('Processing subject: ', subject)
print('Processing session: ', ses)
print('Processing task: ', task)

ica=False # FIXME


if task=='naming':
    from settings_naming import (bandpass_fmin, bandpass_fmax, runs, event_ids,  epoch_tmin, epoch_tmax, conditions, baseline)
elif task=='N400':
    from settings_N400 import (bandpass_fmin, bandpass_fmax, runs, event_ids, epoch_tmin, epoch_tmax, conditions, baseline)
elif task=='envspeech':
    from settings_envspeech import (bandpass_fmin, bandpass_fmax, runs, event_ids, epoch_tmin, epoch_tmax, conditions, baseline)
else:
    print('Undefined task')


raws = list()
for run in runs:
    # Read the raw data
    print('Processing run: ', run)

    raw = mne.io.read_raw_fif(fname.raw(subject=subject,
                                        ses=ses,
                                        task=task,
                                        run=run,
                                        proc=proc, 
                                        preload = True))
    raw.load_data()
    raws.append(raw)

# Concatenate raw data, check that data has been transformed to the same position
raw=mne.concatenate_raws(raws)
#delete to save memory
del raws


#################################################
# Remove artefacts, currently using PCA
if ica==False:
    projs_ecg, _ = mne.preprocessing.compute_proj_ecg(raw, n_grad=2, n_mag=2)
    projs_eog, _ = mne.preprocessing.compute_proj_eog(raw, n_grad=2, n_mag=2)
        
    raw.info['projs'] += projs_ecg
    raw.info['projs'] += projs_eog
    raw.apply_proj()
#################################################

# Filter
raw.filter(bandpass_fmin, bandpass_fmax) # check default filtering
raw.notch_filter(freqs=np.arange(50, 251, 50))

# Get epochs and average
events =  mne.find_events(raw, min_duration = 2/raw.info['sfreq'])
epochs = mne.Epochs(raw, events, event_ids, epoch_tmin, epoch_tmax, baseline = baseline, 
                    preload=True, reject=reject, proj=False) 


# Clean here or separately?
# Apply precomputed ICA
if ica==True:
    ica= mne.preprocessing.read_ica(fname.ica(subject=subject, method='infomax', ses=ses, task=task))

    ica.apply(epochs, exclude=ica.exclude)

# Save epochs
if not os.path.exists(fname.megprocessed_dir(subject=subject, ses=ses) + 'epochs/'):    
    os.makedirs(fname.megprocessed_dir(subject=subject,ses=ses) + 'epochs/')
epochs.save(fname.epochs(subject=subject,ses=ses,task=task),overwrite=True)

# Information about rejections
#epochs.drop_bad() 
#epochs.plot_drop_log(subject=subject)  
#droplog_fname = epochs_directory + subject + '_' + session + '_task-' + task + 'epochs_rejected.png'
#plt.savefig(droplog_fname)
      
# Calculate and save evoked responses
if not os.path.exists(fname.megprocessed_dir(subject=subject, ses=ses) + 'evoked/'):    
    os.makedirs(fname.megprocessed_dir(subject=subject,ses=ses) + 'evoked/')
    
evokeds=list()
for cond in conditions:
    evoked = epochs[cond].average()
    evokeds.append(evoked)

    # Save evoked responses for each condition
    mne.write_evokeds(fname.evoked(subject=subject,ses=ses,task=task),evoked)

#Or save list of evokeds
#evoked_fname = output_directory + subject + '_' + session + '_task-' + task + '-ave.fif'
#mne.write_evokeds(evoked_fname)
