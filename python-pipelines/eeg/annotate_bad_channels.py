"""
This script is used for manually annotation bad EEG channels.
The idea is to run this in an ipython console through:

>>> %run annotate_bad_channels.py <subject number>

The script will load the raw data and make a plot of it. Inside the plot you
can scroll through the data and click on the channel names to mark them as bad.
When you close the plot, the script will print out the channels you have
marked. These channels can then be added to the big `bads` list inside
config_eeg.py.
"""
import argparse
import mne
from mne.preprocessing import annotate_muscle_zscore
import matplotlib.pyplot as plt
import numpy as np

from config_eeg import fname, tasks
from utils import write_bad_channels, read_channels

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', type=str, help='The subject to process')
#parser.add_argument('task', type=str, help='The task to process')
args = parser.parse_args()

task = tasks[1] #FIXME
run = 1

read_known_bad_channels = False # Not implemented
set_known_annotations = False
find_muscle_artefacts = False
set_manual_annotations = False
set_manual_bads = True
save_channels = False # Automatic saving does not work, save manually once finished 
save_annotations = False

#%%
# Read raw data
fname_raw = fname.raw(subject=args.subject, task=task, run=run)
raw = mne.io.read_raw_edf(fname_raw, preload=True)

# File for bad channels
fname_channels = fname.bads(subject=args.subject, task=task, run=run)

# File for annotations
fname_annotations = str(fname.annotations(subject=args.subject, task=task, run=run))


# Open up raw data, mark new bad channels if necessary
raw.pick_types(meg=False, eeg=True, eog=True)

# Filter data, note that example data has sampling frequency 80Hz
raw.filter(1, 70)
raw.notch_filter([50])

# Any existing annotations in the data?
#existing_annotations = raw.annotations

# or use empty existing
existing_annotations = mne.Annotations(onset=[],
                           duration=[],
                           description=[])

#%%

if set_known_annotations:
    
    # Example annotations
    # Add your own, or read from file    
    
    #known_annotations = mne.read_annotations('fname_annotations)

    known_annotations = mne.Annotations(onset=[],
                           duration=[],
                           description=[])
    
    print(known_annotations)

    try:
        # FIXME: does not work if orig_time is not the same
        existing_annotations = existing_annotations + known_annotations
        print('Merging annotations')
    except:
        existing_annotations = known_annotations
        print('Warning: new annotations not compatible with existing annotations, dropping existing')
        
    # Assumes that meas_time and orig_time are the same, works for Neuromag data, check EEG
    raw.set_annotations(existing_annotations) # overwrites existing annotations
    print(raw.annotations)

def annotate_muscles(raw, z_score=5, filter_freq=[110,140]):

    # The threshold is data dependent, check the optimal threshold by plotting
    # ``scores_muscle``.
    threshold_muscle = z_score  # z-score
    # Choose one channel type, if there are axial gradiometers and magnetometers,
    # select magnetometers as they are more sensitive to muscle activity.
    annot_muscle, scores_muscle = annotate_muscle_zscore(
        raw, ch_type="eeg", threshold=threshold_muscle, min_length_good=0.2,
        filter_freq=filter_freq)
    
    return annot_muscle, scores_muscle

if find_muscle_artefacts:
    z_score=10
    annot_muscle, scores_muscle = annotate_muscles(raw,z_score=z_score,filter_freq=[60,70])# high-freq not available in test data

    fig, ax = plt.subplots()
    ax.plot(raw.times, scores_muscle)
    ax.axhline(y=z_score, color='r')
    ax.set(xlabel='time, (s)', ylabel='zscore', title='Muscle activity')
    plt.show()
    
    try:
        # FIXME: does not work if orig_time is not the same
        existing_annotations = existing_annotations + annot_muscle
        print('Merging annotations')
    except:
        existing_annotations = annot_muscle 
        print('Warning: muscle annotations not compatible with existing annotations')
    
    raw.set_annotations(existing_annotations) # Warning! removes other annotations than muscles
    raw.plot(start=5, duration=10,scalings=dict(eog=100E-6, eeg=50E-6))
    
if set_manual_annotations:
    print('Annotate manually, then close')
    fig=raw.plot(scalings=dict(eog=100E-6, eeg=50E-6))
    fig.canvas.key_press_event('a')
    
if read_known_bad_channels:
    # Read bad channels, if they exist
    known_bads = read_channels(fname_channels) # Not implemented

    # Mark existing bad channels
    raw.info['bads'] = known_bads
    
if set_manual_bads:
    print('Choose bad channels manually, then close')
    fig=raw.plot(scalings=dict(eog=100E-6, eeg=50E-6))



#%%

if save_annotations:
    raw.annotations.save(fname_annotations)


# Does not wait for markup of bad channels in interactive plot: need to run manually
if save_channels:
    # All bad channels
    print(raw.info['bads'])
    bads=raw.info['bads']

    # Save bad channels
    write_bad_channels(fname_channels, raw.info, bads)


