"""
===========
Config file
===========

Configuration parameters for the study.
"""

# Task
task = 'envspeech'

# Sessions
sessions = ['01','02']

# Runs
runs = ['01'] 

# Experimental conditions
conditions = {'speech', 'speech-noise','env','env-noise'}
event_ids = {"speech":1, "speech-noise":2, "env":4, "env-noise":8} 
# Time window (relative to stimulus onset)
epoch_tmin, epoch_tmax = -0.2, 1.5

# Baseline window
baseline = (-0.2,0)

# Filter settings
bandpass_fmin, bandpass_fmax = None, 40

# ICA method
ica_method='infomax'

# ICA time window (wrt stimulus) to find artefact from
ica_epoch_tmin, ica_epoch_tmax =  -0.5, 1.5
