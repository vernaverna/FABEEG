#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 15:35:46 2021

@author: heikkiv
"""

import mne
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from collections import defaultdict
import numpy as np
from config_eeg import age_df, fname, subjects
from mne.time_frequency import psd_multitaper, psd_welch
from mne.stats import bootstrap_confidence_interval
from mne.baseline import rescale
from mne import open_report

x = age_df['Age']

# logbins = np.logspace(np.log10(x.min()), np.log10(x.max()), 10) #pretty meh
bins = np.array([0.05, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 18.7])
# determine to which bin/age group subject belongs to
age_groups = np.digitize(x, bins)
fig = plt.figure()
plt.hist(x, bins=bins)
fig.show()
age_df['Age group'] = age_groups

figures = defaultdict(list)

# TODO: download .fif spectras

# Loop through all filenames, then determine the age group from age_df.
# Add the data of the subject to corresponding group holder, from which
# the average is then calculated. OVER ALL CHANNELS OR PER CHANNEL?
subjects = [x.strip('.edf') for x in subjects]
subjects = [x.strip('sub-') for x in subjects]
group_spectra_n1 = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: []}
group_spectra_n2 = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: []}


for subj in subjects:

    try:
        age_group = int(age_df.loc[age_df['File'] == subj]['Age group'])
        evoked_spectra = mne.read_evokeds(fname.psds(subject=subj)) #accidentallly in psds: moved to evoked later

        evoked_N1, evoked_N2 = evoked_spectra[0], evoked_spectra[1]

        group_spectra_n1[age_group].append(evoked_N1)
        group_spectra_n2[age_group].append(evoked_N2)

    except:
        continue



cm = plt.get_cmap('viridis')
colors = [cm(x) for x in np.linspace(0, 1, 8)] 

cm = plt.get_cmap('jet')
colors = [cm(x) for x in np.linspace(0, 1, 19)] 


figs = []
captions = []


for i in range(1, len(group_spectra_n1)+1):

    group_evoked_avg = mne.grand_average(group_spectra_n2[i]) #averaged spectra per channels
    freqs = group_evoked_avg.times
    
    for j in range(len(group_evoked_avg.ch_names)):
        fig = plt.figure(figsize=(10, 7))
        #psds, freqs = psd_multitaper(group_evoked_avg, fmin=1, fmax=40, n_jobs=1)
        #psds = 10. * np.log10(psds)
        psd_avg = np.log10(group_evoked_avg.data[j])
        psds_mean = psd_avg.mean(0) #mean over all channels 
        psds_std = psd_avg.std(0)
    
        plt.plot(freqs, psds_mean, color=colors[i],
                 label='PSD - Channel average')
    
        ci_low, ci_up = bootstrap_confidence_interval(psd_avg, random_state=0)
        plt.fill_between(freqs, ci_low, ci_up,
                         color=colors[i], alpha=.5)
    
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Power Spectral Density (dB)')
    
        agebin = 'Ch: ' + group_evoked_avg.ch_names[j] + 'Age: ' + str(bins[i-1]) + '-' + str(bins[i]) + ' years'
    
        captions.append(agebin)
        figs.append(fig)


# Define regions of interest (ROIs)
# TODO: Confirm these!!! 
roi_viz = ['P3', 'P4', 'O1', 'O2', 'PZ']
roi_temp_R = ['T4', 'T6', 'C4']
roi_temp_L = ['T3', 'T5', 'C3']
roi_front = ['FP1', 'FP2', 'FZ', 'F3', 'F4', 'F7', 'F8']

picks_list = [roi_viz, roi_temp_R, roi_temp_L, roi_front]

roi_figs=[]
captions=[]

for roi in picks_list:
    evokeds = []
    for i in range(1, len(group_spectra_n1)):
        evokeds.append(mne.grand_average(group_spectra_n2[i])) #change comment!
    
    fig = mne.viz.plot_compare_evokeds(evokeds, combine='mean', picks=roi, 
                                       colors=[1,2,3,4,5,6,7,8], cmap='viridis',
                                       show_sensors=True, legend=True, 
                                       title='N2 sleep spectra of age groups')
    roi_figs.append(fig)

roi_figs = sum(roi_figs, []) #hacky


#TODO: make a channel-wise slider!



caps = ['ROI 1', 'ROI 2', 'ROI 3', 'ROI 4']

with open_report('Average_spectra.h5') as report:
    report.add_slider_to_section(figs, captions=captions, title='Averaged N2 PSDs estimates',
                                 section='Channel averages: N2 sleep', replace=True)
    #report.add_figs_to_section(roi_figs, captions=caps, 
    #                           section='Mean spectra of groups', replace=True)
    report.save('Average_spectra.html', overwrite=True, open_browser=True)


# TODO: Compute group-averaged spectras

# TODO: Plot results (distributions; violinplots etc)
