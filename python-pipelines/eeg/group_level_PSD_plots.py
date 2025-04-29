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
from mne.time_frequency import psd_array_welch
from mne.stats import bootstrap_confidence_interval
from mne.baseline import rescale
from mne import open_report

x = age_df['Age']
#bins = np.array([0.05, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 18.7])

#create even-sized bins
bin_labels = ['1','2','3','4','5','6','7','8','9','10']
#bins = pd.qcut(age_df['Age'], q=10, labels=bin_labels)
bins, bin_labs = pd.qcut(age_df['Age'], q=10, retbins=True, labels=bin_labels)

age_df['Age group'] = bins

figures = defaultdict(list)

# TODO: download .fif spectras

# Loop through all filenames, then determine the age group from age_df.
# Add the data of the subject to corresponding group holder, from which
# the average is then calculated. OVER ALL CHANNELS OR PER CHANNEL?
subjects = [x.strip('.edf') for x in subjects]
subjects = [x.strip('sub-') for x in subjects]
group_spectra_n1 = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9:[], 10: []}
group_spectra_n2 = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9:[], 10: []}


for subj in subjects:

    try:
        age_group = int(age_df.loc[age_df['File'] == subj]['Age group'])
        evoked_spectra = mne.read_evokeds(fname.evoked(subject=subj)) #accidentallly in psds: moved to evoked later

        evoked_N1, evoked_N2 = evoked_spectra[1], evoked_spectra[4]

        group_spectra_n1[age_group].append(evoked_N1)
        group_spectra_n2[age_group].append(evoked_N2)

    except:
        continue



#Take 1 example subject for sensor locations
evk_file =  group_spectra_n1[2][1]

# Get sensor locations from layout
layout_from_evk = mne.channels.make_eeg_layout(evk_file.info)
locations = layout_from_evk.pos #4 columns; x, y, width, height
locations = locations[:,0:2] #take only x & y coords
np.savetxt('coords.csv', locations, delimiter=',') #these are needed for later analysis


# Set colormap for visuals 
cm = plt.get_cmap('viridis')
colors = [cm(x) for x in np.linspace(0, 1, 11)] 

#cm = plt.get_cmap('jet')
#colors = [cm(x) for x in np.linspace(0, 1, 19)] 


figs = []
captions = []


for i in range(1, len(group_spectra_n1)+1):
    
    #note: changed to 1 person only
    group_evoked_avg_n1 = group_spectra_n1[i][1] #averaged spectra per channels over "evoked" responses
    group_evoked_avg_n2 = group_spectra_n2[i][1] 
    # group_evoked_avg_n1 = mne.grand_average(group_spectra_n1[i]) #averaged spectra per channels over "evoked" responses
    # group_evoked_avg_n2 = mne.grand_average(group_spectra_n2[i]) 
    freqs = group_evoked_avg_n2.times
    chs = group_evoked_avg_n1.info['ch_names']
    
    fig, ax = plt.subplots(figsize=(6, 4))
    psd_avg_n1 = np.log10(group_evoked_avg_n1.data[8,]) # one channel only
    psd_avg_n2 = np.log10(group_evoked_avg_n2.data[8,])
    
    # psd_avg_n1 = np.log10(group_evoked_avg_n1.data).mean(axis=0) #global average
    # psd_avg_n2 = np.log10(group_evoked_avg_n2.data).mean(axis=0)

    ax.plot(freqs, psd_avg_n1, color=colors[i-1], linestyle='-', label='N1b Sleep')
    ax.plot(freqs, psd_avg_n2, color=colors[10//i], linestyle='-', label='N2c Sleep') 

    agebin = 'Age: ' + str(bin_labs[i-1]) + '-' + str(bin_labs[i]) + ' years'
    #agebin = 'Age: ' + str(bins[i-1]) + ' years'
    
    #plot confidence intervals
    # c1_low, c1_up = bootstrap_confidence_interval(np.log10(group_evoked_avg_n1.data), random_state=9)
    # ax.fill_between(freqs, c1_low, c1_up, color=colors[i-1], alpha=.3)
    # c2_low, c2_up = bootstrap_confidence_interval(np.log10(group_evoked_avg_n2.data), random_state=9)
    # ax.fill_between(freqs, c2_low, c2_up, color=colors[i-1], alpha=.5)
    
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Log Power Spectral Density (dB)')
    ax.legend()
    #plt.title('Global average PSDs '+agebin)
    plt.title('PSD of channel ' + chs[8] )
    
    figs.append(fig)



#%% OTher stuff - ROI plots

# Define regions of interest (ROIs)
# TODO: Confirm these!!! 
roi_viz = ['P3', 'P4', 'O1', 'O2', 'PZ']
roi_temp = ['T3', 'T5', 'C3', 'T4', 'T6', 'C4']
roi_front = ['FP1', 'FP2', 'FZ', 'F3', 'F4', 'F7', 'F8']

picks_list = [roi_viz, roi_temp, roi_front]
caps = ['ROI Visual', 'ROI Temporal', 'ROI Frontal']

roi_figs=[]
log_n1_spect = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9:[], 10: []}

for roi in picks_list[1]:
    evokeds = []
    captions=[]
   
    for i in range(1, len(group_spectra_n1)+1):
        if i < 12:
            spectra_n1 = [evk.copy().pick_channels([roi]) for evk in group_spectra_n1[i]]
            inf = spectra_n1[1].info
            log_group_spectra = [np.log10(evk.data) for evk in spectra_n1]
            log_evkds = [mne.EvokedArray(data, info=inf) for data in log_group_spectra]
            
            log_n1_spect[i] = log_evkds
            
            evkd_group = mne.grand_average(log_evkds)
            freqs = evkd_group.times
            comment= 'Age: '+str(bin_labs[i-1]) + '-' + str(bin_labs[i]) + ' y'
            captions.append(comment)
            evokeds.append(evkd_group.data.mean(0)) 
            
    evk_df = pd.DataFrame(evokeds).T
    evk_df.columns = captions
    evk_df.index = freqs
    
    g = sns.relplot(data=evk_df, kind="line", palette=colors[0:10])
    g.fig.suptitle("N1: channel " + roi)
    

    roi_figs.append(g)

roi_figs = sum(roi_figs, []) #hacky way to get a 1D list of figures

####################################################################
#FIGURE 1A: Global averages 

glob_figs=[]
log_n1_spect = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9:[], 10: []}
log_n2_spect = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9:[], 10: []}


evokeds_n1 = []
evokeds_n2 = []
captions=[]
   
for i in range(1, len(group_spectra_n1)+1):
    if i < 12:
        spectra_n1 = [evk.copy() for evk in group_spectra_n1[i]]
        spectra_n2 = [evk.copy() for evk in group_spectra_n2[i]]
        
        to_discard1= [i for i, e in enumerate(spectra_n1) if np.sum(e.data)==0.0 or (np.isnan(e.data)).all()]
        to_discard2= [i for i, e in enumerate(spectra_n2) if np.sum(e.data)==0.0 or (np.isnan(e.data)).all()]
        
        for index in sorted(to_discard2, reverse=True):
            del spectra_n2[index]
        
        inf = spectra_n1[1].info
        log_group_spectra = [np.log10(evk.data) for evk in spectra_n1]
        log_evkds = [mne.EvokedArray(data, info=inf) for data in log_group_spectra]
        
        
        
        log_group_spectra2 = [np.log10(evk.data) for evk in spectra_n2]
        log_evkds2 = [mne.EvokedArray(data, info=inf) for data in log_group_spectra2]
        
        log_n1_spect[i] = log_evkds
        log_n2_spect[i] = log_evkds2
        
        evkd_group = mne.grand_average(log_evkds)
        evkd_group2 = mne.grand_average(log_evkds2)
        freqs = evkd_group.times
        comment= 'Age: '+str(round(bin_labs[i-1],2)) + '-' + str(round(bin_labs[i],2)) + ' y'
        captions.append(comment)
        evokeds_n1.append(evkd_group.data.mean(0))
        evokeds_n2.append(evkd_group2.data.mean(0))


# https://seaborn.pydata.org/generated/seaborn.relplot.html
evk_df1 = pd.DataFrame(evokeds_n1).T
evk_df1.columns = captions
evk_df1.index = freqs

evk_df2 = pd.DataFrame(evokeds_n2).T
evk_df2.columns = captions
evk_df2.index = freqs

frames=[evk_df1, evk_df2]

evk_df = pd.concat(frames, keys=['N1', 'N2'])
evk_df.reset_index(inplace=True)
evk_df = evk_df.rename(columns={'level_0': 'Sleep', 'level_1': 'Frequency'})

# make the figures separately for both sleep stage labels
g = sns.relplot(data=evk_df1, kind="line", palette=colors[0:10])
g.fig.suptitle("N1 Global average")
glob_figs.append(g)

f = sns.relplot(data=evk_df2, kind="line", palette=colors[0:10])
f.fig.suptitle("N2 Global average")
glob_figs.append(f)

#list(evk_df.columns)
col_palette = sns.color_palette("viridis",10)
cols=[col_palette[3], col_palette[8]]
cols = ['#DCE319FF', '#33638DFF']

h = sns.relplot(data=evk_df, x="Frequency", y="Age: 0.83-1.28 y", hue="Sleep",
                palette=cols, kind='line')
h.fig.suptitle("Group global averages, Age: 0.83-1.28 y")
(h.set_axis_labels("Frequency (Hz)", "log PSD"))
glob_figs.append(h)

# In order to use plot_compare_evokeds, change the group spectra dict keys
# to string

# keys_n_values = group_spectra_n1.items()
# group_spect_n1 = {str(key): value for key, value in keys_n_values}

keys_n_values = log_n1_spect.items()
log_n1_spect = {str(key): value for key, value in keys_n_values}

fig = mne.viz.plot_compare_evokeds(log_n1_spect, combine='mean',  
   cmap='viridis', ylim=dict(eeg=[-1.4e+07,-9.0e+06]), ci=False,
   show_sensors=False, legend=True,
   title='N1 sleep spectra of age groups')


####################################################################

with open_report('Average_spectra.h5') as report:
    report.add_slider_to_section(figs, captions=captions, title='Averaged N2 PSDs estimates',
                                 section='Channel averages: N2 sleep', replace=True)
    #report.add_figs_to_section(roi_figs, captions=caps, 
    #                           section='Mean spectra of groups', replace=True)
    report.save('Average_spectra.html', overwrite=True, open_browser=True)


