"""
Compute the Power Spectral Density (PSD) for each channel.
Can be processed in python cell writing
import subprocess
subprocess.run('/net/theta/fishpool/projects/FABEEG/python-pipelines/eeg/runAll.sh', shell=True)
"""
import argparse

import mne
from mne.io import read_raw_fif
from h5io import write_hdf5


import pandas as pd
import numpy as np
from mne.time_frequency import psd_array_welch
from mne.viz import iter_topography
from mne import open_report, find_layout, pick_info, pick_types
import matplotlib.pyplot as plt

from config_eeg import fname, n_fft, fmax

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', help='The subject to process')
args = parser.parse_args()

age_df = pd.read_csv('ages.csv', index_col=0)

# calculate 1-min estimates or longer segmnets?
long_psd = False

# Initialize PSDS
psds = dict()
raw = read_raw_fif(fname.filt(subject=args.subject), preload=True)
raw.pick_types(eeg=True)

# Add a PSD plot to the report.
info = pick_info(raw.info, pick_types(raw.info, eeg=True))
layout = find_layout(info)
try:
    subj_info = age_df.loc[age_df['File']==args.subject]
    
    #Add cap_info to age_df
    #age_df.at[subj_info.index, 'Cap'] = cap_info
    #age_df.to_csv('ages.csv', index=False)
    
    #Making evoked arrays: comments + new info structure 
    ag = subj_info['Age']
    se = subj_info['Sex']

except:
    print('No metadata from subject {args.subject}!')
    ag = '?'
    se = '?'
    subj_info = pd.DataFrame({'File': [args.subject], 'Sex': ['nan'], 'Age': ['nan'], 'Cap': ['nan']})
    age_df=pd.concat([age_df,subj_info], ignore_index=True)
    
    #update datafrAME
    age_df.to_csv('ages.csv')

sfreq = raw.info['sfreq']
f_freq = n_fft/sfreq #frequency resolution
info1 = mne.create_info(info['ch_names'], ch_types=["eeg"]*19, sfreq=f_freq)
info1.set_montage(raw.get_montage())


# TODO: think this through. would sliding window make more sense?
# maybe some 

# Create events of 30 s 
events_n1 = mne.make_fixed_length_events(raw, id=1, start=0, stop=300.0, duration=30, overlap=0) 
events_n2 = mne.make_fixed_length_events(raw, id=2, start=300.0, stop=900.0, duration=30, overlap=0)
events = np.append(events_n1, events_n2, axis=0) #this is clumsy, but did not come up with anything else
event_dict = {'sleep N1':1, 'sleep N2':2}


# Create epochs from events
epochs = mne.Epochs(raw, events, event_id=event_dict, tmin=0.0, tmax=30, baseline=(0,3)) 

#TODO: want global N1 & N2, as well as 1-min/30s slices? globs for plotting. 

if long_psd:
    time_indices = {'PSD N1' : range(1,10), #match data set size
                    'PSD N2' : range(1,10)}

else:
    time_indices = {'PSD N1a' : range(1,5),
                    'PSD N1b' : range(5,9),
                    'PSD N2a' : range(1,5),
                    'PSD N2b' : range(5,9),
                    'PSD N2c' : range(9,14),
                    'PSD N2d' : range(14,19)}

del raw

# Create evoked responses, but as spectra
evokeds = dict()
for key in time_indices.keys():
    
    if 'N1' in key:
        spectra, freqs = psd_array_welch(epochs['sleep N1'][time_indices[key]].get_data(), 
                                         sfreq=sfreq, fmin=1, fmax=fmax, n_fft=n_fft)
        
        comment = f'Subj: {args.subject}, Age: { ag }, Sex: { se }, Sleep: N1'
        spectra = spectra*1e12 #scaling for EEG_PSD
        
        evokeds[key] = mne.EvokedArray(spectra.mean(axis=0), info=info1, comment=comment)
        psds[key] = 10*np.log10(spectra.mean(axis=0)) #get decibels
        
    else:
        try:
            spectra, freqs = psd_array_welch(epochs['sleep N2'][time_indices[key]].get_data(), 
                                             sfreq=sfreq, fmin=1, fmax=fmax, n_fft=n_fft)
            
            comment = f'Subj: {args.subject}, Age: { ag }, Sex: { se }, Sleep: N2'
            
            spectra = spectra*1e12
            evokeds[key] = mne.EvokedArray(spectra.mean(axis=0), info=info1, comment=comment)
            psds[key] = 10*np.log10(spectra.mean(axis=0)) #get decibels
        
        except IndexError as e:
            print('Not enough data for last segment!')
            continue
    

    
# Add some metadata to the file we are writing
psds['info'] = info1
psds['freqs'] = freqs

if long_psd:
    print("TODO: make a nested dictionary -> dataframe -> pickle")
else:
    write_hdf5(fname.psds(subject=args.subject), psds, overwrite=True)  # save psd

#ELE12390
#FLE131126
#ELE12390
#FLE131126


# # Save resultus to report
# with open_report(fname.report(subject=args.subject)) as report:
#     report.add_figs_to_section(fig, 'PSDs', section='PSDs', replace=True, 
#                                comments='Age: {}'.format( float(subj_info['Age'])) )
#     report.add_slider_to_section(figs2, captions=captions, title='PSDs per channel', section='PSDs', replace=True)
#     report.save(fname.report_html(subject=args.subject),
#                 overwrite=True, open_browser=False)
