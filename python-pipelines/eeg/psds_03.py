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

age_df = pd.read_csv('ages.csv', index_col=1)



# Up to now, we have always skipped bad channels. To make sure we can average
# the PSDs later, each recording needs to have the same channels defined. So,
# we call `interpolate_bads()` to replace all channels marked "bads" by an
# interpolation of the surrounding channels. This ensures that each PSD object
# has a signal designed for all channels, hence they are all compatible.

# Initialize PSDS
psds = dict()
raw = read_raw_fif(fname.filt(subject=args.subject), preload=True)
raw.pick_types(eeg=True)

# Add a PSD plot to the report.
info = pick_info(raw.info, pick_types(raw.info, eeg=True))
layout = find_layout(info)
subj_info = age_df.loc[args.subject]

#Add cap_info to age_df
#age_df.at[subj_info.index, 'Cap'] = cap_info
#age_df.to_csv('ages.csv', index=False)



#Making evoked arrays: comments + new info structure 
ag = subj_info['Age']
se = subj_info['Sex']
sfreq = raw.info['sfreq']

f_freq = n_fft/sfreq #frequency resolution
info1 = mne.create_info(info['ch_names'], ch_types=["eeg"]*19, sfreq=f_freq)
info1.set_montage(raw.get_montage())


# TODO: think this through. would sliding window make more sense?
# maybe some 

# Create events of 30 s 
events_n1 = mne.make_fixed_length_events(raw, id=1, start=0, stop=300.0, duration=10, overlap=5) 
events_n2 = mne.make_fixed_length_events(raw, id=2, start=300.0, stop=900.0, duration=10, overlap=5)
events = np.append(events_n1, events_n2, axis=0) #this is clumsy, but did not come up with anything else
event_dict = {'sleep N1':1, 'sleep N2':2}


# Create epochs from events
epochs = mne.Epochs(raw, events, event_id=event_dict, tmin=0.0, tmax=30, baseline=(0,0)) 


time_indices = {'PSD N1 (1)' : range(6,19), # make into 1 min slice -> needs 12 seqments 
                'PSD N1 (2)' : range(20,33),
                'PSD N2 (1)' : range(2,15),
                'PSD N2 (2)' : range(30,43),
                'PSD N2 (3)' : range(10,12),
                'PSD N2 (4)' : range(12,14)}


# Create evoked responses, but as spectra
evokeds = dict()
for key in time_indices.keys():
    
    if 'N1' in key:
        spectra, freqs = psd_array_welch(epochs['sleep N1'][time_indices[key]].get_data(), 
                                         sfreq=sfreq, fmin=1, fmax=fmax, n_fft=n_fft)
        comment = f'Subj: {args.subject}, Age: { ag }, Sex: { se }, Sleep: N1'
    else:
        spectra, freqs = psd_array_welch(epochs['sleep N2'][time_indices[key]].get_data(), 
                                         sfreq=sfreq, fmin=1, fmax=fmax, n_fft=n_fft)
        comment = f'Subj: {args.subject}, Age: { ag }, Sex: { se }, Sleep: N2'
    
    evokeds[key] = mne.EvokedArray(spectra.mean(axis=0), info=info1, comment=comment)
    psds[key] = np.log10(spectra.mean(axis=0))
    
# Add some metadata to the file we are writing
psds['info'] = info1
psds['freqs'] = freqs

write_hdf5(fname.psds(subject=args.subject), psds, overwrite=True)  # save psd

del raw


def callback(ax, ch_idx):
     """Create a larger PSD plot for when one of the tiny PSD plots is
        clicked."""
     ax.plot(psds['freqs'], psds['PSD N1 (1)'][ch_idx], color='C0',
             label='sleep N1')
     ax.plot(psds['freqs'], psds['PSD N2 (2)'][ch_idx], color='C1',
             label='sleep N2')

     ax.legend()
     ax.set_xlabel('Frequency')
     ax.set_ylabel('PSD')


# Make the big topo figure
# TODO: label naming!!!!!
fig = plt.figure(figsize=(14, 9))

for ax, ch_idx in iter_topography(info1, layout, on_pick=callback, fig=fig,
                       axis_facecolor='white', fig_facecolor='white',
                       axis_spinecolor='white'):

    handles = [
        ax.plot(psds['freqs'], psds['PSD N1 (1)'][ch_idx], color='C0', label='sleep N1'),
        ax.plot(psds['freqs'], psds['PSD N2 (2)'][ch_idx], color='C1', label='sleep N2')
    ]
    
#fig.legend("N1 sleep", "N2 sleep")
#fig.show()


# Added to manually save the channel-wise plots
# NOTE: dropping freqs over 20.5 Hz to have more informative figs, hence the selection [0:x] 
# TODO: remove reference, get 18 figs?

figs2 = []
captions = []
for ch_idx in range(len(info.ch_names)):
    fig2 = plt.figure(figsize=(10,7))
    
    plt.plot(psds['freqs'], psds['sleep N1 (1)'][ch_idx].T, color='C0',
                                            label='sleep N1')
    plt.plot(psds['freqs'], psds['sleep N2 (2)'][ch_idx].T, color='C1',
                                            label='sleep N2')
    plt.yscale('log')
    plt.legend()
    captions.append((info.ch_names[ch_idx]))
    plt.xlabel('Frequency')
    plt.ylabel('PSD')

        
        
    figs2.append(fig2)
    plt.close(fig2)


# # Save resultus to report
# with open_report(fname.report(subject=args.subject)) as report:
#     report.add_figs_to_section(fig, 'PSDs', section='PSDs', replace=True, 
#                                comments='Age: {}'.format( float(subj_info['Age'])) )
#     report.add_slider_to_section(figs2, captions=captions, title='PSDs per channel', section='PSDs', replace=True)
#     report.save(fname.report_html(subject=args.subject),
#                 overwrite=True, open_browser=False)
