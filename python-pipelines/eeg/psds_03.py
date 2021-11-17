"""
Compute the Power Spectral Density (PSD) for each channel.
Can be processed in python cell writing
import subprocess
subprocess.run('/net/theta/fishpool/projects/FABEEG/python-pipelines/eeg/runAll.sh', shell=True)
"""
import argparse

import mne
from mne.io import read_raw_fif
import pandas as pd
import numpy as np
from mne.time_frequency import psd_welch
from mne.externals.h5io import write_hdf5
from mne.viz import iter_topography
from mne import open_report, find_layout, pick_info, pick_types
import matplotlib.pyplot as plt

from config_eeg import fname, n_fft, age_df, fmax

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', help='The subject to process')
args = parser.parse_args()


#Function to change the metadata & pick common channels
def change_metadata(raw): 

    #For testdata only, remove . in channel-names
    raw.rename_channels(lambda x: x.strip('.'))  # remove dots from channel names, FIXME
    raw.rename_channels(lambda x: x.upper())  # capitalize the ch names
    
    # Get reference channel; either Pz or Cz
    ch_info_df = raw.describe(data_frame=True)
    cz_deviation = abs(ch_info_df.iloc[17, 6]-ch_info_df.iloc[17, 4]) #Q3-Q1
    pz_deviation = abs(ch_info_df.iloc[18, 6]-ch_info_df.iloc[18, 4])
    
    if cz_deviation < pz_deviation: #which channel has lower deviation ~> zero signal
        reference = 'CZ'
    else:
        reference = 'PZ'
    
    #Set the reference as well
    #raw.set_eeg_reference(ref_channels=[reference])
    raw.set_eeg_reference('average')

    
    #Drop channels that are not in the shared_channels.txt .file
    allowed_chs = []
    with open('shared_channels.txt') as f:
         allowed_chs = f.read().splitlines()
    
    for ch_name in raw.info['ch_names']:
        if ch_name not in allowed_chs:
            raw.drop_channels(ch_name)

    # Drop photic & pietso (for now)
    raw.drop_channels("PHOTIC")
    raw.drop_channels("PIETSO")
    raw.drop_channels("EKG")
    
    # TODO: Reset sensor types ?
    
    # Sensor locations not provided with data - use standard layout
    ten_twenty_montage = mne.channels.make_standard_montage('standard_1020')
    ten_twenty_montage.ch_names = [CH_NAME.upper() for CH_NAME in ten_twenty_montage.ch_names]

    raw = raw.copy().set_montage(ten_twenty_montage)
    
    return raw



# Up to now, we have always skipped bad channels. To make sure we can average
# the PSDs later, each recording needs to have the same channels defined. So,
# we call `interpolate_bads()` to replace all channels marked "bads" by an
# interpolation of the surrounding channels. This ensures that each PSD object
# has a signal designed for all channels, hence they are all compatible.

# Initialize PSDS
psds = dict()
raw = read_raw_fif(fname.filt(subject=args.subject), preload=True)
#Add metadata to testdata
raw = change_metadata(raw)
#raw.interpolate_bads() # Only works if location is known


# Add a PSD plot to the report.
info = pick_info(raw.info, pick_types(raw.info, eeg=True))
layout = find_layout(info)
subj_info = age_df.loc[age_df['File']==args.subject]


# Making evoked arrays: comments

comment1 = 'Subj: {}, Age: {}, Sex: {}, Sleep: N1'.format(str(subj_info.iloc[0,0]), float(subj_info['Age']),
                                                        str(subj_info.iloc[0,1]))
comment2 = 'Subj: {}, Age: {}, Sex: {}, Sleep: N2'.format(str(subj_info.iloc[0,0]), float(subj_info['Age']),
                                                        str(subj_info.iloc[0,1]))



# Create events of 30 s  
events_n1 = mne.make_fixed_length_events(raw, id=1, start=0, stop=300.0, duration=30, overlap=0) 
events_n2 = mne.make_fixed_length_events(raw, id=2, start=300.0, stop=900.0, duration=30, overlap=0)
events = np.append(events_n1, events_n2, axis=0) #this is clumsy, but did not come up with anything else
event_dict = {'sleep N1':1, 'sleep N2':2}


# Create epochs from events
epochs = mne.Epochs(raw, events, event_id=event_dict, tmin=0.0, tmax=30, baseline=(0,0)) 


# Create evoked responses, but as spectra
evokeds = dict()
#evokeds['N1'] = epochs['sleep N1'].average() #averages over each epoch
#evokeds['N2'] = epochs['sleep N2'].average()



# TODO: make this neater! (now 1 min slices)
info1 = info
n1_spectra, freqs = psd_welch(epochs['sleep N1'][8:10].average(), fmax=fmax, n_fft=n_fft)
n1_spectra2, freqs = psd_welch(epochs['sleep N1'][5:7].average(), fmax=fmax, n_fft=n_fft)
info1['sfreq'] = 1/(freqs[1]-freqs[0]) #change sampling freq so the 'time' axis shows frequencies
evokeds['PSD N1 (1)'] = mne.EvokedArray(n1_spectra, info=info1, comment=comment1)
evokeds['PSD N1 (2)'] = mne.EvokedArray(n1_spectra2, info=info1, comment=comment1)


# TODO: do this inside of a loop!!

# leaving 240 s between n1 & n2
# do the same with n2 spectra; create three 120 s spectras (no overlap between these)
#n2_spectra1, freqs = psd_welch(epochs['sleep N2'][7:12].average(), fmax=fmax, n_fft=n_fft)
#evokeds['PSD N2 (1)'] = mne.EvokedArray(n2_spectra1, info=info1, comment=comment2)
#n2_spectra2, _ = psd_welch(epochs['sleep N2'][14:19].average(), fmax=fmax, n_fft=n_fft)
#evokeds['PSD N2 (2)'] = mne.EvokedArray(n2_spectra1, info=info1, comment=comment2)
#n2_spectra3, _ = psd_welch(epochs['sleep N2'][21:26].average(), fmax=fmax, n_fft=n_fft)
#evokeds['PSD N2 (3)'] = mne.EvokedArray(n2_spectra1, info=info1, comment=comment2)

# create 5 spectras from N2
n2_spectra1, freqs = psd_welch(epochs['sleep N2'][2:4].average(), fmax=fmax, n_fft=n_fft)
evokeds['PSD N2 (1)'] = mne.EvokedArray(n2_spectra1, info=info1, comment=comment2)
n2_spectra2, freqs = psd_welch(epochs['sleep N2'][6:8].average(), fmax=fmax, n_fft=n_fft)
evokeds['PSD N2 (2)'] = mne.EvokedArray(n2_spectra2, info=info1, comment=comment2)
n2_spectra3, freqs = psd_welch(epochs['sleep N2'][10:12].average(), fmax=fmax, n_fft=n_fft)
evokeds['PSD N2 (3)'] = mne.EvokedArray(n2_spectra3, info=info1, comment=comment2)
n2_spectra4, freqs = psd_welch(epochs['sleep N2'][12:14].average(), fmax=fmax, n_fft=n_fft)
evokeds['PSD N2 (4)'] = mne.EvokedArray(n2_spectra4, info=info1, comment=comment2)
n2_spectra5, freqs = psd_welch(epochs['sleep N2'][16:18].average(), fmax=fmax, n_fft=n_fft)
evokeds['PSD N2 (5)'] = mne.EvokedArray(n2_spectra5, info=info1, comment=comment2)


fmin, fmax = freqs[0], freqs[-1]
psds['sleep N1 (1)'] = n1_spectra
psds['sleep N1 (2)'] = n1_spectra2
psds['sleep N2 (1)'] = n2_spectra1
psds['sleep N2 (2)'] = n2_spectra2
psds['sleep N2 (3)'] = n2_spectra3
psds['sleep N2 (4)'] = n2_spectra4
psds['sleep N2 (5)'] = n2_spectra5

# Add some metadata to the file we are writing
psds['info'] = raw.info
psds['freqs'] = freqs
del raw #free up some memory

# save psds and evokeds
write_hdf5(fname.psds(subject=args.subject), psds, overwrite=True) 
mne.write_evokeds(fname.evoked(subject=args.subject), 
                  [evokeds['PSD N1 (1)'], evokeds['PSD N1 (2)'], evokeds['PSD N2 (1)'], 
                   evokeds['PSD N2 (2)'], evokeds['PSD N2 (3)'],
                   evokeds['PSD N2 (4)'], evokeds['PSD N2 (5)']] ) 



def on_pick(ax, ch_idx):
    """Create a larger PSD plot for when one of the tiny PSD plots is
       clicked."""
    ax.plot(psds['freqs'], psds['sleep N1 (1)'][ch_idx], color='C0',
            label='sleep N1')
    ax.plot(psds['freqs'], psds['sleep N2 (2)'][ch_idx], color='C1',
            label='sleep N2')

    ax.legend()
    ax.set_xlabel('Frequency')
    ax.set_ylabel('PSD')


# Make the big topo figure
# TODO: label naming!!!!!
fig = plt.figure(figsize=(14, 9))
axes = iter_topography(info, layout, on_pick=on_pick, fig=fig,
                       axis_facecolor='white', fig_facecolor='white',
                       axis_spinecolor='white')

for ax, ch_idx in axes:
    handles = [
        ax.plot(psds['freqs'], psds['sleep N1 (1)'][ch_idx], color='C0', label='sleep N1'),
        ax.plot(psds['freqs'], psds['sleep N2 (2)'][ch_idx], color='C1', label='sleep N2')
    ]
    
#fig.legend("N1 sleep", "N2 sleep")
#fig.show()


#Added to manually save the channel-wise plots
#NOTE: dropping freqs over 20.5 Hz to have more informative figs, hence the selection [0:x] 
#TODO: remove reference, get 18 figs?

figs2 = []
captions = []
for ch_idx in range(len(info.ch_names)):
    fig2 = plt.figure(figsize=(10,7))
    
    plt.plot(psds['freqs'], psds['sleep N1 (1)'][ch_idx], color='C0',
                                            label='sleep N1')
    plt.plot(psds['freqs'], psds['sleep N2 (2)'][ch_idx], color='C1',
                                            label='sleep N2')
    plt.yscale('log')
    plt.legend()
    captions.append((info.ch_names[ch_idx]))
    plt.xlabel('Frequency')
    plt.ylabel('PSD')

        
        
    figs2.append(fig2)
    plt.close(fig2)


# Save resultus to report
with open_report(fname.report(subject=args.subject)) as report:
    report.add_figs_to_section(fig, 'PSDs', section='PSDs', replace=True, 
                               comments='Age: {}'.format( float(subj_info['Age'])) )
    report.add_slider_to_section(figs2, captions=captions, title='PSDs per channel', section='PSDs', replace=True)
    report.save(fname.report_html(subject=args.subject),
                overwrite=True, open_browser=False)
