# -*- coding: utf-8 -*-
"""
This script is used for plotting the PSDs of a single subject
for quality-control purposes

@heikkiv
"""

import mne
from mne.time_frequency import psd_welch
import matplotlib.pyplot as plt
import numpy as np
from config_eeg import fname, n_fft, age_df
from mne import open_report


def change_metadata(raw): 
    """
    Modify metadata, set correct channel types, and pick common channels 
    from the subject's raw EEG data.

    This function harmonizes raw EEG data by renaming channels, determining 
    the appropriate reference channel (Pz or Cz), and retaining only the 
    channels listed in a shared channels file. It also identifies whether 
    a larger cap was used based on the presence of specific channel names.

    Parameters
    ----------
    raw : mne.io.Raw
        The raw EEG data in FIF format.

    Returns
    -------
    raw : mne.io.Raw
        The processed raw EEG data with modified metadata and selected 
        channels.
    cap_status : str
        Indicates whether the larger cap was used. Possible values:
        - 'FT': Larger cap was used.
        - '-': No larger cap was used.
    """

    #For testdata only, remove . in channel-names
    raw.rename_channels(lambda x: x.strip('.'))  # remove dots from channel names, FIXME
    raw.rename_channels(lambda x: x.upper())  # capitalize the ch names
    
    # Get reference channel; either Pz or Cz
    ch_info_df = raw.describe(data_frame=True)
    cz_deviation = abs(ch_info_df.iloc[17, 6]-ch_info_df.iloc[17, 4]) #Q3-Q1
    pz_deviation = abs(ch_info_df.iloc[18, 6]-ch_info_df.iloc[18, 4])
    
    if cz_deviation < pz_deviation:
        reference = 'CZ'
    else:
        reference = 'PZ'
    
    #Set the reference as well
    raw.set_eeg_reference(ref_channels=[reference])
    
    
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
    
    # TO DO Reset sensor types
    
    # Sensor locations not provided with testdata - use standard layout
    ten_twenty_montage = mne.channels.make_standard_montage('standard_1020')
    ten_twenty_montage.ch_names = [CH_NAME.upper() for CH_NAME in ten_twenty_montage.ch_names]

    raw = raw.copy().set_montage(ten_twenty_montage)
    
    return raw




## PLot age distribution

plt.hist(age_df['Age'], bins=10)

plt.xlabel('Age in years')
plt.ylabel('Frequency')
plt.title('CHILDEEG age distribution')

# result: ~poisson-distribution



filenames = ['FLE12515', 'FLE12431', 'FLE131198', 'FLE141276', 'FLE141752',
             'FLE151099', 'FLE151028']


data_dir='/net/theta/fishpool/projects/FABEEG/childEEG_data/bids/sub-'+filenames[4]+'/eeg/sub-'+filenames[4]+'_task-sleep_eeg.edf'
raw=mne.io.read_raw_edf(data_dir, preload=True)
raw.plot()

# Reading in the data
raw2 = mne.io.read_raw_edf(fname.raw(subject=filenames[-2]), preload=True)
raw2 = change_metadata(raw2)


# Bandpass filtering 
raw2.filter(1,40)
raw2 = raw2.notch_filter(50, picks=['eeg', 'eog'])


# Set reference
raw2.rename_channels(lambda x: x.upper())  # capitalize the ch names
raw2.set_eeg_reference('average')
info = mne.pick_info(raw2.info, mne.pick_types(raw2.info, eeg=True))

raw2.plot(duration=30.0)

plt.savefig("raw_example.svg")

# Define ROIs 
#TODO: Confirm these!!! 
roi_viz = ['P3', 'P4', 'O1', 'O2']
roi_temp = ['T3', 'T4', 'T5', 'T6']
roi_front = ['FP1', 'FP2', 'FZ', 'F3', 'F4', 'F7', 'F8']




# Create events of 30 s  
events_n1 = mne.make_fixed_length_events(raw2, id=1, start=0, stop=300.0, duration=30) #TODO: define overlap?
events_n2 = mne.make_fixed_length_events(raw2, id=2, start=300.0, stop=900.0, duration=30)

events = np.append(events_n1, events_n2, axis=0) #this is clumsy, but did not come up with anything else
event_dict = {'sleep N1':1, 'sleep N2':2}

raw2.plot(events=events, start=270, duration=50, color='gray',
         event_color={1: 'r', 2: 'g'})


# Create epoch from events
epochs = mne.Epochs(raw2, events, event_id=event_dict, tmin=0.0, tmax=30, baseline=(0,0)) 
epochs['sleep N1'].plot_psd(fmax=40, picks='eeg') #visualize
epochs['sleep N1'].plot_psd_topomap() #plots also topomaps



# Create evoked responses, but as spectra
evokeds = dict()
evokeds['N1'] = epochs['sleep N1'].average() #averages over the sensor
evokeds['N2'] = epochs['sleep N2'].average()

n2_spectra, freqs = psd_welch(epochs['sleep N2'].average(), fmax=50, n_fft=n_fft)
info['sfreq']=4.096
evokeds['PSD N2'] = mne.EvokedArray(n2_spectra, info=info, comment='testing')



#some sanity-check visualizations; all ok
mne.viz.plot_compare_evokeds(epochs, cmap=('word length', 'viridis'),
                             picks='eeg', combine='gfp')


evokeds['PSD N2'].plot(picks='eeg', spatial_colors=True, gfp=True)

#try resampling
psdn2_resamp = evokeds['PSD N2'].resample(sfreq=4.096)
psdn2_resamp.plot(picks='eeg', spatial_colors=True, gfp=True)

# Create 30s epochs???
epochs = mne.make_fixed_length_epochs(raw2_ref, duration=30, preload=False) #make events
epochs_n1 = epochs[:10] # first 10 epochs are N1 sleep
epochs_n2 = epochs[10:] # last 20 epochs are N2 sleep


#fig = raw2_ref.plot(duration=10)
#raw2.info['ch_names'][19:]

# TODO: Confirm these!!! 
roi_viz = ['P3', 'P4', 'O1', 'O2']
roi_temp = ['T3', 'T4', 'T5', 'T6']
roi_front = ['FP1', 'FP2', 'FZ', 'F3', 'F4', 'F7', 'F8']


picks = mne.pick_types(raw2_ref, eeg=True, selection=roi_viz) #select a subset of channels



fig = plt.figure()
raw_n2 = raw2.copy().crop(tmin=400.0, tmax=600.0) #rest 600s are N2 sleep 
psds, freqs = psd_welch(raw_n2, fmax=40, n_fft=1024)
plt.plot(freqs, psds[8],
         label='Multitaper PSD - Global average')
plt.yscale('log')

plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectral Density')
fig.show()

#raw2.describe(data_frame=True) #sensor information as a data frame; use these to pick reference



montage = mne.channels.make_standard_montage('standard_1020') #what options are there???
info = mne.create_info(ch_names=montage.ch_names, sfreq=100., ch_types='eeg')
#info.set_montage(montage)

ten_twenty_montage = mne.channels.make_standard_montage('standard_1020')
#raw_1020 = raw2.copy().set_montage(ten_twenty_montage)
fig2 = ten_twenty_montage.plot(kind='topomap', show_names=True)
#fig2.show()


with open_report(fname.report(subject='ELE12632')) as report:
    report.add_figs_to_section(fig, 'Raw data view', section='Raw data', replace=True)
    report.add_figs_to_section(fig2, 'Topomap', section='PSDs', replace=True)
    report.save(fname.report_html(subject='ELE12632'),
                overwrite=True, open_browser=True)



raw2.plot_psd(fmax=50) #estimate


fmax=40
n_fft=2048
psds, freqs = psd_welch(raw2, fmax=fmax, n_fft=n_fft)

plt.plot(psds[15,:])

#store the channel names; includes only EEG





