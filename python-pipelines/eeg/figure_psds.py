"""
Compute the Power Spectral Density (PSD) for each channel.
"""
import mne
from mne.io import read_info
from mne.externals.h5io import read_hdf5
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

from config_eeg import subjects, fname

cond1 = 'restEO'
cond2 = 'restEC'

# Load the PSDs for each subject
psds = [read_hdf5(fname.psds(subject=subject))
        for subject in tqdm(subjects)]
#psds = [read_hdf5(fname.psds(subject='S001'))]

# Create the grand-average (GA) PSD
ga_psds = {
    cond1: np.mean([psd_subj[cond1] for psd_subj in psds], axis=0),
    cond2: np.mean([psd_subj[cond2] for psd_subj in psds], axis=0),
    'freqs': psds[0]['freqs']
}
#data = np.mean([stc.data for stc in stcs], axis=0)
# Function that creates the tiny PSD plots used to create the big topo figure.
def show_func(ax, ch_idx, tmin, tmax, vmin, vmax, ylim):
    ax.plot(ga_psds['freqs'], ga_psds[cond1][ch_idx], color='C0')
    ax.plot(ga_psds['freqs'], ga_psds[cond2][ch_idx], color='C1')

# Function that creates a larger PSD plot for when one of the tiny PSD plots is
# clicked.
def click_func(ax, ch_idx, tmin, tmax, vmin, vmax, ylim, x_label, y_label):
    ax.plot(ga_psds['freqs'], ga_psds[cond1][ch_idx], color='C0', label=cond1)
    ax.plot(ga_psds['freqs'], ga_psds[cond2][ch_idx], color='C1', label=cond2)
    ax.legend()
    ax.set_xlabel('Frequency')
    ax.set_ylabel('PSD')

# Load the channel locations in order to position the tiny PSD plots in the big
# topo figure.
info = read_info(fname.raw(subject="S001", task='restEC', run=1)) #ERROR no fif but edf file instead?
info = mne.pick_info(info, mne.pick_types(info, meg=False, eeg=True, eog=False))
times = [1, 2]
layout = mne.find_layout(info)

# Make the big topo figure
fig = mne.viz.topo._plot_topo(info, times, show_func, click_func, layout,
                              axis_facecolor='white', fig_facecolor='white')
fig.set_size_inches(14, 9)  # Make the figure window larger
plt.savefig(fname.figure_psds)
