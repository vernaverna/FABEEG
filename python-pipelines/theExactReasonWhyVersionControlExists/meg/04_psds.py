"""
Compute the Power Spectral Density (PSD) for each channel.
"""
import argparse

from mne.io import read_raw_fif
from mne.time_frequency import psd_welch
from mne.externals.h5io import write_hdf5
from mne.viz.topo import _plot_topo
from mne import open_report, find_layout, pick_info, pick_types

from config_meg import fname

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', type=int, help='The subject to process')
args = parser.parse_args()

# Compute the PSD for each task
psds = dict()

fmax = 40
n_fft = 512

raw = read_raw_fif(fname.clean(subject=args.subject, task='eyesclosed', run=1),
                   preload=True)
psds['eyesclosed'], freqs = psd_welch(raw, fmax=fmax, n_fft=n_fft)

raw = read_raw_fif(fname.clean(subject=args.subject, task='eyesopen', run=1),
                   preload=True)
psds['eyesopen'], freqs = psd_welch(raw, fmax=fmax, n_fft=n_fft)

raw = read_raw_fif(fname.clean(subject=args.subject, task='pasat', run=1),
                   preload=True)
psds['pasat_run1'], freqs = psd_welch(raw, fmax=fmax, n_fft=n_fft)

raw = read_raw_fif(fname.clean(subject=args.subject, task='pasat', run=2),
                   preload=True)
psds['pasat_run2'], freqs = psd_welch(raw, fmax=fmax, n_fft=n_fft)

psds['info'] = raw.info
psds['freqs'] = freqs

write_hdf5(fname.psds(subject=args.subject), psds, overwrite=True)

# Add a plot contrasting the two pasat runs to the report.
times = [1, 2]
info = pick_info(raw.info, pick_types(raw.info, meg=True))
layout = find_layout(info)


def show_func(ax, ch_idx, tmin, tmax, vmin, vmax, ylim):
    """Create the tiny PSD plots used to create the big figure."""
    ax.plot(psds['freqs'], psds['eyesclosed'][ch_idx], color='C0')
    ax.plot(psds['freqs'], psds['eyesopen'][ch_idx], color='C1')
    ax.plot(psds['freqs'], psds['pasat_run1'][ch_idx], color='C2')
    ax.plot(psds['freqs'], psds['pasat_run2'][ch_idx], color='C3')


# Make the big topo figure
fig = _plot_topo(info, times, show_func, layout=layout,
                 axis_facecolor='white', fig_facecolor='white')
fig.set_size_inches(14, 9)  # Make the figure window larger

with open_report(fname.report(subject=args.subject)) as report:
    report.add_figs_to_section(fig, 'PSDs', section='PSDs', replace=True)
    report.save(fname.report_html(subject=args.subject),
                overwrite=True, open_browser=False)
