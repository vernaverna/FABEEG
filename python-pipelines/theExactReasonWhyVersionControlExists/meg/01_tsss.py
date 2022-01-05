"""
Perform tSSS Maxwell filtering.
"""
import argparse
from collections import defaultdict

from mne.chpi import (compute_chpi_amplitudes, compute_chpi_locs,
                      compute_head_pos)
from mne.io import read_raw_fif
from mne.preprocessing import maxwell_filter
from mne.viz import plot_head_positions
from mne import open_report

import sys
sys.path.append('../')
from config_meg import fname, get_all_fnames, bads

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', type=int, help='The subject to process')
args = parser.parse_args()

# Along the way, we collect figures for quality control
figures = defaultdict(list)

# Apply maxfilter to each raw file with continuous head position tracking.
# Since not every subject has the same raw files (some subject don't have an
# emptyroom recording for example), we use the `get_all_fnames` function to
# figure out which files are actually present.
raw_fnames = get_all_fnames(args.subject, kind='raw')
tsss_fnames = get_all_fnames(args.subject, kind='tsss')
for raw_fname, tsss_fname in zip(raw_fnames, tsss_fnames):
    raw = read_raw_fif(raw_fname, preload=True)
    raw.pick_types(meg=True, eeg=False, eog=True, ecg=True, stim=True)

    # Mark bad channels that were manually annotated earlier.
    raw.info['bads'] = bads[args.subject]

    # Emptyroom recordings don't have a head, so no head position to track.
    if 'emptyroom' in str(raw_fname):
        head_pos = None
        coord_frame = 'meg'
        destination = None
    else:
        # Compute head positions from the HPI coils
        chpi_amplitudes = compute_chpi_amplitudes(raw)
        chpi_locs = compute_chpi_locs(raw.info, chpi_amplitudes)
        head_pos = compute_head_pos(raw.info, chpi_locs)
        coord_frame = 'head'

        # Add a plot of the head positions to the HTML report for quality
        # control
        figures['head_pos'].append(plot_head_positions(head_pos, show=False))

        # We move the head position to a fixed location, so it is the same for
        # all records. This is needed when we want to average recordings later
        # on.
        destination = (0, 0, 0.04)

    # Perform TSSS maxfilter, compensating for head movements
    tsss = maxwell_filter(
        raw,
        calibration=fname.cal(subject=args.subject),
        cross_talk=fname.ctc(subject=args.subject),
        st_duration=min(120, len(raw) / raw.info['sfreq']),
        head_pos=head_pos,
        coord_frame=coord_frame,
        destination=destination,
    )

    # Create figures to quality check the maxfilter result
    figures['raw'].append(
        raw.plot(
            start=60, lowpass=100, show_scrollbars=False, show_scalebars=False
        )
    )
    figures['tsss'].append(
        tsss.plot(
            start=60, lowpass=100, show_scrollbars=False, show_scalebars=False
        )
    )

    # Make sure the output dir exists
    tsss_fname.parent.mkdir(parents=True, exist_ok=True)
    tsss.save(tsss_fname, overwrite=True)

# Write HTML report with the quality control figures
with open_report(fname.report(subject=args.subject)) as report:
    report.add_slider_to_section(
        figures['head_pos'],
        section='maxfilter',
        title='Head position tracking',
        replace=True,
    )
    report.add_slider_to_section(
        figures['raw'],
        section='maxfilter',
        title='Before maxwell TSSS filter',
        replace=True,
    )
    report.add_slider_to_section(
        figures['tsss'],
        section='maxfilter',
        title='After maxwell TSSS filter',
        replace=True,
    )
    report.save(fname.report_html(subject=args.subject),
                overwrite=True, open_browser=False)
