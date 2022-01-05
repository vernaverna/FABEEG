"""
Perform bandpass filtering and notch filtering to get rid of cHPI and powerline frequencies.
"""
import argparse
from collections import defaultdict

from mne.io import read_raw_fif
from mne.chpi import filter_chpi
from mne import open_report

from config_meg import get_all_fnames, fname

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', type=int, help='The subject to process')
args = parser.parse_args()

# Along the way, we collect figures for quality control
figures = defaultdict(list)

# Not all subjects have files for all conditions. These functions grab the
# files that do exist for the subject.
tsss_fnames = get_all_fnames(args.subject, kind='tsss')
filt_fnames = get_all_fnames(args.subject, kind='filt')

for tsss_fname, filt_fname in zip(tsss_fnames, filt_fnames):
    tsss = read_raw_fif(tsss_fname, preload=True)

    # Add a plot of the power spectrum to the list of figures to be placed in
    # the HTML report.
    figures['before_filt'].append(tsss.plot_psd())

    # Emptyroom recordings don't have a head, so no HPI coils. By default,
    # `filter_chpi` will give an error when no HPI coils are present. Setting
    # `allow_line_only` will allow the function to only filter out power line
    # noise.
    filt = filter_chpi(tsss, include_line=True, allow_line_only=True)

    # Highpass filter above 1Hz. This is needed for the ICA to perform well
    # later on.
    filt.filter(1, None)

    # Save the filtered data
    filt.save(filt_fname, overwrite=True)

    # Add a plot of the power spectrum of the filtered data to the list of
    # figures to be placed in the HTML report.
    figures['after_filt'].append(filt.plot_psd())

# Write HTML report with the quality control figures
with open_report(fname.report(subject=args.subject)) as report:
    report.add_slider_to_section(
        figures['before_filt'],
        section='freq filter',
        title='Before frequency filtering',
        replace=True,
    )
    report.add_slider_to_section(
        figures['after_filt'],
        section='freq filter',
        title='After frequency filtering',
        replace=True,
    )
    report.save(fname.report_html(subject=args.subject),
                overwrite=True, open_browser=False)
