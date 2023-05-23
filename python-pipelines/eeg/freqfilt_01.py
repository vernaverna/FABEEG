"""
Perform bandpass filtering and notch filtering to get rid of cHPI and powerline
frequencies.
"""
# Some relevant files are in the parent folder.
#import sys
#sys.path.append('../../python-pipelines')

import argparse
from collections import defaultdict

from mne.io import read_raw_edf
from mne import open_report

from config_eeg import get_all_fnames, fname, bads, fmin, fmax, fnotch

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', help='The subject to process')
args = parser.parse_args()

# Along the way, we collect figures for quality control
figures = []
caps = ['Before', 'After']

# Not all subjects have files for all conditions. These functions grab the
# files that do exist for the subject.
raw_fnames = get_all_fnames(args.subject, kind=None)
filt_fnames = get_all_fnames(args.subject, kind='filt')

for raw_fname, filt_fname in zip(raw_fnames, filt_fnames):
    
    try:
        raw = read_raw_edf(raw_fname, preload=True)
   
        # Remove MEG channels. This is the EEG pipeline after all.
        raw.pick_types(meg=False, eeg=True, eog=True, stim=True)
    
        # Mark bad channels that were manually annotated earlier.
        
        if args.subject in bads.keys(): ### bads is empty!!
            raw.info['bads'] = bads[args.subject]  
        else:
            raw.info['bads'] = []
    
        # Add a plot of the power spectrum to the list of figures to be placed in
        # the HTML report.
        #TODO: Change the fig. y axis scaling! current one is very uninformative 
        #      is prob. due to the dead 'Photic' channel, so take it out! (Non-fatal tho)   
        figures.append(raw.compute_psd(picks=['eeg']).plot(show=False))
    
        # Remove 50Hz power line noise (and the first harmonic: 100Hz)
        filt = raw.notch_filter(fnotch, picks=['eeg', 'eog'])
    
        # Apply bandpass filter 
        filt = raw.filter(fmin, fmax, picks=['eeg', 'eog'])
    
        # Save the filtered data
        filt_fname.parent.mkdir(parents=True, exist_ok=True)
        filt.save(filt_fname, overwrite=True)
    
        # Add a plot of the power spectrum of the filtered data to the list of
        # figures to be placed in the HTML report.
        figures.append(filt.compute_psd(picks=['eeg']).plot(show=False))
        
        #add raw segment figure
        raw1 = raw.pick_types(eeg=True, eog=False, stim=False).crop(tmin=260,tmax=320).load_data()
    
        # Write HTML report with the quality control figures
        #TODO: looks stupid because of reference. consider raw data insrtead?
        with open_report(fname.report(subject=args.subject)) as report: #CHANGED
            report.add_figure( #slider to figs
                fig = figures,
                title='Filtering',
                caption=caps,
                replace=True,
            )
            report.add_raw(raw=raw1, title='Raw', psd=False)  # omit PSD plot

            
            report.save(fname.report_html(subject=args.subject),
                        overwrite=True, open_browser=False)
    
    except FileNotFoundError:
        print(f'could not open raw data file from {args.subject}')
