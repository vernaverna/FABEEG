"""
Remove EOG artifacts through independant component analysis (ICA).
No ECG artifacts are removed at the moment, since they are hard to detect based
on just EEG data.
"""
import argparse
from collections import defaultdict

from mne import Epochs, channels
from mne.io import read_raw_fif
from mne.preprocessing import find_eog_events, ICA, create_eog_epochs, create_ecg_epochs, corrmap
from mne import open_report

from config_eeg import get_all_fnames, task_from_fname, fname, eog_channel

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', help='The subject to process')
parser.add_argument('remove_eog', help='Try to remove also EOG artifacts? True/False', default='False')
args = parser.parse_args()

# Along the way, we collect figures for quality control
figures = defaultdict(list)


# Not all subjects have files for all conditions. These functions grab the
# files that do exist for the subject.

all_fnames = zip(
    get_all_fnames(args.subject, kind='filt'),
    get_all_fnames(args.subject, kind='ica'),
    get_all_fnames(args.subject, kind='clean'),
)

for filt_fname, ica_fname, clean_fname in all_fnames:
    #task = task_from_fname(filt_fname)

    raw_filt = read_raw_fif(filt_fname, preload=True)
    #raw_filt.crop(tmax=300) #only N1 data?
    
    ecg_evoked = create_ecg_epochs(raw_filt).average()
    
    # Perform ICA decomposition
    ica = ICA(n_components=0.99, random_state=0).fit(raw_filt)
    bads_ecg, scores_ecg = ica.find_bads_ecg(raw_filt) 
    
    if args.remove_eog=='True':
        eog_evoked = create_eog_epochs(raw_filt).average()
        eog_evoked.apply_baseline(baseline=(None, -0.2))
        bads_eog, scores_eog = ica.find_bads_eog(eog_evoked)
   
    else:
        bads_eog = []
    # Mark the ECG/EOG components for removal
    bads = bads_ecg + bads_eog
    ica.exclude = bads
    ica.save(ica_fname, overwrite=True)

    # Remove the EC/OG artifact components from the signal.
    raw_ica = ica.apply(raw_filt.copy())
    raw_ica.save(clean_fname, overwrite=True)

    # Put a whole lot of quality control figures in the HTML report.
    
    with open_report(fname.report(subject=args.subject)) as report:
        report.add_ica(
        ica=ica,
        title="ICA cleaning",
        picks=bads,  # only plot the first two components
        inst=raw_filt,
        ecg_evoked=ecg_evoked,
        ecg_scores=scores_ecg,
        )

        report.save(fname.report_html(subject=args.subject),
                    overwrite=True, open_browser=True)
