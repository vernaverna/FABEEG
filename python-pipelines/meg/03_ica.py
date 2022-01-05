"""
Remove EOG and ECG artifacts through independant component analysis (ICA)
"""
import argparse
from collections import defaultdict

from mne import Epochs
from mne.io import read_raw_fif
from mne.preprocessing import find_ecg_events, find_eog_events, ICA
from mne import open_report

from config_meg import get_all_fnames, task_from_fname, fname

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', type=int, help='The subject to process')
args = parser.parse_args()

# Along the way, we collect figures for quality control
figures = defaultdict(list)

# Not all subjects have files for all conditions. These functions grab the
# files that do exist for the subject.
all_fnames = zip(
    # emptyroom doesn't have eye blinks or heart beats
    get_all_fnames(args.subject, kind='raw', exclude='emptyroom'),
    get_all_fnames(args.subject, kind='filt', exclude='emptyroom'),
    get_all_fnames(args.subject, kind='ica', exclude='emptyroom'),
    get_all_fnames(args.subject, kind='clean', exclude='emptyroom'),
)

for raw_fname, filt_fname, ica_fname, clean_fname in all_fnames:
    task = task_from_fname(raw_fname)

    # Run a detection algorithm for the onsets of eye blinks (EOG) and heart
    # beats (ECG). This is best performed on non-TSSS filtered data.
    raw = read_raw_fif(raw_fname, preload=True)
    if task != 'eyesclosed':  # No blinks during eyes closed
        eog_events = find_eog_events(raw)
    ecg_events = find_ecg_events(raw)[0]
    del raw

    # Perform ICA on the filtered data
    raw_filt = read_raw_fif(filt_fname, preload=True)
    ica = ICA(0.95, random_state=0).fit(raw_filt)

    # Find components that are likely capturing EOG and ECG artifacts
    if task == 'eyesclosed' in str(raw_fname):  # No blinks during eyes closed
        bads_eog = []
    else:
        eog_epochs = Epochs(raw_filt, eog_events, tmin=-0.5, tmax=0.5,
                            preload=True)
        bads_eog, scores_eog = ica.find_bads_eog(eog_epochs)
    ecg_epochs = Epochs(raw_filt, ecg_events, tmin=-0.5, tmax=0.5,
                        preload=True)
    bads_ecg, scores_ecg = ica.find_bads_ecg(ecg_epochs)

    ica.exclude = list(set(bads_eog + bads_ecg))
    ica.save(ica_fname)

    # Remove the EOG and ECG artifact components from the signal.
    raw_ica = ica.apply(raw_filt)
    raw_ica.save(clean_fname, overwrite=True)

    # Put a whole lot of quality control figures in the HTML report.
    with open_report(fname.report(subject=args.subject)) as report:
        if 'eyesclosed' not in str(raw_fname):  # No blinks during eyes closed
            report.add_figs_to_section(
                ica.plot_scores(scores_eog, exclude=bads_eog, show=False),
                f'{task}: EOG scores', section='ICA', replace=True)

            report.add_figs_to_section(
                ica.plot_overlay(eog_epochs.average(), show=False),
                f'{task}: EOG overlay', section='ICA', replace=True)

            if len(bads_eog) == 1:
                report.add_figs_to_section(
                    ica.plot_properties(eog_epochs, bads_eog, show=False),
                    [f'{task}: Component {i:02d}' for i in bads_eog],
                    section='ICA', replace=True)
            elif len(bads_eog) > 1:
                report.add_slider_to_section(
                    ica.plot_properties(eog_epochs, bads_eog, show=False),
                    captions=[f'{task}: Component {i:02d}' for i in bads_eog],
                    title=f'{task}: EOG component properties', section='ICA',
                    replace=True)

        report.add_figs_to_section(
            ica.plot_scores(scores_ecg, exclude=bads_ecg, show=False),
            f'{task}: ECG scores', section='ICA', replace=True)

        report.add_figs_to_section(
            ica.plot_overlay(ecg_epochs.average(), show=False),
            f'{task}: ECG overlay', section='ICA', replace=True)

        if len(bads_ecg) == 1:
            report.add_figs_to_section(
                ica.plot_properties(ecg_epochs, bads_ecg, show=False),
                [f'{task}: Component {i:02d}' for i in bads_ecg],
                section='ICA', replace=True)
        elif len(bads_ecg) > 1:
            report.add_slider_to_section(
                ica.plot_properties(ecg_epochs, bads_ecg, show=False),
                captions=[f'Component {i:02d}' for i in bads_ecg],
                title=f'{task}: ECG component properties', section='ICA',
                replace=True)

        report.save(fname.report_html(subject=args.subject),
                    overwrite=True, open_browser=False)
