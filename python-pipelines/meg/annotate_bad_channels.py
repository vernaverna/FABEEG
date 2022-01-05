"""
This script is used for manually annotation bad EEG channels.
The idea is to run this in an ipython console through:

>>> %run annotate_bad_channels.py <subject number>

The script will load the raw data and make a plot of it. Inside the plot you
can scroll through the data and click on the channel names to mark them as bad.
When you close the plot, the script will print out the channels you have
marked. These channels can then be added to the big `bads` list inside
config_meg.py.
"""
import argparse

import mne
from mne.preprocessing import find_bad_channels_maxwell

from config_meg import fname

# Deal with command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', type=int, help='The subject to process')
args = parser.parse_args()

raw = mne.io.read_raw_fif(fname.raw(subject=args.subject, task='pasat', run=1), preload=True)
raw.filter(1, 100)
bads, flats = find_bad_channels_maxwell(raw)
raw.info['bads'] = bads + flats
raw.plot(scalings=dict(eog=100E-6))
print(raw.info['bads'])
