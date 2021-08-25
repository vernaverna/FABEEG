#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 21:46:21 2021

@author: mhusberg
"""
import mne
from mne.datasets import eegbci
from mne.io import concatenate_raws, read_raw_edf


# download data and plot ####################################################
subject = 5  # use data from subject 1. Can be in the range of 1-109 (inclusive).
runs = [1,2]  # use only eynes open and eyes closed data

path = '/home/heikkiv/childEEG_data/bids/'
fnames = eegbci.load_data(subject, runs, path=path)
raws = [read_raw_edf(f, preload=True) for f in fnames]
raw = concatenate_raws(raws)

raw.plot()

#events, _ = mne.events_from_annotations(raw, event_id=dict(T1=2, T2=3))

#picks = mne.pick_channels(raw.info["ch_names"], ["C3", "Cz", "C4"])

