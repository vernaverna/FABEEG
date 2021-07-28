#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script governs the process of EEG filtering & cerates power spectra
of all ~800 subjects
Or should this be done as a bashscript?

The process goes as follows:
    
    0. run config_eeg.py and create_bids.py if not yet done
    1. annotate_bad_channels.py -> HOW to make feasible for all subj?
    2. 01_freqfilt.py -> plotting issue, no coords
    3. 02_ica.py -> same problem
    4. 03_psds.py -> same problem
    5. figure_psds.py
    6. save spectra for each each individual in their own file/folder

  
Created on Mon Jun  7 15:28:23 2021

@author: heikkiv
"""

