Codes for children's sleep-EEG data-analysis.
The preprocessing scripts are in python-pipelines -folder
and further data analysis in BRRR -folder. 


Running things in python-pipelines:

- if data is not in the bids format, run create_bids.py
- type doit in the terminal; this runs tasks 01-04 for all subjects automatically
- spactra_analysis.py is mainly for plotting 

Running things with BRRR:

- you probably need to rearrange the data with file_arranging.R
- runscripts for BRRR: EEG_fingerprint.R and 
