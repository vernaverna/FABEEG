This repository containts analysis scripts used for children's sleep-EEG fingerprinting.
The manuscript related to this work is avilable at \url{}

For preprocessing, refer to python-pipelines -folder and for further data analysis in BRRR -folder. 


Running things in python-pipelines:

- if data is not in the bids format, run create_bids.py
- type doit in the terminal; this runs tasks 01-04 for all subjects automatically
- spactra_analysis.py is mainly for plotting 

Running things with BRRR:

- you probably need to rearrange the data with file_arranging.R
- runscripts for BRRR: EEG_fingerprint.R and brrr.R
- for replicating results in the article: run_stats.R, article_plots.R


The code was developed with python v3.10 and mne v1.5 and R v.4.0
