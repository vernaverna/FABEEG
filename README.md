This repository containts analysis scripts used for children's sleep-EEG fingerprinting.
The manuscript related to this work is avilable at \url{}

For preprocessing, refer to python-pipelines -folder and for further data analysis and inference, see BRRR -folder. 

The code is structured as follows:
```
├── LICENSE
├── README.md
├── BRRR
│   ├── EEG_fingerprint.R
│   ├── article_plots.R
│   ├── brrr.R
│   ├── codes
│   │   ├── ...
│   ├── dataToR
│   │   └── txt2r.R
│   ├── file_arranging.R
│   ├── plot_components.R
│   ├── run_stats.R
│   ├── var
│   │   └── coords
│   │       ├── ch_names.csv
│   │       └── coords.csv
├── python-pipelines
    ├── check_system.py
    ├── config_common.py
    ├── eeg
    │   ├── Makefile
    │   ├── PSD_04b.py
    │   ├── annotate_bad_channels.py
    │   ├── bandpower_04a.py
    │   ├── config_eeg.py
    │   ├── create_bids.py
    │   ├── dodo.py
    │   ├── download_data.py
    │   ├── freqfilt_01.py
    │   ├── group_level_PSD_plots.py
    │   ├── ica_02.py
    │   ├── plot_singlesubj_spectras.py
    │   ├── psds_03.py
    │   ├── psds_stability_figures.py
    │   ├── runAll.sh
    ├── fnames.py
    └── requirements.txt
```

Running things in python-pipelines:

- if data is not in the bids format, run create_bids.py
- type doit in the terminal; this runs tasks 01-04 for all subjects automatically
- spactra_analysis.py is mainly for plotting 

Running things with BRRR:

- you probably need to rearrange the data with file_arranging.R
- runscripts for BRRR: EEG_fingerprint.R and brrr.R
- for replicating results in the article: run_stats.R, article_plots.R


The code was developed with python v3.10 and mne v1.5 and R v.4.0
