This repository containts analysis scripts used for children's sleep-EEG fingerprinting.
For preprocessing, refer to python-pipelines -folder and for further data analysis and inference, see BRRR -folder. 

The manuscript related to this work is avilable at [url]


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

**python-pipelines**

- scripts with suffixes *01.py--*04.py are related to preprocessing
- `plot_singlesubj_spectras.py` is mainly for plotting; `psds_stability_figures.py` contains statistical analysis and within-sleep stage stability estimates 

**BRRR**

- `EEG_fingerprint.R` and `brrr.R` are used for running the BRRR analysis with the data
- functions related to probabilistic inference are sored under the `codes` -folder
- for replicating statistical results and plots, see `run_stats.R` and `article_plots.R`


The code was developed with python v3.10 and mne v1.5 and R v.4.4
