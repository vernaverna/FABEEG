"""
===========
Config file
===========

Configuration parameters for the study.
"""

import os
from getpass import getuser
from socket import gethostname

###############################################################################
# Determine which user is running the scripts on which machine and set the path
# where the data is stored and how many CPU cores to use.

user = getuser()  # Username of the user running the scripts
host = gethostname()  # Hostname of the machine running the scripts

# You want to add your machine to this list
#if host == 'nbe-024' and user == 'vanvlm1':
if host.startswith('vdiubuntu') and user == 'mhusberg':
    # VDI machine
    raw_data_dir = '/m/nbe/scratch/rubberboot/test/bids/'
    processed_data_dir = '/m/nbe/scratch/rubberboot/test/bids/derivatives/'
    reports_dir = os.path.join('/m/nbe/scratch/rubberboot/test/reports',user)
    figures_dir = os.path.join('/m/nbe/scratch/rubberboot/test/figures',user)
    n_jobs = 2  # VDI gives you 2 cores
    matplotlib_backend = 'Qt5Agg'
elif host.endswith('triton.aalto.fi'):
    # Triton cluster
    raw_data_dir = '/m/nbe/scratch/rubberboot/test/bids/'
    processed_data_dir = '/m/nbe/scratch/rubberboot/test/bids/derivatives/'
    reports_dir = os.path.join('/m/nbe/scratch/rubberboot/test/reports',user)
    figures_dir = os.path.join('/m/nbe/scratch/rubberboot/test/figures',user)
    n_jobs = 1
    matplotlib_backend = 'Agg'  # No graphics on triton
elif host == 'nbe-065' and user == 'hkoivikk':
    # Hanna's workstation
    raw_data_dir = '/m/nbe/scratch/rubberboot/test/bids/'
    processed_data_dir = '/m/nbe/scratch/rubberboot/test/bids/derivatives/'
    reports_dir = os.path.join('/m/nbe/scratch/rubberboot/test/reports',user)
    figures_dir = os.path.join('/m/nbe/scratch/rubberboot/test/figures',user)
    n_jobs = 4
    matplotlib_backend = 'Qt5Agg'
elif host == 'altair' and user == 'heikkiv' :
    # Verna's workstation in BioMag
    raw_data_dir = '/net/theta/fishpool/projects/FABEEG/childEEG_data/bids/'
    processed_data_dir = '/net/theta/fishpool/projects/FABEEG/childEEG_data/bids/derivatives/'
    reports_dir = os.path.join('/net/theta/fishpool/projects/FABEEG/childEEG_data/reports',user)
    figures_dir = os.path.join('/net/theta/fishpool/projects/FABEEG/childEEG_data/figures',user)
    n_jobs = 4
    matplotlib_backend = 'Qt5Agg'
else:
    raise ValueError(f'Please enter the details of your system ({user}@{host}) in config_common.py')

# For BLAS to use the right amount of cores
os.environ['OMP_NUM_THREADS'] = str(n_jobs)

# Configure the graphics backend
import matplotlib
matplotlib.use(matplotlib_backend)
