"""
Examine inter-individual variation and intra-individual stability in PSD.
Plot the power spectral densities on individual level.

@author: heikkiv
"""
import mne
from mne.io import read_info
from h5io import read_hdf5
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm

from config_eeg import subjects, fname

#assuming shorter data seqments
conditions = ['PSD N1a', 'PSD N1b', 'PSD N2a', 'PSD N2b', 'PSD N2c', 'PSD N2d']
subjects = [x.strip('sub-') for x in subjects]

# Load the PSDs for each subject
psds = [read_hdf5(fname.psds(subject=subject))
        for subject in tqdm(subjects)]

# also the supplimentary data 
age_df = pd.read_csv('ages.csv', index_col=1)
age_df = age_df.drop(columns=['Unnamed: 0'])

# Create holder dict 
PSD_dict = {}

for i in range(len(subjects)):
    subj=subjects[i]
    
    try:
        subj_info = age_df.loc[subj].to_dict()
        for cond in conditions: #loop thru conditions
            psds_dict=psds[i]      
            subj_info[cond] = psds_dict[cond]
            
        PSD_dict[subj] = subj_info
        freqs = psds_dict['freqs']
    except:
        print(f'New PSD data missing from subject {subj}')

    
PSD_df = pd.DataFrame.from_dict(PSD_dict, orient='index')

# save the big dataframe to pickle for later use
PSD_df.to_pickle('PSD_dataframe_.pkl')


#%%

def plot_glob_intra_individual(PSD_df, subject, freqs):
    """
    Calculates and plots the global mean PSD of all data segments
    of one subject in a same plot

    Parameters
    ----------
    PSD_df : pandas DataFrame
        PSD data frame
    subject : str
        which subject to plot?
    freqs : array-like
        determines the x-axis

    Returns
    -------
    fig : plt figure

    """
    
    subj_data = PSD_df.loc[subject]
    data = subj_data[3:] #choose only the psd 
    metadata=subj_data['Sex']+ ', '+ str(subj_data['Age'])+'y'
    
    means = [np.mean(datum, axis=0) for datum in data]
    conditions = data.index
    
    means = [np.mean(datum, axis=0) for datum in data]
    conditions = data.index  
    
    plot_df = pd.DataFrame(means).T
    plot_df.columns=conditions
    plot_df['Freq. (Hz)'] = freqs
    plot_df = pd.melt(plot_df, ['Freq. (Hz)'], var_name='segment', value_name='Power')
    plot_df['sleep'] = ['N1' if '1' in plot_df['segment'][i] else 'N2' for i in range(len(plot_df))]
    #plot_df = plot_df.set_index('freq')
    
    
    
    # set up figure aesthetics
    cm = plt.get_cmap('viridis')
    colors = [cm(x) for x in np.linspace(0, 1, len(conditions))] 
    palette = {conditions[i]:colors[i] for i in range(len(colors))}
    dashes = [(conditions[i],[]) if '2' in conditions[i] else (conditions[i],[6,2]) for i in range(len(colors)) ]
    dashes=dict(dashes)
    
    fig = sns.lineplot(x='Freq. (Hz)', y='Power', hue='segment', style='sleep', 
                       palette='viridis', data=plot_df)
    
    return fig, metadata


subj=subjects[44]

fig, metadata = plot_glob_intra_individual(PSD_df, subj, freqs)
plt.title(f'Global average PSDs, {subj} ({metadata})')


#data = np.mean([stc.data for stc in stcs], axis=0)
# Function that creates the tiny PSD plots used to create the big topo figure.
def show_func(ax, ch_idx, tmin, tmax, vmin, vmax, ylim):
    ax.plot(ga_psds['freqs'], ga_psds[cond1][ch_idx], color='C0')
    ax.plot(ga_psds['freqs'], ga_psds[cond2][ch_idx], color='C1')

# Function that creates a larger PSD plot for when one of the tiny PSD plots is
# clicked.
def click_func(ax, ch_idx, tmin, tmax, vmin, vmax, ylim, x_label, y_label):
    ax.plot(ga_psds['freqs'], ga_psds[cond1][ch_idx], color='C0', label=cond1)
    ax.plot(ga_psds['freqs'], ga_psds[cond2][ch_idx], color='C1', label=cond2)
    ax.legend()
    ax.set_xlabel('Frequency')
    ax.set_ylabel('PSD')

# Load the channel locations in order to position the tiny PSD plots in the big
# topo figure.
info = read_info(fname.raw(subject="S001", task='restEC', run=1)) #ERROR no fif but edf file instead?
info = mne.pick_info(info, mne.pick_types(info, meg=False, eeg=True, eog=False))
times = [1, 2]
layout = mne.find_layout(info)

# Make the big topo figure
fig = mne.viz.topo._plot_topo(info, times, show_func, click_func, layout,
                              axis_facecolor='white', fig_facecolor='white')
fig.set_size_inches(14, 9)  # Make the figure window larger
plt.savefig(fname.figure_psds)
