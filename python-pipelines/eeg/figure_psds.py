"""
Examine inter-individual variation and intra-individual stability in PSD.
Plot the power spectral densities on individual level.

@author: heikkiv
"""
import mne
import os
from mne.io import read_info
from mne import find_layout
from h5io import read_hdf5
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from mne.viz import iter_topography
from config_eeg import fname

data_dir = '/net/theta/fishpool/projects/FABEEG/childEEG_data/bids/derivatives/'

#assuming shorter data seqments
conditions = ['PSD N1a', 'PSD N1b', 'PSD N2a', 'PSD N2b', 'PSD N2c', 'PSD N2d']

all_subjects = os.listdir(data_dir)
subjects = [x.strip('sub-') for x in all_subjects]

# Load the PSDs for each subject
psds = [read_hdf5(fname.psds(subject=subject))
        for subject in tqdm(subjects)]

# also the supplimentary data 
age_df = pd.read_csv('ages.csv', index_col=1)
age_df = age_df.drop(columns=['Unnamed: 0'])
age_df = age_df[~age_df.index.duplicated(keep='first')]

# Create holder dict 
PSD_dict = {}

for i in range(len(subjects)):
    subj=subjects[i]  
    try:
        subj_info = age_df.loc[subj].to_dict()
        for cond in conditions: #loop thru conditions
            psds_dict=psds[i]
            
            try:
                subj_info[cond] = psds_dict[cond]
            except KeyError:
                subj_info[cond] = np.nan #some have empty N2 sequences
                
        PSD_dict[subj] = subj_info
        freqs = psds_dict['freqs']
        info = psds_dict['info']
    
    except KeyError:
        print(f'PSD data missing from subject {subj}')

PSD_df = pd.DataFrame.from_dict(PSD_dict, orient='index')

# save the big dataframe to pickle for later use
PSD_df.to_pickle('PSD_dataframe_.pkl')
chs = info['ch_names']

#%% PLOTTING FUNCTIONS

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
    
    N1_df = plot_df.loc[plot_df['sleep']=='N1']
    N2_df = plot_df.loc[plot_df['sleep']=='N2']
    #plot_df = plot_df.set_index('freq')
    
    
    
    # set up figure aesthetics - use two palettes rather than one?
    #cm = sns.color_palette("crest", 6)
    cm_1 = sns.color_palette("dark:orangered", 2)
    cm_2 = sns.color_palette("dark:#5A9_r", 4)
    #colors = [cm(x) for x in np.linspace(0, 1, len(conditions))] 
    #palette = {conditions[i]:colors[i] for i in range(len(colors))}
    #dashes = [(conditions[i],[]) if '2' in conditions[i] else (conditions[i],[6,2]) for i in range(len(colors)) ]
    #dashes=dict(dashes)
    
    fig, ax = plt.subplots(figsize=(8,6)) 
    sns.lineplot(x='Freq. (Hz)', y='Power', hue='segment', #style='sleep', 
            palette=cm_1, data=N1_df)
    sns.lineplot(x='Freq. (Hz)', y='Power', hue='segment', #style='sleep', 
            palette=cm_2, data=N2_df, linestyle='dashed')
    
    # fig = sns.lineplot(data=plot_df, x='Freq. (Hz)', y='Power', hue='segment', 
    #                    style='sleep',  palette=cm)
    
    return fig, metadata



def sd_mean_inter_age_groups(PSD_df, freqs, sleep='PSD N1a'):
    """
    
    Parameters
    ----------
    PSD_df : pandas DataFrmae
        PSD data frame
    freqs : frequencies 
        needed fot x-axis
    sleep : str. optional
        which sleep stage to plot?. The default is 'PSD N1a'.

    Returns
    -------
    cohort_n_mean: list of mean and SD data per age group
    """
    plot_df = PSD_df[['Sex', 'Age', 'Cap', sleep]]
    
    #create age groups: even-sized bins #TODO: move this to outer layer?
    bin_labels = ['1','2','3','4','5','6','7','8','9','10']
    bins, bin_labs = pd.qcut(plot_df['Age'], q=10, retbins=True, labels=bin_labels)
    bin_names = [str( round(bin_labs[i-1],1)) + ' – ' + str(round(bin_labs[i],1)) for i in range(1,len(bin_labs))]
    
    bins, bin_labs = pd.qcut(plot_df['Age'], q=10, retbins=True, labels=bin_names)
    
    plot_df['Age group'] = bins
    
    
    cohrt_groups = plot_df.groupby('Age group')
    
    group_ch_means = []
    group_ch_sd = []
    names = [] 
    for name, cohort_df in cohrt_groups:
        cohort_data = cohort_df[sleep]
        cohort_data.dropna(inplace=True) #get rid of nan -values
        Arr = np.array([i for i in cohort_data]) #nsubj x n_chs x n_freq
        mean_data = np.nanmean(Arr, axis=0) #get mean data over all subjects within cohort
        sd_data= np.nanstd(Arr, axis=0)
        group_ch_means.append(mean_data)
        group_ch_sd.append(sd_data)
        names.append(name)
    
    cohort_n_mean1 = zip(names, group_ch_means, group_ch_sd)
    cohort_n_mean = [(name, psd_data, sd_data) for name, psd_data, sd_data in cohort_n_mean1]
    
    return cohort_n_mean
    



#%%% Plots
sleep='PSD N1b'
cohort_n_mean =  sd_mean_inter_age_groups(PSD_df, freqs, sleep=sleep)

def my_callback1(ax, ch_idx):
    """
    This block of code is executed once you click on one of the channel axes
    in the plot. To work with the viz internals, this function should only take
    two parameters, the axis and the channel or data index.
    """
    
    # Set colormap for visuals 
    cm = plt.get_cmap('viridis')
    colors = [cm(x) for x in np.linspace(0, 1, len(cohort_n_mean))] 

    i=0
    for name, mean_data, sd_data in cohort_n_mean:
        ax.plot(freqs, mean_data[ch_idx], color=colors[i], label=name) #for all the group averages
        #lower =  mean_data[ch_idx]-sd_data[ch_idx]
        #upper =  mean_data[ch_idx]+sd_data[ch_idx]
        #ax.fill_between(freqs, lower, upper, color=colors[i], alpha=.15)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Power (dB)')
        ax.legend(loc="upper right", title='Age (in years)')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #ax.set_title(sleep) #overrides channel name :(
        i=i+1


template = mne.io.read_raw_fif(fname.filt(subject=subj), preload=True)
#loop through all channels and create an axis for them
for ax, idx in iter_topography(template.info, 
                               fig_facecolor='white',
                               axis_facecolor='white',
                               axis_spinecolor='white',
                               on_pick=my_callback1):
    ax.plot(cohort_n_mean[0][1][idx], color='green') #just to show some general output for the big figure
    
plt.gcf().suptitle(f'Power spectral densities in {sleep}')
plt.show()


### single-subject plots

subj=subjects[5]

fig, metadata = plot_glob_intra_individual(PSD_df, subj, freqs)
plt.title(f'Global average PSDs, {subj} ({metadata})')


