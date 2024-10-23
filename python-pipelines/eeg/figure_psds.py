"""
Examine inter-individual variation and intra-individual stability in PSD.
Plot the power spectral densities on individual level.

@author: heikkiv
"""
import mne
import os
import scipy
from mne.io import read_info
from mne import find_layout
from h5io import read_hdf5
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from mne.viz import iter_topography
from mne.stats import permutation_cluster_test
from mne.channels import find_ch_adjacency
from mne.stats import combine_adjacency, spatio_temporal_cluster_test
from config_eeg import fname, bandpower, f_bands, read_template_subject

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
            freqs = psds_dict['freqs']

            cond_bpw = cond.replace('PSD', 'BPW')
            try:
                subj_info[cond] = psds_dict[cond]
                
                data_bandpower = []
                
                for band in f_bands:
                     fmin, fmax = band[0], band[1]
                     bandpwr = bandpower(psd=psds_dict[cond], f=freqs, fmin=fmin, fmax=fmax)
                     data_bandpower.append(bandpwr)
                
                data_bandpower = np.array(data_bandpower)
                #sums over the absolute power of all channels
                abs_power = np.sum(data_bandpower, axis=0)
                #normalizes data
                #if all (abs_power < -1e-5):
                data_bandpower = data_bandpower/abs_power #did not scale with 100 yet
                
                subj_info[cond_bpw] = np.array(data_bandpower).T
            
            
            except KeyError:
                subj_info[cond] = np.nan #some have empty N2 sequences
                
        PSD_dict[subj] = subj_info
        info = psds_dict['info']
    
    except KeyError:
        print(f'PSD data missing from subject {subj}')

PSD_df = pd.DataFrame.from_dict(PSD_dict, orient='index')

# save the big dataframe to pickle for later use
PSD_df.to_pickle('PSD_dataframe_with_bandpower.pkl')
chs = info['ch_names']

#%% PLOTTING FUNCTIONS AND STATS   

def plot_glob_intra_individual(PSD_df, subject, freqs, perform_clustering=False):
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
    perform_clustering : boolean
        perform clustering test? Defaults to false.

    Returns
    -------
    fig : plt figure

    """
    
    subj_data = PSD_df.loc[subject]
    data = subj_data[3:] #choose only the psd 
    metadata=subj_data['Sex']+ ', '+ str(subj_data['Age'])+'y'
    
    means = [np.nanmean(datum, axis=0) for datum in data]
    conditions = data.index  
    
    plot_df = pd.DataFrame(means).T
    plot_df.columns=conditions
    plot_df['Freq. (Hz)'] = freqs
    plot_df = pd.melt(plot_df, ['Freq. (Hz)'], var_name='segment', value_name='Power')
    plot_df['sleep'] = ['N1' if '1' in plot_df['segment'][i] else 'N2' for i in range(len(plot_df))]
    
    N1_df = plot_df.loc[plot_df['sleep']=='N1']
    N2_df = plot_df.loc[plot_df['sleep']=='N2']
    #plot_df = plot_df.set_index('freq')
    
    if perform_clustering:
        # clustering based on mean data over N1 vs N2 segments? or f test?
        N1_data, N2_data = np.array(means[0:2]), np.array(means[2:])
        m,n=N1_data.shape 
        k=N2_data.shape[0]
        data_array = [N1_data.reshape(m,n,1), N2_data.reshape(k,n,1)]
        pval = 0.05  # arbitrary
        dfn = 1 #two conditions - 1
        dfd = len(means) - 2  # degrees of freedom denom.
        thresh = scipy.stats.f.ppf(1 - pval, dfn=dfn, dfd=dfd)  # F distribution
        
        X=data_array #np.array([*data_array]) #unpack data
        F_obs, clusters, cluster_pv, H0 = permutation_cluster_test(
            X,
            n_permutations=1000,
            threshold=thresh,
            tail=0,
            )
    
    # set up figure aesthetics - use two palettes rather than one?
    #cm = sns.color_palette("crest", 6)
    cm = plt.get_cmap('viridis')
    colors = [cm(x) for x in np.linspace(0, 1, 10)] 

    cm_1 = [colors[0], colors[2]]
    cm_2 = [colors[3], colors[4], colors[6], colors[8]]
    #colors = [cm(x) for x in np.linspace(0, 1, len(conditions))] 
    #palette = {conditions[i]:colors[i] for i in range(len(colors))}
    #dashes = [(conditions[i],[]) if '2' in conditions[i] else (conditions[i],[6,2]) for i in range(len(colors)) ]
    #dashes=dict(dashes)
    
    fig, ax = plt.subplots(figsize=(8,6)) 
    sns.lineplot(x='Freq. (Hz)', y='Power', hue='segment', #label='N1', 
            palette=cm_1, data=N1_df)
    sns.lineplot(x='Freq. (Hz)', y='Power', hue='segment', #label='N2', 
            palette=cm_2, data=N2_df, linestyle='dashed')
    sns.despine(fig, bottom=False, left=False)
    fig.suptitle(f'Global average PSDs, {subj} ({metadata})')
    plt.ylabel('log-PSD')
    plt.close()

    fig2, ax = plt.subplots(figsize=(8,6)) 
    sns.lineplot(data=plot_df, x='Freq. (Hz)', y='Power', hue='segment', 
                        style='sleep',  palette=colors[1:10])
    if perform_clustering:
        for i_c, c in enumerate(clusters):
            c=c[0]
            if cluster_pv[i_c] <= 0.05:
                h = ax.axvspan(freqs[c[0]], freqs[c[-1]], color="grey", alpha=0.2)
    else:
        clusters, cluster_pv = None, None
          
    sns.despine(fig2, bottom=False, left=False)
    fig2.suptitle(f'Global average PSDs, {subj} ({metadata})')
    plt.ylabel('log-PSD')
    plt.close()
    
    return fig, fig2, metadata, clusters, cluster_pv



def sd_mean_inter_age_groups(PSD_df, freqs, sleep='PSD N1a', do_bandpower=False):
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
    stat_df : pd DataFrame containing the statics about the overall power betwwn groups
    """
    plot_df = PSD_df.loc[:,['Sex', 'Age', 'Cap', sleep]]
    
    #create age groups: even-sized bins #TODO: move this to outer layer?
    bin_labels = ['1','2','3','4','5','6','7','8','9','10']
    bins, bin_labs = pd.qcut(plot_df['Age'], q=10, retbins=True, labels=bin_labels)
    bin_names = [str( round(bin_labs[i-1],1)) + ' – ' + str(round(bin_labs[i],1)) for i in range(1,len(bin_labs))]
    
    bins, bin_labs = pd.qcut(plot_df['Age'], q=10, retbins=True, labels=bin_names)
    
    plot_df['Age group'] = bins
    AUC_recs = []
    
    cohrt_groups = plot_df.groupby('Age group', observed=False)
    
    group_ch_means = []
    group_ch_sd = []
    names = [] 
    
    
    
    for name, cohort_df in cohrt_groups:
        cohort_data = cohort_df[sleep]
        cohort_data.dropna(inplace=True) #get rid of nan -values
        Arr = np.array([i for i in cohort_data]) #nsubj x n_chs x n_freq
        
        avg_Arr = np.nanmean(Arr, axis=1) #average over channel dimension
        if do_bandpower:
            
            for i in range(avg_Arr.shape[1]):
                bandpowers = avg_Arr[:,i]
                bandnames= [f'band {i}']*avg_Arr.shape[0]
                namevec= [name]*avg_Arr.shape[0]
                
                AUC_recs.append(  list( zip(namevec, bandnames, bandpowers, cohort_df['Sex'].values, cohort_df.index.values, ) ))
                
        else:
            #compute AUC values ("total power overall (from grand average)")
            AUCs = np.trapz(avg_Arr, x=freqs, axis=-1)         
            namevec = [name]*AUCs.shape[0]
            AUC_recs.append(  list( zip(namevec, AUCs, cohort_df['Sex'].values, cohort_df.index.values, ) ))
            
        mean_data = np.nanmean(Arr, axis=0) #get mean data over all subjects within cohort
        sd_data= np.nanstd(Arr, axis=0)
        group_ch_means.append(mean_data)
        group_ch_sd.append(sd_data)
        names.append(name)
    
    cohort_n_mean1 = zip(names, group_ch_means, group_ch_sd)
    cohort_n_mean = [(name, psd_data, sd_data) for name, psd_data, sd_data in cohort_n_mean1]
    
    
    # make AUC df for later analysis
    AUC_records = sum(AUC_recs, []) #unlist
    
    if do_bandpower:
        AUC_df = pd.DataFrame.from_records(AUC_records, columns=['Age group', 'Band', 'AUC', 'Sex', 'Subject'])
    else:
        AUC_df = pd.DataFrame.from_records(AUC_records, columns=['Age group', 'AUC', 'Sex', 'Subject'])
    
    
    
    return cohort_n_mean, AUC_df
    

def plot_glob_age_groups(cohort_n_mean, freqs, sleep):
    
    cm = plt.get_cmap('viridis')
    colors = [cm(x) for x in np.linspace(0, 1, len(cohort_n_mean))] 
    
    fig, ax = plt.subplots(figsize=(8,6)) 
    
    for ii, (name, mean_data, sd_data) in enumerate(cohort_n_mean):
        glob_mean = np.mean(mean_data, axis=0)
        ax.plot(freqs, glob_mean, color=colors[ii], label=name)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Power $log_{10}$(dB $\mu V^{2}$/ Hz )')
        ax.legend(loc="upper right", title='Age (in years)')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        
    return fig


def plot_individual_psd_pasta(subj, PSD_df, freqs, chs, segment='PSD N2a'):
    """
    

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    PSD_df : TYPE
        DESCRIPTION.

    Returns
    -------
    fig : a matplotlib figure

    """
    
    
    subj_data=PSD_df.loc[subj][segment]
    n_chs=len(subj_data)
    
    # Set colormap for visuals 
    cm = plt.get_cmap('viridis')
    colors = [cm(x) for x in np.linspace(0, 1, n_chs)] 
    
    fig, ax = plt.subplots(figsize=(8,6)) 
    
    for ii in range(n_chs):
        ax.plot(freqs, subj_data[ii,:], color=colors[ii], alpha=0.6, label=chs[ii])
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Power $log_{10}$(dB $\mu V^{2}$/ Hz )')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    
    plt.legend()
    
    return fig

def auc_over_freqrange(PSD_df, freqs, freqrange=(12,15), relative=True):
    
    
    # average sata from PSD segmetns within one sleep stage
    N1_average = np.array(PSD_df[['PSD N1a', 'PSD N1b']]).mean(axis=1)
    N2_average = np.array(PSD_df[['PSD N2a', 'PSD N2b', 'PSD N2c', 'PSD N2d']]).mean(axis=1)
    
    data_df = PSD_df.iloc[:,0:2]
    
    data_df['N1 data'] = N1_average
    data_df['N2 data'] = N2_average

    # Group the data according to age groups
    bin_labels = ['1','2','3','4','5','6','7','8','9','10']
    bins, bin_labs = pd.qcut(data_df['Age'], q=10, retbins=True, labels=bin_labels)
    bin_names = [str( round(bin_labs[i-1],1)) + ' – ' + str(round(bin_labs[i],1)) for i in range(1,len(bin_labs))]
    
    bins, bin_labs = pd.qcut(data_df['Age'], q=10, retbins=True, labels=bin_names)
    
    data_df['Age group'] = bins
    AUC_recs = []
    
    cohrt_groups = data_df.groupby('Age group', observed=False)
    
    
    # determine the freq. indices of interest
    f_idxs =  np.where((freqs >= freqrange[0]) & (freqs <= freqrange[1]))[0]
    
    
    for name, cohort_df in cohrt_groups:
        cohort_data = cohort_df[['N1 data', 'N2 data']]
        cohort_data = cohort_data.dropna() #get rid of nan -values
        # average over channels: nsubj x n_chs x n_freq -> nsubj x n_freq
        N1_Arr = np.array(np.stack(cohort_data['N1 data'].values)).mean(axis=1) 
        N2_Arr = np.array(np.stack(cohort_data['N2 data'].values)).mean(axis=1)
       
        #compute AUC values ("total power overall (from grand average)")
        AUCs_n1 = np.trapz(N1_Arr[:,f_idxs], x=freqs[f_idxs], axis=-1)
        AUCs_n2 = np.trapz(N2_Arr[:,f_idxs], x=freqs[f_idxs], axis=-1) 
        
        # devide these by total power to get relative aucs 
        if relative:
            tot_n1 = np.trapz(N1_Arr, x=freqs, axis=-1)
            tot_n2 = np.trapz(N2_Arr, x=freqs, axis=-1)
            
            AUCs_n1 = AUCs_n1/tot_n1
            AUCs_n2 = AUCs_n2/tot_n2

            
        namevec = [name]*AUCs_n1.shape[0]
        AUC_recs.append(  list( zip(namevec, AUCs_n1, AUCs_n2, cohort_df['Sex'].values, cohort_df.index.values, ) ))            

    
    # make AUC df for later analysis
    AUC_records = sum(AUC_recs, []) #unlist  
    AUC_df = pd.DataFrame.from_records(AUC_records, columns=['Age group', 'AUC N1', 'AUC N2', 'Sex', 'Subject'])
    
    AUC_df = AUC_df.dropna()
    
    # return the data in a long format 
    AUC_df_long = pd.melt(AUC_df, id_vars=['Subject', 'Sex', 'Age group'],
                          value_vars=['AUC N1', 'AUC N2'], var_name='Sleep', value_name='AUC')
    AUC_df_long['Sleep'] = [name.strip('AUC ') for name in AUC_df_long['Sleep'].values]
    
    return AUC_df_long, AUC_df
    

def within_sleep_stability(PSD_df, stage='N1'):

    raw = read_template_subject()
    adjacency, ch_names = find_ch_adjacency(raw.info, ch_type="eeg")
  
    subj_X = [] 
    for subject in PSD_df.index:
        subj_data = PSD_df.loc[subject]
        data = subj_data[3:] #choose only the psd 
        #metadata=subj_data['Sex']+ ', '+ str(subj_data['Age'])+'y'
        
        #means = [np.nanmean(datum, axis=0) for datum in data]
        means = [datum for datum in data]
        # conditions = data.index  
        
        # plot_df = pd.DataFrame(means).T
        # plot_df.columns=conditions
        # plot_df['Freq. (Hz)'] = freqs
        # plot_df = pd.melt(plot_df, ['Freq. (Hz)'], var_name='segment', value_name='Power')
        # plot_df['sleep'] = ['N1' if '1' in plot_df['segment'][i] else 'N2' for i in range(len(plot_df))]
     
        # # clustering based on mean segments? or f test?
        try:
            N1_data, N2_data = np.array(means[0:2]), np.array(means[2:])
            
            if stage == 'N1':
                subj_X.append(N1_data.transpose(0,2,1))
            else:    
                subj_X.append(N2_data.transpose(0,2,1))


        except ValueError:
          print('Not enough N2 segments, skipping!')
          
     
    subj_X = np.array(subj_X)
    X =[ subj_X[:, i, :] for i in range(subj_X.shape[1])]
    
    alpha_cluster_forming = 0.001
    pval = 0.05
 
    n_conditions = len(X)
    n_observations = len(X[0])
    dfn = n_conditions - 1
    dfd = n_observations - n_conditions
    
    thresh = scipy.stats.f.ppf(1 - alpha_cluster_forming, dfn=dfn, dfd=dfd)  # F distribution
        
        # X=data_array #np.array([*data_array]) #unpack data
    F_obs, clusters, cluster_pv, H0 = spatio_temporal_cluster_test(
        X,
        n_permutations=1000,
        threshold=thresh,
        tail=0,
        )
        
    #print(f'F-statistic: {F_obs}')
    if np.any(cluster_pv < pval):
        good_cluster_inds = np.where(cluster_pv < pval)[0]
    
        return clusters, good_cluster_inds



#%%% Plots
PSD_df = pd.read_pickle('PSD_dataframe_with_bandpower.pkl')
sleep= 'BPW N1b' #'PSD N1a'
cohort_n_mean, AUC_df = sd_mean_inter_age_groups(PSD_df, freqs, sleep=sleep, do_bandpower=True)

PSD_df = PSD_df.drop(['BPW N1a', 'BPW N1b', 'BPW N2a', 'BPW N2b', 'BPW N2c', 'BPW N2d'], axis=1)

#get absolute value from aucs -> makes comparisons more intuitive
#AUC_df['AUC'] = np.abs(AUC_df['AUC'])

fig = plot_glob_age_groups(cohort_n_mean, freqs, sleep)
plt.gcf().suptitle(f'Power spectral densities in {sleep}', fontsize=14)
plt.show()



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

subj=subjects[15] #365, 399, 179, 48, 449

fig, fig2, metadata, clusters, cluster_pv = plot_glob_intra_individual(PSD_df, subj, freqs, perform_clustering=True)
fig2.show()


#### also plot all the channels in one plot - for pipeline fig.

ind_fig = plot_individual_psd_pasta(subj, PSD_df, freqs, chs, segment='PSD N1b')

#%% Stats
import scipy.stats as stats
import pingouin as pg
from statsmodels.formula.api import mixedlm

#stability stats
within_sleep_stability(PSD_df, stage='N1')




# examine the difference of sleep spindle area over age groups
# using mixed linear model with subjects as random effects
AUC_df_long, AUC_df = auc_over_freqrange(PSD_df, freqs, freqrange=(12,15), relative=False)
AUC_df_long = AUC_df_long.rename(columns={'Age group': 'Group'})
AUC_df_long['AUC'] = np.log10(30+AUC_df_long['AUC'])
model = mixedlm('AUC ~ Sleep', data=AUC_df_long,
                groups="Subject", vc_formula={"Group": "0 + Group"}, re_formula="~Sleep")
result = model.fit()
print(result.summary())

# diagnostic plots
residuals = result.resid
stats.probplot(residuals, dist="norm", plot=plt)
plt.show()

plt.scatter(result.fittedvalues, residuals)
plt.xlabel("Fitted values")
plt.ylabel("Residuals")
plt.title("Residuals vs Fitted Values")
plt.show()

# plot the data as well


#Other option: difference between measurements (N2-N1), does that differ between groups?


#Another option: check the difference in groups separately 
aucs_grouped = AUC_df.groupby('Age group')

for name, group_df in aucs_grouped:
    stat, p_val = stats.wilcoxon(group_df['AUC N1'], group_df['AUC N2'])
    print(name)
    print(f"Wilcoxon statistic: {stat}, p-value: {p_val}")


###########################################
# AUC stats
aucs_grouped = AUC_df.groupby(['Age group']) #for relative bandpower
group_stats = aucs_grouped.describe()
group_stats = group_stats.iloc[[0,1,2,3,4,6,7,8,9,5],:] 

# plot auc distribution boxplots as per age group
cm = plt.get_cmap('viridis')
colors = [cm(x) for x in np.linspace(0, 1, 10)] 
vd_palette = sns.color_palette("viridis", n_colors=10, as_cmap=True)


band_index=6
#AUC_df2 = AUC_df[AUC_df['Band']==f'band {band_index}']
#band_freqs = f_bands[band_index]
AUC_df2 = AUC_df
AUC_df2['AUC'] = AUC_df2['AUC N2']

g = sns.FacetGrid(AUC_df2, row='Age group', hue='Age group', aspect=15, height=0.75, palette=colors)
# add the densities kdeplots for age groups
g.map(sns.kdeplot, 'AUC',
      bw_adjust=1, clip_on=False,
      fill=True, alpha=0.5, linewidth=1, legend=False)
# here we add a white line that represents the contour of each kdeplot
g.map(sns.kdeplot, 'AUC', 
      bw_adjust=1, clip_on=False, 
      color="w", lw=0.2)
# here we add a horizontal line for each plot
g.map(plt.axhline, y=0,
      lw=1, clip_on=False)
# we loop over the FacetGrid figure axes (g.axes.flat)
# ax.lines[-1].get_color() enables to access the last line's color in each matplotlib.Axes
label_dict = dict(zip([1,2,3,4,5,6,7,8,9,10],
                      group_stats.iloc[:,1].values) )
                      #list(np.unique(AUC_df['Age group']))) )

for i, ax in enumerate(g.axes.flat):
    kdeline = ax.lines[0]

    xs, ys, = kdeline.get_xdata(), kdeline.get_ydata()
    #middle = xs.mean()
    middle=label_dict[i+1]
    #left, middle, right = np.percentile(xs, [25, 50, 75])
    ax.vlines(middle, 0, np.interp(middle, xs, ys), color=ax.lines[-1].get_color(), ls='--')
    #ax.fill_between(xs, 0, ys, where=(left <= xs) & (xs <= right), interpolate=True, 
    #                facecolor=ax.lines[-1].get_color(), alpha=0.2)

    
    ax.text(220, 0.001, f'Mean: {middle:.2f}',
            fontweight='bold', fontsize=10,
            color=ax.lines[-1].get_color())
    
# eventually we remove axes titles, yticks and spines
g.set_titles("")
g.set(yticks=[])
g.set_ylabels()
g.despine(bottom=True, left=True)

# we use matplotlib.Figure.subplots_adjust() function to get the subplots to overlap
g.fig.subplots_adjust(hspace=-0.01)

plt.setp(ax.get_xticklabels(), fontsize=12)
plt.xlabel('Area under curve', fontsize=15)
g.fig.suptitle('AUC in Grand Average N2 PSDs per Age Group',
               #ha='right',
               fontsize=18,
               fontweight=20)
g.add_legend(markerscale=14, title_fontsize=18, fancybox=True)
plt.show()

# ok now stats for real for real

age_group1 = AUC_df.loc[AUC_df['Age group']=='0.0 – 0.5']['AUC']
age_group2 = AUC_df.loc[AUC_df['Age group']=='0.5 – 0.8']['AUC']
age_group3 = AUC_df.loc[AUC_df['Age group']=='0.8 – 1.3']['AUC']
age_group4 = AUC_df.loc[AUC_df['Age group']=='1.3 – 1.9']['AUC']
age_group5 = AUC_df.loc[AUC_df['Age group']=='1.9 – 2.8']['AUC']
age_group6 = AUC_df.loc[AUC_df['Age group']=='2.8 – 4.4']['AUC']
age_group7 = AUC_df.loc[AUC_df['Age group']=='4.4 – 6.2']['AUC']
age_group8 = AUC_df.loc[AUC_df['Age group']=='6.2 – 8.1']['AUC']
age_group9 = AUC_df.loc[AUC_df['Age group']=='8.1 – 11.3']['AUC']
age_group10 = AUC_df.loc[AUC_df['Age group']=='11.3 – 18.6']['AUC']


stats.f_oneway(age_group1,age_group2,age_group3,age_group4,age_group5,
               age_group6,age_group7,age_group8,age_group9,age_group10)
age_group_labs=['0.0 – 0.5','0.5 – 0.8','0.8 – 1.3','1.3 – 1.9','1.9 – 2.8',
                '2.8 – 4.4','4.4 – 6.2','6.2 – 8.1','8.1 – 11.3','11.3 – 18.6']
# test pairwise correlations?
age_pairs=[]
for ag1 in range(9):
    for ag2 in range(ag1+1, 10):
        age_pairs.append( (age_group_labs[ag1], age_group_labs[ag2]) )

corrected_p = 0.05/len(age_pairs)
print(f'Bonferroni-corrected p-value: {corrected_p}')
for group1, group2 in age_pairs:
    res=stats.ttest_ind(AUC_df.loc[AUC_df['Age group']==group1]['AUC'],
                          AUC_df.loc[AUC_df['Age group']==group2]['AUC'])
    if res.pvalue < corrected_p:
        print(f'{group1} vs {group2}')
        print(res)













