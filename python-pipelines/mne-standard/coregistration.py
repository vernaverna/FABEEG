# -*- coding: utf-8 -*-
"""
Coordinate frame alignment.
Coregistration works properly only with anaconda2!
Works with new mne (0.20.0) and new mayavi (4.7.1), still freezes sometimes and needs restarting
Checking coregistration works only with anaconda3!

In terminal, before starting e.g. spyder
ml freesurfer
export FREESURFER_HOME=/work/modules/Ubuntu/14.04/amd64/t314/freesurfer/6.0.0

"""
import mne


subject = 'sub-s01y1'
data_path = '/m/nbe/project/aging/recordings/'
subjects_dir = data_path + 'BIDS/derivatives/freesurfer/'

## Coregistration 

mne.gui.coregistration(subject=subject, subjects_dir=subjects_dir)

# 1. Load digitalization source (head shape source) from any trans_tsss file
# 2. Edit Fiducaials (left side panel)
# 3. Save fiducials (default)
# 4. Fit fiducials
# 5. Fit ICP (head shape)
# 6. Omit outliers and fit ICP again
# 7. When happy with co-registation, save (sub-p02_ses-01_coreg-transfile.fif)
# 8. Run separately for each HPI digitalization (ses-01, ses-02)


## Check coregistration

#test-file:
coreg='%smeg/aging_s03_y3/200710/sub-s03y3-ses-01-raw-trans.fif' % (data_path)
test_data='%smeg/aging_s03_y3/200710/sub_s03y3_ses_01_task_envspeech_run_01.fif' % (data_path)

#another test-file
#coreg_ses = 'ses-01'
#other_ses = 'ses-02'
#other_digitalization_source = 'task-N400_run-01_proc-tsssmctrans_meg'
#
#coreg = '%sBIDS/%s/%s/meg/%s_coreg-transfile.fif' % (data_path, subject, coreg_ses, subject)
#test_data = '%sBIDS/%s/%s/meg/%s_%s_%s.fif' % (data_path, subject, other_ses, subject, other_ses, other_digitalization_source)

info = mne.io.read_info(test_data)
## surfaces: 'head-dense' = high resolution, 'head' = low reso, 'brain' = pial 
mne.viz.plot_alignment(info, coreg, subject=subject, dig=True, meg=['helmet', 'sensors'], subjects_dir=subjects_dir, surfaces='head-dense')
