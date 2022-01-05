import subprocess
import argparse
import mne
from config import fname
import os


# Be verbose
mne.set_log_level('INFO')

# Handle command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', metavar='sub###', help='The subject to process')
parser.add_argument('-ses',  help='The session to process')
parser.add_argument('-task', help='The task to process')
parser.add_argument('-run', help='The run to process')
args = parser.parse_args()
subject = args.subject
ses = args.ses
task = args.task
run = args.run
print('Run tSSS on subject:', subject)
print('Running task:', task)
                    
input_file = fname.raw(subject=subject, task=task, run=run, ses=ses, proc='raw')
output_file = fname.raw(subject=subject, task=task, run=run, ses=ses, proc='tsssmctrans')       
trans_file =  fname.raw(subject=subject, task=task, run='01', ses=ses, proc='raw')
pos_file = fname.pos(subject=subject, task=task, run=run, ses=ses)
    
# extract bad channels of the current subject
#current_bad = bad_channels.query(f'subject=={subject} and session=={ses} and run=={run}').channels.iloc[0]

# args to be used in maxfilter
#args = ['maxfilter', '-f', input_file, '-o', output_file, '-trans', \
#            trans_file, '-movecomp', 'inter', '-st', '60', \
#            '-v','-autobad','off','-hpicons','-origin','fit',\
#            '-frame','head', '-bad'] +  current_bad.split()

args = ['maxfilter', '-f', input_file, '-o', output_file, '-trans', \
                trans_file, '-movecomp', 'inter', '-st', '10', \
                '-v','-autobad','on','-hpicons','-origin','fit',\
                '-in', '8', '-out', '3', '-frame','head', \
                 '-hp', pos_file ]   

# save the log and errors
log_output = open(fname.tsss_log(subject=subject, run=run, ses=ses, task=task), 
                          "w") 
        
# run maxfilter
subprocess.run(args=args, stdout=log_output,stderr=log_output)
