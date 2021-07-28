#!/bin/bash


while read LINE; do
    #echo ${LINE}
    python3 inspect_psds.py $LINE 
done < eeg_subjects.txt


