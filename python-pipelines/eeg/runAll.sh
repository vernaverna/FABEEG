#!/bin/bash

while read LINE; do
    #echo ${LINE}
    python3 freqfilt_01.py $LINE
    python3 ica_02.py $LINE True
    python3 psds_03.py $LINE
    
done < EEG_subjects.txt


