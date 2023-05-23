#!/bin/bash


while read LINE; do
    #echo ${LINE}
    python3 freqfilt_01.py $LINE 
done < EEG_subjects.txt


