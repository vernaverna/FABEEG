#!/bin/bash


while read LINE; do
    #echo ${LINE}
    python3 psds_03.py $LINE 
done < EEG_subjects.txt


