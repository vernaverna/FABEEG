#!/bin/bash


while read LINE; do
    #echo ${LINE}
    #python3 psds_03.py $LINE 
    python3 PSD_04b.py $LINE 
done < EEG_subjects.txt


