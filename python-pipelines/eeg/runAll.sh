#!/bin/bash


while read LINE; do
    #echo ${LINE}
    python3 ica_02.py $LINE False 
done < EEG_subjects.txt


