############################################
# This script is for plotting the spectra  #
############################################

library("ggplot2")
library("viridis")

setwd("/projects/FABEEG/BRRR")
load(paste0("data/",ex,"spectrum.RData"))
ages = read.csv('data/age_df.csv')
ages <- ages[,-1]


#ja sit taas vähän pseudukkaa

# valihte ne tyypit joilla ages[which(ages$Age>10),]
# .. ja valitse ne Y:stä
# laske jäljellejääneille tyypeille Global Average spectra [14x19] -> 14x1
# ..ja tallenna se dataframeen
# plottaa dataframe lineplottina
# 
# lisäksi
# laske keskiarvospektri
# .. ja yksilöistä johtuva SD
# ja plottaa tämä myös