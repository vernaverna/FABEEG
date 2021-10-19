############################################
# This script is for plotting the spectra  #
############################################

library("ggplot2")
library("viridis")
library("reshape2")
library("dplyr")

ex="N2A"

setwd("/projects/FABEEG/BRRR")
load(paste0("data/",ex,"spectrum.RData"))
ages = read.csv('data/age_df.csv')
ages <- ages[,-1]

subj = ages[which(ages$Age>10),]$File #choose the age group
Y=Y[names(Y) %in% subj == TRUE] #TODO: there are some in ages that are not in Y
subjects=subj

#first try with one subject only 
ind_data <- log10(Y[[1]])
rownames(ind_data) <- read.csv("var/ch_names.csv", header = F)$V1
avg <- apply(ind_data, 2, mean) #average over all channels

freq_bins <- apply(freq, MARGIN=1, FUN=function(x){paste(x[1], x[2], sep="-")}) #create freq.bins

df <- data.frame(cbind(cbind(avg, freq_bins), t(ind_data)))
df2 <- melt(df, id.vars=c('freq_bins', 'avg')) #change to long format (ggplot likes dis)

df2 %>%
  filter(variable == "avg") %>%
  ggplot(aes(freq_bins, value)) +
  labs(y = "log power", x="freq") +
  geom_line(data = df2, aes(group = variable), color = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data=df2, aes(y = avg), color = "blue")

# Trying bins 
library("ggridges")
df2 %>%
  ggplot(aes(x=freq_bins, y=value, fill=variable)) +
  geom_density_ridges(alpha=0.6, stat="identity") +
  theme_ridges() 

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
  
  
  
  
  
  
  
  
  