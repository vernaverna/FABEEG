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
ages = read.csv('data/new_age_df.csv')
ages <- ages[,-1]


## PLOT AGE HISTOGRAM ##

#TODO: think about the colour 
p <- ggplot(data=ages, aes(x=Age)) +
  geom_histogram(binwidth = 1, fill='yellowgreen', color="#e9ecef", alpha=0.9) + 
  ylab('Count') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) #no grid
p
ggsave('figures/agehist.pdf', plot=p)

subj = ages[which(ages$Age<5),]$File #choose the age group
Y=Y[names(Y) %in% subj == TRUE] #TODO: there are some in ages that are not in Y
subjects=subj

#first try with one subject only 
ind_data <- log10(Y[[34]])
rownames(ind_data) <- read.csv("var/ch_names.csv", header = F)$V1
avg <- apply(ind_data, 2, mean) #average over all channels

freq_bins <- apply(freq, MARGIN=1, FUN=function(x){paste(x[1], x[2], sep="-")}) #create freq.bins

df <- data.frame(cbind(cbind(avg, freq_bins), t(ind_data)))
df2 <- melt(df, id.vars=c('freq_bins', 'avg')) #change to long format (ggplot likes dis)
df2$avg <- as.numeric(df2$avg)
df2$value <- as.numeric(df2$value)
df2$freq_bins <- factor(df2$freq_bins, levels=freq_bins) #reorder freq bins as factor

#line plot
#library("ggridges")
df2 %>%
  ggplot(aes(x=freq_bins, y=value, color=variable, group=variable)) +
  geom_line(alpha=0.6) + 
  geom_line(aes(x=freq_bins, y=avg), color="black") + theme_minimal()


# Trying small multiples

df2 %>% 
  ggplot(aes(x=freq_bins, y=value, group=variable, colour=variable)) +
  geom_line() + ggtitle("Log-bandpowers of 1 subject") + 
  facet_wrap(~variable, scale="free_y") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))



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
  
  
  
  
  
  
  
  
  