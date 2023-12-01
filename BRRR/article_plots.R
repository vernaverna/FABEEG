# collection of plots for teh article.

setwd("/projects/FABEEG/BRRR")
library("penalizedLDA")
library("ggplot2")
library("dplyr")
library("cvms")
library("RColorBrewer")
library("corrplot")


load("results/full/NEW_alldata_2N1_BRRR_K12.RData")
Y <- res$data$phenotypes
X <- res$data$genotypes
subj <- row.names(X)
ages = read.csv('data/new_age_df.csv')
ages <- ages[,-1]
ages[ages==" "] <- NA #replace empty strings with NA-values
ages = ages[ages$File %in% unique(subj)==TRUE, ] #to get rid of some extra subjects that should not be there
ages['Sex'] <- as.factor(ages$Sex)


# FIGURE 6
# make a stacked histogram with subject info
age_df <- ages %>% mutate(age_group = cut(Age, breaks=14))
age_df <- na.omit(age_df)
prop_df <- age_df %>% group_by(age_group, Sex) %>% summarise(n = n()) %>% mutate(prop=n/sum(n) )
prop_df['prop'] <- round(prop_df$prop, 2)
prop_df$prop[prop_df$Sex=='M'] <- NA
prop_df <- prop_df %>% mutate(total=sum(n) )
age_df <- merge(x=age_df,y=prop_df, by=c('age_group', 'Sex'), all=TRUE)

pg <- ggplot(age_df, aes(x=age_group, fill=Sex)) + 
  geom_bar(position='stack') + theme_bw() + ylab("Frequency") + xlab('Age') +
  geom_text(aes(y=total, label=prop), vjust=-0.4, color='grey', size=5) + 
  scale_fill_hue(l=70, c=90) +
  scale_x_discrete(labels=c("0","","","4","","","8","","","12","","","16","")) +
  theme(legend.key = element_rect(linewidth = 16), 
        legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold", size = 18),
        axis.text.x = element_text(size=13), axis.text.y = element_text(size=13),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

 
# NONNIH.



