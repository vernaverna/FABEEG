##########################################
# Linear regression on the BRRR solution #
#                                        #
# by Verna Heikkinen 2022                #
##########################################

setwd("/projects/FABEEG/BRRR")
library("ggplot2")
library("ggpubr")
library("stats")
library("rstatix")
library("lme4")
library("car")
library("dplyr")
library("factoextra")

load("results/full/all_2N2_BRRR_K12.RData") 
comp <- res$scaling
Y <- res$data$phenotypes
X <- res$data$genotypes
subj <- row.names(X)
lat_space <- Y%*%comp

#get the ages out
ages = read.csv('data/new_age_df.csv')
ages <- ages[,-1]
ages = ages[ages$File %in% unique(subj)==TRUE, ] #to get rid of some extra subjects that should not be there


#build a dataframe 
nsubj <- length(unique(subj))
# adding together N2 mappings  plus some covariates
lat_map = as.data.frame(rbind(lat_space))
lat_map['spectra'] = c(rep('N2A', nsubj), rep('N2B', nsubj))
lat_map['age'] = rep(ages[which(ages$File%in%subj),]$Age, 2) #get age data
lat_map['group'] = rep(round(ages[which(ages$File%in%subj),]$Age, 0), 2) #get age data
lat_map['sex'] = rep(ages[which(ages$File%in%subj),]$Sex, 2)
lat_map['cap'] = rep(ages[which(ages$File%in%subj),]$Cap, 2)
lat_map['subject'] = subj[1:nsubj]


#Regression
fit <- lm(age~V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12,data=lat_map)
summary(fit)
#all are significant... kinda makes sense.

coefficients(fit) # model coefficients
confint(fit, level=0.95) # CIs for model parameters
fitted(fit) # predicted values
residuals(fit) # residuals
anova(fit) # anova table
vcov(fit) # covariance matrix for model parameters
influence(fit) # regression diagnostics


#Plots... or nah?
ggplot(data=lat_map,aes(V1, age)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_minimal() + #ylim(7,20) +
  labs(x='V1', y='Age (years)', title='Linear Regression Plot') 


#TODO: glmer with gender as a random effect?


