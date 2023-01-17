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


#Check some correlations
cor.test(lat_map$V1, as.numeric(lat_map$sex=='F')) #nyce

cor.test(lat_map$V1, as.numeric(lat_map$cap=='FT'))
cor.test(lat_map$age, as.numeric(lat_map$cap=='FT')) #kinda the point.. eh?

#Regression
fit <- lm(age~V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12,data=lat_map)
fit2 <- lm(age~V1, data=lat_map)
summary(fit)
#all are significant... kinda makes sense.

coefficients(fit) # model coefficients
confint(fit, level=0.95) # CIs for model parameters
fitted(fit) # predicted values
residuals(fit) # residuals
anova(fit) # anova table
vcov(fit) # covariance matrix for model parameters
influence(fit) # regression diagnostics
plot(fitted(fit), resid(fit))

#Plots... or nah?
ggplot(data=lat_map,aes(V1, age)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_minimal() + #ylim(7,20) +
  labs(x='V1', y='Age (years)', title='Linear Regression Plot') 


#TODO: glmer with gender as a random effect?

#TODO: RECONSTRUCTION ERROR and its correlates? inspired by mr. Kimmo Alakulju

dist_func <- function(dis='L2', x, y){
  if(dis=='L2'){
    return( sqrt(sum(x-y)**2 ))
  } else if(dis=='cos'){
    return( sum(x*y) / (sqrt(sum(x**2))*sqrt(sum(y**2))) )
  } else {
    print("dumb")
  }
}

Yhat <- lat_space %*% t(comp)
I=nrow(Yhat)
M=nsubj

rec_error <- vector(mode='list', length=I)
names(rec_error) <- dimnames(lat_space)[[1]]#[1:M]

for(i in 1:I){ #M
  #v1 <- dist_func(dis="L2",Y[i,], Yhat[i,])
  #v2 <- dist_func(dis="L2",Y[(i+M),], Yhat[(i+M),])
  #rec_error[[i]] <- mean(c(v1,v2))
  rec_error[[i]] <- dist_func(dis="L2",Y[i,], Yhat[i,])
} 

cor.test(unlist(rec_error), lat_map$age)
cor.test(unlist(rec_error), as.numeric(lat_map$sex=='F'))
