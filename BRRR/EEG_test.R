setwd("/projects/FABEEG/BRRR")

library("ggplot2")
library("cvms")
#library(plotly)
#library(htmlwidgets)


#every value of the variance parameter of Ω can
#be immediately interpreted as the percentage of variance explained by the noise model as compared to the covariates

#############################################
# By Susa, testing just to see what this is #
#############################################

# So apparently some kind of permutation test

Longitudinal_set_correlation<- function(brrr_model, oos_target, oos_structure, significance_level){
  # copied from BRRR utilities to do a tryout with correlation instead of manhattan distance
  #' function to compute out-of-sample validation for a given model
  #' where the model and the out of sample data come from the same subjects.
  #' This aims to use an inverted validation process,
  #' i.e. it does not predict based on the out of sample covariates
  #' as these would likely be the same for both sample and
  #' validation set
  #' this only considers the case of two observations being paired
  #' (i.e. at most 2 siblings from a family, for the HCP MEG data that is the case, more compicated structures
  #' will need additional complexity in the code, which should only be implemented, once an actual need arises)
  #' Parameters:
  #' brrr_model: the model to be validated as output by brrr
  #' oos_target: the out of sample data to use for validation
  #' oos_structure: the structure the model is based on, i.e family structure
  #'   this needs to be of the one confound level per column - form, with pseudo-boolean
  #'   entries (1/0)
  
  res<-list()
  # average gamma because Bayesian! (all the gamma in the trace are samples of the "real" gamma)
  res$scaling <- ginv(averageGamma(brrr_model))
  res$oosproj <- oos_target %*% res$scaling
  res$isproj <- brrr_model$data$target %*% res$scaling
  
  corr_to_inSample <- cor(t(res$oosproj),t(res$isproj))
  max_corr <- apply(corr_to_inSample, 2, max)
  test_class <- diag(corr_to_inSample)
  classification<- test_class == max_corr
  #print(classification)
  res$class_rate = sum(classification)/length(classification)
  projection_size <- dim(corr_to_inSample)[1]
  perm_rates <- vector(mode="numeric",length = significance_level)
  for (i in 1:length(perm_rates)){
    permutation <- sample(projection_size,size=projection_size)
    perm_class <- vector(mode = "numeric",length = projection_size)
    k <-1
    for (j in permutation){
      perm_class[k] <- corr_to_inSample[j, k]
      k <- k+1
    }
    perm_rates[i] <- sum(perm_class == max_corr)/length(perm_class)
  }
  res$permutedClassificationRates<-perm_rates
  
  
  return(res)
}



create_plot<-function(validation_result){
  permuted_rates = validation_result$permutedClassificationRates
  plot_data=as.data.frame(permuted_rates)
  pp<-ggplot(data=plot_data ,aes(x=permuted_rates))+
    geom_histogram(binwidth = 0.005,color="black", fill="white")+
    geom_vline(xintercept = validation_result$class_rate,color="red")
  show(pp)
  #return(ggplotly(pp))
}


brrr_model = res
brrr_model$data$target=brrr_model$data$phenotypes

validation_result = Longitudinal_set_correlation(brrr_model=brrr_model, oos_target=Y2, 
                                                oos_structure=X[1:591,],10000)

create_plot(validation_result)

