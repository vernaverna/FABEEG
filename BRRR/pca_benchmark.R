##################################
# Principal component regression #
#       by Verna Heikkinen       #
##################################
library(corrplot)

pca_regress <- function(X=Y, response){
  pca_data <- prcomp(X) #data matrix is already centered by data preparation function
  scaled_response = scale(response,center=T,scale=T) #the response vector has to be normalized as well
  
  comp = pca_data$x #take the principal components
  plot(summary(pca_data)$importance[3,]) #TODO: save into figs folder?
  
  #Do regression
  regr_data = as.data.frame(cbind(scaled_response, comp))
  colnames(regr_data)[1] = "response"
  model <- lm(response ~., data = regr_data) #regress against all PC:s
  
  NAs = which(is.na(model$coefficients)) #take NA indexes
  model$coefficients[NAs] <- 1
  #TODO: drop those indexes from the analysis
  
  coeff.estimates <- as.matrix(model$coefficients[2:length(model$coefficients)]) #ignore intercept
  eigen.matrix <- as.matrix(pca_data$rotation) #extract eigenvalues
  # Model coefficients = eigenvalues %*% PCR coefficients
  coeff <- eigen.matrix %*% coeff.estimates
  
  pred <- X%*%coeff
}

