##################################
# Principal component regression #
#       by Verna Heikkinen       #
##################################
library(corrplot)
library(pls)

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

pls_regress <-function(X, y){
  scaled_resposnse <- scale(y,center=T,scale=T) #the response vector has to be normalized as well
  
  #split data to test & training sets
  df = as.data.frame(cbind(scaled_response, X))
  train <- df[1:150,]
  test <- df[151:nrow(df),2:ncol(df)]
  y_test <-df[151:nrow(df),1]
  
  #fit the model
  model <- pcr(V1 ~., data=train, scale=TRUE, validation="CV")
  pred <- predict(model, test, ncomp=20)
  
  #RMSE
  sqrt(mean((pred - y_test)^2))
  #TODO: move these too
  #validationplot(model, val.type = "MSEP")
}


