# operate the process here...

# Determine the following parameters:
K <- 6 # number of components (reduced rank) to use (max number of classes - 1)
model <- 'brrr' # model type ('penlda'=penalized LDA, 'brrr'=BRRR)
CVfolds <- 10 # how many cross-validation folds to use
omitMag <- TRUE  # logical indicating whether to omit magnetometers or not
omitFreq <- c(0, 90) # the lowest and highest frequency value to take into account
plotname <- "siblings_allfreq_nomag"
network <- TRUE # logical value for producing network figures. NOTE: this only works with R version 3.5 or higher 
n.iter <- 500


load("data/ecspectrum.RData")
source("var.R")

## Cross-validation for class prediction
D <- acc(CVfolds=CVfolds,model=model,discTop=K, omitMag = omitMag, omitFreq = omitFreq, network=network, n.iter=n.iter, ex='ec')

if(!network){
  M <- length(classes)
  groups <- names(classes)
  conf <- matrix(NA,M,M,dimnames=list(paste("True",groups),paste("Pred",groups)))
  labs <- paste0("Distance to ",colnames(D))
  R <- max(abs(D))
  
  # The following code is used for plotting the results
  if(M == 2) {
    model.name <- if(model == 'penlda') '_penLDA' else '_brrr'
    pdf(file=paste0("figs/",plotname,model.name,".pdf"),width=7,height=7)
    plot(NULL, xlab=labs[1], ylab=labs[2], xlim=c(0,R), ylim=c(0,R))
    abline(a=0,b=1,col="#444444")
  }
  
  # A loop creating the classification matrix:  
  for(m in 1:M) {
    #new sub-matrix for each data group
    C <- D[classes[[m]], 1:M, K]
    c <- dim(C)[1]
    if(M == 2) {
      points(D[classes[[m]],1,K],D[classes[[m]],2,K],pch=m)
    }
    for (n in 1:M) { # for each column in C
      s <- c()
      for (d in 1:c) {
        # list of initials that are the smallest in each row
        a <- sum(C[d,n] == min(C[d, ]))  
        s = append(s, a)
        # the results are stored in the conf -matrix
        conf[m, n] <- sum(s)
      }
      
    }
    
  }
  
  # The score shows how many subjects were classified correctly 
  score <- sum(diag(conf)) / sum(conf)
  print("the precentage of classification accuracy: ")
  print(score)
  
  if(M == 2) {
    for(d in 1:2) {
      text(R/2,R/2,paste("Predicted as",groups[d]),pos=c(3,1)[d],srt=45)
      legend(c("topleft","bottomright")[d],paste0(groups,"(",conf[,d],")"),pch=1:M)
    }
    dev.off()
  }
  print("the sensitivity and specificity of the classification: ")
  sensitivity <- conf[1,1]/(conf[1,1]+conf[1,2])
  specificity <- conf[2,2]/(conf[2,1]+conf[2,2]) 
  ss <- c(sensitivity, specificity)
  ## might be vice verca, but does not necessarily matter?
  print(ss)
  print(conf)
}


