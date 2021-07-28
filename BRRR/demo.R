setwd("/projects/FABEEG/BRRR")

## Create artificial data
N <- 500  #500 samples
classes <- sample(c(-1,0,1), N, replace=TRUE) #Three classes
W <- c(.2, .6, -.3) #Mapping from classes to observations
x <- outer(classes, W) #The observations; outer product 
x <- x[,] + rnorm(prod(dim(x))) #Noise


# Learn BRRR model
source("brrr.R")
res <- brrr(matrix(classes, ncol=1), x, K=1,n.iter=100,init="penLDA")

# The true and inferred mapping (note: scaling is different due to Psi in BRRR)
mapping <- data.frame(actual=W, final_posterior_sample=c(res$model$brr$context$Gamma))
print(mapping)

#Posterior samples for Gamma. n.iter should be long enough so that these are converged
# -- initialization will have effect as well
GammaPosterior <- simplify2array(res$traces$Gamma)


# !!!!!!!!
#every value of the variance parameter of Ω can
#be immediately interpreted as the percentage of variance explained by the noise model as compared to the covariates


load("/projects/FABEEG/BRRR/data/N2spectrum.RData")
ages = read.csv('data/age_df.csv')
ages <- ages[,-1]

ages = subset(ages, ages$File %in% names(Y))

# Searching the names/identifiers of the subjects
subjects <- unlist(individuals)
S <- length(subjects)
A <- vector("list",length=S)
names(A) <- subjects
chs <- 1:19 #channels


#Omits frequencies; not in this case
omitFreq = c(0,45)
frequencies <- which(freq[,2] >= omitFreq[1] & freq[,2] <= omitFreq[2])
freq <- freq[frequencies, ]


# Saving the data into matrix A
for(s in names(Y)) {
  
  if(s %in% ages$File){ #temporal solution to get over the filename hassle
    tmp <- t(Y[[s]])
    tmp <- log10(tmp)
    if(any(is.na(tmp))) browser()
    A[[s]] <- tmp[frequencies,chs]
    
    corrupted=NULL
    
  } else {
    corrupted = s
  } 
}

if(!is.null(corrupted)){
  obs <- subjects[-which(subjects==corrupted)]
  A = A[names(A) %in% corrupted == FALSE]  
} else { 
  obs <- subjects
}


print("Data dimension per subject:")
print(dim(A[[1]]))
keepFeat <- NA

# Making MEG matrix X and response vec y and filling them
X <- matrix(NA,length(A),prod(dim(A[[1]])),dimnames=list(names(A),c()))
y = matrix(NA, length(A), dimnames=list(names(A)))
colnames(X) <- c(outer(paste0("s",1:nrow(A[[1]]),"."),1:ncol(A[[1]]),paste0))

for(i in obs) { 
  tmp <- A[[i]]
  if(nrow(tmp)==ncol(tmp)) tmp[lower.tri(tmp)] <- 0
  X[i,] <- c(tmp)
  y[i,] <- ages[ages$File==i,]$Age #get the age of the obs
}

keepFeat <- which(apply(X,2,var,na.rm=T)>0)
X <- X[,keepFeat]

print("LDA matrix dimension:")
print(dim(X))

# Fit BRRR model
res <- brrr(X, t(y), K=6, n.iter=500, init="penLDA")

final_posterior_gamma = res$model$brr$context$Gamma
final_posterior_psi = res$model$brr$context$Psi

#actual fit; weird?
psi = t(final_posterior_gamma)
gamma = t(final_posterior_psi)
omega = res$model$brr$context$Omega

Y_pred = (X%*%psi + omega)%*%gamma
print(Y_pred)
