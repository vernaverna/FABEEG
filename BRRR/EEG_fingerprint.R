# operate the process here...
setwd("/projects/FABEEG/BRRR")
library("penalizedLDA")
library("ggplot2")
library("cvms")


## TODO: use less individuals! try with old and young kids as well

#' Function for preparing the EEG spectra for analysis
#' 
#' @param spectra list, tells which spectra are used for the model: e.g. c('N1C', 'N2B')
#' @param validation_set which spectra to use in validation, e.g. 'N2C'
#' 
# ex = N1 or N2, depending which data is read in
prepare_data <- function(spectra, validation_set){
  
  
  # Choosing the channels and frequencies
  chs <- 1:19
  omitFreq = c(0,60)
  
  #extract the age groups
  ages = read.csv('data/new_age_df.csv')
  ages <- ages[,-1]
  
  use_all=T #should we use all subjects in training?
  age_gap=c(1,19) #exclude some of the younger children?
  #Cap='-'
  
  data_Y = vector(mode='list',length=length(spectra)) #containers for targets Y and covariates X
  data_X = vector(mode='list',length=length(spectra))
  
  names(data_Y) <- spectra
    
  #Go through each data file
  for(i in 1:length(spectra)){
    # Reading data into workspace
    spec = spectra[i]
    loadfile <- paste0("data/",spec,"spectrum.RData")
    
    # Check if the file exists
    if(file.exists(loadfile)) {
      load(loadfile)
      
      if(use_all==F){
        set.seed(11)
        #load("/projects/FABEEG/BRRR/ref_subjects.RData")
        #individuals=obs
        individuals = sample(individuals, 180)
        Y = Y[names(Y) %in% individuals]
        
      }
      
      # Omits frequencies
      frequencies <- which(freq[,2] >= omitFreq[1] & freq[,2] <= omitFreq[2])
      freq <- freq[frequencies, ]
      
      
      
      # Searching the names/identifiers of the subjects
      subjects <- unlist(individuals)
      S <- length(subjects)
      A <- vector("list",length=S)
      names(A) <- subjects
      
      corrupted = c() #container for files that are dropped
      
      # Saving the data into matrix A
      for(s in names(Y)) {
        
        if(s %in% ages$File){ #temporal solution to get over the filename hassle
          age_data = ages[ages$File==s,]
          if(age_data$Age >= age_gap[1] && age_data$Age <= age_gap[2]){
            #if(age_data$Cap == Cap){
            tmp <- t(Y[[s]]) #transposed
            tmp <- log10(tmp) 
            if(any(is.na(tmp))){ #browser()
              corrupted = c(corrupted, s)
            } 
            A[[s]] <- tmp[frequencies,chs] 
          } else {
            corrupted = c(corrupted, s)
          }
        } else {
          corrupted = c(corrupted, s)
        } 
      }
      
      if(length(corrupted)>0){
        obs <- subjects[subjects %in% corrupted == FALSE]
        A = A[names(A) %in% corrupted == FALSE]
      } else { 
        obs <- subjects
      }
      
    } else {
      print(paste0("File '",loadfile,"' does not exist, returning!"))
      return(0)
    }
    
    
    # The number of classes and storing the subjects
    M <- length(obs) #nclasses
    S <- length(obs) #nsubj
    
    print("Data dimension per subject:")
    print(dim(A[[1]]))
    keepFeat <- NA
    
    # Making MEG matrix X and response vec y and filling them
    Y <- matrix(NA,length(A),prod(dim(A[[1]])),dimnames=list(names(A),c()))
    x = matrix(NA, length(A), dimnames=list(names(A)))
    colnames(Y) <- c(outer(paste0("s",1:nrow(A[[1]]),"."),1:ncol(A[[1]]),paste0))
    
    
    for(j in obs) { 
      tmp <- A[[j]]
      if(nrow(tmp)==ncol(tmp)) tmp[lower.tri(tmp)] <- 0
      Y[j,] <- c(tmp)
    }
    
    keepFeat <- which(apply(Y,2,var,na.rm=T)>0)
    Y <- Y[,keepFeat]
    
    print("LDA matrix dimension:")
    print(dim(Y))
    
    data_Y[[i]] <- Y
    
  }

  # Combine data into 1 big matrix

  val=which(names(data_Y)==validation_set)
  validation_data = data_Y[[val]] #extract validation data
  data_Y[[val]] <- NULL
  
  data_Y <- do.call(rbind, data_Y) #make into matrices
  # TODO: make sure that data_Y has the same elements!!!!
  
  #choose those subjects who are both in the validation set and train set
  common_subjects = intersect(rownames(data_Y), rownames(validation_data)) 
  to_keep = which(rownames(data_Y)%in%common_subjects)
  data_Y <- data_Y[to_keep,]
  
  #Center and scale Y & validation data
  tmp <- scale(data_Y[!is.na(data_Y[,1]),],center=T,scale=T)
  Y1 <- scale(data_Y, attr(tmp,"scaled:center"), attr(tmp,"scaled:scale"))
  
  tmp1 <- scale(validation_data[!is.na(validation_data[,1]),],center=T,scale=T)
  Y_val <- scale(validation_data, attr(tmp1,"scaled:center"), attr(tmp1,"scaled:scale"))
  
  
  #Construct identifier matrix (covariates)
  subj = dimnames(Y1)[[1]]
  M = length(unique(subj)) #number of groups (=subjects)
  S = length(subj) #number of observations
  
  x <- rep(NA,M); names(x) <- unique(subj)
  for(s in 1:M){
    x[subj[s]] <- s
  }
  x=c(x,x)
  
  X <- matrix(0,S,M,dimnames=list(subj,paste0("subj_",unique(subj)))) 
  for(i in 1:length(x)) if(!is.na(x[i])) X[i,x[i]] <- 1
  
  #Clean up the subject info dataframe
  ages = ages[ages$File %in% unique(subj)==TRUE, ] #to get rid of some extra subjects that should not be there
  
  
  #Trying with added Z
  # Z1 = ages[ages$File %in% unique(subj),]
  # reorder_idx <- match(rownames(Y1), Z1$File)
  # Z1 <- Z1[reorder_idx,]
  # 
  # Z = matrix(c(Z1$Age, Z1$Age)) 
  

  return(list(Y1, x, X, ages, Y_val, Z=NA))
}


n2_data <- prepare_data(spectra = c("N2A","N2B","N2D"), validation_set = "N2D")
Y = n2_data[[1]]
X = n2_data[[3]]
x = n2_data[[2]]
ages = n2_data[[4]]

Y2 = n2_data[[5]]

### TRAINING ###

# The model is trained using two sets of N2 data, and the within-sample performance is evaluated using
# MSE, PTVE and accuracy (L1 distances in the projection)

#TODO: out-of sample prediction!

source("brrr.R")
#pred <- X*NA
res <- brrr(X=X,Y=Y,K=20,Z=NA,n.iter=1000,thin=5,init="LDA",fam=x) #fit the model
res$scaling <- ginv(averagePsi(res)%*%averageGamma(res)) #others have used projection (Y*(psi*gamma)-1)
W <- res$scaling

#save(res, file = "results/full/over5_indN2_BRRR_K15.RData")
lat_map <- Y%*%ginv(averageGamma(res))
lat_map_n2 <- Y2%*%W #mapping to latent space with N2_D data!# 

D <- matrix(NA, ncol(X), ncol(X), dimnames=list(unique(names(x)), paste0("mean",unique(names(x))) ) ) #distance matrix# # # 

for(testidx in 1:nrow(lat_map)){ #calculates the distances between individual and group mean in lat.space   
  for(m in 1:ncol(X)){     
    group_members <- rownames(lat_map)[m]   
    idxs = which(row.names(lat_map) %in% group_members) #     
    group_mean <- colMeans(lat_map[idxs,]) #mean over all individuals, vector of length K
    D[testidx,m] <- sum(abs(lat_map[testidx,]-group_mean)) #L1 distance #TODO: CHECK this idiot
    #D[testidx,m] <- sqrt(sum( (lat_map_n2[testidx,]-group_mean)**2 ))#L2 distance  
  } 
}





#### VALIDATION ####

PROJ <- D*0 #initialize 'projection' matrix

for(r in 1:nrow(D)){ #assign age groups based on lat. space distances
  index <- which.min(D[r,])
  PROJ[r,index] <- 1+PROJ[r,index]
}
colnames(PROJ) <- c(paste0("class",rownames(D)))

accuracy = sum(diag(PROJ))/nrow(D)
print(paste("Model accuracy:", accuracy)) #alright this is how I like it ;)


nsubj <- length(unique(subj))
# adding together N2 mappings  plus some covariates
lat_map = as.data.frame(rbind(lat_map_n2[1:nsubj,], lat_map))
lat_map['condition'] = c(rep('test', nsubj), rep('train', 2*nsubj))
lat_map['age'] = rep(ages[which(ages$File%in%subj),]$Age, 3) #get age data
lat_map['group'] = rep(round(ages[which(ages$File%in%subj),]$Age, 0), 3) #get age data
lat_map['sex'] = rep(ages[which(ages$File%in%subj),]$Sex, 3)

subj_number <- factor(c(x,seq(1,nsubj)))

ggplot(data=as.data.frame(lat_map), aes(lat_map[,1], lat_map[,2], shape=condition, col=factor(sex))) + 
  geom_point(aes(size=age)) + ggtitle("Subjects in latent mapping ") + 
  xlab("Component #1") + ylab("Component #2")



# trying 3D plots

library("scatterplot3d")
library("viridis")

shapes <- c(16, 17)
shapes <- shapes[factor(lat_map$condition)]
colors <- viridis_pal(option = "D")(length(unique(lat_map$group)))
colors <- colors[factor(lat_map$group)]

scatterplot3d(lat_map[,3:5], pch=shapes, color=colors, angle=20)


#trying t-SNE 
#TODO: test for larger data set


library("Rtsne")
D_ages <- ages[ages$File %in% rownames(D),]
#match the ordering that is mixed due to random sampling
reorder_indexes <- match(rownames(D), D_ages$File)
D_ages <- D_ages[reorder_indexes,]


tsne <- Rtsne(D, dims=2, is_distance = T,
              perplexity=5, verbose=TRUE, max_iter = 500, check_duplicates = F)

#plot(tsne$Y,col=colors, asp=1)
#scatterplot3d(tsne$Y, pch=4, color=D_ages$Age, angle=60)



#trying distance matrix with GGplot
# TODO: something is wrong
D_ages$Y1 <- tsne$Y[,1]
D_ages$Y2 <- tsne$Y[,2]


ggplot(D_ages, aes(x=Y1, y=Y2, color=Age)) +
  geom_point(size=4) + theme_minimal()




# Sammon mapping

#https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/sammon.html



#how distance correlates with age?
# -> make a matrix with shape identical to D, fill with age diffs and 
#    find correlation between these matrices (or each row vector I guess?)


G <- D*0
D_ages <- ages[ages$File %in% rownames(D),]
#match the ordering that is mixed due to random sampling
reorder_indexes <- match(rownames(D), D_ages$File)
D_ages <- D_ages[reorder_indexes,]

for(r in rownames(D)){
  age_diffs = abs(D_ages[D_ages$File==r,]$Age - D_ages$Age) #age difference
  G[r,] <- age_diffs
}
C=cor(t(D), t(G)) # calculates correlation between rows of G and D
heatmap(C)

cors=c()
for(r in 1:nrow(G)){cors=c(cors, cor(G[r,], D[r,]))}
D_ages$cors <- cors


#try hierarchical clustering with dissimilarity matrix

#create another distance matrix from test data
#TODO: calculate per component and do the clustering analysis
library("dendextend")

for(k in 1:1){
  S <- D*0 #S is dissimilarity matrix
  colnames(S) <- rownames(S)
  for(testidx in 1:nrow(D)){ #calculates the distances between individuals   
    for(m in 1:M){     
      group_members <- rownames(lat_map_n2)[m]   
      idxs = which(row.names(lat_map_n2) %in% group_members) #     
      #group_data <- lat_map_n2[idxs,][1:k]
      group_data <- lat_map_n2[c(idxs), 1:k]
      S[testidx,m] <- sum(abs(lat_map_n2[testidx,1:k]-group_data)) #L1 distance#  
    } 
  }
  
  colnames(S) <- D_ages$Cap
  rownames(S) <- D_ages$Cap
  
  hc <- hclust(as.dist(S), method="ward.D")
  plot(hc, main=paste0("Clustering on component #", 1, "-", k))
  
}


heatmap(S)

dend <- as.dendrogram(hc)
dend %>% set("branches_k_color", k = 5) %>% plot(main = "Nice defaults")


















