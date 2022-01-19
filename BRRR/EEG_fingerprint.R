# operate the process here...
setwd("/projects/FABEEG/BRRR")
library("penalizedLDA")
library("ggplot2")
library("cvms")


## TODO: use less individuals! try with old and young kids as well

# Function for preparing the data for analysis
# ex = N1 or N2, depending which data is read in
prepare_data <- function(ex){
  
  # Choosing the channels and frequencies
  chs <- 1:19
  omitFreq = c(0,60)
  
  # Reading data into workspace
  loadfile <- paste0("data/",ex,"spectrum.RData")
  
  #extract the age groups
  ages = read.csv('data/new_age_df.csv')
  ages <- ages[,-1]
  
  use_all=F #should we use all subjects in training?
  #age_gap=c(1,19) #exclude some of the younger children?
  Cap='FT'
  #TODO: change to range?
  
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
    
    corrupted = c()
    
    # Saving the data into matrix A
    for(s in names(Y)) {
      
      if(s %in% ages$File){ #temporal solution to get over the filename hassle
        age_data = ages[ages$File==s,]
        #if(age_data$Age >= age_gap[1] && age_data$Age <= age_gap[2]){
        if(age_data$Cap == Cap){
          tmp <- t(Y[[s]]) #transposed
          tmp <- log10(tmp) #TODO: is this necessary???
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
    #TODO: even age groups!!!!!
    
    if(length(corrupted)>0){
      obs <- subjects[subjects %in% corrupted == FALSE]
      A = A[names(A) %in% corrupted == FALSE]
      ages = ages[ages$File %in% corrupted == FALSE,]
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
  
  
  for(i in obs) { 
    tmp <- A[[i]]
    if(nrow(tmp)==ncol(tmp)) tmp[lower.tri(tmp)] <- 0
    Y[i,] <- c(tmp)
  }
  
  keepFeat <- which(apply(Y,2,var,na.rm=T)>0)
  Y <- Y[,keepFeat]
  
  print("LDA matrix dimension:")
  print(dim(Y))
  
  # Scaling & centering the Y for BRRR 
  tmp <- scale(Y[!is.na(Y[,1]),],center=T,scale=T)
  Y <- scale(Y, attr(tmp,"scaled:center"), attr(tmp,"scaled:scale"))
  
  
  # Saving classes to x
  x <- rep(NA,S); names(x) <- obs
  for(s in 1:S){
    x[obs[s]] <- s
  }
  
  # Making design matrix out of classes
  X <- matrix(0,S,M,dimnames=list(obs,paste0("class",1:M))) 
  for(i in 1:length(x)) if(!is.na(x[i])) X[i,x[i]] <- 1
  if(M==2) {
    print("Condensing two classes into one +-1 covariate.")
    X <- X[,1,drop=FALSE] 
    X[X==0] <- -1
  }
  
  return(list(Y, x, X, M, ages))
}


n2_data <- prepare_data(ex="N2B")
n2b_data <- prepare_data(ex="N2A")
n2c_data <- prepare_data(ex="N2D")


validation_data <- prepare_data(ex="N2C")
Y4 <- validation_data[[1]]

Y1 = n2_data[[1]]
Y2 = n2b_data[[1]]
Y3 = n2c_data[[1]]

if(length(rownames(Y3)) != length(rownames(Y2))){
  to_del= which(!rownames(Y2) %in% rownames(Y3))
  Y2 <- Y2[-to_del,]
  Y1 <- Y1[-to_del,]
  Y4 <- Y4[-to_del,]
}


Y = rbind(Y1, Y2, Y3)

x = c(n2c_data[[2]], n2c_data[[2]], n2c_data[[2]])
X = rbind(n2c_data[[3]], n2c_data[[3]], n2c_data[[3]]) #row binding
M = n2c_data[[4]]
subj = dimnames(Y)[[1]]
ages = validation_data[[5]]

#Trying with added Z
#TODO: put inside prepare-data!
Z1 = ages[ages$File %in% rownames(Y1),]
reorder_idx <- match(rownames(Y1), Z1$File)
Z1 <- Z1[reorder_idx,]

Z = matrix(c(Z1$Age, Z1$Age, Z1$Age)) 

# TODO: add CV here, var.R does not work at all
### TRAINING ###

# The model is trained using two sets of N2 data, and the performance is evaluated using third set

source("brrr.R")
#pred <- X*NA
res <- brrr(X=X,Y=Y,K=6,Z=NA,n.iter=500,thin=5,init="LDA", fam =x) #fit the model
res$scaling <- ginv(averageGamma(res))
W <- res$scaling

#save(res, file = "results/full/293_indN2_BRRR_K6.RData")
lat_map <- Y%*%W
lat_map_n2 <- Y3%*%W #mapping to latent space with N2_C data!# 

D <- matrix(NA, nrow(lat_map_n2), ncol(X), dimnames=list(unique(names(x)), paste0("mean",unique(names(x))) ) ) #distance matrix# # # 

for(testidx in 1:nrow(lat_map_n2)){ #calculates the distances between individual and group mean in lat.space   
  for(m in 1:M){     
    group_members <- rownames(lat_map)[m]   
    idxs = which(row.names(lat_map) %in% group_members) #     
    group_mean <- colMeans(lat_map[idxs,]) #mean over all individuals, vector of length K
    D[testidx,m] <- sum(abs(lat_map_n2[testidx,]-group_mean)) #L1 distance #TODO: CHECK this idiot
    #D[testidx,m] <- sqrt(sum( (lat_map_n2[testidx,]-group_mean)**2 ))#L2 distance  
  } 
}





#### VALIDATION ####

PROJ <- D*0 #initialize 'projection' matrix

for(r in 1:nrow(D)){ #assign age groups based on lat. space distances
  index <- which.min(D[r,])
  PROJ[r,index] <- 1+PROJ[r,index]
}
colnames(PROJ) <- c(paste0("class",1:M))

accuracy = sum(diag(PROJ))/M
print(paste("Model accuracy:", accuracy))


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

library("Rtsne")
age_col = round(ages[which(ages$File%in%subj),]$Age, 0)
sex_col = ages[which(ages$File%in%subj),]$Sex

colors = rainbow(length(unique(age_col)))
names(colors) = unique(age_col)

tsne <- Rtsne(D, dims=2, is_distance = T,
              perplexity=4, verbose=TRUE, max_iter = 500, check_duplicates = F)
#plot(tsne$Y, t='n', main="tsne")
#text(tsne$Y, labels=x, col=colors[x])

plot(tsne$Y,col=colors, asp=1)



#trying distance matrix with GGplot
# TODO: something is wrong
data_df <- data.frame(tsne$Y)
rownames(data_df) <- rownames(D)
data_df$age <- ages[which(ages$File%in%subj),]$Age
data_df$sex <- ages[which(ages$File%in%subj),]$Sex
data_df$group <- round(ages[which(ages$File%in%subj),]$Age, 0)

ggplot(data_df, aes(x=X1, y=X2, color=age, shape=sex)) +
  geom_point(size=4) + theme_minimal()



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


















