# This script reads in the data from csv files


set.seed(121)
savepath <- paste0(getwd(), "/data/n2C_sleep/") #folder that contains spectra
types <- list.files(savepath) #gives the number of groups/families/individuals

#for testing purposes, let us work only with a subset of individuals
#types <- sample(types, 320)

# Generate empty lists for classes and Y
individuals <- list()
Y <- list()

for(t in 1:length(types)) {

  type <- types[t]
  individuals[[type]] <- c()
  datapath <- savepath
  files <- type #pick individual file

  Reps <- matrix(NA,0,2)
  samples <- rep(NA,length(files))
  
  for(j in 1:length(files)) {
    file <- files[j]
    #takes MEG matrices per subject for temporary storage
    X <- t(read.table(paste0(datapath,files[j]),sep=",",header=F))
    subj <- file
    subj <- sub("n2_3.csv", "", subj)
    Y[[subj]] <- X #*1e11
    #Y[[length(Y)+1]] <- X*1e11
    individuals[[type]] <- c(individuals[[type]], subj) #Old naming
  }
  
}
subjects <- unlist(individuals)
names(Y) <- subjects
freq <- read.csv(paste0(getwd(),"/data/f.txt"),header=F)
fname <- paste0(getwd(), "/data/N2Cspectrum.RData") #name whichever spectra you may use 
save(Y, subjects, individuals, freq, file=fname) #saving data


