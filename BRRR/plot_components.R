######################################
# Based on script by Joonas Haakana  #
#                                    #
######################################

#TODO: add miniheads? (nope) + same individual, different data??

data <- "EEG_ind" # datatype

if(data == "EEG_ind") { #TODO: fix repetition
  #fname <- paste0("results/full/over7_2N1_BRRR_K12.RData")
  file = paste0("results/full/NEW_alldata_2N1_BRRR_K12.RData")
  load(file)
  datafile <- paste0("data/N1Aspectrum.RData") #only to get frequencies
  load(datafile)


  keepFeat <- which(apply(res$data$phenotypes,2,var,na.rm=T)>0)
  #Y <- cbind(res$scaling[,1:5], t(res$model$A)) #to get confounder weight matrix too
  Y <- res$scaling 
  #Y <- t(exp (res$data$phenotypes) )[,1:12] #take 12 individuals, N2A data 
  net <- list(omitMag=F,Y=Y,keepFeat=keepFeat, #create the net object only now
              penLDA=4,omitN=2000,l1d=200)
  net$lambda <- 0
  
  
  net$freq <- freq[c(-1),]
  
  coords <- read.table("var/coords.csv",sep=",")
  ch_names <- read.table("var/ch_names.csv", sep=",")
  
} else if(data == 'raw'){
  

  spectra <- c('1A','1B','2A','2B','2C','2D')
  Ys <- matrix(ncol=length(spectra), nrow = 266) 
  #Ys <- vector(mode = 'list', length=length(spectra))
  
  for(i in 1:length(spectra)){
    spectrum=spectra[i]
    datafile <- paste0("data/N",spectrum,"spectrum.RData")
    load(datafile)
    #Ys[,i] <- exp( log10(c(Y$FLE151021)) ) #scaling in order to make differences more visible
    Ys[,i] <- c(Y$FLE151021)
  }
  
  Y=Ys
  keepFeat <- which(apply(Y,1,var,na.rm=T)>0)
  net <- list(omitMag=F,Y=Ys,keepFeat=keepFeat, #create the net object only now
              penLDA=4,omitN=2000,l1d=200)
  net$lambda <- 0
  
  
  net$freq <- freq
  
  coords <- read.table("var/coords.csv",sep=",")
  ch_names <- read.table("var/ch_names.csv", sep=",")
  
  
} else {
  fname <- paste0("results/full/age_over4_N2_BRRR_K6.RData")
  load(fname)
  datafile <- paste0("data/N2Cspectrum.RData")
  load(datafile)
  
  keepFeat <- which(apply(res$data$phenotypes,2,var,na.rm=T)>0)
  
  net <- list(omitMag=F,Y=res2$discrim,keepFeat=keepFeat, #create the net object only now
                        penLDA=4,omitN=2000,l1d=200)
  #net <- list(omitMag=F,Y=res$scaling,keepFeat=keepFeat, #create the net object only now
  #            penLDA=4,omitN=2000,l1d=200)
  net$lambda <- 0
  

  net$freq <- freq
  
  coords <- read.table("var/coords.csv",sep=",")
  ch_names <- read.table("var/ch_names.csv", sep=",")
}

coords = coords - 0.5
coords = coords * 200

####
####

filename <- paste0("figures/NEW_all_2N1_", data, ".pdf") # pdf file for saving plots
pdf(file=filename,width=20,height=30)
plotLabels <- paste0("K",1:ncol(net$Y)) # plot lables
if(data=='raw'){
  plotLabels <- paste0("N", spectra)
}

maxK <- ncol(net$Y) # number of latent components
n_bands <- nrow(net$freq) # number of frequency bands
n_channels <- nrow(net$Y)/n_bands # number of meg channels


# extract data 
X <- list()
b <- array(0,c(n_bands, n_channels))
for(k in 1:maxK) {
  X[[k]] <- array(0,c(n_bands, n_channels))
  X[[k]][net$keepFeat] <- net$Y[,k]
  a <- X[[k]]
  if(data=="hcp") {
    for(i in 1:n_channels){
      ch <- strtoi(gsub("A","",ch_names[i,1]))
      b[, ch] <- a[,i]
    }
    X[[k]][net$keepFeat] <- b 
  }
}

# define colors for plotting
scol <- colorRampPalette(c("#CCCC00","orange","#FF5555","gray"))(floor(nrow(net$freq)/2))
scol <- c(scol, colorRampPalette(c("#888888","#5555FF","#33CC33"))(ceiling(nrow(net$freq)/2)))


par(mfrow=c(3,2),mar=rep(1,4),oma=c(0,0,0,12)) # create subplot
for(K in 1:maxK) {

  xlim <- c(-100, 95) 
  ylim <- c(-110, 85)
  plot(NULL, xlim=xlim, ylim=ylim, xlab="", ylab="", bty="n", axes=F)
  
  text(-100,80,plotLabels[K],cex=4,pos=4,xpd=NA) # component number
  
  Y <- X[[K]]

  for(ch in 1:n_channels) {
    x0 <- coords[ch,1] #coordinate locations
    y0 <- coords[ch,2]
    x1 <- 9 #size parameters
    y1 <- 6
    
    y <- Y[,ch]
    y <- y/max(abs(Y))
    
    xseq <- x0+seq(-x1,x1,length=n_bands)
    
    xs <- c(xseq-diff(xseq[1:2])/2, xseq[length(xseq)+diff(xseq[1:2])/2])
    
    # plot components on scalp
    for(i in 1:n_bands){
      polygon(c(xs[i],xs[i+1],xs[i+1],xs[i]), c(0,0,y[i],y[i])*y1+y0,col=scol[i], border=NA)
    } 
    lines(xseq, y*y1+y0, type="l", pch=20, lwd=1, cex=0.2)
    lines(x0+c(-x1,x1)*1.15, rep(y0, 2), col="black")
    
  }
  
  # plot frequency colorbar
  if(K==6) {

    par(xpd=NA)
    xc <- rep((1:6)*0.8-(7*0.8)/2, 4)
    yc <- rep((4:1)*0.8-(7*0.8)/2, each=6)
    
    cba <- 110; 
    cbw <- 5; 
    yl <- 180*par()$mfrow[1]
    y0 <- -60
    
    # plot colorbar
    for(i in 1:n_bands){
      polygon(c(0,cbw,cbw,0)+cba, c(0,0,yl/n_bands,yl/n_bands)+yl*(i-1)/n_bands+y0, col=scol[i], border=NA)
    }
    
    lines(rep(cba+cbw,2), c(0,1)*yl+y0)
    freq <- c(net$freq[1,1], net$freq[,2])

    id50 <- which.max(freq[freq<50])
    freq[id50] <- paste0(net$freq[id50,1],"\n",freq[id50])

    
    for(f in 1:length(freq)) {
      h <- yl*(f-1)/n_bands+y0
      dd <- c(-cbw,cbw)/5
      lines(cba+cbw+dd, rep(h,2))
      if(nrow(net$freq)<50 | f%%5==0) text(cba+cbw*1.2,h,freq[f],pos=4,cex=sqrt(6))
    }
    text(cba+cbw+dd[1], y0-yl/n_bands*c(2,1,1)[3], "Freq.\n(Hz)", pos=3, cex=sqrt(6))
    par(xpd=F)
    
  }
  
}


dev.off()
