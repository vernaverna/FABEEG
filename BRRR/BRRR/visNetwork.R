visNetwork <- function(net,maxK=98,onlyPdf=TRUE,graph=TRUE,changePar=T,norm="comp",condense=FALSE,
                       ltype="l",levelplot=FALSE,consistent=FALSE,lineCols=NULL,diffOrder=F,
                       plotLabels=NULL,matchSign=FALSE,topLeg=NULL,tightMargin=TRUE,lty=NULL) {
  library(ellipse)
  if(all(is.na(net$Y)) && !is.null(net$X)) net$Y <- t(net$X)
  if(is.null(plotLabels)) plotLabels <- paste0("K",1:ncol(net$Y))
  maxK <- min(maxK,ncol(net$Y))
  lda <- TRUE; l <- 2
  D <- c(nrow(net$freq),nrow(net$Y)/nrow(net$freq))
  pcaK <- NA #net$pcaK  
  X <- list()
  maxK <- min(ncol(net$Y),maxK)
  P <- if(length(dim(net$Y))==2) 0 else dim(net$Y)[3]
  print(P)
  for(k in 1:maxK) {
    if(is.na(pcaK)) {
      X[[k]] <- array(0,c(D[1],D[2],max(P,1)))
      if(P==0) {
        X[[k]][net$keepFeat] <- net$Y[,k]
      } else {
        for(p in 1:P) {
          s <- if(matchSign && p>1) sign(cor(net$Y[,k,p], net$Y[,k,1])) else 1 
          X[[k]][,,p] <- s*net$Y[,k,p]
        }
      }
      
    } else {
      stop("TODO")
      W <- matrix(NA,pcaK,pcaK); diag(W) <- 0
      W[is.na(W)] <- net$Y[,k]
      X[[k]] <- matrix(0,D,D)
      for(p1 in 1:pcaK) {
        for(p2 in 1:pcaK) {
          X[[k]] <- X[[k]] + W[p1,p2]*outer(net$pca[,p1],net$pca[,p2])
          X[[k]] <- X[[k]] + W[p2,p1]*outer(net$pca[,p2],net$pca[,p1])
        }
      }
      diag(X[[k]]) <- 0
    }
  }
  D <- D[2]
  Y <- X[[1]]
  #}
  coords <- read.table("var/coords/coords.csv",sep=",")
  chNames <- read.table("var/coords/ch_names.csv",sep=",",stringsAsFactors=FALSE)[1:306]
  chNames <- sub("MEG","",chNames)
  if(D==204) {
    coords <- coords[-seq(3,306,by=3),]
    chNames <- chNames[-seq(3,306,by=3)]
  } else if(D%in%c(16,24,68,102)) {
    if(D>30) {load("data/MEGclust3.RData"); names(clustId) <- 1:length(clustId)}
    else load("data/MEGclusters.RData")
    chNames <- rep(names(clustId), each=c(2,3,2,3)[D==c(16,24,68,102)])
    cc <- matrix(NA,D,2)
    id <- 1:2; if(D%in%c(24,102)) id <- c(id,0)
    for(cl in 1:length(clustId))
      for(i in id) {
        cc[(cl-1)*length(id)+which(id==i),] <- colMeans(coords[clustId[[cl]][clustId[[cl]]%%3==i],])
      } 
    coords <- cc
    browser()
  }
  
  columns <- 2
  rows <- min(3, ceiling(maxK/columns))
  if(condense) {columns <- condense; rows <- 1}
  
  filename <- paste0("figs/",sub(".RData",".pdf",sub("results/full/","",net$fname)))
  if(levelplot) filename <- sub(".pdf","level.pdf",filename)
  print(filename)
  width <- if(levelplot) 8 else 20
  height <- if(levelplot) 4 else 10*rows
  if(condense) {height <- 3; width <- height*columns}
  if(onlyPdf) pdf(file=filename,width=width,height=height)
  
  if(lda) {
    ends <- 1:maxK
    if(changePar) {
      side <- if(condense) 0 else 12
      par(mfrow=c(rows, columns),mar=rep(1,4),oma=c(0,0,0,side))
    }
    if(P>1 && length(X)>200) { #--Most varying subjects first--
      #X[[k]] <- array(0,c(D[1],D[2],P))
      subjSd <- sapply(X,function(x) mean(apply(x,1:2,sd)))
      
      if(diffOrder) ends <- order(subjSd,decreasing=TRUE)[1:maxK] 
      else ends <- order(net$fam,decreasing=FALSE)[1:maxK] #Family order
    }
    
  } else {
    ends <- c(100, 500, 1500, -100, -1000)
    if(ncol(Y)==16) ends <- c(1, 5, 10, 20, 50)
    if(changePar) par(mfrow=c(rows,columns),mar=rep(1,4))
  }
  
  if(!lda) {
    hist(Y,100,xlab="",ylab="",yaxt="n",main=paste("Distance to siblings (1 to unrelated):",
                                                   paste0(net$exp,collapse=",")))
    if(is.null(net$freq)) diag(Y) <- 1
    id <- order(Y)
    id <- arrayInd(id,dim(Y))
    for(end in ends) {
      if(end>0) abline(v=Y[id[end,,drop=F]],col="blue")
      else abline(v=Y[id[nrow(id)+end,,drop=F]],col="red")
    }
  }
  if(!is.null(net$freq)) {
    #if(lda) {
    #  plot(NULL,xlim=c(.5,1.5),ylim=c(0,1),xlab="",ylab="",axes=F)
    #} else {
    #  ends <- c(1, 5, 10, 50, 2000)
    #}
    if(nrow(Y)<=24) {by = 0.8; d <- c(6,4); ptcex <- 1.8}
    else {by = 0.3; d <- c(16,13); ptcex <- 0.5}
    xc <- rep((1:d[1])*by-((d[1]+1)*by)/2, d[2])
    yc <- rep((d[2]:1)*by-((d[2]+1)*by)/2, each=d[1])
    scol <- colorRampPalette(c("#CCCC00","orange","#FF5555","gray"))(floor(nrow(net$freq)/2))
    scol <- c(scol, colorRampPalette(c("#888888","#5555FF","#33CC33"))(ceiling(nrow(net$freq)/2)))
    #Colorbar
    n <- nrow(net$freq)#+1
    cba <- 1; cbw <- 1/50; yl <- par("usr")[4]*.9
  }
  
  cols <- colorRampPalette(c("blue","white","red"))(201)
  
  if(levelplot) {
    load("data/MEGclusters.RData")
    clust <- clust[c(1:4,8:5)]
    Ks <- maxK #6
    waves <- c("Delta","Theta","Alpha","Beta","Gamma")
    #matrix(c(0,3,4,7,8,13,16,31,32,100),nrow=2)
    freqs <- list(1,2:3,4:5,7:10,11:21); names(freqs) <- waves
    Y <- array(NA,dim=c(length(clust),length(waves),Ks),
               dimnames=list(names(clust),waves,paste0("K",1:Ks)))
    if(!is.null(plotLabels)) dimnames(Y)[[3]] <- plotLabels
    for(k in 1:Ks) {
      for(area in names(clust)) {
        for(wave in waves) {
          chs <- which(as.numeric(chNames)%in%clust[[area]])
          f <- freqs[[wave]]
          Y[area,wave,k] <- mean(X[[k]][f, chs, 1])
        }
      }
    }
    cols <- colorRampPalette(rev(c("orange","red","white","blue","cyan")))
    ax <- list(alternating=1,tck=c(1,0))
    kk <- c(4:6,1:3)
    while(min(kk)<=maxK) {
      kk <- kk[kk <= dim(Y)[3]]
      M <- max(abs(Y[,,kk]),na.rm=T)
      print(levelplot(Y[,ncol(Y):1,kk],col.regions=cols,at=seq(-M,M,length=100),
                      xlab="Region",ylab="Frequency band",layout=c(min(Ks,3),min(2,ceiling(Ks/3))),
                      scales=list(x=ax,y=ax), panel = function(...){
                        panel.levelplot(...); panel.abline(v=length(clust)/2+.5)}))
      kk <- kk+6
    }
    
    
  } else {
    for(end in ends) {
      xlim <- range(coords[,1])
      xlim <- xlim+c(-10,10)
      ylim <- range(coords[,2]); ylim[2] <- ylim[2]+5
      plot(NULL,xlim=xlim,ylim=ylim,xlab="",ylab="",bty="n",axes=F)
      #draw.circle(0, -10, 80, nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)
      el <- ellipse(diag(1,2), scale=c(38.5,31.5), centre=c(-3,-6), npoints=100)
      lines(el[,1], el[,2]) #Scalp
      if(condense) lines(c(-10,0,10)-3, c(0,10,0)+4.7+max(coords[,2]),xpd=NA) #Nose
      else lines(c(-5,0,5)-3, c(0,5,0)+5+max(coords[,2])) #Nose
      ff <- if(is.null(net$fam)) "" else paste0("F",net$fam[end])
      lSize <- 2.8
      lpos <- if(condense) c(-110,67) else c(-97,65)
      if(nchar(plotLabels[end])>10) {
        lSize <- 1.7; lpos <- c(-97, 65) #-95, 62
      }
      if(is.null(plotLabels)) text(-70,60,paste(ff,c("Top","K")[l],end),cex=3)
      else text(lpos[1],lpos[2],plotLabels[end],cex=columns*sqrt(rows/3)*lSize,pos=4,xpd=NA)
      if(!graph & !is.null(net$freq))
        for(ch in 1:nrow(coords)) text(coords[ch,1], coords[ch,2], chNames[ch],cex=0.7)
      if(end>0) A <- 1:end
      else A <- nrow(id):(nrow(id)+end)
      if(!is.null(topLeg)) legend("topright",topLeg$text,col=topLeg$col,lty=1)
      
      if(lda) {
        Y <- X[[end]]
        #id <- order(abs(Y),decreasing=T)
        #id <- arrayInd(id,dim(Y))
        #A <- 1:nrow(id) #min(1000,nrow(id))
        cc <- 0
      } else {
        cc <- 1
      }
      M <- max(abs(Y-cc),na.rm=T) #Y[id[A,,drop=F]]
      f <- function(x,p=c(4,1)[l]) {y <- abs((x-cc)/M)^p; y[y<0.2] <- 0; return(y)}
      g <- function(x,p=c(3,1)[l]) {cols[round(((x-cc)/M)^p*100)+101]}
      
      if(graph & !is.null(net$freq)) {
        chs <- if(condense) c(9,123,157,185,95,59,49) else 1:D
        for(ch in chs) {
          x0 <- coords[ch,1]; y0 <- coords[ch,2]
          x1 <- 5; y1 <- 4
          if(condense) {x1 <- 20; y1 <- 20}
          if(!tightMargin) y1 <- 3
          #lines(x0+c(-x1,x1,x1,-x1,-x1), y0+c(-y1,-y1,y1,y1,-y1), col="#EEEEEE")
          
          y <- Y[,ch,]
          if(norm=="comp") y <- y/max(abs(Y))
          if(norm=="all") y <- y/max(abs(unlist(X)),na.rm=T)
          #Colorbar
          xseq <- x0+seq(-x1,x1,length=nrow(Y))
          if(P<=1) {
            xs <- c(xseq-diff(xseq[1:2])/2, xseq[length(xseq)+diff(xseq[1:2])/2])
            cba <- 1 ;yl <- par("usr")[4]*.9
            if(FALSE) { #Channel numbering
              text(mean(xseq),y0,ch,cex=0.7)
            } else {
              lwd <- 0.2
              if(n<50) lwd <- 1
              if(is.null(net$noColor)) {
                for(i in 1:n) #Coloring
                  polygon(c(xs[i],xs[i+1],xs[i+1],xs[i]), c(0,0,y[i],y[i])*y1+y0,col=scol[i], border=NA)
              } else {
                if(n<50) lwd <- 4
              }
              
              lines(xseq, y*y1+y0, type=ltype, pch=20, lwd=lwd, cex=0.2)
              #lines(rep(coords[ch,1],2), coords[ch,2] + range(yc))
              if(ltype!="p" & is.null(net$no0)) lines(x0+c(-x1,x1)*1.15, rep(y0, 2), col="black") #w=0
              if(ltype=="p")
                for(i in seq(5,nrow(Y),by=4))
                  lines(rep(mean(xseq[(i-1):i]),2), c(-y1,y1)/2+y0, lwd=0.1)
            }
            
          } else if(!all(is.na(y))) {
            if(is.null(lty)) lty <- rep(1,P)
            lwd <- if(P>10) 0.2 else 1
            if(is.null(lineCols)) lineCols <- c("black",rep("gray",P-1))
            if(consistent | TRUE) lines(x0+c(-x1,x1)*1.15, rep(y0, 2), col="black") #w=0
            for(p in P:1) { #Reversed order to get the first to top
              lines(xseq, y[,p]*y1+y0, type=ltype, pch=20, lwd=lwd, col=lineCols[p], lty=lty[p])
              if(all(y[,p] < min(y[,-p])) && !consistent)
                text(max(xseq),y[nrow(y),p]*y1+y0, p-1, col="red")
            }
          }
        }
        
      } else {
        for(i in rev(A)) {
          chs <- id[i,]
          x <- Y[id[i,,drop=F]]
          if(is.null(net$freq))
            lines(coords[chs,1], coords[chs,2], lwd=f(x), col=g(x))
          else {
            points(coords[chs[2],1]+xc[chs[1]], coords[chs[2],2]+yc[chs[1]], pch=c(15,NA,20)[sign(x-cc)+2],
                   col=scol[chs[1]], cex=abs(f(x,1/2))*ptcex*c(0.6,NA,1)[sign(x-cc)+2])
            #sqrt(abs(1-x))^2*c(1,9)[l]*c(1,NA,3)[sign(x-1)+2]
          }
        }
      }
      
      
      if(is.null(net$freq)) {
        if(lda) {
          tmp <- rowMeans(Y)+colMeans(Y)
          col <- g(tmp/max(abs(tmp))*M)
          p <- rowMeans(abs(Y))+colMeans(abs(Y))
          p <- p/max(p)*3
        } else {
          col <- g(1-sign(1-Y[id[A[1],,drop=F]])*M*0.9)
          p <- tabulate(id[A,])
          p <- p/max(p)
        }
        points(coords[1:length(p),1], coords[1:length(p),2], pch=20,
               cex=sqrt(p)*c(8,2)[l], col=col)
      }
      
      
      ## COLORBAR
      if(!condense && !is.null(net$freq) && which(ends==end)%%(rows*columns)==0) {
        par(xpd=NA)
        xc <- rep((1:d[1])*by-((d[1]+1)*by)/2, d[2])
        yc <- rep((d[2]:1)*by-((d[2]+1)*by)/2, each=d[1])
        cba <- 100; cbw <- 5; yl <- 140*par()$mfrow[1]
        y0 <- -60
        for(i in 1:n){
          polygon(c(0,cbw,cbw,0)+cba, c(0,0,yl/n,yl/n)+yl*(i-1)/n+y0,
                  col=scol[i], border=NA)
        }
        #Colorbar: axis
        lines(rep(cba+cbw,2),c(0,1)*yl+y0)
        freq <- c(net$freq[1,1],net$freq[,2])
        if(nrow(net$freq)==21) {
          id50 <- which.max(freq[freq<50])
          freq[id50] <- paste0(net$freq[id50,1],"\n",freq[id50])
        }
        for(f in 1:length(freq)) {
          h <- yl*(f-1)/n+y0
          dd <- c(-cbw,cbw)/5
          lines(cba+cbw+dd, rep(h,2))
          if(nrow(net$freq)<50 | f%%5==0) text(cba+cbw*1.2,h,freq[f],pos=4,cex=sqrt(rows*3))
        }
        text(cba+cbw+dd[1], y0-yl/n*c(2,1,1)[rows], "Freq.\n(Hz)",pos=3,cex=sqrt(rows*3)) #h+yl/n/2
        par(xpd=F)
        
        if(TRUE) { #miniheads
          library("plotrix")
          for(i in 1:2) {
            xm <- c(100,100)[i]
            ym <- c(375, 398)[i]
            el <- ellipse(diag(1,2), scale=c(38.5,31.5)/9, centre=c(xm,ym), npoints=100)
            lines(el[,1], el[,2],xpd=NA)
            lines(c(-10,0,10)/5+xm, c(0,10,0)/5+ym+8.4,xpd=NA) #Nose
            l <- 8
            arLwd <- 4
            if(i==1) {
              ang <- 25
              draw.arc(x=xm, y=ym, radius=l, deg1=90+ang, deg2=360+90-ang, lwd=arLwd, xpd=NA)
              pos <- function(phi, l=6.6) { #Why does this need a separate l?
                x0 <- xm+l*cos(phi/360*2*pi)
                y0 <- ym+l*sin(phi/360*2*pi)
                return(c(x0, y0))
              }
              arrows(pos(90+ang)[1], pos(90+ang)[2], pos(90+ang-1)[1], pos(90+ang-1)[2], length=0.1, lwd=arLwd, xpd=NA)
              
            } else if(i==2) {
              for(j in 1:3) {
                lx <- if(j==1) l else sign(j-2.5)*l/2
                ly <- if(j==1) 0 else l*2/3
                arrows(xm-lx, ym-ly, xm+lx, ym+ly, code=3, length=0.1, lwd=arLwd, xpd=NA)
              }
            }
          }
        }
      }
      
      
      if(end==1 && FALSE) {
        print("Shown Y range:")
        print(range(Y[id[A,]]))
        print("Resulting in lwd range:")
        print(range(f(Y[id[A,]])))
      }
    }
  }
  
  if(!onlyPdf) dev.copy(pdf,file=filename,width=width,height=height)
  dev.off()
  print(filename)
  
  #MOve
  #if(returnData) return(list(X=X,Y=Y,fam=fam,freq=freq,keepFeat=keepFeat))
  return(Y)
}
