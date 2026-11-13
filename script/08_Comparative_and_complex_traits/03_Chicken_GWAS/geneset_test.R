#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Geneset test.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

################################################################################
#  gfm functions
################################################################################

  computeG <- function(W=NULL,pdf=TRUE) {
    SS <- tcrossprod(W)        # compute crossproduct all SNPs
    N <- tcrossprod(!W==0)     # compute number of observation all SNPs
    G <- SS/N
    if (pdf) G <- makepdf(G)
    return(list(G=G,SS=SS,N=N))
  }  
  
  
  cvgfm <- function(fm=NULL, vfm=NULL, data=NULL, Glist=NULL, validate=NULL, tol=1e-6, control=NULL) {
    yvar <- unlist( strsplit( as.character(fm),split="~", fixed=TRUE))[2]
    cvar <- sapply(vfm, function(x) { unlist( strsplit( as.character(x),split="~", fixed=TRUE))[2] })
    ncvar <- length(cvar)
    # full data set
    fit <- gfm(fm=fm, vfm=vfm, data=data, Glist=Glist, tol=tol, control=control)
    blup <- vector(mode="list",length=ncvar)
    for (j in 1:ncvar) {
      blup[[j]] <- cbind(blup[[j]],fit$blup[fit$blup[,2]==j,5])        
    }
    # training data set
    cvblup <- vector(mode="list",length=ncvar)
    for ( i in 1:ncol(validate)){
      cvdata <- data
      miss <- (as.character(cvdata[,cvar[1]]))%in%validate[,i]
      cvdata[miss,yvar] <- -999999
      cvfit <- gfm(fm=fm, vfm=vfm, data=cvdata, Glist=Glist, tol=tol, control=control)
      for (j in 1:ncvar) {
        cvblup[[j]] <- cbind(cvblup[[j]],cvfit$blup[cvfit$blup[,2]==j,5])        
      }
    }
    if (ncvar==1) blupt <- blup[[1]]
    if (ncvar==2) blupt <- blup[[1]] + blup[[2]]
    if (ncvar==1) cvblupt <- cvblup[[1]]
    if (ncvar==2) cvblupt <- cvblup[[1]] + cvblup[[2]]
    cors <- beta <- NULL
    for ( i in 1:ncol(validate)) {
      rws <- validate[,i]
      cors <- c(cors,cor(blupt[rws,1], cvblupt[rws,i] ))
      beta <- c(cors,lm(blupt[rws,1]~ cvblupt[rws,i] )$coefficients)
    }
    beta <- matrix(beta,ncol=2, byrow=TRUE)
    return(list(cvblup=cvblup,blup=blup, cors=cors, beta=beta, fit=fit ))
  }

   ginv <- function(G=NULL, tol=NULL) {
     rn <- rownames(G)
     e <- eigen(G)                                    # eigen value decomposition of the matrix G
     U <- e$vectors                                   # eigen vectors
     e <- e$values                                    # eigen values
     ie <- e
     ie[e>tol] <- 1/e[e>tol]
     #ie[e<tol] <- tol
     ie[e<tol] <- 0
     D <- diag(ie)                                   # set inverse D to 1/e
     G <- U%*%D%*%t(U) 				                        # compute inverse
     ldet <- sum(log(e[e>tol])) 
     #ldet <- sum(log(e)) 
     colnames(G) <- rownames(G) <- rn
     return(list(G=G,ldet=ldet))                      # log determinant 
   } 
     
   vec2mat <- function(vec=NULL,n=NULL, rcnames=NULL) {
   	 X <- diag(n)
   	 X[lower.tri(X, diag=TRUE)] <- vec
   	 X <- X + t(X) - diag(diag(X))  
   	 if(!is.null(rcnames)) { rownames(X) <- colnames(X) <- rcnames}
     X
   }
    
   writeG <- function(G=NULL, filename=NULL, clear=TRUE, ldet=NULL) {
     if (clear) {file.remove(filename)}
     nr <- nrow(G) 
     if(!is.null(ldet)) { write.table( t(c(0,0,ldet )),filename,quote=F,sep=" ",row.names=F, col.names=F, append=TRUE)}
     for (i in 1:nr) {
      out <- data.frame(rownames(G)[i], rownames(G)[i:nr], G[i:nr,i])
      write.table( out,filename,quote=F,sep=" ",row.names=F, col.names=F, append=TRUE)
     }
   }
    
   gfm <- function(fm=NULL, vfm=NULL, data=NULL, Glist=NULL, tol=1e-6, control=NULL) {
     #library(Rdmu)
     corlist <- NULL
       fmsplit <- unlist(strsplit(as.character(fm)[3],split="(", fixed=TRUE))
       rws <- grep(")", fmsplit,fixed=TRUE)
       fmsplit <- fmsplit[rws]
       fmsplit <- gsub("1 | ", "", fmsplit, fixed=TRUE)
       fmsplit <- gsub("+", "",fmsplit, fixed=TRUE)
       fmsplit <- gsub(")", "",fmsplit, fixed=TRUE)
       vvar <- gsub(" ", "",fmsplit, fixed=TRUE)
     
     if (!is.null(vfm)) {
       cvar <- NULL
       cvar <- sapply(vfm, function(x) {
          x <- unlist( strsplit( as.character(x),split="~", fixed=TRUE))[-1]
       })
       cvar <- t(cvar)
       rownames(cvar) <- cvar[,1]
       cvar <- cvar[,-1]

       corG <- NULL
       for ( i in 1:length(cvar) ) {
        iG <- ginv(Glist[[cvar[i]]], tol=tol)
        fn <- paste(names(cvar)[i],".txt",sep="")
        writeG(G=iG$G, filename=fn, ldet=iG$ldet)
        corG[[i]] <- COR(data[,names(cvar)[i]], file=fn)
        attr(corG[[i]], "varname") <- names(cvar)[i]

       }
       names(corG) <- names(cvar) 
       corlist <- corG
     }
  
     fit <- Rdmuai( formula=fm, corlist=corlist, data=data, control=c(control,list("SOL"=TRUE, "residuals"=TRUE)) )
     #fit <- dmuai( formula=fm, corlist=corlist, data=data, control=c(control,list("SOL"=TRUE, "residuals"=TRUE)) )

     # added this 16/2 2015  
     if (is.logical(fit$SOL)) {
       fit$SOL <- as.matrix(read.table("test.SOL", as.is = TRUE)[-1, , drop = FALSE])
       fit$PAROUT <- as.matrix(read.table("test.PAROUT"))
       fit$LLIK <- scan("test.LLIK", what = character(0), quiet = TRUE)
     }
     # end added this 16/2 2015

     blue <- fit$SOL[fit$SOL[,1]==2,c(2,4:6,8:9)]
     blup <- fit$SOL[fit$SOL[,1]==3,c(2,4:6,8:9)]
     residuals <- fit$residuals
     fitted <- fit$fitted
     rownames(blup) <- blup[,3]
     
     if (!is.matrix(blue)) names(blue) <- c("Trait","Fixed Effect-ID", "Level-ID", "N-obs", "E-BLUE", "SE") 
     if (is.matrix(blue)) colnames(blue) <- c("Trait","Fixed Effect-ID", "Level-ID", "N-obs", "E-BLUE", "SE") 
     colnames(blup) <- c("Trait","Random Effect-ID", "Level-ID", "N-obs", "E-BLUP", "SEP") 
 
     sigma <- fit$PAROUT[,4]
       names(sigma) <- c(vvar,"E")
     llik <- as.numeric(fit$LLIK[12])
     res <- list(sigma=sigma, blup=blup, blue=blue, llik=llik, residuals=residuals,fitted=fitted)
     return(res)
   }

    
################################################################################ 
#  Functions needed for weighted G-Blup 
################################################################################
  
      wblup <- function(fm=NULL,vfm=NULL,W=NULL, miss=NULL, sets=NULL, data=NULL, niter=3, tol=0.00001) {

       fmsplit <- unlist(strsplit(as.character(fm)[3],split="(", fixed=TRUE))
         rws <- grep(")", fmsplit,fixed=TRUE)
         fmsplit <- fmsplit[rws]
         fmsplit <- gsub("1 | ", "", fmsplit, fixed=TRUE)
         fmsplit <- gsub("+", "",fmsplit, fixed=TRUE)
         fmsplit <- gsub(")", "",fmsplit, fixed=TRUE)
       vvar <- gsub(" ", "",fmsplit, fixed=TRUE)
       fmr <- paste("~ 0 + ", vvar)
  
       fmf <- unlist(strsplit(as.character(fm)[3],split="+ (", fixed=TRUE))[1]
       fmf <- paste("~ ", fmf)

       vfm <- list(as.formula(paste(vvar[1], "~ G")))
       
       rvar <- as.character(fm)[2]
       data <- droplevels(data)
       
       y <- data[,rvar]   
        n <- length(y)
        m <- ncol(W)      
        d <- rep(1,m)                   	# initial weights
        b <- vb <- NULL
        g <- NULL
        sigma <- NULL
        llik <- NULL
        trH <- trT <- NULL

       X <- model.matrix(as.formula(fmf), data=data)	# design matrix for fixed effects
       Z <- NULL    				                          # design matrix for random effects
    	  nr <- length(fmr)
    	  nlr <- NULL 
    	  for ( i in 1:nr ){
     	   Zf <- model.matrix(as.formula(fmr[i]),data=data) 
       	  if (i==1) colnames(Zf) <- gsub(vvar[i],"",colnames(Zf))
       	  if (i==1) Zf <- Zf[,rownames(W)] 
       	  nlr[i] <- ncol(Zf) 
     	   Z <- cbind(Z,Zf)	
    	  }
        for ( i in 1:niter) {
          WD <- t(t(W)*sqrt(as.vector(d)))
          WW <- WD%*%t(WD)
          G <- WW/sum(d)
          fit  <- gfm(  fm=fm, vfm=vfm, data=data, Glist=list(G=G), tol=tol )
            Sg <- fit$sigma[1]
            Se <- fit$sigma[nr+1]
            print(fit$sigma)
            ghat <- fit$blup[colnames(G),5]
          Gi <- ginv((WW/sum(d))*Sg,tol=tol)$G
          iG <- ginv((WW/sum(d)),tol=tol)$G
          Gci <- diag( rep(1/fit$sigma[1:nr], times=nlr) )    
          Gci[1:nlr[1],1:nlr[1]] <- Gi
    	    mme <- setMME( X=X, Z=Z, Gci=Gci, Ri=1/Se, y=y, G=G*Sg, tol=tol)
            vghat <- mme$vghat
            trH <- rbind(trH,mme$trH)
            remove(mme)  
          WWi <- ginv(WW,tol=tol)$G
          bhat <- crossprod(WD,WWi%*%ghat)    
          WWWi <- crossprod(WD,WWi)
            remove(WWi)
            remove(WD)
          WWWivg <- tcrossprod(WWWi,vghat)
    	    vbhat <- rowSums(WWWivg*WWWi)
    	      remove(WWWi)
    	      remove(vghat)
          d <- abs(bhat)/sqrt(vbhat)
            d[bhat==0] <- 0            # check this
          d <- (m*d)/sum(d)
            #plot(d)
            print(sum(d))
          llik <- c(llik,fit$llik)
          sigma <- rbind(sigma,fit$sigma)
          b <- cbind(b,bhat)
          vb <- cbind(vb,vbhat)
    	    g <- cbind(g,ghat)
        }
        return(list(g=g,b=b,vb=vb, sigma=sigma,llik=llik, trH=trH, n=n))
    }
    
    
    
   makepdf <- function(G=NULL, tol=0.0001) {
     rn <- rownames(G)
     e <- eigen(G)                                    # eigen value decomposition of the matrix G
     U <- e$vectors                                   # eigen vectors
     e <- e$values                                    # eigen values
     e[e<tol] <- tol
     D <- diag(e)                                   
     G <- U%*%D%*%t(U)   			                        # compute pdf
     colnames(G) <- rownames(G) <- rn
     return(G)                      
   } 
    
  
     setMME <- function(X=NULL,Z=NULL,Gci=NULL,Ri=NULL,y=NULL,G=NULL, tol=tol) {
          Cxx <- t(X)%*%X*Ri
          Cxz <- t(X)%*%Z*Ri
          Czx <- t(Cxz)
          Czz <- t(Z)%*%Z*Ri + Gci
          #trT <- sum(diag(solve(Czz)))
          #print(trT)
          Xy <- t(X)%*%y*Ri
          Zy <- t(Z)%*%y*Ri
          Cx <- cbind(Cxx,Cxz)
          Cz <- cbind(Czx,Czz)
          C <- rbind(Cx,Cz)
          rhs <- rbind(Xy,Zy)
          iC <- solve(C)
          #iC <- ginv(C, tol=0.0001)$G
          iCxx <- iC[1:nrow(Cxx),1:nrow(Cxx)]
          iCxz <- iC[1:nrow(Cxx),(nrow(Cxx)+1):nrow(C)]
          iCzx <- t(iCxz)
          iCzz <- iC[(nrow(Cxx)+1):nrow(C),(nrow(Cxx)+1):nrow(C)]
          sol <- iC%*%rhs 
          sse  <- t(sol)%*%rhs
          bhat <- sol[1:nrow(Cxx)] 
          ghat <- sol[(nrow(Cxx)+1):nrow(C)] 
          #vbhat <- iCxx
          #trHb <- sum(rowSums(X*t(cbind(iCxx,iCxz)%*%t(cbind(X,Z))*Ri))) 
          #trHg <- sum(rowSums(Z*t(cbind(iCzx,iCzz)%*%t(cbind(X,Z))*Ri))) 
          trH <- sum(rowSums(cbind(X,Z)*t(iC%*%t(cbind(X,Z))*Ri))) 
          #print(c(trHb,trHg,trH))
          #vghat <- makepdf(G-iCzz)
          #vghat <- G-iCzz[1:ncol(G),1:ncol(G)]
          #vghat <- makepdf(G-iCzz[1:ncol(G),1:ncol(G)])
          vghat <- makepdf(G, tol=tol)-iCzz[1:ncol(G),1:ncol(G)]
          #vghat <- iCzz[1:ncol(G),1:ncol(G)]
          #pevg <- makepdf(iCzz)
          #return(list(iC=iC,sse=sse,sol=sol,vghat=vghat,pevg=pevg,trH=c(trH=trH,trHg=trHg,trHb=trHb)))
          return(list(iC=iC,sse=sse,sol=sol,vghat=vghat,trH=trH))
        }
    
    

    aireml <- function(G=NULL, y=NULL, X=NULL, Z=NULL, theta=NULL, restricted=NULL, maxit=100)
      {
        np <- length(theta)
        if (is.null(restricted)) restricted <- rep(FALSE,np)
      	n <- length(y)
      	if (is.null(Z)) { for ( i in 1:np){ Z <- c(Z, list(diag(n))) } }
        for ( i in 1:np){ G[[i]] <- Z[[i]]%*%G[[i]]%*%t(Z[[i]]) }
        AI <- matrix(0, ncol=np, nrow=np)
      	s <- matrix(0, ncol=1, nrow=np)
        tol <- 100
        i <- 0
        while ( max(tol)>10e-7 ) {
        V <- matrix(0,n,n)
    		i <- i + 1
        for ( j in 1:np) { V <- V + G[[j]]*theta[j] }
    		Vinv <- solve(V)
    		P <- Vinv - Vinv%*%X%*%solve(t(X)%*%Vinv%*%X)%*%t(X)%*%Vinv
        for ( j in 1:np) {
         for ( k in j:np) {
        	AI[j,k] <- t(y)%*%P%*%G[[j]]%*%P%*%G[[k]]%*%P%*%y
      		AI[k,j] <- AI[j,k]
        }
          s[j,1] <- sum(diag(G[[j]]%*%P))-t(y)%*%P%*%G[[j]]%*%P%*%y
        }
        theta0 <- theta
    		theta0[!restricted] <- theta[!restricted] - solve(AI[!restricted,!restricted])%*%s[!restricted]
    		isNeg <- theta0<0
        if (any(isNeg)) {
          fixed <- !restricted
          fixed[isNeg] <- FALSE
          theta0[fixed] <- theta[fixed] - solve(AI[fixed,fixed])%*%s[fixed]
          theta0[isNeg] <- 0.000000001
        }
        tol <- abs(theta - theta0)
    		theta <- theta0
    	  print(paste(c("Iteration:",i,"Theta:",round(theta,5)), sep=""))
        if (i==maxit) break
        }

        for ( j in 1:np) { V <- V + G[[j]]*theta[j] }
	      eigV <- eigen(V)
	      ldV <- sum(log(eigV$values))
    		Vinv <- solve(V)
        eigXVX <- eigen(t(X)%*%Vinv%*%X)
	      ldXVX<-sum(log(eigXVX$values))
    		P <- Vinv - Vinv%*%X%*%solve(t(X)%*%Vinv%*%X)%*%t(X)%*%Vinv
	      llik <- -0.5*(ldV+ldXVX+t(y)%*%P%*%y)

      	return(list( "theta"=theta, llik=llik, niter=i   ))
      }

    ttest <- function(y=NULL,X=NULL) {
      # y is centered vector of phenotypes
      # X is a centered and scaled matrix of genotypes
      xy <- colSums(X*y)
      xx <- colSums(X**2) 
      b <- xy/xx
      yhat <- t(t(X)*b)
      sse <- colSums((y-yhat)**2)
      dfe <- nrow(X)-2 
      tt <- b*sqrt(dfe)/sqrt(sse/xx)
      ptt <- 2*pt(-abs(tt),df=dfe)
      return(list(b=b,t=tt,p=ptt))
    }

   hssvs <- function( y=NULL, X=NULL, set=NULL, p1=NULL, g0=NULL, nsamp=100, hgprior=list(sce0=0.001, scg0=0.001, dfe0=4, dfg0=4)) {
     
       n <- length(y)                           # number of observations
       p <- ncol(X)                             # number of regression variables
       dxx <- colSums(X**2)                     # diagonal elements of the X'X matrix
       p0 <- 1-p1                               # prior proportion of elements in gamma to be set to one

       b <- rep(0,p)                            # initialize b
       g <- rep(0,p)                            # intialize g
       
       nset <- length(set)                      # number of sets
       pset <- sapply(set,length)
       for ( i in 1:nset) {
         n1 <- as.integer(pset[i]*p1) 
         cset <- abs(cor(y, X[,set[[i]]]))
         o <- order(cset,decreasing=TRUE)[1:n1]
         o <- set[[i]][o]
         g[o] <- 1
       }

       bset <- NULL                             # beta's within sets
       gset <- NULL                             # g's within sets
       for ( i in 1:nset ) {
         bset[[i]] <- b[set[[i]]]
         gset[[i]] <- g[set[[i]]]
       }

       mu <- 0                                  # initilaize mu
       sigma2 <- var(y)/2                       # initialize sigma2
       g1 <- 1                                  # initialize slab variance
       if (is.null(g0)) g0 <- 1e-09             # initialize spike variance

       g0 <- rep(g0,nset)                       # spike variance for each set
       g1 <- rep(g1,nset)                       # slab variance for each set
       
       #Xs <- matrix(0,nrow=n,ncol=nset)         # matrix (n*nset) of genomic values one for each set

       e <- as.vector(y)                        # initialize residuals
 
       for ( i in 1:nsamp ) {

         for ( j in 1:nset) {
           cls <- set[[j]]
           samp <- ssvs(e=e, X=X[,cls], b=bset[[j]], dxx=dxx[cls], mu=mu, g=gset[[j]], sigma2=sigma2, p0=p0, p1=p1, g0=g0[j], g1=g1[j], hgprior=hgprior)
             e <- samp$e
             bset[[j]] <- samp$b
             mu <- samp$mu
             sigma2 <- samp$sigma2
             gset[[j]] <- samp$g
             g0[j] <- samp$g0
             g1[j] <- samp$g1
            # Xs[,j] <- X[,cls]%*%bset[[j]]
         }
         #  h2 <- apply(Xs,2,var)
         #  h2 <- h2/var(y)
         #  barplot(h2)
           print("")
           print(c("round",i, "gset", sapply(gset,sum)))
         #  print(c("round",i, "mu, sigma2, cor(y,Xs)", c(mu,sigma2,cor(y,rowSums(Xs)))))
         #  print(c("round",i, "g0",g0 ))
         #  print(c("round",i, "g1",g1 ))
       }
  
  
     }
     
     ssvs <- function(e=e, X=X, b=b, dxx=dxx, mu=mu, g=g, sigma2=sigma2, p0=p0, p1=p1, g0=g0, g1=g1, hgprior=hgprior) {
  
         n <- length(e)                           # number of observations
         p <- ncol(X)                             # number of regression variables
  
         sce0 <- hgprior$sce0                     # prior residual sums of squares
         scg0 <- hgprior$scg0                     # prior slab residual sums of squares
         dfe0 <- hgprior$dfe0                     # prior residual degrees of freedom
         dfg0 <- hgprior$dfg0                     # prior slab residual degrees of freedom
  
     
       # sample b and update residuals 
         bC <- b
         shrinkage <- rep(sigma2/g1,times=p)
         shrinkage[g==0] <- sigma2/g0 
         xxs <- dxx + shrinkage
         for ( j in 1:p){
           rhs <- sum(X[,j]*e) + dxx[j]*bC[j]
           b[j] <- rnorm(1,mean=rhs/xxs[j],sd=sigma2/xxs[j])
           e  <- e - X[,j]*(b[j] - bC[j]) 
         }
     
       # sample g0 
         scg00 <- sum((b**2)[g==0]) + scg0
         dfg00 <- sum(1-g) + dfg0 - 2 
         g0 <- scg00/rchisq(n=1, df=dfg00, ncp=0) 
                       
       # sample g1 
         scg1 <- sum((b**2)[g==1]) + scg0
         dfg1 <- sum(g) + dfg0 - 2 
         for (k in 1:1000) {                      # scaling needs to be revised
          g1 <- scg1/rchisq(n=1, df=dfg1, ncp=0) 
          if(g1>g0) break                         # scaling needs to be revised
          if(k==1000) g1 <- g0                    # scaling needs to be revised
         }                                        # scaling needs to be revised
                       
                       
       # sample g 
         ratio1 <- p1*(exp(-(b**2)/g1)/sqrt(g1))
         ratio0 <- p0*(exp(-(b**2)/g0)/sqrt(g0))
         ratio <- ratio1/ratio0
         u <- 1/(1-runif(p,0,1))
         g[1:p] <- 0
         g[ratio>u] <- 1
     
  
       # sample mu and update residuals 
         muC <- mu
         mu <- rnorm(1,mean=mean(e+muC),sd=sigma2/n)
         e <- e - mu + muC 
  
       # sample sigma 
         sce <- sum(e**2) + sce0
         dfe <- n + dfe0 - 2 
         sigma2 <- sce/rchisq(n=1, df=dfe, ncp=0) 
         
         
       return(list(e=e, mu=mu,b=b,sigma2=sigma2,g=g, g0=g0, g1=g1))
       
     }

################################################################################
# some set based functions
################################################################################
  
    bluptt <- function(fit=NULL,n=NULL) {
      tv <- fit$b/sqrt(fit$vb)
      dfe <- fit$n-fit$trH[,1]
      p <- 2*pt(-abs(tv),df=dfe)
      return(list(value=tv,p=p))
    }

    hyperGSet <- function(p=NULL,sets=NULL,threshold=0.05) {
      N <- length(p)
      Na <- sum(p<threshold)
      Nna <- N-Na
      Nf <- sapply(sets,length)
      Naf <- sapply(sets, function(x) { sum(p[x]<threshold) })
      Nnaf <- Nf-Naf
      Nanf <- Na-Naf
      Nnanf <- Nna-Nnaf
      phyperg <- 1-phyper(Naf-1, Nf, N-Nf, Na)
      phyperg
    }

    blupGSet <- function(b=NULL,g=NULL, W=NULL, sets=NULL,nperm=100) {
     m <- length(m)
     nsets <- length(sets)
     Wb <- t(t(W)*b)
     setT <- matrix(0,nrow=nsets,ncol=4)
     psetTC <- matrix(0,nrow=nsets,ncol=3)
     psetTS <- matrix(0,nrow=nsets,ncol=3)
     for ( i in 1:nsets) {
      cls <- gsets[[i]]
      gS <- rowSums(Wb[,cls])
      nS <- length(cls)
      tstat <- c(cor(g,gS), sum((g*gS)**2)/nS, var(gS))
      setT[i,2:4] <- tstat
        tstatO <- tstat[-1]
        for ( j in 1:nperm) {
          n0 <- sample(1:m,1)
          cls0 <- (c(n0:m,1:(n0-1)))[cls]
          if(n0==1) cls0 <- cls 
          if(n0==m) cls0 <- cls 
          gS0 <- rowSums(Wb[,cls0])
          gS00 <- sample(gS0)
          tstatC <- c( cor(g,gS0), sum((g*gS0)**2)/nS, var(gS0))
          tstatS <- c( cor(g,gS00), sum((g*gS00)**2)/nS, var(gS00))
          psetTC[i,] <- psetTC[i,] + as.numeric(tstat<tstatC)
          psetTS[i,] <- psetTS[i,] + as.numeric(tstat<tstatS)
          print(c(i,j))
        }
      }
      colnames(setT) <- c("nSet","Corr(g,gS)","MSSggS","Var(gS)")
      colnames(psetTC) <- c("Corr(g,gS)","MSSggS","Var(gS)")
      colnames(psetTS) <- c("Corr(g,gS)","MSSggS","Var(gS)")
      rownames(setT) <- rownames(psetTS) <- rownames(psetTC) <- names(sets)
      return(list(setT=setT, pC=psetTC,pS=psetTS))
    }

    blupGFM <- function(b=NULL,g=NULL, W=NULL, sets=NULL) {
     m <- length(m)
     nsets <- length(sets)
     Wb <- t(t(W)*b)
     gsets <- matrix(0,nrow=nrow(W),ncol=nsets)
     for ( i in 1:nsets) {
      cls <- gsets[[i]]
      gsets[,i] <- rowSums(Wb[,cls])
      }
      colnames(gsets) <- names(sets)
      rownames(gsets) <- rownames(W)
      gsets
    }

    setTest <- function(stat=NULL, sets=NULL, nperm=NULL, method="sum") {
      if(method=="mean") setT <- sapply( sets,function(x){ mean(stat[x]) })
      if(method=="sum") setT <- sapply( sets,function(x){ sum(stat[x]) })
      if(method=="max") setT <- sapply( sets,function(x){ max(stat[x]) })
      if (!is.null(nperm)) {
       p <- rep(0,length(sets)) 
       n <- length(stat)
       nset <- sapply(sets,length)
       rws <- 1:n
       names(rws) <- names(stat)
       sets <- lapply(sets, function(x) { rws[x] })
       for ( i in 1:nperm ) {
        rws <- sample(1:n,1)
        o <- c(rws:n,1:(rws-1))
        pstat <- stat[o]
        if(method=="mean") setTP <- sapply(sets,function(x){ mean(pstat[x]) })
        if(method=="sum") setTP <- sapply(sets,function(x){ sum(pstat[x]) })
        if(method=="max") setTP <- sapply(sets,function(x){ max(pstat[x]) })
        p <- p + as.numeric(setT>setTP) 
      }
      p <- 1-p/nperm
      setT <- data.frame(setT,nset,p)
      }
      return(setT)
    }

  cvsTest <- function(b=NULL,g=NULL, W=NULL, sets=NULL, nperm=100){
    Wb <-   t(t(W)*as.vector(b))
    cvs <- colSums(as.vector(g)*Wb)
    setT <- setTest( stat=cvs, sets=sets, nperm=nperm, method="sum")$p
    names(setT) <- names(sets)
    #setT <- setTest( stat=cvs, sets=sets, nperm=nperm, method="sum")
    #rownames(setT) <- names(sets)
    return(setT)
  }

  scoreTest <- function(e=NULL,W=NULL, sets=NULL, nperm=100){
    we2 <- as.vector((t(W)%*%e)**2)   
    names(we2) <- colnames(W)       
    setT <- setTest( stat=we2, sets=sets, nperm=nperm, method="sum")$p
    return(setT)
  }

    blupGFM <- function(b=NULL,g=NULL, W=NULL, sets=NULL) {
     m <- length(m)
     nsets <- length(sets)
     Wb <- t(t(W)*b)
     gv <- matrix(0,nrow=nrow(W),ncol=nsets)
     for ( i in 1:nsets) {
      cls <- sets[[i]]
      gv[,i] <- rowSums(Wb[,cls])
      }
      colnames(gv) <- names(sets)
      rownames(gv) <- rownames(W)
      gv
    }

################################################################################
