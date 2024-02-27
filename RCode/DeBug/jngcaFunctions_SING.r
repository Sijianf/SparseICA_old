
#----------------------------------
# Benjamin Risk and Irina Gaynanova
# Contact: brisk@emory.edu
# Functions supporting SING
#------------------------------


# This function for easy parallelization needs to be revised to work better on other systems. Right now it involves my user-specific directory structure;
# some issues with passing variables to cluster
mlcaFP_parallelized <- function(xData, n.comp = ncol(xData), ncores=20, W.list = NULL, whiten = c('eigenvec','none'), restarts=NULL, distribution=c('tiltedgaussian','logistic','JB'), df=0, initFOBI = FALSE, scratch="~/temp",keepall = FALSE,sourcefile='~/Dropbox/JINGCA/Programs/Functions/jngcaFunctions.R',...) {
  warning("If running mlcaFP_parallelized, you need to manually change the path that contains jngca_functions in order to load the functions onto the workers. Sorry for the inconvenience. See comments in the program.")
  require(snow)
  require(JADE)
# require(parallel)  
#  if (is.null(ncores)) {
#    ncores = detectCores()
#  }
  system(paste('mkdir',scratch))
  
  cl = makeSOCKcluster(rep("localhost",ncores))
  
  if(is.null(restarts)) restarts=ncores
  
  source(sourcefile)
  #clusterExport(cl,"sourcefile")
  #clusterEvalQ(cl,sourcefile)
  
  # can't get this to work with path name variable!!!
  # Need to edit this manually to get it to work:
  invisible(clusterEvalQ(cl,source('~/Dropbox/JINGCA/Programs/Functions/jngcaFunctions.R')))
 
  invisible(clusterEvalQ(cl,.libPaths('~/Rlibs')))
  .libPaths('~/Rlibs')
 
  n.Voxels = nrow(xData)
  xData <- scale(xData, center=TRUE, scale=FALSE)
  whiten=match.arg(whiten)
  
  distribution = match.arg(distribution)
  # perform whitening once:
  if (whiten=='eigenvec') {
    zdata = svd(xData)$u*sqrt(n.Voxels) 
  } else {
    zdata = xData
  }
  
  n.Time = ncol(zdata)
  
  #Generate initializations from the principal subspace:
  W.temp = gen.inits(p=n.comp,d=n.comp,runs=ceiling(restarts/2)-1,orth.method='svd')
  zeros = matrix(0,n.Time-n.comp,n.comp)
  W.temp = lapply(W.temp,FUN = function(x) rbind(x,zeros)) #pad with zeros
  W.temp[[1]] = diag(n.Time)[,1:n.comp] #one initialization from PCA solution.
  
  #Generate initializations from the column space of the data:
  W.list=gen.inits(p=n.Time,d=n.comp,runs=ceiling(restarts/2),orth.method='svd')
  W.list = c(W.temp,W.list)
  
  #Generate an initialization from the FOBI solution:
  if (initFOBI) {
    require(JADE)
    estFOBI = FOBI(zdata)
    W.temp = list(t(estFOBI[['W']])[,1:n.comp])
    W.list = c(W.temp,W.list)
    rm(estFOBI)
  }
  
  newruns = length(W.list)
  mymlcaFP <- function(x,W.list,xData,n.comp,maxit=maxit,distribution,df) {
    est = mlcaFP(xData = xData, n.comp = n.comp, W.list=list(W.list[[x]]), whiten="none", maxit=300, distribution=distribution, df=df)
    
    save(est,file=paste0(scratch,'/est',x,'.RData'))
  }

  # estimate on whitened data:
  clusterApply(cl=cl,x=as.list(1:newruns),fun=mymlcaFP,W.list=W.list,xData=zdata,n.comp=n.comp,distribution=distribution,df=df)
  stopCluster(cl)
  
  ###################
  # 2. Find argmax
  ##Get list of estimates:
  setwd(scratch)
  list.est = NULL
  for (x in 1:newruns) {
    list.est = c(list.est,paste0("est",x,'.RData'))
  }

  loglik.v = numeric(newruns)
  
  ##Find best likelihood:
  for(i in 1:newruns) {
    load(list.est[i])
    loglik.v[i] = est[['loglik']]
  }

  load(list.est[which.max(loglik.v)])
  # # clean up intermediate files:
   if (!keepall) {
     for (i in 1:newruns) {
     tmpfile = paste0("est",i,".RData")
     file.remove(file.path(scratch,tmpfile))
     }
   }
  #   
 
  # equivalent to est.M.ols for since Xdata has been centered above
   Mhat <- t(est$S)%*%xData

  # 19 April 2019: changed named of Shat to S and Mhat to M
  # 16 July 2019: Deleted ordering since this was done in mlcaFP
#  return(list(S=Shat,M=Mhat,loglik.maxS=sort(loglik.maxS,decreasing=TRUE),nongaussianity=nongaussianity,loglik.allinits=loglik.v))
  return(list(S=est$S,M=Mhat,Ws=est$Ws,nongaussianity=est$nongaussianity,loglik.allinits=loglik.v,whitener=NULL))
  
  }

  
  
# IGAY: edited on June 5th, 2019 to correct orderingof Ws  
# BRisk: edited on 17 July 2019 to estimate M for xData, not for whitened data
mlcaFP <- function(xData, n.comp = ncol(xData), W.list = NULL, whiten = c('eigenvec','sqrtprec','none'), maxit = 1000, eps = 1e-06, verbose = FALSE, restarts.pbyd = 0, restarts.dbyd = 0, distribution=c('tiltedgaussian','logistic','JB'), density=FALSE, out.all=FALSE, orth.method=c('svd','givens'), max.comp = FALSE, reinit.max.comp=FALSE, df=0, irlba=FALSE,...) {

    #note: small changes from mlcaFP from the JASA paper: 
      # 1) output Mhat. 
      # 2) order by skewness, with option for n.comp=1
  
    #former option:
    #alg.typ = c('symmetric','deflation'),
    #alg.typ = match.arg(alg.typ)
    alg.typ = 'symmetric'

  
    distribution = match.arg(distribution)
    whiten=match.arg(whiten)

    if(restarts.dbyd>0 && whiten!='eigenvec') stop('Use whiten=eigenvec with restarts.dbyd')
    ## whiten:
    
    if(irlba) require(irlba)
    if(max.comp) { #if statement evaluates to true for all max.comp!=0
      s.comp = n.comp
      n.comp = max.comp
    }
    if(max.comp=='TRUE') stop('max.comp should be an integer or FALSE')
    if(reinit.max.comp && max.comp==FALSE) stop('Can not reinitialize from max.comp solution if max.comp==FALSE')
    if(reinit.max.comp && alg.typ=='deflation') stop('reinit.max.comp not yet written for deflation algorithm')
    
    #require(multidcov)
    if(distribution=='tiltedgaussian') {
      Gfunc = tiltedgaussian
      require(ProDenICA)  
    }
    if(distribution=='tiltedgaussian' && df==0) stop('df must be greater than 0 for tiltedgaussian')
    if(distribution=='logistic'  && df>0) stop('df should be set to zero when using logistic')
    if(distribution=='logistic') Gfunc = logistic
    if(distribution=='JB'  && df>0) stop('df should be set to zero when using JB')
    if(distribution=='JB') Gfunc = jb.stat
    if(!is.null(W.list) & class(W.list)!='list') stop('W.list must be a list')
    if(length(W.list) && (restarts.pbyd || restarts.dbyd)) stop('restarts.pbyd and restarts.dbyd must be equal to zero when supplying W.list')

    orth.method= match.arg(orth.method)
    p = ncol(xData)
    nRow = nrow(xData)
    d = n.comp
    
    # center xData such that ones%*%xData = 0
    xData <- scale(xData, center=TRUE, scale=FALSE)
    if (d > p) stop('d must be less than or equal to p')
    if (whiten=='eigenvec') {
      # Use whitener=='eigenvec' so that restarts.dbyd initiates from the
      # span of the first d eigenvectors.
      temp = whitener(X = xData,n.comp = p,irlba=irlba)
      zData = temp$Z
      whitener = temp$whitener
      rm(temp)
      }  else if (whiten=='sqrtprec') {
         est.sigma = cov(xData)  ## Use eigenvalue decomposition rather than SVD.
         evd.sigma = svd(est.sigma)
         whitener = evd.sigma$u%*%diag(evd.sigma$d^(-1/2))%*%t(evd.sigma$u)
         zData = xData%*%whitener
        }
    else {
      # no whitening:
      zData = xData
      whitener = diag(p)
    }
    # warning('TO DO: Look at whitening methods and check for inconsistent options')
  if (is.null(W.list)) {
    if(restarts.pbyd) W.list = gen.inits(p=p,d=d,runs=restarts.pbyd,orth.method=orth.method)
    if(restarts.dbyd) {
      W.temp = gen.inits(p=d,d=d,runs=restarts.dbyd,orth.method=orth.method)
      #pad with zeros:
      zeros = matrix(0,p-d,d)
      W.temp = lapply(W.temp,FUN = function(x) rbind(x,zeros))
      W.list = c(W.list,W.temp)
    }
  }
  ## If restarts.pbyd and restarts.dbyd both equal zero:
  if (is.null(W.list)) W.list = gen.inits(p=p,d=d,runs=1,orth.method=orth.method)
  runs = length(W.list)
  out.list = NULL
  loglik.v = numeric(runs)
  for(k in 1:runs) {
    W0 = as.matrix(W.list[[k]])
    if(alg.typ == 'symmetric') {
      out.list[[k]] = lca.par(xData=zData,W0=W0,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=n.comp,df=df, ...)
    out.list[[k]]$df = df
    }
    if(alg.typ == 'deflation') {
      out.list[[k]] = lca.def(xData=zData,W0=W0,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=n.comp,df=df,...)
    out.list[[k]]$df = df
    }
    if(max.comp) {
      flist0 = list()
      for (j in 1:d) flist0[[j]] <- Gfunc(out.list[[k]]$S[, j], ...)
      loglik.S <- apply(sapply(flist0, "[[", "Gs"),2,sum)
      orderedLL = order(loglik.S,decreasing=TRUE)
      out.list[[k]]$S = out.list[[k]]$S[,orderedLL[1:s.comp]]
      out.list[[k]]$Ws = out.list[[k]]$Ws[,orderedLL[1:s.comp]]
      out.list[[k]]$loglik = sum(loglik.S[orderedLL[1:s.comp]])
      loglik.v[k] = out.list[[k]]$loglik
      } else {
      loglik.v[k] = out.list[[k]]$loglik
    }
  }
  for(i in 1:runs){
    out.list[[i]]$distribution=distribution
    out.list[[i]]$whitener = whitener
  }
  out = out.list[[which.max(loglik.v)]]
  
  maxS = out$S
  loglik.maxS = numeric(n.comp)
  
  #order by negentropy:
  densities.maxS = apply(maxS,2,Gfunc)
  
  # note: sum here leads to correct jb statistic
  # IGAY: this is the first place where ordering is applied to maxS (out$S)
  for(i in 1:n.comp) loglik.maxS[i] = sum(densities.maxS[[i]][['Gs']])
  o.maxS = maxS[,order(loglik.maxS,decreasing=TRUE)]
  rm(maxS)
  
  #force skewness to be positive:
  # added option for one component on 6 May 2019
  if(n.comp>1) {
    Shat = rightskew(S = o.maxS,order.skew = FALSE)
  } else {
    Shat = sign(mean(o.maxS^3))*o.maxS
  }  
  # equivalent to est.M.ols for since 0 has been centered above
  # BRisk note: This gives the mixing matrix for the singular vectors, not for the data
  # so it is M t(V) D^{-1}, 
  
  Mhat <- t(Shat)%*%xData
  
  # IGAY: this is where out is reassigned for S and M, but not Ws with right ordering
  out$S = Shat
  out$M = Mhat
  # IGAY: fix the ordering of Ws based on loglik as well, fix on June 5th, 2019
  out$Ws = out$Ws[, order(loglik.maxS,decreasing=TRUE)]
  
  out$nongaussianity = sort(loglik.maxS,decreasing = TRUE)
  
  
  
  if(reinit.max.comp) {
    out = lca.par(xData=zData,W0=out$Ws,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=s.comp,df=df, ...)
    out$df = df
    out.list[[k+1]] = out
  }
  if(out.all==TRUE) {
    out.list[[which.max(loglik.v)]] = out
    out.list 
  }  else {
    out
}
}
#---------------------------------

#-------------------------
gen.inits <- function(p,d,runs,orth.method=c('svd','givens')) {
  orth.method=match.arg(orth.method)
  W.list = list()
  for(i in 1:runs) {
    if(orth.method=='givens') {
      W.list[[i]] <- as.matrix(theta2W(runif(n=choose(p,2),min=0,max=2*pi)))[,1:d]
      } else {
        temp = matrix(rnorm(p*d),p,d)
        W.list[[i]] <- svd(temp)$u
        }
  }
  W.list
}

#-------------------------------------

#-----------------------
# logistic <- function(xData, scale=sqrt(3)/pi, df=0) {
#   #maximizes likelihood given s then calculates gradient w.r.t. w.hat
#   #df is not used
#   xData = as.vector(xData)
#   list(Gs = -xData/scale - log(scale) - 2*log(1+exp(-xData/scale)), gs = -1/scale + 2*exp(-xData/scale)/(scale*(1+exp(-xData/scale))), gps = (2*exp(-2*xData/scale) - 2*exp(-xData/scale)*(1+exp(-xData/scale))) / (scale^2*(1+exp(-xData/scale))^2))
# }

logistic <- function(xData, scale=sqrt(3)/pi, df=0) {
  #maximizes likelihood given s then calculates gradient w.r.t. w.hat
  #df is not used
  xData = as.vector(xData)
  list(Gs = -xData/scale - log(scale) - 2*log(1+exp(-xData/scale)), gs = -1/scale + 2*exp(-xData/scale)/(scale*(1+exp(-xData/scale))), gps = - 2*exp(-xData/scale) / (scale^2*(1+exp(-xData/scale))^2))
}

jb.stat <- function(x, df=0) {
  n <- length(x)
  s <- sum(x^3)
  k <- sum(x^4)
  Gs <- x^3 * s / n^2 + (x^4 * k / n^2 + 9 / n - 6 * x^4 / n) / 4
  gs <- 6 * x^2 * s / n^2 + (8 * x^3 * (k / n - 3) / n) / 4
  gps <- 6 * (3 * x^4 + 2 * x * s) / n^2 + (24 * x^2 * (k / n - 3) / n + 32 * x^6 / n^2) / 4
  list(Gs = Gs, gs = gs, gps = gps)
}
  
#-----------------------
orthogonalize = function (W) {
  ##For arbitrary W, this is equivalent to (WW^T)^{-1/2} W
  temp <- svd(W)
  tcrossprod(temp$u,temp$v)
}

#-------------
orthonormalize <- function(xk,X,k) {
  #Gram-Schmidt orthogonalization
  #assumes columns of X have length equal to 1
  if(k!=1) {
    t <- numeric(length(xk))
    for (u in 1:(k-1)) {
      a <- sum(xk * X[,u])
      t <- t + a * X[,u]
    }
    xk <- xk - t
    }
  xk / sqrt(crossprod(xk))
}

genWfromWs <- function(Ws) {
  d = ncol(Ws)
  p = nrow(Ws)
  tempW = cbind(Ws,diag(p)[,(d+1):p])
  for(k in (d+1):p) {
    oldWk = tempW[,k]
    tempWk = tempW[,k]
    for(j in 1:(k-1)) {
      tempWj = tempW[,j]
      tempWk = tempWk - tempWj * crossprod(tempWj,oldWk)/crossprod(tempWj,tempWj)
  }
  tempW[,k] = tempWk/sqrt(crossprod(tempWk))
  }
  tempW
}

temp.orthogonalize <- function(V,W) {
  #orthogonalizes the vector V to all columns in W
  #and returns cbind(W,orthV)
  oldWk = V
  tempWk = V
  tempW=cbind(W,V)
  k=ncol(W)+1
  for(j in 1:(k-1)) {
      tempWj = tempW[,j]
      tempWk = tempWk - tempWj * crossprod(tempWj,oldWk)/crossprod(tempWj)
    }
  tempW[,k] = tempWk/sqrt(crossprod(tempWk))
  tempW
}


#------------------------------------------------
# symmetric algorithm:
lca.par <- function(xData,W0,Gfunc,maxit,verbose,density,eps,n.comp,df,...) {
  W0 = as.matrix(W0)
  d = ncol(W0)
  if(n.comp!=d) stop('W0 needs to be p x d')
  p = ncol(xData)
  nRow = nrow(xData)
  s <- xData %*% W0
  flist <- as.list(1:d)
  ##  Densities of likelihood components:
  for (j in 1:d) flist[[j]] <- Gfunc(s[, j], df=df,...)
  flist0 <- flist
  crit0 <- mean(sapply(flist0, "[[", "Gs"))
  nit <- 0
  nw <- 10
  repeat {
    nit <- nit + 1
    gS <- sapply(flist0, "[[", "gs")
    gpS <- sapply(flist0, "[[", "gps")
    #t1 <- t(xData) %*% gS/nRow
    t1 <- crossprod(xData,gS)/nRow
    t2 <- apply(gpS, 2, mean)
    if(d>1) W1 <- t1 - W0%*%diag(t2) else W1 <- t1 - W0*t2
    W1 <- orthogonalize(W1)
    if(d>1) nw <- frobICA(t(W0), t(W1))^2 else nw <- mean((W0-W1)^2) #Uses a measure that works for non-square matrices -- MSE. The measure is defined for M so here we use transpose of W.
    W0 <- W1
    s <- xData %*% W0
    for (j in 1:d) flist0[[j]] <- Gfunc(s[, j], df=df, ...)
    crit0 <- mean(sapply(flist0, "[[", "Gs"))
    if (verbose) cat("Iter", nit, "G", crit0, "Delta", nw, "\n")
    if ((nit >= maxit)) {
      warning('Max iter reached')
      break
    }
    if (nw < eps) break
  }
  out = list(Ws = W0, loglik = d*nRow*crit0, S = s)
  if(density) out$density = lapply(flist0, "[[", "density")
  out
}

#--------------------------------------
myMixmat <-  function (p = 2) {
  a <- matrix(rnorm(p * p), p, p)
  sa <- svd(a)
  d <- sort(runif(p,min=1,max=10))
  mat <- sa$u %*% (sa$v * d)
  attr(mat, "condition") <- d[p]/d[1]
  mat
}

#------------------------
#---------------------------------------
standardizeM <- function(M) {
  #M is d x p
  diag(diag(M%*%t(M))^(-1/2))%*%M
}

#--------------------
tiltedgaussian = function (xData, df = 8, B = 100, ...) {
  #This function is based on ProDenICA::GPois by Trevor Hastie
  #NOTE: Assumes data are zero mean.
  require(gam)
  n <- length(xData)
  sd.x = sd(xData)
  rx <- c(min(xData)-0.1*sd.x, max(xData)+0.1*sd.x) 
  xg <- seq(from = rx[1], to = rx[2], length = B)
  gaps <- diff(rx)/(B - 1)
  xcuts <- c(rx[1] - gaps/2, xg[-B] + gaps/2, rx[2] + gaps/2)
  #NOTE: I use the response variable that corresponds to the LCA paper.
  #This differs from the GPois algorithm in ProDenICA
  ys <- as.vector(table(cut(xData, xcuts)))/(gaps*n)
  pois.fit <- suppressWarnings(gam(ys ~ s(xg, df)+offset(dnorm(xg,log=TRUE)), family = poisson, ...))
  Gs <- predict(pois.fit) #log tilt function predicted at grid locations (note: predict on gam object can not be used to obtain derivatives)
  # the gam object with the predict function can not be used directly to obtain the derivatives
  # of the smoothing spline.
  # Here, we refit another iteration of the IRWLS algorithm used in gam:
  # Note: working residuals = (y - mu0)/mu0
  # weights = mu0
  # fitted(pois.fit) = mu0
  # predict(pois.fit) = eta0 = log(mu0)
  sGs = Gs #+ log(sum(dnorm(xg))/sum(fitted(pois.fit)))
  z0 <- sGs + residuals(pois.fit, type='working')
  pois.refit <- smooth.spline(x=xg, y=z0, w=fitted(pois.fit),df=df) #obtain the log tilt function in an object that can be used to obtain derivatives
  Gs <- predict(pois.refit, xData, deriv = 0)$y
  gs <- predict(pois.refit, xData, deriv = 1)$y
  gps <- predict(pois.refit, xData, deriv = 2)$y
  fGs <- function(x) predict(pois.refit,x,deriv=0)$y
  fgs <- function(x) predict(pois.refit,x,deriv=1)$y
  fgps <- function(x) predict(pois.refit,x,deriv=2)$y
  list(Gs = Gs, gs = gs, gps = gps, fGs = fGs, fgs=fgs, fgps=fgps)
}
#---------------------------------------------

#----------------------------
# estimate mixing matrix from estimates of components:
est.M.ols <- function(sData,xData,intercept=TRUE) {
  if(intercept) coef(lm(xData~sData))[-1,] else coef(lm(xData~sData-1))
}
# NOTE: for centered X, equivalent to t(sData)%*%xData/(V-1)

#-----------------------------------
# order by likelihood
# option for positive skewness
order.likelihood <- function(S,positive.skew=TRUE,distribution=c('logistic','tiltedgaussian','logcosh','JB'),out.loglik=FALSE,...) {
  distribution = match.arg(distribution)
  nObs = nrow(S)
  d = ncol(S)
  if(distribution=='tiltedgaussian') Gfunc = tiltedgaussian
  if(distribution=='logistic') Gfunc = logistic
  if(distribution=='logcosh') Gfunc = ProDenICA::G1
  if(distribution=='JB') Gfunc = jb.stat
  if(positive.skew) {
    skewness <- function(x, n = nObs) (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
    skew = apply(S, 2, skewness)
    sign = -1 * (skew < 0) + 1 * (skew > 0)
    S = S %*% diag(sign)
  }
  flist0 = list()
  for (j in 1:d) flist0[[j]] <- Gfunc(S[, j], ...)
  loglik.S <- apply(sapply(flist0, "[[", "Gs"),2,sum)
  orderedLL = order(loglik.S,decreasing=TRUE)
  S = S[,orderedLL]
  if(out.loglik) return(list(S=S,loglik=sort(loglik.S,decreasing=TRUE))) else S
}

marginal.likelihoods <- function(S,distribution=c('logistic','tiltedgaussian','logcosh','GPois','JB'),...)
{
  distribution = match.arg(distribution)
  if(distribution=='tiltedgaussian') Gfunc = tiltedgaussian
  if(distribution=='logistic') Gfunc = logistic
  if(distribution=='logcosh') Gfunc = ProDenICA::G1
  if(distribution=='GPois') Gfunc = ProDenICA::GPois
  if(distribution=='JB') Gfunc = jb.stat
  d = ncol(S)
  flist0 = list()
  for (j in 1:d) flist0[[j]] <- Gfunc(S[, j], ...)
  apply(sapply(flist0, "[[", "Gs"),2,sum)
}



#-------------------------
#Match based on L2 distances and Hungarian
matchICA=function(S,template,M=NULL) {
  require(clue)
  n.comp=ncol(S)
  n.obs=nrow(S)
  if(n.comp>n.obs) warning('Input should be n x d')
  if(!all(dim(template)==dim(S))) warning('Template should be n x d')
  S = t(S)
  template = t(template)
  l2.mat1=matrix(NA,nrow=n.comp,ncol=n.comp)
  l2.mat2=l2.mat1
  for (j in 1:n.comp) {
    for (i in 1:n.comp) {
      l2.mat1[i,j]=sum((template[i,]-S[j,])^2)/n.obs
      l2.mat2[i,j]=sum((template[i,]+S[j,])^2)/n.obs
    }
  }
  l2.mat1=sqrt(l2.mat1)
  l2.mat2=sqrt(l2.mat2)
  l2.mat=l2.mat1*(l2.mat1<=l2.mat2)+l2.mat2*(l2.mat2<l2.mat1)
  map=as.vector(solve_LSAP(l2.mat))
  l2.1=diag(l2.mat1[,map])
  l2.2=diag(l2.mat2[,map])
  sign.change=-1*(l2.2<l2.1)+1*(l2.1<=l2.2)
  perm=diag(n.comp)[,map]%*%diag(sign.change)

  s.perm=t(perm)%*%S
  if(!is.null(M)) {
    M.perm=t(M)%*%perm
    return(list(S=t(s.perm),M=t(M.perm)))
  }  else {
    t(s.perm)
  }
}

# IGAY: edited on June 4th, 2019 to return perm matrix and omangles ordering
# IGAY: edited on July 15th, 2019 to avoid exiting with error when rx>ry
# IGAY: edited on Aug 12th, 2019 to avoid dropping to vector type from matrix when original rank is 1
angleMatchICA=function(Mx,My,Sx=NULL,Sy=NULL) {
  # match the colums of Mx and My, using the 
  # n x p parameterization of the JIN decomposition
  # assumes 
  require(clue)
  rx = ncol(Mx)
  ry = ncol(My)
  #if(rx>ry) stop('rx must be less than or equal to ry')
  n.comp = max(ry, rx) # IGAY: adjusted this on Jyly 19th, 2019 to be maximal of the two
  
  n.obs=nrow(Mx)
  if(n.comp>=n.obs) warning('Input should be n x r')
  if(n.obs!=nrow(My)) warning('Mx and My need to have the same number of rows')
  Mx = t(Mx)
  Mx = Mx / sqrt(apply(Mx^2,1,sum))
  
  
  My = t(My)
  My = My / sqrt(apply(My^2,1,sum))
  
  if(rx<ry) {
    Mx = rbind(Mx,matrix(0,ry-rx,n.obs))
  }
  if(ry<rx) {
    My = rbind(My,matrix(0,rx-ry,n.obs))
  }
  angle.mat1=acos(Mx%*%t(My))
  angle.mat2=acos(-1*Mx%*%t(My))
  angle.mat=angle.mat1*(angle.mat1<=angle.mat2)+angle.mat2*(angle.mat2<angle.mat1)
  map=as.vector(solve_LSAP(angle.mat))
  angle.mat1.perm = angle.mat1[,map]
  angle.mat2.perm = angle.mat2[,map]
  angle1=diag(angle.mat1.perm)
  angle2=diag(angle.mat2.perm)
  matchedangles = apply(cbind(angle1,angle2),1,min)
  allangles = angle.mat1.perm*(angle.mat1.perm<=angle.mat2.perm)+angle.mat2.perm*(angle.mat2.perm<angle.mat1.perm)
  sign.change=-1*(angle2<angle1)+1*(angle1<=angle2)
  perm=diag(n.comp)[,map]%*%diag(sign.change)
  
  My.perm=t(perm)%*%My

  # reorder components by their matched angles 
  smatchedangles = sort(matchedangles)
  omangles = order(matchedangles)
  sallangles = allangles[omangles[1:rx],omangles[1:ry]]
  sMx = Mx[omangles[1:rx],, drop = F]
  sMy.perm = My.perm[omangles[1:ry],, drop = F]

  if(!is.null(Sy)) {
  Sy.perm=t(perm)%*%Sy
  sSy.perm = Sy.perm[omangles[1:ry],, drop = F]
  sSx = Sx[omangles[1:rx],, drop = F]
  return(list(Mx=t(sMx),My = t(sMy.perm),matchedangles = smatchedangles,allangles = sallangles,Sx = sSx, Sy = sSy.perm, perm = perm, omangles = omangles))
  }
  else {
    return(list(Mx=t(sMx),My = t(sMy.perm),matchedangles = smatchedangles,allangles = sallangles, perm = perm, omangles = omangles))
  }
}


#-------------------------------------
# Match mixing matrices:
# This function does not require M to be square:
frobICA<-function(M1=NULL,M2=NULL,S1=NULL,S2=NULL,standardize=FALSE) {
  #MODEL: X = S M + E, so M is d x p
  #standardize: if standardize==TRUE, then standardizes rows of M1 and M2
  #to have unit norm; if using S1 and S2, standardizes columns to have unit variance.
  #standardize=TRUE makes the measure scale invariant.

  require(clue)
  tfun = function(x) all(x==0)
  if(is.null(M1) && is.null(M2) && is.null(S1) && is.null(S2)) stop("need to supply either M1 and M2 or S1 and S2")
  if(!is.null(M1) && !is.null(M2) && !is.null(S1) && !is.null(S2)) {
    stop("provide either (M1 and M2) or (S1 and S2) but not both (M1,M2) and (S1,S2)")
  }
  if(!is.null(M1) && nrow(M1) > ncol(M1)) stop("The input appears to be S1 and S2, but the arguments were not specified; re-run with S1=<object> and S2=<object>")

  if(is.null(M1)) {
    nS = nrow(S1)
    if(nS!=nrow(S2)) stop('S1 and S2 must have the same number of rows')
    if(sum(apply(S1,2,tfun)) + sum(apply(S2,2,tfun))) stop('frobICA not defined when S1 or S2 has a column of all zeros')
    if(standardize) {
      S1 = scale(S1)
      S2 = scale(S2)
    }
    p = ncol(S1)
    q = ncol(S2)
    if(p < q) {
      S1 = cbind(S1,matrix(0,nS,(q-p)))
    }
    if(q < p) {
      S2 = cbind(S2,matrix(0,nS,(p-q)))
    }
    Stemp = matchICA(S=S1,template=S2)
    n.comp = max(q,p)
    indices = c(1:n.comp)[!(apply(Stemp,2,tfun) | apply(S2,2,tfun))]
    return(sqrt(sum((Stemp[,indices] - S2[,indices])^2))/sqrt(nS*min(p,q)))
  }

  else {
    if(sum(apply(M1,1,tfun)) + sum(apply(M2,1,tfun))) stop('frobICA not defined when M1 or M2 has a row of all zeros')
    if(standardize) {
      temp = diag((diag(M1%*%t(M1)))^(-1/2))
      M1 = temp%*%M1
      temp = diag((diag(M2%*%t(M2)))^(-1/2))
      M2 = temp%*%M2
    }
    p = ncol(M1)
    if(p!=ncol(M2)) stop("M1 and M2 must have the same number of columns")
    d = nrow(M1)
    q = nrow(M2)
    n.comp=max(d,q)
    if(n.comp > p) warning("M should be d x p")
    if(d<q) {
      M1 = rbind(M1,matrix(0,(q-d),p))
    }
    if(q<d) {
      M2 = rbind(M2,matrix(0,(d-q),p))
    }
    l2.mat1=l2.mat2=matrix(NA,nrow=n.comp,ncol=n.comp)
    for (j in 1:n.comp) {
      for (i in 1:n.comp) {
        #since signs are arbitrary, take min of plus and minus:
        l2.mat1[i,j]=sum((M2[i,]-M1[j,])^2)
        l2.mat2[i,j]=sum((M2[i,]+M1[j,])^2)
      }
    }
    l2.mat1=sqrt(l2.mat1)
    l2.mat2=sqrt(l2.mat2)
    #take the min of plus/min l2 distances. This is okay because solve_LSAP is one to one
    l2.mat=l2.mat1*(l2.mat1<=l2.mat2)+l2.mat2*(l2.mat2<l2.mat1)
    map=as.vector(solve_LSAP(l2.mat))
    #retain relevant l2 distances:
    l2.1=diag(l2.mat1[,map])
    l2.2=diag(l2.mat2[,map])
    #sign.change is for re-ordered matrix 2
    sign.change=-1*(l2.2<l2.1)+1*(l2.1<=l2.2)
    perm=diag(n.comp)[,map]%*%diag(sign.change)
    M.perm=t(perm)%*%M1
    indices = c(1:n.comp)[!(apply(M.perm,1,tfun) | apply(M2,1,tfun))]
    return(sqrt(sum((M.perm[indices,]-M2[indices,])^2))/sqrt(p*min(d,q)))
  }
}


#----------------
#----------------------------------------
#Function to make most extreme values for the skewed tail positive, i.e., force all distributions to be right skewed, and order ICs by skewness.
rightskew=function(S,M=NULL,order.skew=TRUE) {
  #S: n x d matrix
  #A: d x d` corresponding to X = S A
  #If order = TRUE, then ICs are organized from HIGHEST to LOWEST skewness where skewness is forced to be positive for all ICs.
  nObs <- nrow(S)
  if(ncol(S)>nObs) stop('S must be n x d')
  skewness<-function(x,n = nObs) (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  skew=apply(S,2,skewness)
  sign=-1*(skew<0)+1*(skew>0)
  S.new=S%*%diag(sign)
  if(!is.null(M)) M.new=t(diag(sign))%*%M
  if(order.skew==TRUE) {
    skew.new=apply(S.new,2,skewness)
    perm=order(skew.new,decreasing=TRUE)
    S.new=S.new[,perm]
    if(!is.null(M)) M.new=M.new[perm,]
  }
  if(is.null(M)) {
    S.new
  } else
    return(list(S=S.new,M=M.new))
}
#---------------------------------------------

mProDenICA <- function(X, n.comp = ncol(X), restarts=0, tol=1e-07,maxit=100,G = c('GPois','G0','G1'),verbose=FALSE,whiten=FALSE,...) {
  ##NOTE: the restarts in ProDenICA evaluate the likelihood at a sample of orthogonal matrices, identifies the random matrix associated with highest likelihood, and then estimates ICs for this single initialization. Here, I initiate from the entire set of random matrices.
  ##NOTE: Restarts defined differently here than in ProDenICA. ProDenICA is initiatialized from restarts+1 initial values.
  ##NOTE: G defined differently from ProDenICA's Gfunc; here it is a string
  require(ProDenICA)
  G = match.arg(G)
  if(G=='G0') Gfunc=G0
  if(G=='G1') Gfunc=G1
  if(G=='GPois') Gfunc=GPois
  est.list = list()
  runs = restarts+1
  obj.v = numeric(runs)
  theta.list = lapply(rep(choose(n.comp, 2), runs), runif, min = 0, max = 2 * pi)
  W.list = lapply(theta.list, theta2W)
  if(whiten) {
    a<- whitener(X=X,n.comp=n.comp)
    zData <- a$Z
  } else {
    zData <- X[,1:n.comp]
  }

  for(i in 1:runs) {
    est.list[[i]] = ProDenICA(x=zData, k=n.comp, W0=W.list[[i]], whiten=FALSE, maxit = maxit, thresh = tol, trace=verbose, restarts=0, Gfunc=Gfunc,...)
    if (G=='G1') obj.v[i] = calc.negent.hyvarinen(s=est.list[[i]]$s) else obj.v[i] = est.list[[i]]$negentropy
  }
  out = est.list[[which.max(obj.v)]]
  if(G=='G1') out$negentropy=obj.v[which.max(obj.v)]
  if(whiten) out$whitener=a$whitener
  out
}

theta2W = function(theta)
{
  #<author hidden>
  # For a vector of angles theta, returns W, a d x d Givens rotation matrix:
  # W = Q.1,d %*% ... %*% Q.d-1,d %*% Q.1,d-1 %*% ... %*% Q.1,3 %*% Q.2,3 %*% Q.1,2
  ##  if(theta < 0  || pi < theta){stop("theta must be in the interval [0,pi]")}
  d = (sqrt(8*length(theta)+1)+1)/2
  if(d - floor(d) != 0){stop("theta must have length: d(d-1)/2")}
  W = diag(d)
  index = 1
  for(j in 1:(d-1)){
    for(i in (j+1):d){
      Q.ij = givens.rotation(theta[index], d, c(i,j))
      W = Q.ij %*% W
      index = index + 1
    }
  }
  W
}

givens.rotation <- function(theta=0, d=2, which=c(1,2))
{
  # For a given angle theta, returns a d x d Givens rotation matrix
  #
  # Ex: for i < j , d = 2:  (c -s)
  #                         (s  c)
  c = cos(theta)
  s = sin(theta)
  M = diag(d)
  a = which[1]
  b = which[2]
  M[a,a] =  c
  M[b,b] =  c
  M[a,b] = -s
  M[b,a] =  s
  M
}

# Whitening Function:
whitener <- function(X,n.comp=ncol(X),center.row=FALSE,irlba=FALSE) {
  require(MASS)
  #X must be n x d
  if(ncol(X)>nrow(X)) warning('X is whitened with respect to columns')
  #Corresponds to model of the form X.center = S A, where S are orthogonal with covariance = identity.
  x.center=scale(X,center=TRUE,scale=FALSE)
  if(center.row==TRUE) x.center = x.center - rowMeans(x.center)
  n.rep=dim(x.center)[1]
  if(irlba==FALSE) svd.x=svd(x.center,nu=n.comp,nv=n.comp)
  if(irlba==TRUE) {
    require(irlba)
    svd.x=irlba(x.center,nu=n.comp,nv=n.comp)
  }
  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)
  return(list(whitener=t(ginv(svd.x$v%*%diag(svd.x$d[1:n.comp])/sqrt(n.rep-1))),Z=sqrt(n.rep-1)*svd.x$u,mean=apply(X,2,mean)))
}



# Returns square root of the precision matrix for whitening:
covwhitener <- function(X,n.comp=ncol(X),center.row=FALSE) {
  require(MASS)
  #X must be n x d
  if(ncol(X)>nrow(X)) warning('X is whitened with respect to columns')
  #Creates model of the form X.center = S A, where S are orthogonal with covariance = identity.
  x.center=scale(X,center=TRUE,scale=FALSE)
  if(center.row==TRUE) x.center = x.center - rowMeans(x.center)
  n.rep=dim(x.center)[1]
  covmat = cov(x.center)
  evdcov = eigen(covmat,symmetric = TRUE)
  whitener = evdcov$vectors%*%diag(1/sqrt(evdcov$values))%*%t(evdcov$vectors)
  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)
  return(list(whitener=whitener,Z=x.center%*%whitener,mean=apply(X,2,mean)))
}



rtwonorm <- function(n, mean=c(0,5), sd=c(2/3,1), prob=0.5) {
  k <- rbinom(n,1,prob=prob)
  k*rnorm(n,mean[1],sd[1])+(1-k)*rnorm(n,mean[2],sd[2])
}

rmixnorm <- function(n, pars = list(mean=c(0,5), sd = c(2/3,1), prob=c(0.25,0.75))) {
  probs = pars[['prob']]
  means = pars[['mean']]
  sigmas = pars[['sd']]
  if(sum(probs)!=1) stop('Probabilities must sum to one')
  z = rmultinom(n=n, size=1, prob = probs)
  k = length(probs)
  #use rnorm recycling:
  x = rnorm(k*n,means,sigmas)
  dim(x) = c(k,n)
  apply(x*z,2,sum)
}

SimTwoNorms <- function(n.samples, distribution=c('mix-sub','mix-super'),snr,noisyICA=FALSE) {
  distribution = match.arg(distribution)
  if(distribution=='mix-sub') {
    mean=c(-1.7,1.7); sd=c(1,1); prob=0.75
  }
  if(distribution=='mix-super') {
    mean=c(0,5); sd=c(2/3,1); prob=0.95
  }
  sim.S <- rtwonorm(n=2*n.samples, mean=mean, sd=sd, prob=prob)
  dim(sim.S) <- c(n.samples,2)
  sim.M = myMixmat(5)
  sim.W = solve(sim.M)
  if(noisyICA) {
    sim.N <- matrix(rnorm(n=5*n.samples,mean=0,sd=1),nrow=n.samples,ncol=5)
  } else {
    sim.N <- matrix(rnorm(n=3*n.samples,mean=0,sd=1),nrow=n.samples,ncol=3)
  }

  sim.Ms = sim.M[1:2,]
  sim.Xs = sim.S%*%sim.Ms
  if(noisyICA) {
    sim.Mn = NULL
    sim.Xn <- sim.N
  } else {
    sim.Mn <- sim.M[3:5,]
    sim.Xn <- sim.N%*%sim.Mn
  }
  #alpha = 1/sqrt(mean(sim.Xs^2))
  alpha = 1/sd(as.vector(sim.Xs))
  sim.Xs = sim.Xs*alpha
  mean.S = apply(sim.S,2,mean)
  temp.S = scale(sim.S,center=TRUE,scale=FALSE)
  scalingMat = diag(apply(temp.S,2,sd))
  scaled.sim.S = temp.S%*%solve(scalingMat)
  scaled.sim.Ms = sqrt(snr)*alpha*scalingMat%*%sim.Ms
  #sim.Xn = sim.Xn/sqrt(mean(sim.Xn^2))
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn))
  sim.Xs = sqrt(snr)*sim.Xs #equivalent to scaled.sim.S%*%(alpha*sqrt(snr)*scalingMat%*%sim.Ms)+alpha*sqrt(snr)*matrix(mean.S,nrow=n.samples,ncol=2,byrow=TRUE)%*%sim.Ms
  #since we only recover scaled.sim.S, "true Ms" for, e.g., IFA, is defined as in scaled.sim.Ms
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)
  return(list(sim.S = sim.S, sim.N = sim.N, sim.Ms = sim.Ms, sim.Mn = sim.Mn, sim.X=sim.X, scaled.sim.S = scale(sim.S),scaled.sim.Ms = scaled.sim.Ms,scaled.sim.X = scale(sim.X), whitened.sim.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
}


SimLCA <- function(n.samples, distribution = c('logistic','t','gumbel'), nu.vector = c(sqrt(3)/pi,sqrt(3)/pi), snr, k = NULL, noisyICA=FALSE) {
  ##k = number of noise components
  require(MASS)
  if(is.null(k)) {
    if (noisyICA) k = 5 else k = 3
  }
  d = length(nu.vector)
  p = ifelse(noisyICA==TRUE,k,k+d)
  distribution=match.arg(distribution)
  sim.S = NULL
  if(distribution == 'logistic') for (i in 1:d) sim.S = cbind(sim.S,rlogis(n=n.samples,scale=nu.vector[i]))
  if(distribution == 't') for (i in 1:d) sim.S = cbind(sim.S,rt(n=n.samples,df=nu.vector[i]))
  if(distribution == 'gumbel') {
    require(evd)
    for (i in 1:d) sim.S = cbind(sim.S,rgumbel(n=n.samples,scale=nu.vector[i]))
  }
  sim.N = matrix(rnorm(n.samples*k,sd=1),nrow=n.samples)
  sim.M = myMixmat(p)
  sim.Ms = sim.M[1:d,]
  sim.Xs = sim.S%*%sim.Ms
  if(noisyICA) {
    sim.Mn = NULL
    sim.Xn = sim.N
  } else {
    sim.Mn = sim.M[(d+1):p,]
    sim.Xn = sim.N%*%sim.Mn
  }
  #alpha = 1/sqrt(mean(sim.Xs^2))
  alpha = 1/sd(as.vector(sim.Xs))
  sim.Xs = sim.Xs*alpha
  mean.S = apply(sim.S,2,mean)
  temp.S = scale(sim.S,center=TRUE,scale=FALSE)
  scalingMat = diag(apply(temp.S,2,sd))
  scaled.sim.S = temp.S%*%solve(scalingMat)
  scaled.sim.Ms = sqrt(snr)*alpha*scalingMat%*%sim.Ms
  #sim.Xn = sim.Xn/sqrt(mean(sim.Xn^2))
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn)) #eigenvalues sum to p
  sim.Xs = sqrt(snr)*sim.Xs #equivalent to scaled.sim.S%*%(alpha*sqrt(snr)*scalingMat%*%sim.Ms)+alpha*sqrt(snr)*matrix(mean.S,nrow=n.samples,ncol=2,byrow=TRUE)%*%sim.Ms
  #since we only recover scaled.sim.S, "true Ms" for, e.g., IFA, is defined as in scaled.sim.Ms
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)

  return(list(sim.S = sim.S, sim.N = sim.N, sim.Ms = sim.Ms, scaled.sim.Ms = scaled.sim.Ms, sim.Mn = sim.Mn, sim.X=sim.X, scaled.sim.S = scale(sim.S), scaled.sim.X = scale(sim.X), whitened.sim.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
}



################################################
SimFMRI123 = function(snr = 1, noisyICA=FALSE, nTR=50, nImages=1, phi=0.5, dim.data=c(33,33), var.inactive=0.0001) {
  ##ASSUME 1,000 samples
  require(neuRosim)
  require(steadyICA)
  m = nImages
  #Latent components are fixed for each simulation:
  x1 = rep(3,5)
  y1 = c(3:7)
  s1.coords = cbind(x1,y1)
  s1 = specifyregion(dim = dim.data, coord = s1.coords, form = "manual")
  s1[s1!=0] = seq(0.5,1,length=length(x1))
  x2 = c(8,8,8,9,10,9,10,10,10,9,8)
  y2 = c(15,14,13,13,13,15,15,16,17,17,17)
  
  s2.coords = cbind(c(x2,x2+7),c(y2,y2))
  s2 = specifyregion(dim=dim.data, coord = s2.coords, form = 'manual')
  s2[s2!=0] = seq(0.5,1,length=2*length(x2))
  
  x3 = c(13,14,15,15,15,14,13,15,15,14,13)
  y3 = c(19,19,19,20,21,21,21,22,23,23,23)
  
  s3.coords = cbind(c(x3,x3+7,x3+14),c(y3,y3,y3))
  s3 = specifyregion(dim=dim.data, coord = s3.coords, form = 'manual')
  s3[s3!=0] = seq(0.5,1,length=3*length(x3))
  
  sim.S = cbind(as.vector(s1),as.vector(s2),as.vector(s3))
  
  if(m>1) {
    t.sim.S = sim.S
    for(i in 1:(m-1)) t.sim.S = rbind(t.sim.S,sim.S)
    sim.S = t.sim.S
    rm(t.sim.S)
  }
  
  ## Add small amount of Gaussian noise to inactive voxels
  nInactive = sum(sim.S == 0)
  baseline = rnorm(nInactive,mean=0,sd=sqrt(var.inactive))
  sim.S[sim.S==0] = baseline
  
  ##For noise, simulate Gaussian random field. Unique for each simulation:
  if(noisyICA)  nscan = nTR else nscan = nTR-3
  sim.GRF = NULL
  for(k in 1:m) {
    t.sim.GRF <- spatialnoise(dim = dim.data, sigma=1, nscan = nscan, method = "gaussRF", FWHM = 6)
    dim(t.sim.GRF) <- c(prod(dim.data),nscan)
    sim.GRF = rbind(sim.GRF,t.sim.GRF)
  }
  
  ##Mixmat:
  #create timecourses for latent components:
  totaltime <- nTR
  nOnsets = 5+1
  onsets <- seq(from=1, to=totaltime, length=nOnsets)
  dur <- totaltime/10
  #s <- stimfunction(totaltime = totaltime, onsets = onsets, durations = dur, accuracy = 1)
  row1 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(1,3)]), durations = list(dur), effectsize = 1, TR = 1, conv = "gamma")
  row2 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,5)]), durations = list(dur), effectsize = 1, TR=1, conv='gamma')
  #NOTE: Time courses can not be identical.
  row3 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,4)]), durations=list(dur), effectsize=1, TR=1, conv='gamma')
  
  sim.Ms = matrix(c(row1,row2,row3),nrow=3,byrow=TRUE)
  sim.Xs = sim.S%*%sim.Ms
  
  if(noisyICA)  {
    sim.Mn = NULL
    sim.Xn = sim.GRF
    for(t in 2:nTR) sim.Xn[,t] = phi*sim.Xn[,t-1]+sim.Xn[,t]
  }  else {
    sim.Mn = matrix(rnorm(nscan*nTR,0,1),nrow=nscan,ncol=nTR)
    for(t in 2:nTR) sim.Mn[,t] = phi*sim.Mn[,t-1] + sim.Mn[,t]
    sim.Xn = sim.GRF%*%sim.Mn
  }
  #sim.Xs = sim.Xs/sqrt(mean(sim.Xs^2)) 
  #sim.Xn = sim.Xn/sqrt(mean(sim.Xn^2)) 
  sim.Xs = sim.Xs/sd(as.vector(sim.Xs)) #standardize so we can control SNR
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn)) 
  sim.Xs = sqrt(snr)*sim.Xs
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)
  
  if(noisyICA) { 
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.Xn, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  } else {
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.GRF, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  }
}
#--------------------------------------------

# function to create network matrices from vectorized lower diagonals:

vec2net = function(invector,make.diag=1) {
  #invector: choose(p,2) x 1, where p is the number of nodes
  nNode = (1 + sqrt(1+8*length(invector)))/2
  outNet = matrix(0,nNode,nNode)
  outNet[lower.tri(outNet)] = invector
  dim(outNet) = c(nNode,nNode)
  outNet = outNet + t(outNet)
  diag(outNet) = make.diag
  outNet
}


##############



whitener.evd = function(xData) {
  est.sigma = cov(xData)  ## Use eigenvalue decomposition rather than SVD.
  evd.sigma = svd(est.sigma)
  whitener = evd.sigma$u%*%diag(evd.sigma$d^(-1/2))%*%t(evd.sigma$u)
  list(whitener=whitener,zData = xData%*%whitener)
}

# here subjects are by rows, columns correspond to components
aveM = function(mjX,mjY) {
  
  mjX = t(t(mjX) / sqrt(apply(mjX^2,2,sum))) # each column has eucledean norm 1
  mjY = t(t(mjY) / sqrt(apply(mjY^2,2,sum)))
  n = nrow(mjX) # number of subjects
  rj = ncol(mjX) # number of components
  aveMj = matrix(0,n,rj)
  for (j in 1:rj) {
    # Need to take sign into account
    signXY = sign(sum(mjX[, j] * mjY[, j])) # sign
    temp = (mjX[,j] + mjY[,j] * signXY)/2
    aveMj[,j] = temp/sqrt(sum(temp^2))
  }
  aveMj
}


# BRisk: 30 March 2020: Note: uses Hungarian. This is now deprecated.
jin_match = function(invLx,invLy,Ux,Uy,nperm,multicore=0,alpha=0.01){
  warning('This function is now deprecated')
  # takes arguments as in updateUboth
  # invLx: inverse of whitener for X n x n
  # Ux: px x n 
  # Conducts permutation test for joint rank
  
  #multicore: an integer; if =0 or FALSE, then uses a single-core
  # else uses doParallel
  Mx = invLx%*%t(Ux)
  My = invLy%*%t(Uy)
  nSubject = ncol(Ux)
  px = nrow(Ux)
  py = nrow(Uy)
  if(px>py) stop('px must be less than or equal to py -- switch order of X and Y')
  
  #  if (corr) {
  # scaling is not necessary, done in angleMatchICA
  # however, centering is
  # to ensure Mx has zero row means, one needs to make sure the
  # data input to the ngca function is double centered
  
  #    Mx = t(scale(Mx))
  #    My = t(scale(My))
  #  }
  
  a = all(apply(Mx,2,mean)<1e-10) & all(apply(My,2,mean)<1e-10)
  if (!a) {
    stop('Columns of Mx are not mean 0 -- correlations and angles are not equivalent! re-fit ngca with double centering')
  }
  
  matchedResults = angleMatchICA(Mx,My,Ux,Uy)
  # trick to use for Ux and Uy: treat like Sx and Sy
  
  
  # equivalent to correlations:
  corrmatched = cos(matchedResults$matchedangles)
  
  # fisher z-transform: (eliminate variance dependence on size of correlation)
  #zcorrmatched = tanh(corrmatched) # does not impact results, monotonic transformation

  
  #if(!is.null(seed)) set.seed(seed)
  # get sample of max z-corr from random permutations:
  if (multicore>0) {
    require(doParallel)
    require(doRNG)
    registerDoParallel(multicore)
    corrperm = foreach (b = 1:nperm, .combine=rbind) %dorng% {
      new.My = My[sample(1:nSubject),]
      max(abs(cor(Mx,new.My)))
    }
  }
  
  pperm = NULL
  maxrj = min(px,py)
  for (j in 1:maxrj) pperm[j] = mean(corrmatched[j]<corrperm)
  
  
  rj = sum(pperm<alpha)

  Ujx = matchedResults$Sx[1:rj,]
  Ujy = matchedResults$Sy[1:rj,]
  # should be equivalent to matchedResults, but is a little off
  #cMjx = invLx%*%t(Ujx)
  #cMjy = invLy%*%t(Ujy)
  #diag(cor(cMjx,cMjy))
  
  #matchedResults$Mx[1:5,1:5]
  #cMjx[1:5,1:5]
  #a = matchedResults$Mx[,1:rj]
  #sum((cMjx-a)^2)
  
  #a = cor(matchedResults$Mx,matchedResults$My)
  #diag(a[1:10,1:10])
  # # correlations are the same
  
  
  Mjx = matchedResults$Mx[,1:rj]
  Mjy = matchedResults$My[,1:rj]
  
  
  return(list(Ujx=Ujx,Ujy=Ujy,invLx=invLx,invLy=invLy,Mjx=Mjx,Mjy=Mjy,fwer_alpha = quantile(corrperm,1-alpha)))
}

# Use greedy algorithm:
greedymatch = function(Mx,My,Ux,Uy) {
    # input:
    # Mx: n x px
    # My: n x py
    # For this matching to coincide with angle matching, the columns must have zero mean.
    # NOTE: returns permuted columns of Mx and My. This function does not do any scaling or sign flipping.
  
    # check the means are approximately zero:
    checkX = all(colMeans(Mx)<1e-05)
    checkY = all(colMeans(My)<1e-05)
    if(!checkX) warning('Columns of Mx do not appear to be centered -- correlations and angle match give different results')
    if(!checkY) warning('Columns of My do not appear to be centered -- correlations and angle match give different results')
    
    allCorr = abs(cor(Mx,My))
    px = ncol(Mx)
    py = ncol(My)
    allCorr = abs(cor(Mx,My))
    minpxpy = min(px,py)
    mapX = numeric(minpxpy)
    mapY = numeric(minpxpy)
    selcorr = numeric(minpxpy)
    for (t in 1:minpxpy) {
      selcorr[t] = max(allCorr)
      tindices=which(allCorr==selcorr[t],arr.ind=TRUE)
      mapX[t]=tindices[1]
      mapY[t]=tindices[2]
      allCorr[tindices[1],] = 0
      allCorr[,tindices[2]] = 0
    }
    
    # reorder Mx and My
    # grab all columns (applies when px ne py), non-matched columns
    # can still impact the selection of the penalty:
    notinX = c(1:px)[!c(1:px)%in%mapX]
    mapX = c(mapX,notinX)
    notinY = c(1:py)[!c(1:py)%in%mapY]
    mapY = c(mapY,notinY)
    reorderMx = Mx[,mapX]
    reorderMy = My[,mapY]
    reorderUx = Ux[mapX,]
    reorderUy = Uy[mapY,]
    return(list('Mx' = reorderMx, 'My' = reorderMy,'Ux' = reorderUx, 'Uy' = reorderUy, 'correlations' = selcorr, 'mapX' = mapX, 'mapY' = mapY))
}

# BRisk: 30 March 2020
permTestJointRank = function(MatchedMx,MatchedMy,nperm=1000,alpha=0.01,multicore=0) {
  # Mx: nSubject x px
  # multicore: if = 0 or 1, uses a single processor; multicore > 1 requires packages doParallel and doRNG
  nSubject = nrow(MatchedMx)
  nSubject2 = nrow(MatchedMy)
  px = ncol(MatchedMx)
  py = ncol(MatchedMy)
  if (nSubject!=nSubject2) {
    stop('Mx and My have different numbers of rows. The number of subjects must be equal in X and Y.')
  }
    if (multicore>1) {
    require(doParallel)
    require(doRNG)
    registerDoParallel(multicore)

    corrperm = foreach (b = 1:nperm, .combine=rbind) %dorng% {
      new.My = MatchedMy[sample(1:nSubject),]
      max(abs(cor(MatchedMx,new.My)))
    }
  } else {
    corrperm = numeric(nperm)
    for (k in 1:nperm) {
      new.My = MatchedMy[sample(1:nSubject),]
      corrperm[k] = max(abs(cor(MatchedMx,new.My)))
    }
  }
  
  pperm = NULL
  maxrj = min(px,py)
  corrmatched = numeric(maxrj)
  for (i in 1:maxrj) {
    corrmatched[i] = abs(cor(MatchedMx[,i],MatchedMy[,i]))
  }
  for (j in 1:maxrj) pperm[j] = mean(corrmatched[j]<corrperm)
  rj = sum(pperm<alpha)

  return(list(rj=rj,pvalues = pperm,fwer_alpha = quantile(corrperm,1-alpha)))
}


# BRisk 27 January 2020: added comments, did not change any calculations
computeRho <- function(JBvalX, JBvalY, Mxall, Myall, rjoint){
  # Compute the minimal value of rho so that the joint components have lower starting objective than individual
  # JBvalX, JBvalY - JB values for all components, first are assumed to be joint
  # Mxall, Myall - corresponding mixing matrices 
  # rjoint - number of joint components
  rho = 0
  nx = length(JBvalX)
  ny = length(JBvalY)
  if ((nx == rjoint)|(ny == rjoint)) {return(0)}
  
  Mxscaled = Mxall %*% diag(sqrt(1/colSums(Mxall^2)))
  Myscaled = Myall %*% diag(sqrt(1/colSums(Myall^2)))
  
  Mxy = crossprod(Mxscaled, Myscaled)
  
  chordalJ = rep(0, rjoint)
  for (k in 1:rjoint){
    chordalJ[k] = 2 - 2 * Mxy[k, k]^2
  }
  
  for (i in (rjoint + 1):nx){ # ith individual from X
    for (j in (rjoint + 1):ny){ # jth individual from Y
      
      ratio = rep(NA, rjoint) # rhs/lhs for each joint
      #print(paste('pairs',i,j))
      for (k in 1:rjoint){ #kth joint
        rhs = max(JBvalX[i] - JBvalX[k] + JBvalY[j] - JBvalY[k], 0) 
        # Brisk comment: 2 - 2*Mxy[i,j]^2 is chordal norm between individual components in X and Y
        lhs = 2 - 2 * Mxy[i, j]^2 - chordalJ[k]
        if (lhs <= 0){
          ratio[k] = NA
        }else{
          ratio[k] = rhs/lhs
        }
      }
      rho = max(rho, max(ratio, na.rm = T)) 
    }
  }
  return(rho)
}




#IGAY: adjust function so that the path is not hard-coded
plotNetwork = function(component,title='',qmin=0.005, qmax=0.995, path = '~/Dropbox/JINGCA/Data/community_affiliation_mmpplus.csv',make.diag=NA) {
  # component:
  # vectorized network of length choose(n,2)
  require(ggplot2)
  require(grid)
  require(scales)
  
  # load communities for plotting:
  mmp_modules = read.csv(path)
  mmp_order = order(mmp_modules$Community_Vector)
  
  #check community labels:
  #table(mmp_modules$Community_Label)
  #table(mmp_modules$Community_Label,mmp_modules$Community_Vector)
  
  labels = c('VI','SM','DS','VS','DM','CE','SC')
  coords = c(0,70.5,124.5,148.5,197.5,293.5,360.5)

  
  zmin = quantile(component,qmin)
  zmax = quantile(component,qmax)
  
  netmat = vec2net(component,make.diag)
  
  meltsub = create.graph.long(netmat,mmp_order)
  #g2 = ggplot(meltsub, aes(X1,X2,fill=value))+ geom_tile()+ scale_fill_gradient2(low = "blue",  high = "red",limits=c(zmin,zmax),oob=squish)+labs(title = paste0("Component ",component), x = "Node 1", y = "Node 2")+coord_cartesian(clip='off',xlim=c(-0,390))
  
  g2 = ggplot(meltsub, aes(X1,X2,fill=value))+ geom_tile()+ scale_fill_gradient2(low = "blue",  high = "red",limits=c(zmin,zmax),oob=squish)+labs(title = title, x = "Node 1", y = "Node 2")+coord_cartesian(clip='off',xlim=c(-0,390))
  
  for (i in 1:7) {
    if (i!=3) {
      g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+10), xmax = (coords[i]+10), ymin = -7, ymax = -7)
    } else{
      g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+1), xmax = (coords[i]+1), ymin = -7, ymax = -7)
    }
  }
  # which nodes are prominent:
  loadingsummary = apply(abs(netmat),1,sum,na.rm=TRUE)
  loadingsum2 = loadingsummary[mmp_order]
  
  Community = factor(mmp_modules$Community_Label)[mmp_order]
  
  g3 = qplot(c(1:379),loadingsum2,col=Community,size=I(3))+xlab('MMP Index')+ylab('L1 Norm of the Rows')
  
  return(list(netmatfig = g2, loadingsfig = g3, netmat=netmat, loadingsummary = loadingsummary))
}

############
# Function for plotting networks with ggplot
create.graph.long = function(gmatrix,sort_indices=NULL) {
  nnode = nrow(gmatrix)
  X1 = c(1:nnode)%x%rep(1,nnode)
  X2 =  rep(1,nnode)%x%c(1:nnode)
  if (!is.null(sort_indices)) {
    gmatrix = gmatrix[sort_indices,sort_indices]
  }
  value = as.vector(as.matrix(gmatrix))
  data.frame(X1,X2,value)
}


###############################################################
# Auxillary functions taken from other people
##############################################################

# Original function by Yunfeng Zhang - calculate the power of a matrix
"%^%" <- function(S, power){
  out = eigen(S)
  nonzero = which(out$values > 1e-8)
  out$vectors[,nonzero] %*% (out$values[nonzero]^power * t(out$vectors[,nonzero]))
}

###############################################################
# Irina's auxillary functions for updates
###############################################################
# Calculate subspace distance using two semi-orthogonal matrices of the same size
# Here U1, U2 are both r times n with orthonormal rows
subspace_dist <- function(U1, U2){
  # Use frobICA
  frobICA(U1, U2)
}

###############################################################
# Update functions
##############################################################

# Function that calculates T(U), JB gradient with respect to U - DONE, but need to test more
calculateT <- function(U, X, alpha = 0.8, adjusted = F){
  # Calculate u^{top}X for each u_l and Xj, this is UX
  UX = U %*% X # r times p
  p = ncol(X)
  # Calculate all individual components
  gamma = rowMeans(UX^3) # length r
  kappa = rowMeans(UX^4 - 3) # length r
  prod1 = tcrossprod(UX^2, X)/p # r by n
  prod2 = tcrossprod(UX^3, X)/p # r by n
  # TU must be r by n
  TU = 6*alpha*diag(gamma)%*%prod1 + 8*(1-alpha)*diag(kappa)%*%prod2
  if (adjusted){
    # Use gradient adjustment as in Virta et al
    TU = TU - 24 * (1-alpha) * diag(kappa)%*%U
  }
  return(TU)
}

# Funciton that calculates full gradient with respect to current U - one function for X and Y
calculateG <- function(U, DataW, invL, A, rho, alpha = 0.8, r0 = nrow(U)){
  TU <- t(calculateT(U, X = DataW, alpha = alpha, adjusted = F)) # this makes it n by r
  GU <- -TU
  invLU <- tcrossprod(invL, U)
  normsLU2 <- colSums(invLU^2)
  for (i in 1:r0){
    uLa <- as.numeric(crossprod(invLU[ , i], A[i, ]))
    GU[ , i] <- GU[ , i] - 4 * rho * uLa * (invL %*% A[i, ]/normsLU2[i] - uLa * invL %*% invLU[ , i, drop = F]/(normsLU2[i]^2)) 
  }
  return(GU)
}


# Curvlinear algorithm with r0 joint scores
updateUboth_r0 <- function(Ux, Uy, xData, yData, invLx, invLy, rho, tau = 0.01, alpha = 0.8, maxiter = 1000, tol = 1e-6, r0){
  
  tau1 = tau
  # Form standardized UX
  normLX = tcrossprod(invLx, Ux[1:r0, , drop = F])
  normLX = diag(sqrt(1/diag(crossprod(normLX))))%*% t(normLX)
  
  # Form standardized UY
  normLY = tcrossprod(invLy, Uy[1:r0, , drop = F])
  normLY = diag(sqrt(1/diag(crossprod(normLY))))%*% t(normLY)
  
  # Calculate objective value at current point
  obj1 = objectiveJoint(Ux, Uy, xData, yData, normLX, normLY, rho, alpha = alpha)
  print(obj1)
  
  error = 100
  iter = 0
  obj <- rep(0, maxiter)
  r = nrow(Ux)
  n = ncol(Ux)
  while ((error > tol)&(iter < maxiter)){
    iter <- iter + 1
    # Calculate gradient-type function at current point
    GUx <- calculateG(Ux, xData, invLx, A = normLY, rho, alpha = alpha, r0 = r0) # this is n by r
    GUy <- calculateG(Uy, yData, invLy, A = normLX, rho, alpha = alpha, r0 = r0) # this is n by r
    # Calculate skew-symmetric W
    Wx = GUx%*%Ux - crossprod(Ux, t(GUx))
    Wy = GUy%*%Uy - crossprod(Uy, t(GUy))
    
    # Generate new point Y
    Vtaux = tcrossprod(Ux, solve(diag(n) + tau*Wx/2, diag(n) - tau*Wx/2))
    Vtauy = tcrossprod(Uy, solve(diag(n) + tau*Wy/2, diag(n) - tau*Wy/2))
    # Form standardized UX
    normLXtau = tcrossprod(invLx, Vtaux[1:r0, ])
    normLXtau = diag(sqrt(1/diag(crossprod(normLXtau)))) %*% t(normLXtau)
    
    # Form standardized UY
    normLYtau = tcrossprod(invLy, Vtauy[1:r0,])
    normLYtau = diag(sqrt(1/diag(crossprod(normLYtau))))%*% t(normLYtau)
    
    # Check the value of objective at new Ytau
    obj2 = objectiveJoint(Vtaux, Vtauy, xData, yData, normLXtau, normLYtau, rho, alpha = alpha)
    # Make sure going along descent direction
    while ((obj2 > obj1)&(tau > 1e-14)){  
      # Overshoot
      tau = 0.8*tau
      # New Vtau
      Vtaux = tcrossprod(Ux, solve(diag(n) + tau*Wx/2, diag(n) - tau*Wx/2))
      Vtauy = tcrossprod(Uy, solve(diag(n) + tau*Wy/2, diag(n) - tau*Wy/2))
      
      # Form standardized UX
      normLXtau = tcrossprod(invLx, Vtaux[1:r0,])
      normLXtau = diag(sqrt(1/diag(crossprod(normLXtau)))) %*% t(normLXtau)
      
      # Form standardized UY
      normLYtau = tcrossprod(invLy, Vtauy[1:r0,])
      normLYtau = diag(sqrt(1/diag(crossprod(normLYtau)))) %*% t(normLYtau)
      
      # New objective
      obj2 = objectiveJoint(Vtaux, Vtauy, xData, yData, normLXtau, normLYtau, rho, alpha = alpha)
    }
    if (obj2 < obj1){
      obj1 = obj2
      obj[iter] = obj2
      print(paste("Objective function value", obj2))
      Unewx = Vtaux
      Unewy = Vtauy
      normLY = normLYtau
      normLX = normLXtau
    }else{
      obj[iter] = obj1
      print(paste("Objective function value", obj1))
      Unewx = Ux
      Unewy = Uy
    }
    # How large is the change in objective
    error <- subspace_dist(Ux, Unewx) + subspace_dist(Uy, Unewy)
    Ux <- Unewx
    Uy <- Unewy
  }
  return(list(Ux = Unewx, Uy=Unewy, tau = tau, error = error, iter = iter, obj = obj[1:iter]))
}

# One coordinate update at a time, can again use the same function for X and Y
# Do curvlinear update for now
updateUboth <- function(Ux, Uy, xData, yData, invLx, invLy, rho, tau = 0.001, alpha = 0.8, maxiter = 1000, tol = 1e-6, r0 = nrow(Ux)){
  
  # Check whether only joint, or both joint and individual
  if ((r0 < nrow(Ux))|(r0 < nrow(Uy))){
    return(updateUboth_r0(Ux, Uy, xData, yData, invLx, invLy, rho, tau = tau, alpha = alpha, maxiter = maxiter, tol = tol, r0 = r0))
  }
  
  # Only joint
  tau1 = tau
  
  # Form standardized UX
  normLX = tcrossprod(invLx, Ux)
  normLX = diag(sqrt(1/diag(crossprod(normLX))))%*% t(normLX)
  
  # Form standardized UY
  normLY = tcrossprod(invLy, Uy)
  normLY = diag(sqrt(1/diag(crossprod(normLY))))%*% t(normLY)
  
  
  # Calculate objective value at current point
  obj1 = objectiveJoint(Ux, Uy, xData, yData, normLX, normLY, rho, alpha = 0.8)
  print(paste("Objective function value", obj1))
  
  error = 100
  iter = 0
  obj <- rep(0, maxiter)
  r = nrow(Ux)
  n = ncol(Ux)
  while ((error > tol)&(iter <= maxiter)){
    iter <- iter + 1
    # Calculate gradient-type function at current point
    GUx <- calculateG(Ux, xData, invLx, A = normLY, rho, alpha = alpha, r0 = r) # this is n by r
    GUy <- calculateG(Uy, yData, invLy, A = normLX, rho, alpha = alpha, r0 = r) # this is n by r
    # Calculate skew-symmetric W
    Wx = GUx%*%Ux - crossprod(Ux, t(GUx))
    Wy = GUy%*%Uy - crossprod(Uy, t(GUy))
    
    # Generate new point Y, same tau for both
    Vtaux = tcrossprod(Ux, solve(diag(n) + tau*Wx/2, diag(n) - tau*Wx/2))
    Vtauy = tcrossprod(Uy, solve(diag(n) + tau*Wy/2, diag(n) - tau*Wy/2))
    # Form standardized UX
    normLXtau = tcrossprod(invLx, Vtaux)
    normLXtau = diag(sqrt(1/diag(crossprod(normLXtau))))%*% t(normLXtau)
    
    # Form standardized UY
    normLYtau = tcrossprod(invLy, Vtauy)
    normLYtau = diag(sqrt(1/diag(crossprod(normLYtau))))%*% t(normLYtau)
    
    # Check the value of objective at new Ytau
    obj2 = objectiveJoint(Vtaux, Vtauy, xData, yData, normLXtau, normLYtau, rho, alpha = alpha)
    # Make sure going along descent direction
    while ((obj2 > obj1)&(tau > 1e-14)){  
      #while ((obj2 > obj1)&(tau > 1e-7)){
      # Overshoot
      tau = 0.8*tau
      # New Vtau
      Vtaux = tcrossprod(Ux, solve(diag(n) + tau*Wx/2,diag(n)-tau*Wx/2))
      Vtauy = tcrossprod(Uy, solve(diag(n) + tau*Wy/2,diag(n)-tau*Wy/2))
      # Form standardized UX
      normLXtau = tcrossprod(invLx, Vtaux)
      normLXtau = diag(sqrt(1/diag(crossprod(normLXtau))))%*% t(normLXtau)
      
      # Form standardized UY
      normLYtau = tcrossprod(invLy, Vtauy)
      normLYtau = diag(sqrt(1/diag(crossprod(normLYtau))))%*% t(normLYtau)
      # New objective
      obj2 = objectiveJoint(Vtaux, Vtauy, xData, yData, normLXtau, normLYtau, rho, alpha = alpha)
      
    }
    # print(tau)
    if (obj2 < obj1){
      obj1 = obj2
      obj[iter] = obj2
      print(paste("Objective function value", obj2))
      Unewx = Vtaux
      Unewy = Vtauy
      normLY = normLYtau
      normLX = normLXtau
    }else{
      obj[iter] = obj1
      print(paste("Objective function value", obj1))
      Unewx = Ux
      Unewy = Uy
    }
    # How large is the change in objective
    error <- subspace_dist(Ux, Unewx) + subspace_dist(Uy, Unewy)
    Ux <- Unewx
    Uy <- Unewy
  }
  return(list(Ux = Unewx, Uy = Unewy, tau = tau, error = error, iter = iter, obj = obj[1:iter]))
}

#######################################################
# Objective functions - to monitor the convergence
####################################################
# Function that calculates JB
calculateJB <- function(U, X, alpha = 0.8){
  # Calculate u^{top}X for each u_l and Xj, this is UX
  UX = U %*% X # r times p
  p = ncol(X)
  # Calculate all individual components
  gamma = rowMeans(UX^3) # length r
  kappa = rowMeans(UX^4 - 3) # length r
  
  # TU must be r by n
  JB = sum(alpha*gamma^2 + (1-alpha)*kappa^2)
  
  return(JB)
}
# Function that calculates joint objective value
objectiveJoint <- function(Ux, Uy, X, Y, normLX, normLY, rho, alpha = 0.8){
  JBpartX = calculateJB(Uy, Y, alpha = alpha)
  JBpartY = calculateJB(Ux, X, alpha = alpha)
  Innerproduct = sum(rowSums(normLX*normLY)^2)
  return(-JBpartX - JBpartY- 2*rho*Innerproduct)
}

# Separate function for chordal distance of two vectors
chordal <- function(x, y){
  sum((tcrossprod(x)/sum(x^2) - tcrossprod(y)/sum(y^2))^2)
}

############
# BRisk: returns 4/5 of Jin's function
calculateJBofS <- function(S, alpha = 0.8){
  p = ncol(S)
  if (p < nrow(S)) warning("S should be r x p")
  # Calculate all individual components
  gamma = rowMeans(S^3) # length r
  kappa = rowMeans(S^4 - 3) # length r
  
  JB = sum(alpha*gamma^2 + (1-alpha)*kappa^2)
  
  return(JB)
}


signchange = function(S,M=NULL) {
  # S: 59,412 x r_J
  signskew = sign(apply(S,2,function(x) mean(x^3)))
  newS = t(signskew*t(S))
  if(!is.null(M)) newM = t(signskew*t(M))
  ifelse(is.null(M),return(newS),return(list(S=newS,M=newM)))
}

###############
###############
# BRISK: generateData_v2 alters individual component in second dataset to be a little more sparse, 
# which makes it more realistic. (The original scenario was a more pathological example where logis fail
# but JB succeeds...)
# IGAY: added centering to mj so already column-centered approximately
generateData_v2 <- function(nsubject = 48, snr = c(0.2, 0.2), vars = c(0.01,0.01)){
  # Generate mixing matrices
  n1 = round(nsubject/2)
  mj1 = c(rep( 1, n1), rep(-1, nsubject - n1)) + rnorm(nsubject)
  mj2 = c(rep(-1, n1), rep( 1, nsubject - n1)) + rnorm(nsubject)
  mj = cbind(mj1, mj2)
  # mj = mj - matrix(colMeans(mj), nsubject, 2, byrow = T)
  
  # Create X components:
  # grab the 1, 2, 3 components used in the LCA paper; snr doesn't matter here as just grab the components in S
  simData = SimFMRI123(var.inactive = vars[1]) #simulates LCA model with 3 LCs 
  # joint and individual signal components:
  px = nrow(simData$S)
  simS = scale(simData$S)
  
  # Create joint structure for X
  sjX = t(simS[,2:3])
  djX = mj%*%sjX
  
  # Create individual structure for X
  siX = t(simS[,1])
  n4 = round(nsubject/4)
  miX = c(rep(-1,n4),rep(1,n4),rep(-1,n4),rep(1,nsubject-3*n4))+rnorm(nsubject)
  # miX = miX - mean(miX)
  diX = miX%*%siX
  
  # Calculate Frobenius norm of the signal
  signalXF2 = sum((djX + diX)^2)
  
  # Generate noise
  nX = t(scale(matrix(rnorm((nsubject-3)*px),px)))
  mnX = matrix(rnorm((nsubject-3)*nsubject),nsubject)
  # mnX = mnX - matrix(colMeans(mnX), nsubject, nsubject - 3, byrow = T)
  dnX = mnX%*%nX
  
  # Adjust the noise with snr ratio
  # Wrt to Frobenius norm
  dnX = dnX * sqrt(signalXF2/(sum(dnX^2)*snr[1]))
  
  
  # Create data matrix X
  dX = djX + diX + dnX
  
  # Calculate R^2 values for X joint
  R2x = sum(djX^2)/sum(dX^2)
  
  
  # Create Y components:
  # components that will represent network communities:
  # use a block structure:
  temp1 = c(rep(1,10),numeric(90))
  temp1 = temp1%*%t(temp1)
  
  temp2 = c(numeric(10),rep(1,20),numeric(70))
  temp2 = temp2%*%t(temp2)
  
  # BRISK: edited this component to be more sparse
  temp3 = c(numeric(40),rep(1,30),numeric(30))
  temp3 = temp3%*%t(temp3)
  
  temp4 = c(numeric(80),rep(1,20))
  temp4 = temp4%*%t(temp4)
  
  # Add small noise within the block structure
  var.noise = vars[2]
  
  # Create joint structure for Y
  sjY = cbind(temp1[lower.tri(temp1)],temp2[lower.tri(temp2)])
  inactive = sum(sjY==0)
  sjY[sjY==0] = rnorm(inactive, mean = 0, sd = sqrt(var.noise))
  sjY = t(scale(sjY))
  scalemj = t(t(mj)*c(-5,2))
  djY = scalemj%*%sjY
  py = ncol(sjY)
  
  # Create individual structure for Y
  siY = cbind(temp3[lower.tri(temp3)],temp4[lower.tri(temp4)])
  inactive = sum(siY==0)
  siY[siY==0] = rnorm(inactive, mean = 0, sd = sqrt(var.noise))
  siY = t(scale(siY))
  n8 = round(nsubject/8)
  miY = cbind(c(rep(1,n8),rep(-1,n8),rep(1,n8),rep(-1,n8),rep(1,n8),rep(-1,n8),rep(1,n8),rep(-1,nsubject-7*n8))+rnorm(nsubject),c(rep(1,n1),rep(-1,nsubject-n1))+rnorm(nsubject))
  # miY = miY - matrix(colMeans(miY), nsubject, 2, byrow = T)
  diY = miY%*%siY
  
  # Calculate Frobenius norm of the signal
  signalYF2 = sum((djY + diY)^2)
  
  # Generate noise for Y
  nY = t(scale(matrix(rnorm((nsubject-4)*py),py)))
  mnY = matrix(rnorm((nsubject-4)*nsubject),nsubject)
  # mnY = mnY - matrix(colMeans(mnY), nsubject, nsubject - 4, byrow = T)
  dnY = mnY%*%nY
  
  # Adjust the noise with snr ratio
  # Wrt to Frobenius norm
  dnY = dnY * sqrt(signalYF2/(sum(dnY^2)*snr[2]))
  
  # Create data matrix Y
  dY = djY + diY + dnY
  
  # Calculate R^2 values for X joint
  R2y = sum(djY^2)/sum(dY^2)
  
  return(list(dX = dX, dY = dY, mj = mj, sjX = sjX, sjY = sjY, siX = siX, siY = siY, snr = snr, R2x = R2x, R2y = R2y))
}




# Wrapper functions for rank estimation using permutatation approach
#######################################################################

permmatRank_sequential_JB = function(xdata,maxcompdata=ncol(xdata),ncomplist,nperms,ninitialization_data=10,ninitialization_perm=5) {
  #xdata: p x n subjects 
  #maxcompdata: number of components to estimate from data. This should be the maximum number of possible non-Gaussian components. By default equal to n
  #ncomplist is a vector corresponding to the components to be tested for Gaussianity. 
  # In the simplest case, this, will be 1:maxcompdata 
  #for each ncomp in ncomplist, ncomp refers to a test of whether the ncomp$th$ component is Gaussian
  #the maximum of this should be less than or equal to maxcompdata
  #nperms: number of samples for each permutation test
  #ncores: number of cores to use via registerDoParallel
  
  # subfunctions used in this function:
  permmatRank = function(xdata, ncomp, nperms,
                         ninitialization_perm) {
    nsubjects = ncol(xdata)
    px = nrow(xdata)
    ngauss = nsubjects - ncomp + 1
    permmatJB = rep(0, nperms)
    for(k in 1:nperms){
      # sample n - r subjects
      tempX = xdata[ , sample(1:nsubjects, ngauss)]
      newX = matrix(0, px, ngauss)
      for (j in 1:ngauss) {
        # permute jth column
        pmat = sample(1:px)
        newX[,j] = tempX[pmat, j]
      } 
      # estimate non-gaussianity of the component
      permmatJB[k] = mlcaFP(newX, n.comp = 1, restarts.pbyd=ninitialization_perm,distribution='JB')$nongaussianity
    }
    permmatJB
  } 
  
  
  # Estimate model with components=maxcompdata :
  estX = mlcaFP(xdata,restarts.pbyd=round(ninitialization_data/2),restarts.dbyd=round(ninitialization_data/2),distribution='JB', n.comp=maxcompdata)
  
  # Construct p-values for reach component
  nc = length(ncomplist)
  permmatJB_bigmat = matrix(0, nrow = nperms, ncol = nc)
  pvalues = rep(NA, nc)
  
  for (i in 1:length(ncomplist)) {
    # Estimate residuals from the first r-1 components
    if (ncomplist[i] > 1){
      newxdata = lm(xdata~estX$S[,1:(ncomplist[i]-1)])$residuals
    }else{
      newxdata = xdata
    }
    # Find values of non-gaussianity from all permutations
    permmatJB_bigmat[,i] = permmatRank(newxdata, ncomplist[i], nperms, ninitialization_perm)
    # Calculate corresponding p-value
    pvalues[i] = mean(estX$nongaussianity[ncomplist[i]]<permmatJB_bigmat[,i])
  }
  colnames(permmatJB_bigmat) = ncomplist
  return(list(pvalues = pvalues,sample=permmatJB_bigmat))
}


permmatRank_joint = function(matchedResults, nperms = 100){
  
  # Calcualte correlations based on original Ms via angle-matching
  Mx = t(matchedResults$Mx)
  My = t(matchedResults$My)
  corrmatched = cos(matchedResults$matchedangles)
  n = ncol(Mx)
  
  # Calculate maximal correlations based on permutations
  corrperm = rep(NA, nperms)
  for(b in 1:nperms){
    # Resample the subjects in the 2nd mixing matrix
    new.My= My[,sample(1:n), drop = F]
    # Calculate maximal correlation
    corrperm[b] = max(abs(atanh(cor(t(Mx),t(new.My)))))
  }
  
  # Calculate p-values from permutations
  rx = nrow(Mx)
  ry = nrow(My)
  maxrj = min(rx,ry) #maximum possible number of joint components
  
  pperm = rep(NA, maxrj)
  for (j in 1:maxrj) {
    pperm[j] = mean(atanh(corrmatched[j]) < corrperm)
  }
  
  return(list(pvalues = pperm, corrperm = corrperm, corrmatched = corrmatched))
}



###############
###############

makecifti <- function(table, filename, template = '~/Dropbox/MyHCP/Data/100307/MNINonLinear/tfMRI_MOTOR_LR_Atlas.dtseries.nii', toolboxloc1 = '~/Dropbox/Applications2/hcp_ciftimatlab', toolboxloc2 = '~/Dropbox/mfunctions/robertoostenveld-cifti-matlab-27383b8',wbloc = '~/Applications2/workbench/bin_macosx64/wb_command',scratchloc = '~',matlabloc = '/Applications/MATLAB_R2019b.app/bin/matlab'){
  
  # Adapted from https://mandymejia.wordpress.com/2016/10/28/r-function-to-write-cifti-files/ 
  # table = Vxp matrix
  # table can be a vector if only a single CIFTI is needed
  # V = number of locations in both or one hemisphere
  # p = number of time points in .dtseries.nii
  # filename = vector length p, containing file names to be written (no extension)
  # template = name of existing CIFTI file that can be used as a template
  # toolboxloc = location and name of folder containing cifti-matlab or fieldtrip toolbox
  # scratchloc = location where CSV files will be written then deleted
  
  # this function writes a MATLAB script and executes it
  # must have MATLAB installed
  
  # Write table to CSV in scratch location
  fname <- file.path(scratchloc, 'tmp.csv')
  write.table(table, file=fname, row.names=FALSE, col.names=FALSE, sep=',')
  
  # Write MATLAB script
  line1 <- paste0("addpath '", toolboxloc1, "'")
  line2 <- paste0("addpath '", toolboxloc2, "'")
  line3 <- paste0("cifti = ciftiopen('",template,"','",wbloc,"');") 
  line4 <- paste0("data = csvread('",scratchloc,"/tmp.csv')",";")
  line5 <- "[nrow,ncol]=size(data);"
  
  # Zero out entires, allowing for files without subcortical voxels:
  line6 <- "cifti.cdata = NaN(91282,ncol);"
  line7 <- "cifti.cdata(1:nrow,:)=data;"
  line8 <- paste0("ciftisave(cifti,'",filename,"','",wbloc,"');")
  matlab_lines <- c(line1, line2, line3, line4, line5, line6, line7, line8)
  writeLines(matlab_lines, con=file.path(scratchloc, 'myscript.m'))
  
  system(paste0(matlabloc," -nodisplay -r \"run('~/myscript.m'); exit\""))
  
  file.remove(fname)
  file.remove(file.path(scratchloc, 'myscript.m'))
}





##############
makeciftimmp <- function(v360, filename, template = '~/Dropbox/MyHCP/Data/100307/MNINonLinear/tfMRI_MOTOR_LR_Atlas.dtseries.nii', toolboxloc1 = '~/Dropbox/Applications2/hcp_ciftimatlab', toolboxloc2 = '~/Dropbox/mfunctions/robertoostenveld-cifti-matlab-27383b8',wbloc = '~/Applications2/workbench/bin_macosx64/wb_command',scratchloc = '~',matlabloc = '/Applications/MATLAB_R2019b.app/bin/matlab',mmpatlasloc="~/Dropbox/HCP_rsfmri_mmp/Data/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii"){
  
  # This function modifies makecifti to create a cifti file from 360 values corresponding to the glasser mmp atlas
  # makecifti is adapted from https://mandymejia.wordpress.com/2016/10/28/r-function-to-write-cifti-files/ 
  
  # Write table to CSV in scratch location
  fname <- file.path(scratchloc, 'tmp.csv')
  write.table(v360, file=fname, row.names=FALSE, col.names=FALSE, sep=',')
  
  # Write MATLAB script
  line1 <- paste0("addpath '", toolboxloc1, "'")
  line2 <- paste0("addpath '", toolboxloc2, "'")
  line3 <- paste0("cifti = ciftiopen('",template,"','",wbloc,"');") 
  line4 = paste0("mmpatlas = ciftiopen('",mmpatlasloc,"','",wbloc,"');")
  line5 = paste0("ncii = size(cifti.cdata,1);")
  line6 = paste0("ncortex = 59412;")
  line7 <- paste0("data = csvread('",scratchloc,"/tmp.csv')",";")
  line8 <- "[~,ncol]=size(data);"
  
  # Zero out entires, allowing for files without subcortical voxels:
  line9 <- "cifti.cdata = NaN(91282,ncol);"
  # TO DO: Expand to subcortical
  line10 <- "mmpatlas_plus = [mmpatlas.cdata;zeros(ncii-ncortex,1)];"
  line11 <- "for r = 1:360" 
  line12 <- "indices = mmpatlas_plus==r;"
  line13 <- "cifti.cdata(indices,:) = repmat(data(r,:),[sum(indices),1]);"
  line14 <- "end"
  line15 <- paste0("ciftisave(cifti,'",filename,"','",wbloc,"');")
  matlab_lines <- c(line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15)
  writeLines(matlab_lines, con=file.path(scratchloc, 'myscript.m'))
  
  system(paste0(matlabloc," -nodisplay -r \"run('~/myscript.m'); exit\""))
  
  file.remove(fname)
  file.remove(file.path(scratchloc, 'myscript.m'))
}