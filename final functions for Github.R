
library(mvtnorm)

# helper functions
kernelGP <- function(X1, X2, xi=1.0, tau=1.0){
  
  X1mat <- matrix(rep(X1,length(X2)),nrow=length(X1))
  X2mat <- matrix(rep(X2,length(X1)),nrow=length(X1),byrow=TRUE)
  
  tau**2 * exp(- xi**2 * (X1mat - X2mat) * (X1mat - X2mat) / 2 )
  
}
plotGP <- function(grid,samples,ylim=NA,xlim=NA, main="main",xlab="x",ylab="y",lwd=1){
  if(is.na(ylim[1])){
    ylim <- c(min(samples),max(samples))
  }
  if(is.na(xlim[1])){
    xlim <- c(min(grid),max(grid))
  }
  plot(grid,samples[1,],type="l",ylim=ylim,xlim=xlim,main=main,xlab=xlab,ylab=ylab,lwd=lwd)
  if(nrow(samples) > 1){
    for(p in 2:nrow(samples)){
      lines(grid,samples[p,],lwd=lwd)
    }
  }
}
student <- function(x,mu=0,s=1,df=1){
  exp(
    lgamma((df+1)/2) - lgamma(df/2) -.5 * log(pi * s**2 * df) - (df+1)/2 *
      log(1 + 1/df * (x - mu)**2/s**2)
  )
}
halfstudent <- function(x,s=1,df=1){
  2 * exp(
    lgamma((df+1)/2) - lgamma(df/2) - .5 * log(pi * s**2 * df) - (df+1)/2 *
      log(1 + 1/df * x**2/s**2)
  )
}
loghalfstudent <- function(x,s=1,df=1,log=FALSE){
  # density of parameter for which the exponential has a half t distribution, i.e.,
  # x = log(y), so that y ~ half-t(s,df)
  ifelse(log==FALSE,
         exp(
           log(2) + lgamma((df+1)/2) - lgamma(df/2) - .5 * log(pi * s**2 * df) - (df+1)/2 *
             log(1 + 1/df * exp(2 * x)/s**2) + x
         ),
         log(2) + lgamma((df+1)/2) - lgamma(df/2) - .5 * log(pi * s**2 * df) - (df+1)/2 *
           log(1 + 1/df * exp(2 * x)/s**2) + x
  )
}
loghalfnormal <- function(x,s){
  -.5 * log(2*pi*s**2) - exp(2*x)/(2*s**2) + x
}

dFkern <- function(x,nu=1,delta=1,b=1){
  exp(
    (nu/2-1)*log(x) -(nu+delta)/2*log(1+x/b)
  )
}
logIG <- function(x,alpha,beta,log=FALSE){
  ifelse(log==FALSE,
         exp(
           alpha*log(beta) - lgamma(alpha) - x*alpha - beta / exp(x)
         ),
         alpha*log(beta) - lgamma(alpha) - x*alpha - beta / exp(x)
  )
}
logIGkern <- function(x,alpha,beta,log=FALSE){
  ifelse(log==FALSE,
         exp(
           - x*alpha - beta / exp(x)
         ),
         - x*alpha - beta / exp(x)
  )
}

# Computation of marginal likelihoods
GPPtest4 <- function(yobs,xobs,Zobs=NULL,xi_hyper=c(exp(-1),1),sigma2_hyper=c(0,0),
                     samsizeMCMC=2e3,samsizeIS=1e4){
  #yobs is a vector of observations of dependent variable
  #xobs is a vector of observations of the key predictor
  #NOTE. it is recommended to have xobs to have a mean of 0 (see paper). 
  #Zobs is a matrix with covariates (e.g., column of ones for an intercept)
  #xi has half Student prior with hyperparameters xi_hyper=c(scale=exp(-1),df=1)
  #sigma2_hyper=c(0,0) to place a Jeffreys prior on sigma2
  #samsizeMCMC is the sample size to approximate the posterior
  #samsizeIS is the sample size for the importance sampling estimate
  #the outcome is a list where the first element contains the posterior draws,
  #  the second element is the marginal likelihood under M1,
  #  the third element is the marginal likelihood under M0.
  
  # properties data
  n <- length(yobs)
  
  if(is.null(Zobs)){
    k <- 0
  }else{
    Zobs <- as.matrix(Zobs)
    k <- ncol(Zobs)
  }
  
  #initial values
  logxi <- log(xi_hyper[1])
  sigma2 <- 1
  eta <- rep(0,length=k)
  
  xDiffMmat <- (matrix(rep(xobs,each=n),nrow=n) - 
                  matrix(rep(xobs,each=n),nrow=n,byrow=T))**2
  xxobs <- xobs %*% t(xobs)
  sumxtx <- sum(xobs**2)
  g <- n
  
  RWsd <- 1
  #initial unstandardized kernel
  covm_xi <- exp(-exp(2*logxi)/2*xDiffMmat) * xxobs
  covm <- sigma2 * g / sumxtx * covm_xi + sigma2 * diag(n)
  if(k==0){
    mean1 <- matrix(rep(0,n),ncol=1)
  }else{
    mean1 <- Zobs%*%eta
  }
  logpost <- dmvnorm(yobs,mean=mean1,sigma=covm,log=TRUE)
  
  #start sampling
  samsize <- samsizeMCMC
  check1 <- 30 #check every 'check1' draws whether random walk sd needs to be increased/decreased
  drawsMat <- matrix(0,nrow=samsize,ncol=2)
  acceptMat <- matrix(0,samsize,1)
  colnames(drawsMat) <- c("logxi","sigma2") 
  colnames(acceptMat) <- "logxi"
  #drawsMat <- cbind(drawsMat,matrix(0,nrow=samsize,ncol=n))
  #colnames(drawsMat)[-(1:3)] <- paste0("beta",1:n)
  if(k > 0){
    drawsMat <- cbind(drawsMat,matrix(0,nrow=samsize,ncol=k))
    colnames(drawsMat)[-(1:(2))] <- paste0("eta",1:k)
  }
  #wer <- Sys.time()
  pb = txtProgressBar(min = 0, max = samsize, initial = 0)
  
  for(s in 1:samsize){
    
    # draw logxi
    logxi_can <- rnorm(1,mean=logxi,sd=RWsd[1])
    covm_xi_can <- exp( - exp(logxi_can)**2/2 * xDiffMmat) * xxobs
    covm_can <- sigma2 * g / sumxtx * covm_xi_can + sigma2 * diag(n)
    logpost_can <- dmvnorm(yobs,mean=mean1,sigma=covm_can,log=TRUE)
    MHprob2 <- exp(
      logpost_can - logpost +
        loghalfstudent(logxi_can,s=xi_hyper[1],df=xi_hyper[2],log=TRUE) - 
        loghalfstudent(logxi,s=xi_hyper[1],df=xi_hyper[2],log=TRUE)
    )
    if(runif(1) < MHprob2){
      logxi <- logxi_can
      covm_xi <- covm_xi_can
      covm <- covm_can
      logpost <- logpost_can
      acceptMat[s,1] <- 1 
    }
    
    # draw sigma2
    chl <- chol(solve(covm / sigma2))
    sigma2 <- extraDistr::rinvgamma(1,alpha = sigma2_hyper[1]+n/2, 
                                    beta = sum((chl %*% (yobs - mean1))**2)/2+sigma2_hyper[2])
    covm <- sigma2 * g / sumxtx * covm_xi + sigma2 * diag(n)
    
    # draw eta
    if(k > 0){
      covmEta <- solve(t(Zobs) %*% solve(covm) %*% Zobs)
      meanEta <- c(covmEta %*% t(Zobs) %*% solve(covm) %*% yobs)
      eta <- c(rmvnorm(1,mean=meanEta,sigma=covmEta))
      mean1 <- Zobs %*% eta
    }
    
    #store draws
    drawsMat[s,] <- c(logxi,sigma2,eta)  
    
    if(ceiling(s/check1)==s/check1){ #adjust sd of random walk based on acceptance proportions of last 'check1' draws
      probs <- mean(acceptMat[(s-check1+1):s,])
      upper1 <- .5
      lower1 <- .15
      RWsd[probs>upper1] <- RWsd[probs>upper1] * ( (probs[probs>upper1]-upper1)/(1-upper1) + 1)
      RWsd[probs<lower1] <- RWsd[probs<lower1] * 1 / ( 2 - (probs[probs<lower1])/lower1 )
    }
    setTxtProgressBar(pb,s)
  }
  cat("\n","posterior kernel approximated","\n")
  #Sys.time() - wer
  discard <- .5 #discard half for burn-in
  drawsMat <- as.matrix(drawsMat[(discard*samsize+1):samsize,]) #discard with 20% for burn-in
  
  drawsMat_dummy <- drawsMat
  drawsMat_dummy[,2] <- log(drawsMat_dummy[,2])
  meanIS <- apply(drawsMat_dummy,2,mean)
  covIS <- cov(drawsMat_dummy) * 4 #create more spread for importance sample estimate
  
  nIS <- samsizeIS
  sampleIS <- rmvnorm(nIS,mean=meanIS,sigma=covIS)
  
  logint1 <- unlist(lapply(1:nIS,function(s){
    logxi <- sampleIS[s,1]
    logsigma2 <- sampleIS[s,2] 
    eta <- sampleIS[s,-(1:2)]
    
    if(k==0){
      mean1 <- matrix(rep(0,n),ncol=1)
    }else{
      mean1 <- Zobs%*%eta
    }
    
    covm <- exp(logsigma2) * g / sumxtx * exp(-exp(2*logxi)/2*xDiffMmat) * xxobs + 
      exp(logsigma2) * diag(n)
    
    # compute integrand
    dmvnorm(yobs,mean=mean1,sigma=covm,log=TRUE) + 
      loghalfstudent(logxi,xi_hyper[1],xi_hyper[2],log=TRUE) +
      logIGkern(logsigma2,alpha=sigma2_hyper[1],beta=sigma2_hyper[2],log=TRUE) -
      dmvnorm(sampleIS[s,],mean=meanIS,sigma=covIS,log=TRUE)
  }))
  #marginal likelihood H1:nonlinear
  m1 <- log(mean(exp(logint1 - max(logint1)))) + max(logint1)
  
  # marginal likelihood computation for H0: linear effect
  Xmat <- cbind(xobs,Zobs)
  XXi <- solve(t(Xmat)%*%Xmat)
  lm0 <- lm(yobs ~ -1 + Xmat)
  
  samsize0 <- 1e3
  sigma2Draws0 <- extraDistr::rinvgamma(samsize0,sigma2_hyper[1] + n/2,
                                        sigma2_hyper[2] + sum((lm0$residuals)**2)/2) #approximate posterior
  betagammaDraws0 <- matrix(unlist(lapply(1:samsize0,function(s){
    covm1 <- XXi*sigma2Draws0[s]
    rmvnorm(1,mean=lm0$coefficients,sigma=covm1)
  })),byrow=TRUE,nrow=samsize0)
  
  drawsMat0_dummy <- cbind(log(sigma2Draws0),betagammaDraws0[,-1])
  
  #mean and covm of importance sample distr
  meanIS0 <- apply(drawsMat0_dummy,2,mean)
  covmIS0 <- cov(drawsMat0_dummy) * 4
  sampleIS0 <- rmvnorm(samsizeIS,mean=meanIS0,sigma=covmIS0)
  logint0 <- unlist(lapply(1:nrow(sampleIS0),function(s){
    
    logsigma2_0 <- sampleIS0[s,1]
    eta_0 <- sampleIS0[s,-1]
    
    mean1 <- if(k==0){
      mean1 <- matrix(rep(0,n),ncol=1)
    }else{
      mean1 <- Zobs%*%eta_0
    }
    covm0 <- exp(logsigma2_0) * g / sumxtx * xxobs + 
      exp(logsigma2_0) * diag(n)
    
    # compute integrand
    dmvnorm(yobs,mean=mean1,sigma=covm0,log=TRUE) +
      logIGkern(logsigma2_0,alpha=sigma2_hyper[1],beta=sigma2_hyper[2],log=TRUE) -
      dmvnorm(sampleIS0[s,],mean=meanIS0,sigma=covmIS0,log=TRUE)
  }))
  m0 <- log(mean(exp(logint0 - max(logint0)))) + max(logint0)
  
  cat("marginal likelihoods computed","\n")
  
  cat("log(B01) = ",m0-m1,"\n")
  return(list(drawsGPP=drawsMat,m1=m1,m0=m0))
}
# Posterior prediction 
GPPpred4 <- function(X_train,Y_train, Zobs, eta=0, xi=1.0, sigma_y=1,draws=1){
  
  Y_train_eta <- c(Y_train - Zobs %*% eta)
  n <- length(Y_train)
  g <- n
  sumxx <- sum(X_train**2)
  XX_train <- X_train %*% t(X_train)
  
  kxx <- kernelGP(X1=X_train, X2=X_train, xi=xi, tau=sqrt(sigma_y**2 * g / sumxx))
  
  covm_prior <- kxx * XX_train
  covm_data <- sigma_y**2 * diag(n)
  covm2i <- solve(covm_prior + covm_data)
  # compute mean and covariance matrix 
  mean_betax <- c(covm_prior %*% covm2i %*% Y_train_eta)
  covm_betax <- covm_prior - covm_prior %*% covm2i %*% covm_prior
  # draw random "predictive" value
  draw_betax <- rmvnorm(draws,mean=mean_betax,sigma=covm_betax)
  # only the use of the mean (mean_betax) makes sense
  return(list(mean=mean_betax,covm=covm_betax,draw=draw_betax))
}

# Example analysis
set.seed(135)
n <- 50
xobs <- sort(runif(n,min=-3,max=3))
mu = rep(0,length(xobs))
sigma2 <- .1**2
error <- rnorm(n,sd=sqrt(sigma2))
yobs <- .6*dnorm(xobs)*xobs + error
plot(xobs,yobs)
GPPresults <- GPPtest4(yobs,xobs,Zobs=matrix(1,n,1))
# Draw from posterior predictive distribution
preddraw1 <- 50
predGP4 <- do.call(rbind,lapply(1:preddraw1,function(dr){
  rndraw <- sample(nrow(GPPresults[[1]]))[1]
  c(GPPpred4(X_train=xobs,
             Y_train=yobs,
             Zobs=matrix(1,length(xobs),1),
             eta=GPPresults[[1]][rndraw,3],
             xi=exp(GPPresults[[1]][rndraw,1]),
             sigma_y=sqrt(GPPresults[[1]][rndraw,2]),
             draws=1)[[3]]) + GPPresults[[1]][rndraw,3]
}))
#par(mfrow = c(2,1))
plot(xobs,yobs,lwd=2)
for(dr in 1:nrow(predGP4)){
  lines(xobs,predGP4[dr,],col=4)
}
points(xobs,yobs,lwd=2)





