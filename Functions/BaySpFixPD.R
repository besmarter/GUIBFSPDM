#######################################
#      Posterior Model Probability    #
#          Contiguity Matrix          #
#          September 9, 2016          #
#   Professor: Andres Ramirez Hassan  #
#        aramir21@eafit.edu.co        #
#        www.besmarter-team.org       #
#######################################



#rm(list=ls())
library(Matrix)
#crossprod and tcrossprod
library(splm)
#splm
# library(plm)
#plm: model without explicit spatial effects 
library(MASS)
#mvrnorm
library(MCMCpack)
#rinvgamma
library(utils)
library(tcltk)
#Progress bar

BaySpFixPD<- function(formula, cmatrix=NULL, data, N, T1, lambda=TRUE, rho=TRUE, durbin=TRUE, prior, Mcmc)
  {
  pandterm = function(message) {
    stop(message, call. = FALSE)
  }
  if (missing(formula)) {
    pandterm('Requires formula argument')
  }
  if (lambda==TRUE & missing(cmatrix) | rho==TRUE & missing(cmatrix) | durbin==TRUE & missing(cmatrix)) {
    pandterm('Requires contiguity matrix')
  }
  if (lambda==TRUE | rho==TRUE | durbin==TRUE){
    Ws<-as.matrix(cmatrix)
  }
  #Contiguity matrix
  if (missing(data)) {
    pandterm('Requires data argument')
  }
  mf <- model.frame(formula=formula, data=data)
  Z <- as.matrix(mf[,-1])
  y <- as.matrix(mf[,1])
  
  if (missing(N)) {
    pandterm('Requires number of cross sectional units')
  }
  if (lambda==FALSE & rho==FALSE & durbin==FALSE) {
    Ws<-diag(N)
  }
  if (missing(T1)) {
    pandterm('Requires number of periods')
  }
  if (T1*N != nrow(Z)) {
    pandterm('Requires a balanced panel data set')
  }
  if (N != nrow(Ws) | N != ncol(Ws)) {
    pandterm('Wrong dimension: check contiguity matrix')
  }
  if (missing(lambda)) {
    pandterm('Requires lambda=TRUE or lambda=FALSE')
  }
  if (missing(rho)) {
    pandterm('Requires rho=TRUE or rho=FALSE')
  }
  if (missing(durbin)) {
    pandterm('Requires durbin=TRUE or durbin=FALSE')
  }
  it<-rep(1,T1) 
  #Vector of ones
  i<-rep(1:N,T1)
  #Cross sections units
  It<-diag(T1)
  In<-diag(N)
  It1<-diag(T1-1)
  M<-It-tcrossprod(it)/T1       
  #Mean deviation matrix
  V<-eigen(M)$vectors      
  #Orthonormal eigen vectors matrix
  F<-V[,-T1]                
  #First T-1 eigenvectors
  Q<-kronecker(t(F),In)    
  #Transformation matrix
  Qy<-Q%*%y                #y* in paper
  QZ<-Q%*%Z                #Z* in paper
  if (durbin == TRUE){  
    QZW<- kronecker(It1,Ws)%*%QZ
    QXWX<- cbind(QZ,QZW)
  }
  else {QXWX<-QZ}               #X* in paper
  k<- ncol(QXWX)
  formula1<-Qy~QXWX
  #########################################
  # Hyperparameters Prior Distributions   #
  #########################################
  #Prior parameters
  if (missing(prior)) {
    bp <- c(rep(0, k))                 #Prior betas vector
    Bp <- 1000*diag(k)                 #Prior covariance prior matrix
    Bpi <- 0.001*diag(k)               #Inverse prior covariance prior matrix
    ap <- 0.001                        #Prior Alpha
    tp <- 0.001                        #Prior tau
  }
  else {
    if (is.null(prior$bp)) {
      bp <- c(rep(0, k))
    }
    else {
      bp <- prior$bp
    }
    if (is.null(prior$Bp)) {
      Bp <- 1000*diag(k)
      Bpi <- solve(Bp)
    }
    else {
      Bp <- prior$Bp
      Bpi <- solve(Bp)
    }
    
    if (is.null(prior$ap)) {
      ap <- 0.001
    }
    else {
      ap <- prior$ap
    }
    if (is.null(prior$tp)) {
      tp <- 0.001
    }
    else {
      tp = prior$tp
    }
  }
  if (ncol(Bp) != nrow(Bp) | ncol(Bp) != k | nrow(Bp) != k) {
    pandterm('Bad dimensions for Bp')
  }
  if (length(bp) != k) {
    pandterm(paste('bp wrong length, length= ', length(bp)))
  }
  if (missing(Mcmc)) {
    pandterm('Requires Mcmc argument')
  }
  else {
    if (is.null(Mcmc$R)) {
      pandterm('Requires Mcmc element R')
    }
    else {
      R = Mcmc$R
    }
    if (is.null(Mcmc$keep)) {
      keep = 1
    }
    else {
      keep = Mcmc$keep
    }
    if (is.null(Mcmc$burnin)) {
      burnin = 1
    }
    else {
      burnin = Mcmc$burnin
    }
  }

  if (lambda == TRUE) {
    lag<-1
  }
  else {lag<-0}
  if (rho == TRUE) {
    spef<- 'kkp'
  }
  else {spef<-'none'}
  
  Wf <- mat2listw(Ws, row.names = NULL)

  if (lambda == TRUE | rho == TRUE){
  sararMV<- spml(formula=formula1, data=as.data.frame(cbind(rep(1:N,T1-1),rep(1:(T1-1),each=N),Qy,QXWX)), index=NULL, listw=Wf, lag=lag, spatial.error=spef, model='within',
                     effect='individual', method='eigen', na.action=na.fail, quiet=TRUE, zero.policy=NULL, tol.solve=1e-10,
                     control=list(), legacy=FALSE)
  }

  if (lambda == TRUE & rho == TRUE){
    lambda1<-sararMV$coefficients[1]
    rho1<-sararMV$coefficients[2]
    b1<-sararMV$coefficients[-c(1,2)]
    s1<- sararMV$sigma2
    s.lambda<- sararMV$vcov[1,1]^0.5           #Standard error of lambda from MV
    s.rho<- sararMV$vcov[2,2]^0.5              #Standard error of lambda from MV
  }
  if (lambda == TRUE & rho == FALSE){
    lambda1<-sararMV$coefficients[1]
    rho1<-0
    b1<-sararMV$coefficients[-c(1)]
    s1<- sararMV$sigma2
    s.lambda<- sararMV$vcov[1,1]^0.5           #Standard error of lambda from MV
    s.rho<- NULL                               #Standard error of lambda from MV
  }
  if (lambda == FALSE & rho == TRUE){
    lambda1<-0
    rho1<-sararMV$coefficients[1]
    b1<-sararMV$coefficients[-c(1)]
    s1<- sararMV$sigma2
    s.lambda<- NULL                            #Standard error of lambda from MV
    s.rho<- sararMV$vcov[1,1]^0.5              #Standard error of lambda from MV
  }
  if (lambda == FALSE & rho == FALSE){
    formula1<-Qy~QXWX-1
    sararMV<- lm(formula=formula1, data=as.data.frame(Qy,QXWX))
    b1<-sararMV$coefficients
    s1<- sum(sararMV$residuals^2)/sararMV$df.residual
    lambda1<-0
    rho1<-0
    s.lambda<- NULL                            #Standard error of lambda from MV
    s.rho<- NULL                               #Standard error of lambda from MV
  }

  if (is.null(Mcmc$tuning.lambda)) {
    s.lambda <- s.lambda
  }
  else {
    s.lambda <- Mcmc$tuning.lambda
  }

  if (is.null(Mcmc$tuning.rho)) {
    s.rho <- s.rho
  }
  else {
    s.rho <- Mcmc$tuning.rho
  }

  b.post <- matrix(double(R*k), ncol = k)
  s_v.post <- double(R)
  lambda.post <- double(R)
  rho.post <- double(R)
  lambda.rate <- double(R)
  rho.rate <- double(R)
  # pb <- winProgressBar(title = "progress bar", min = 0,
  #                      max = R, width = 100)
  withProgress(message = 'Making calculations', value = 0, {
  for(l in 1:R){
    #######Gibbs sampler Betas#########
    ###################################
    A.post<-In-lambda1*Ws                                              #Matrix A in paper
    Ai.post<-solve(A.post)                                             #A inverse
    B.post<-In-rho1*Ws                                                 #Matrix B in paper
    Bi.post<-solve(B.post)                                             #B inverse
    S.post<-kronecker(It1,Bi.post)%*%kronecker(It1,t(Bi.post))         #Sigma matrix
    Si.post<-solve(S.post)                                             #Inverse Sigma
    O.post<-kronecker(It1,Ai.post)%*%kronecker(It1,t(Ai.post))         #Omega matrix
    Oi.post<-solve(O.post)                                             #Inverse Omega
    GG.post<-t(QXWX)%*%Si.post%*%QXWX+Bpi                              #GG matrix in paper
    GGi.post<-solve(GG.post)                                           #Inverse GG (Covariance matrix B post)
    Qy_bar.post<-kronecker(It1,A.post)%*%Qy                            #y~* in paper
    bMV<-solve(t(QXWX)%*%Si.post%*%QXWX)%*%t(QXWX)%*%Si.post%*%Qy_bar.post  #Max Likehood Beta
    bmean.post<-GGi.post%*%(Bpi%*%bp+(t(QXWX)%*%Si.post%*%QXWX)%*%bMV) #B* in paper (Mean vector B post)
    b1<-mvrnorm(1,bmean.post,s1*GGi.post)           #Posterior distribution Betas

    #######Gibbs sampler Sigma#########
    ###################################
    sum.post1<- t(Qy_bar.post-QXWX%*%bmean.post)%*%Si.post%*%(Qy_bar.post-QXWX%*%bmean.post)+t(bp-bmean.post)%*%Bpi%*%(bp-bmean.post)
    a.post<-0.5*(N*(T1-1)+2*ap)                                         #Posterior inverse gamma shape parameter
    t.post<-0.5*(sum.post1+2*tp)                                       #Posterior inverse gamma sclae parameter
    s1<-rinvgamma(1,a.post,t.post)                            #Posterior distribution sigma_v

    ###Metropilis-Hastings lambda#####
    ###################################
    if (lambda == TRUE){
      repeat{
        e.lambda<-runif(1, min=-1.96*s.lambda, max=1.96*s.lambda)
        lambda.pot<-lambda1+e.lambda
        if (lambda.pot<1 & lambda.pot>-1){break}
      }
      A.pot<-In-lambda.pot*Ws                                    #Matrix A in paper using lambda.pot
      Qy_bar.pot<-kronecker(It1,A.pot)%*%Qy                      #y~* in paper
      v<-kronecker(It1,B.post)%*%(Qy_bar.pot-QXWX%*%b1)          #v~* in paper
      sum.pot<-t(v)%*%v                                          #Squared sum using lambda.pot
      v<-kronecker(It1,B.post)%*%(Qy_bar.post-QXWX%*%b1)         #v~* in paper
      sum.post<-t(v)%*%v
      alpha.lambda<-min(exp(-(1/(2*s1))*(sum.pot-sum.post))*(det(A.pot)/det(A.post))^(T1-1),1) #Alpha in paper
      v1<-runif(1,0,1)                                           #Uniform random number to select lambda
      if (v1<=alpha.lambda) {lambda1<-lambda.pot ; lambda.rate[l]<-1} else {lambda1<-lambda1 ; lambda.rate[l]<-0} #Posterior distribution lambda and non rejection rate
    }
    else{lambda1<-0}
    #####Metropilis-Hastings rho######
    ###################################
    if (rho == TRUE){
    repeat{
      e.rho<-runif(1, min=-1.96*s.rho, max=1.96*s.rho)
      rho.pot<-rho1+e.rho
      if (rho.pot<1 & rho.pot>-1){break}
    }
    A.post<-In-lambda1*Ws                                       #Matrix A in paper using new lambda
    Qy_bar.post<-kronecker(It1,A.post)%*%Qy                     #y~* in paper
    B.pot<-In-rho.pot*Ws                                        #Matrix B in paper using rho.pot
    v<-kronecker(It1,B.pot)%*%(Qy_bar.post-QXWX%*%b1)           #v~* in paper
    sum.pot<-t(v)%*%v                                           #Squared sum using rho.pot
    v<-kronecker(It1,B.post)%*%(Qy_bar.post-QXWX%*%b1)          #v~* in paper
    sum.post<-t(v)%*%v                                          #Squared sum using rho[l-1]
    alpha.rho<-min(exp(-(1/(2*s1))*(sum.pot-sum.post))*(det(B.pot)/det(B.post))^(T1-1),1) #Alpha in paper
    v1<-runif(1,0,1)                                           #Uniform random number to select lambda
    if (v1<=alpha.rho) {rho1<-rho.pot ; rho.rate[l]<-1} else {rho1<-rho1 ; rho.rate[l]<-0} #Posterior distribution lambda and non rejection rate
    }
    else{rho1<-0}
    lambda.post[l]<-lambda1
    rho.post[l]<-rho1
    s_v.post[l]<-s1
    b.post[l,]<-b1
    lambda.rat<-sum(lambda.rate)/R
    rho.rat<-sum(rho.rate)/R
    # setWinProgressBar(pb, l, title=paste( round(l/R*100, 0),
    #                                         "% done"))
    incProgress(1/l, detail = paste('Doing iteration', l))
    }
})
  lambda.post<-lambda.post[seq(burnin,R,keep)]
  rho.post<-rho.post[seq(burnin,R,keep)]
  s_v.post<-s_v.post[seq(burnin,R,keep)]
  b.post<-b.post[seq(burnin,R,keep),]

  attributes(b.post)$class = c('mcmc')
  attributes(b.post)$mcpar = c(burnin, R, keep)
  attributes(s_v.post)$class = c('mcmc')
  attributes(s_v.post)$mcpar = c(burnin, R, keep)
  attributes(lambda.post)$class = c('mcmc')
  attributes(lambda.post)$mcpar = c(burnin, R, keep)
  attributes(rho.post)$class = c('mcmc')
  attributes(rho.post)$mcpar = c(burnin, R, keep)
  # close(pb)
  return(list(b.post = b.post, s_v.post = s_v.post, lambda.post = lambda.post, rho.post = rho.post, lambda.rat = lambda.rat, rho.rat = rho.rat))
}


