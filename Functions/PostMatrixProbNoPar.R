#######################################
#      Posterior Model Probability    #
#          Contiguity Matrix          #
#          September 9, 2016          #
#   Professor: Andres Ramirez Hassan  #
#        aramir21@eafit.edu.co        #
#        www.besmarter-team.org       #
#######################################



#rm(list=ls())
library(cubature)
#adaptIntegrate function
library(Matrix)
#crossprod and tcrossprod
library(pbapply)
#Progress bar

PostMatProb<- function(formula, matrices, data, n, t, lambda=TRUE, rho=TRUE){
  pandterm = function(message) {
    stop(message, call. = FALSE)
  }
  if (missing(formula)) {
    pandterm("Requires formula argument")
  }
  if (missing(matrices)) {
    pandterm("Requires list of contiguity matrices")
  }
  if (missing(data)) {
    pandterm("Requires data argument")
  }
  mf <- model.frame(formula=formula, data=data)
  X <- as.matrix(mf[,-1])
  y <- as.matrix(mf[,1])
  
  if (missing(n)) {
    pandterm("Requires number of cross sectional units")
  }
  if (missing(t)) {
    pandterm("Requires number of periods")
  }
  N<<- n
  T1<<- t
  if (T1*N != nrow(X)) {
    pandterm("Requires a balanced panel data set")
  }
  if (N != nrow(matrices)^0.5) {
    pandterm("Wrong dimension: check contiguity matrices")
  }
  m<- ncol(matrices) 
  Ws<<-lapply(1:m, function(x) matrix(matrices[,x],N,N)) 
  #Create a list with these matrices
  if (missing(lambda)) {
    pandterm("Requires lambda=TRUE or lambda=FALSE")
  }
  if (missing(rho)) {
    pandterm("Requires rho=TRUE or rho=FALSE")
  }
  if (lambda == FALSE & rho == FALSE) {
    pandterm("Requires SAR and/or SEM effects")
  }
  k = ncol(X)
  #Number of controls
  it<-rep(1,T1) 
  #Vector of ones
  i<-rep(1:N,T1)
  #Cross sections units
  t<-rep(1:T1,each=N) 
  #Time units
  It<-diag(T1)
  In<<-diag(N)
  It1<<-diag(T1-1)
  M<-It-tcrossprod(it)/T1       
  #Mean deviation matrix
  V<-eigen(M)$vectors      
  #Orthonormal eigen vectors matrix
  F<-V[,-T1]                
  #First T-1 eigenvectors
  Q<-kronecker(t(F),In)    
  #Transformation matrix
  Qy<<-Q%*%y
  #y* in paper
  QX<<-Q%*%X
  #X* in paper
  
  #########################################
  # Hyperparameters Prior Distributions   #
  #########################################
  Bp<-1000*diag(k)                    
  #Covariance prior matrix
  KK<-chol(Bp)
  KKI<- solve(KK)
  Bpi<<-tcrossprod(KKI)
  #Inverse Covariance prior matrix
  bp<<-matrix(rep(0,k),k,1)    
  #Prior mean vector
  ap<<-0.001                                   
  #Prior Alpha
  tp<-0.001                                   
  #Prior tau
  #########################################
  #     Posterior Matrix Probability      #
  #########################################
  #              Marginal                 #
  #########################################
  
  if (lambda==TRUE & rho == TRUE){
  Mat<<-function(Wx) {
    Int<-function(z){
      lambda<-z[1]
      rho<-z[2]
      A<-In-lambda*Wx                                     #Matrix A in paper
      B<-In-rho*Wx                                        #Matrix B in paper
      Si<-kronecker(It1,t(B)%*%B)                         #Inverse Sigma
      GG<-t(QX)%*%Si%*%QX+Bpi                             #GG matrix in paper
      KK1<-chol(GG)
      KK1I<- solve(KK1)
      GGi<-tcrossprod(KK1I)                               #Inverse GG in paper
      Qy_bar<-kronecker(It1,A)%*%Qy                       #y~* in paper
      KK2<-chol(t(QX)%*%Si%*%QX)
      KK2I<- solve(KK2)
      bMV<-tcrossprod(KK2I)%*%t(QX)%*%Si%*%Qy_bar   #Max Likehood Beta
      bmean<-GGi%*%(Bpi%*%bp+(t(QX)%*%Si%*%QX)%*%bMV)     #B* in paper (Mean vector B post)
      v<-kronecker(It1,B)%*%(kronecker(It1,A)%*%Qy-QX%*%bmean) #v~* in paper using bmean.post
      suma<-(crossprod(v)+t(bp-bmean)%*%Bpi%*%(bp-bmean)+ap/2)^(-(N*(T1-1)/2+ap)) #Delta in paper
      MarLik<-(det(A))^(T1-1)*(det(B))^(T1-1)*suma*(det(GG))^(-1/2) #Marginal Likelihood in paper
      return(MarLik)
    }
  }
  }

  if (lambda==TRUE & rho == FALSE){
    Mat<<-function(Wx) {
      Int<-function(z){
        lambda<-z
        A<-In-lambda*Wx                                     #Matrix A in paper
        B<-In                                               #Matrix B in paper
        Si<-kronecker(It1,t(B)%*%B)                         #Inverse Sigma
        GG<-t(QX)%*%Si%*%QX+Bpi                             #GG matrix in paper
        KK1<-chol(GG)
        KK1I<- solve(KK1)
        GGi<-tcrossprod(KK1I)                               #Inverse GG in paper
        Qy_bar<-kronecker(It1,A)%*%Qy                       #y~* in paper
        KK2<-chol(t(QX)%*%Si%*%QX)
        KK2I<- solve(KK2)
        bMV<-tcrossprod(KK2I)%*%t(QX)%*%Si%*%Qy_bar   #Max Likehood Beta
        bmean<-GGi%*%(Bpi%*%bp+(t(QX)%*%Si%*%QX)%*%bMV)     #B* in paper (Mean vector B post)
        v<-kronecker(It1,B)%*%(kronecker(It1,A)%*%Qy-QX%*%bmean) #v~* in paper using bmean.post
        suma<-(crossprod(v)+t(bp-bmean)%*%Bpi%*%(bp-bmean)+ap/2)^(-(N*(T1-1)/2+ap)) #Delta in paper
        MarLik<-(det(A))^(T1-1)*suma*(det(GG))^(-1/2)        #Marginal Likelihood in paper
        return(MarLik)
      }
    }
  }

  if (lambda==FALSE & rho == TRUE){
    Mat<<-function(Wx) {
      Int<-function(z){
        rho<-z
        A<-In                                               #Matrix A in paper
        B<-In-rho*Wx                                        #Matrix B in paper
        Si<-kronecker(It1,t(B)%*%B)                         #Inverse Sigma
        GG<-t(QX)%*%Si%*%QX+Bpi                             #GG matrix in paper
        KK1<-chol(GG)
        KK1I<- solve(KK1)
        GGi<-tcrossprod(KK1I)                               #Inverse GG in paper
        Qy_bar<-kronecker(It1,A)%*%Qy                       #y~* in paper
        KK2<-chol(t(QX)%*%Si%*%QX)
        KK2I<- solve(KK2)
        bMV<-tcrossprod(KK2I)%*%t(QX)%*%Si%*%Qy_bar   #Max Likehood Beta
        bmean<-GGi%*%(Bpi%*%bp+(t(QX)%*%Si%*%QX)%*%bMV)     #B* in paper (Mean vector B post)
        v<-kronecker(It1,B)%*%(kronecker(It1,A)%*%Qy-QX%*%bmean) #v~* in paper using bmean.post
        suma<-(crossprod(v)+t(bp-bmean)%*%Bpi%*%(bp-bmean)+ap/2)^(-(N*(T1-1)/2+ap)) #Delta in paper
        MarLik<-(det(B))^(T1-1)*suma*(det(GG))^(-1/2)        #Marginal Likelihood in paper
        return(MarLik)
      }
    }
  }
  
  if (lambda == TRUE & rho == TRUE){
    ll<<- c(-1,-1)
    ul<<- c(1,1)
  }
  else {
    ll<<- -1
    ul<<- 1
  }
  # pb <- winProgressBar(title = "progress bar", min = 0, max = 1, width = 100)
  pboptions(type='win')
  PropMat<-lapply(Ws, function(x) adaptIntegrate(Mat(x), lowerLimit = ll, upperLimit = ul))
  PropMatIntNew<<-c(sapply(1:m, function(l) unlist(PropMat[[l]]$integral)))
  ########Posterior Model Probaility#############################
  PMP<- sapply(PropMatIntNew,function(x) x/sum(PropMatIntNew)) #The highest PMP should be in 5th position!!!
  # closepb(pb)
  # return(list(PropMat=PropMat,PropMatIntNew=PropMatIntNew,PMP=PMP))
  return(PMP)
}

  
