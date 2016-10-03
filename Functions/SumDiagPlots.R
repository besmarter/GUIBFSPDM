#############################
#####Auxiliary functions#####
#############################

model.formula<- function(formula, data=list(), ...)
{
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  return(list(X=X,y=y))
}


SumDiag.m211<- function(betadraw,sigmasqdraw,lambdadraw,rhodraw,lambdarat,rhorat){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryScale<- summary(sigmasqdraw)
  GewekeTestScale<- geweke.diag(sigmasqdraw)
  RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestScale<- heidel.diag(sigmasqdraw)
  SummaryLambda<- summary(lambdadraw)
  GewekeTestLambda<- geweke.diag(lambdadraw)
  RafteryTestLambda<- raftery.diag(lambdadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLambda<- heidel.diag(lambdadraw)
  SummaryRho<- summary(rhodraw)
  GewekeTestRho<- geweke.diag(rhodraw)
  RafteryTestRho<- raftery.diag(rhodraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestRho<- heidel.diag(rhodraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale,
              SummaryLambda = SummaryLambda, GewekeTestLambda = GewekeTestLambda, RafteryTestLambda = RafteryTestLambda, HeidelTestLocation = HeidelTestLambda, SummaryRho = SummaryRho, GewekeTestRho = GewekeTestRho, RafteryTestRho = RafteryTestRho, HeidelTestRho = HeidelTestRho,
              lambdarat,rhorat))
}

SumDiag.m212<- function(betadraw,sigmasqdraw,lambdadraw,rhodraw,lambdarat,rhorat){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryScale<- summary(sigmasqdraw)
  GewekeTestScale<- geweke.diag(sigmasqdraw)
  RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestScale<- heidel.diag(sigmasqdraw)
  SummaryLambda<- summary(lambdadraw)
  GewekeTestLambda<- geweke.diag(lambdadraw)
  RafteryTestLambda<- raftery.diag(lambdadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLambda<- heidel.diag(lambdadraw)
  SummaryRho<- summary(rhodraw)
  GewekeTestRho<- geweke.diag(rhodraw)
  RafteryTestRho<- raftery.diag(rhodraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestRho<- heidel.diag(rhodraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale,
              SummaryLambda = SummaryLambda, GewekeTestLambda = GewekeTestLambda, RafteryTestLambda = RafteryTestLambda, HeidelTestLocation = HeidelTestLambda, SummaryRho = SummaryRho, GewekeTestRho = GewekeTestRho, RafteryTestRho = RafteryTestRho, HeidelTestRho = HeidelTestRho,
              lambdarat,rhorat))
}

SumDiag.m213<- function(betadraw,sigmasqdraw,lambdadraw,lambdarat){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryScale<- summary(sigmasqdraw)
  GewekeTestScale<- geweke.diag(sigmasqdraw)
  RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestScale<- heidel.diag(sigmasqdraw)
  SummaryLambda<- summary(lambdadraw)
  GewekeTestLambda<- geweke.diag(lambdadraw)
  RafteryTestLambda<- raftery.diag(lambdadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLambda<- heidel.diag(lambdadraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale,
              SummaryLambda = SummaryLambda, GewekeTestLambda = GewekeTestLambda, RafteryTestLambda = RafteryTestLambda, HeidelTestLocation = HeidelTestLambda,lambdarat))
}

SumDiag.m214<- function(betadraw,sigmasqdraw,rhodraw,rhorat){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryScale<- summary(sigmasqdraw)
  GewekeTestScale<- geweke.diag(sigmasqdraw)
  RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestScale<- heidel.diag(sigmasqdraw)
  SummaryRho<- summary(rhodraw)
  GewekeTestRho<- geweke.diag(rhodraw)
  RafteryTestRho<- raftery.diag(rhodraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestRho<- heidel.diag(rhodraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale,
              SummaryRho = SummaryRho, GewekeTestRho = GewekeTestRho, RafteryTestRho = RafteryTestRho, HeidelTestRho = HeidelTestRho,rhorat))
}

SumDiag.m215<- function(betadraw,sigmasqdraw,lambdadraw,lambdarat){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryScale<- summary(sigmasqdraw)
  GewekeTestScale<- geweke.diag(sigmasqdraw)
  RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestScale<- heidel.diag(sigmasqdraw)
  SummaryLambda<- summary(lambdadraw)
  GewekeTestLambda<- geweke.diag(lambdadraw)
  RafteryTestLambda<- raftery.diag(lambdadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLambda<- heidel.diag(lambdadraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale,
              SummaryLambda = SummaryLambda, GewekeTestLambda = GewekeTestLambda, RafteryTestLambda = RafteryTestLambda, HeidelTestLocation = HeidelTestLambda,lambdarat))
}

SumDiag.m216<- function(betadraw,sigmasqdraw,rhodraw,rhorat){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryScale<- summary(sigmasqdraw)
  GewekeTestScale<- geweke.diag(sigmasqdraw)
  RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestScale<- heidel.diag(sigmasqdraw)
  SummaryRho<- summary(rhodraw)
  GewekeTestRho<- geweke.diag(rhodraw)
  RafteryTestRho<- raftery.diag(rhodraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestRho<- heidel.diag(rhodraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale,
              SummaryRho = SummaryRho, GewekeTestRho = GewekeTestRho, RafteryTestRho = RafteryTestRho, HeidelTestRho = HeidelTestRho,rhorat))
}

SumDiag.m217<- function(betadraw,sigmasqdraw,lambdadraw,lambdarat){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryScale<- summary(sigmasqdraw)
  GewekeTestScale<- geweke.diag(sigmasqdraw)
  RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestScale<- heidel.diag(sigmasqdraw)
  SummaryLambda<- summary(lambdadraw)
  GewekeTestLambda<- geweke.diag(lambdadraw)
  RafteryTestLambda<- raftery.diag(lambdadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLambda<- heidel.diag(lambdadraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale,
              SummaryLambda = SummaryLambda, GewekeTestLambda = GewekeTestLambda, RafteryTestLambda = RafteryTestLambda, HeidelTestLocation = HeidelTestLambda,lambdarat))
}

SumDiag.m218<- function(betadraw,sigmasqdraw){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryScale<- summary(sigmasqdraw)
  GewekeTestScale<- geweke.diag(sigmasqdraw)
  RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestScale<- heidel.diag(sigmasqdraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale))
}

Plot<- function(betadraw){
  hist(betadraw, breaks=20,freq=FALSE, xlab="Parameter", main="", col="lightgreen")
  lines(density(betadraw,na.rm = TRUE), col="red", lwd=2)
  abline(h = NULL, v = c(quantile(betadraw,c(0.025, 0.975))), col = "purple", lwd=2)
  text(quantile(betadraw,c(0.025)),y=1, "Quantile 2.5%", col = "black", adj = c(0,-0.5), cex=0.75)
  text(quantile(betadraw,c(0.975)),y=1, "Quantile 97.5%", col = "black", adj = c(0,-0.5), cex=0.75)
  abline(h = NULL, v = c(quantile(betadraw,c(0.5))), col = "red", lwd=3)
  abline(h = NULL, v = mean(betadraw), col = "blue", lwd=2)
  legend("topleft",inset=.05,cex = 0.75,c("Median","Mean"),horiz=TRUE,lty=c(1,1),lwd=c(1,1),col=c("red","blue"),bg="grey96")
}

Plot.trace<- function(betadraw){traceplot(mcmc(betadraw), main = "Trace Plot", xlab="Iteration", ylab= "Parameter", col= "blue")}
Plot.corr<- function(betadraw){autocorr.plot(mcmc(betadraw), main = "Autocorrelation", col="blue")}
