rm(list=ls())

library(shiny)
library(coda)
library(utils)
library(tcltk)

###### Functions ########
source('Functions/SumDiagPlots.R')
source('Functions/PostMatrixProbNoPar.R')
source('Functions/BaySpFixPD.R')


path<-getwd()
unlink(file.path(path,'Posterior Graphs'),recursive=TRUE)


############################################################################
############################################################################

shinyServer(function(input, output) {
  ###### Data NavBar 1. PMP #########
  dataInput1a <- reactive({
    inFile1a <- input$file1a
    if (is.null(inFile1a))
      return(NULL)
    read.csv(inFile1a$datapath, header=input$header1a, sep=input$sep1a)
  })
  
  dataInput1b <- reactive({
    inFile1b <- input$file1b
    if (is.null(inFile1b))
      return(NULL)
    read.csv(inFile1b$datapath, header=input$header1b, sep=input$sep1b)
  })
  
  ######Data NavBar 2. Models #########
  dataInput2a <- reactive({
    inFile2a <- input$file2a
    if (is.null(inFile2a))
      return(NULL)
    read.csv(inFile2a$datapath, header=input$header2a, sep=input$sep2a)
  })
  
  dataInput2b <- reactive({
    inFile2b <- input$file2b
    if (is.null(inFile2b))
      return(NULL)
    read.csv(inFile2b$datapath, header=input$header2b, sep=input$sep2b)
  })
  
  
  ######Formulas NavBar 2. Models ######### 
  
  sumtextM2a <- reactive({
    model.formula(input$Formula2a,dataInput2b())
  })
  
  ######## 1. PMP #########
  PMP<- eventReactive(input$goButton11, {
    withProgress(message = 'Making calculations', value = 0, {
      args <- switch(input$M11,
                     'm111' = list(as.formula(input$Formula1a),dataInput1a(),dataInput1b(),as.numeric(input$Units1),as.numeric(input$Time1),TRUE,TRUE),
                     'm112' = list(as.formula(input$Formula1a),dataInput1a(),dataInput1b(),as.numeric(input$Units1),as.numeric(input$Time1),TRUE,FALSE),
                     'm113' = list(as.formula(input$Formula1a),dataInput1a(),dataInput1b(),as.numeric(input$Units1),as.numeric(input$Time1),FALSE,TRUE)
      )
      do.call(PostMatProb, args)
    })
  })
  
  ####### 2. Models #############
  Posteriors2 <- eventReactive(input$goButton21, {
    #c(100,0,0,0,0,100,0,0,0,0,100,0,0,0,0,100)
    #c(100,0,0,0,0,0,100,0,0,0,0,0,100,0,0,0,0,0,100,0,0,0,0,0,100)
    #c(100,0,0,0,100,0,0,0,100)
    #c(100,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,100)
    #c(100,0,0,0,0,0,0,100,0,0,0,0,0,0,100,0,0,0,0,0,0,100,0,0,0,0,0,0,100,0,0,0,0,0,0,100)
    #c(0,0,0,0,0,0,0)

      if(input$PMeanL==''){bp<-NULL}
      else{
        tas<- sub('[c]',',',isolate(input$PMeanL))
        tas1<- sub('[(]',',',tas)
        tas2<- sub('[)]',',',tas1)
        s3<- strsplit(tas2,',')
        s4<- c(sapply(s3, as.numeric))
        bp<- s4[!is.na(s4)]
      }

      if(input$PVarL==''){Bp<-NULL}
      else{
        tas3<- sub('[c]',',',isolate(input$PVarL))
        tas4<- sub('[(]',',',tas3)
        tas5<- sub('[)]',',',tas4)
        s5<- strsplit(tas5,',')
        s6<- c(sapply(s5, as.numeric))
        s6<- s6[!is.na(s6)]
        if(input$M21=='m211' | input$M21=='m213' | input$M21=='m214' | input$M21=='m217'){
        n1<- (ncol(sumtextM2a()$X)-1)*2
        }
        else{n1<- ncol(sumtextM2a()$X)-1}
        Bp<- matrix(s6,byrow=TRUE,n1,n1)
      }

      if(input$PShL==''){ap<-NULL}
      else{
        ap<- isolate(as.numeric(input$PShL))
      }

      if(input$PScL==''){tp<-NULL}
      else{
        tp<- isolate(as.numeric(input$PScL))
      }
    
    if(input$TLambda1==''){tuning.lambda<-NULL}
    else{
      tuning.lambda<- isolate(as.numeric(input$TLambda1))
    }
    
    if(input$TRho1==''){tuning.rho<-NULL}
    else{
      tuning.rho<- isolate(as.numeric(input$TRho1))
    }
    MCMC<- list(R=input$it,keep=as.numeric(input$keep),burnin=input$burnin,tuning.lambda=tuning.lambda,tuning.rho=tuning.rho)
    Prior<- list(bp=bp,Bp=Bp,ap=ap,tp=tp,tuning.lambda=tuning.lambda,tuning.rho=tuning.rho)

    if(input$M21=='m210')
      return()
    else {
      args <- switch(input$M21,
                     'm211' = list(formula=input$Formula2a,cmatrix=dataInput2a(),data=dataInput2b(),N=as.numeric(input$Units2a),T1=as.numeric(input$Time2a),lambda=TRUE, rho=TRUE, durbin=TRUE, prior=Prior, Mcmc=MCMC),
                     'm212' = list(formula=input$Formula2a,cmatrix=dataInput2a(),data=dataInput2b(),N=as.numeric(input$Units2a),T1=as.numeric(input$Time2a),lambda=TRUE, rho=TRUE, durbin=FALSE, prior=Prior, Mcmc=MCMC),
                     'm213' = list(formula=input$Formula2a,cmatrix=dataInput2a(),data=dataInput2b(),N=as.numeric(input$Units2a),T1=as.numeric(input$Time2a),lambda=TRUE, rho=FALSE, durbin=TRUE, prior=Prior, Mcmc=MCMC),
                     'm214' = list(formula=input$Formula2a,cmatrix=dataInput2a(),data=dataInput2b(),N=as.numeric(input$Units2a),T1=as.numeric(input$Time2a),lambda=FALSE, rho=TRUE, durbin=TRUE, prior=Prior, Mcmc=MCMC),
                     'm215' = list(formula=input$Formula2a,cmatrix=dataInput2a(),data=dataInput2b(),N=as.numeric(input$Units2a),T1=as.numeric(input$Time2a),lambda=TRUE, rho=FALSE, durbin=FALSE, prior=Prior, Mcmc=MCMC),
                     'm216' = list(formula=input$Formula2a,cmatrix=dataInput2a(),data=dataInput2b(),N=as.numeric(input$Units2a),T1=as.numeric(input$Time2a),lambda=FALSE, rho=TRUE, durbin=FALSE, prior=Prior, Mcmc=MCMC),
                     'm217' = list(formula=input$Formula2a,cmatrix=dataInput2a(),data=dataInput2b(),N=as.numeric(input$Units2a),T1=as.numeric(input$Time2a),lambda=FALSE, rho=FALSE, durbin=TRUE, prior=Prior, Mcmc=MCMC),
                     'm218' = list(formula=input$Formula2a,cmatrix=dataInput2a(),data=dataInput2b(),N=as.numeric(input$Units2a),T1=as.numeric(input$Time2a),lambda=FALSE, rho=FALSE, durbin=FALSE, prior=Prior, Mcmc=MCMC)
      )}
    do.call(BaySpFixPD,args)
  })
  
  ####### 2.1 Post Matrix Prob: Download Posterior Chains##########
  
  output$download11 <- downloadHandler(
    filename = function() { 
      paste('Posterior Chains', '.csv', sep=',') 
    },
    
    content = function(file) {
      
      if(input$M11=='m110')
        content<- return()
      
      else{write.csv(PMP(), file)}
    }
  )
  
  ####### 2.1 Models: Download Posterior Chains##########
  
  output$download21 <- downloadHandler(
    filename = function() { 
      paste('Posterior Chains', '.csv', sep=',') 
    },
    
    content = function(file) {
      
      if(input$M21=='m210')
        content<- return()
      
      switch(input$M21,
             'm211' = post21<- cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[],Posteriors2()$rho.post[]),
             'm212' = post21<- cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[],Posteriors2()$rho.post[]),
             'm213' = post21<- cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[]),
             'm214' = post21<- cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$rho.post[]),
             'm215' = post21<- cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[]),
             'm216' = post21<- cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$rho.post[]),
             'm217' = post21<- cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[]),
             'm218' = post21<- cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[])
      )     
      write.csv(post21, file)
    }
  )
  

  ####### 1.1 PMP: Print ##########
  
  output$summary11 <- renderPrint({
    if(input$M11=='m110'){
      return('SELECT A MODEL')}
    
    else{print(PMP())}
  })
  
  ####### 2.1 Bayesian Panel Data Models: Print ##########
    output$summary21 <- renderPrint({
    if(input$M21=='m210'){
      return('SELECT A MODEL')}

    else{switch(input$M21,
                  'm211' = SumDiag.m211(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[],Posteriors2()$rho.post[],Posteriors2()$lambda.rat[],Posteriors2()$rho.rat[]),
                  'm212' = SumDiag.m212(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[],Posteriors2()$rho.post[],Posteriors2()$lambda.rat[],Posteriors2()$rho.rat[]),
                  'm213' = SumDiag.m213(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[],Posteriors2()$lambda.rat[]),
                  'm214' = SumDiag.m214(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$rho.post[],Posteriors2()$rho.rat[]),
                  'm215' = SumDiag.m215(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[],Posteriors2()$lambda.rat[]),
                  'm216' = SumDiag.m216(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$rho.post[],Posteriors2()$rho.rat[]),
                  'm217' = SumDiag.m217(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[],Posteriors2()$lambda.rat[]),
                  'm218' = SumDiag.m218(Posteriors2()$b.post[,],Posteriors2()$s_v.post[])
                    )}
  })
  
  ####### 2.1 Models: Graphs Posterior Chains##########  
  output$plot21 <- renderPlot({
    unlink(file.path(path,'Posterior Graphs'),recursive=TRUE)
    dir.create(file.path(path,'Posterior Graphs'),showWarnings = FALSE)
    setwd(file.path(path,'Posterior Graphs'))
    
    graphs21<- function(post21){
      nc<-ncol(post21)
      for (i in 1:nc) {
        pdf(paste('Density Plot',paste(i,'.pdf', sep = '', collapse = NULL)))
        Plot(post21[,i])
        dev.off()
        setEPS()
        postscript(paste('Density Plot',paste(i,'.eps', sep = '', collapse = NULL)))
        Plot(post21[,i])
        dev.off()
        pdf(paste('Trace Plot',paste(i,'.pdf', sep = '', collapse = NULL)))
        Plot.trace(post21[,i])
        dev.off()
        setEPS()
        postscript(paste('Trace Plot',paste(i,'.eps', sep = '', collapse = NULL)))
        Plot.trace(post21[,i])
        dev.off()
        pdf(paste('Autocorrelation Plot',paste(i,'.pdf', sep = '', collapse = NULL)))
        Plot.corr(post21[,i])
        dev.off()
        setEPS()
        postscript(paste('Autocorrelation Plot',paste(i,'.eps', sep = '', collapse = NULL)))
        Plot.corr(post21[,i])
        dev.off()
      }
    }
    switch(input$M21,
           'm211' = graphs21(cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[],Posteriors2()$rho.post[])),
           'm212' = graphs21(cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[],Posteriors2()$rho.post[])),
           'm213' = graphs21(cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[])),
           'm214' = graphs21(cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$rho.post[])),
           'm215' = graphs21(cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[])),
           'm216' = graphs21(cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$rho.post[])),
           'm217' = graphs21(cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[],Posteriors2()$lambda.post[])),
           'm218' = graphs21(cbind(Posteriors2()$b.post[,],Posteriors2()$s_v.post[]))
    )
    setwd('..')
  })
  output$multiDownload21 <- downloadHandler(
    filename = function() {
      paste("Posterior Graphs", "zip", sep=".")
    },
    
    content = function(file) {
      zip(zipfile=file, files='Posterior Graphs')
    },
    contentType = "application/zip"
  )
  
  
  
})
