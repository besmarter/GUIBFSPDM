####################################################
#           Posterior Model Probability            #
# Bayesian General Nesting Spatial Panel Data Model#
#                Contiguity Matrix                 #
#                   Simulation                     #
#                September 12, 2016                #
#         Professor: Andres Ramirez Hassan         #
#             aramir21@eafit.edu.co                #
#            besmarter.team@gmail.com              #
#             www.besmarter-team.org               #
####################################################
setwd('C:/ANDRES/Portatil/GNCV/UI/UI GNSPD/Data')
data<-read.csv(file='DataProductionUSA.csv', header=TRUE, sep=',')
attach(data)
str(data)
Ws<-as.matrix(read.csv("Matrix.csv",header=TRUE,sep=",")) #This is the 19th matrix in file Matrices
#################SIMULATING y##################
N<-48                       #Number cross sectional units
T1<-17                      #Number time periods
b<-c(0.20,0.35,0.55,-0.007) #Population location parameters 
rho<-0.5                    #Spatial Autocorrelation disturbances
lambda<-0.4                 #Spatial Autocorrelation dependent variable
s_v<-0.0075                 #Standard deviation model
In<-diag(N)
A<-In-lambda*Ws             #Matrix A in paper
Ai<-solve(A)                #A inverse
B<-In-rho*Ws                #Matrix B in paper
Bi<-solve(B)                #B inverse
v<-rnorm(N*T1,0,s_v)        #iid disturbances
It<-diag(T1)                #Identity matrix dim T
e<-kronecker(It,Bi)%*%v     #Disturbances spatially correlated
aux<-kronecker(It,Ai)       #Auxiliar matrix to generate y
set.seed(12345)             #Seed
u<-colMeans(matrix(0.4*log(pcap)+rnorm(N*T1),T1,N)) #Fixed Effects
it<-rep(1,T1)               #Vector of ones
formula<-log(gsp)~log(pcap)+log(pc)+log(emp)+unemp #Formula 
mf <- model.frame(formula=formula, data=data)      #Data from formula
X<- as.matrix(mf[,-1])                             #Regressors
y<-aux%*%X%*%b+aux%*%kronecker(it,In)%*%u+aux%*%e  #Generate y
colnames(y)<-c("y")                                #Name 'y'
t<-rep(1:T1,each=N)                                #Time units
dat<-as.data.frame(cbind(id,t,y,pcap,pc,emp,unemp))                 #Simulated data frame
write.table(dat,file='Experimental Data Set.csv',sep=',',col.names=TRUE)  #Final simulated data set

