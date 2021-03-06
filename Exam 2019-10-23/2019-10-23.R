#########################################Assignment 1 ####################################################
### Graphical Models

library(kernlab)
library(gRain)
library(bnlearn)

monty_network = model2network("[D][P][M|D:P]")
cptD = matrix(c(1/3, 1/3, 1/3), ncol = 3, dimnames = list(NULL, c("D1", "D2", "D3")))
cptP = matrix(c(1/3, 1/3, 1/3), ncol = 3, dimnames = list(NULL, c("P1", "P2", "P3")))
cptM = c(
  0,.5,.5,
  0,0,1,
  0,1,0,
  0,0,1,
  .5,0,.5,
  1,0,0,
  0,1,0,
  1,0,0,
  .5,.5,0)
dim(cptM) = c(3,3,3)
dimnames(cptM) = list("M" = c("M1", "M2", "M3"), "D" =  c("D1", "D2", "D3"), "P" = c("P1", "P2", "P3"))
MHfit = custom.fit(monty_network, list(D = cptD, P = cptP, M = cptM))
MHcom = compile(as.grain(MHfit))


### Approximate Inference
#Applied on the bn.fit object created from custom.fit. Will generate random samples conditional on the evidence using the method specified in the 
#method argument. Can then create a table in order to verify what conditional probability is the highest
table(cpdist(MHfit, nodes = "P", evidence = TRUE))
table(cpdist(MHfit, nodes = "P", (D=="D1" & M == "M2")))
table(cpdist(MHfit, nodes = "P", (D=="D1" & M == "M3")))

### Exact inference
MHfitEv = setFinding(MHcom, nodes = c(""), states = c(""))
querygrain(MHfitEv, c("P"))
MHfitEv<-setFinding(MHcom,nodes=c("D","M"),states=c("D1","M2"))
querygrain(MHfitEv,c("P"))
MHfitEv<-setFinding(MHcom,nodes=c("D","M"),states=c("D1","M3"))
querygrain(MHfitEv,c("P"))


### Part 2

XOR_BN = model2network("[A][B][C|A:B]")
cptA = matrix(c(1/2, 1/2), ncol = 2, dimnames = list(NULL, c("0", "1")))
cptB = matrix(c(1/2, 1/2), ncol = 2, dimnames = list(NULL, c("0", "1")))
#cptC = matrix(c(.5,.5,.5,.5,.5,.5,.5,.5), ncol = 4, dimnames = list(NULL, c("0", "1", "1", "0")))
cptC = c(
  1, 0,
  0, 1,
  0, 1,
  1, 0)
  
dim(cptC) = c(2,2,2)
dimnames(cptC) = list("C" = c("0", "1"), "A" = c("0", "1"), "B" = c("0", "1"))
XOR_fit = custom.fit(XOR_BN, list(A = cptA, B = cptB, C = cptC))
XOR_comp = compile(as.grain(XOR_fit))
XOR_samples = rbn(XOR_fit, 1000)


for(i in 1:10){
  # samples = sample(1:1000, 100, replace = FALSE)
  XOR_samples = rbn(XOR_fit, 1000)
  hc1 = hc(XOR_samples, restart = 0)
  plot(hc)
}

# Hc algorithms is not asymtotically correct under faithfullness => There is no guarantee that we
# find the true structure even with unlimited data.  In this case the HC algotihm faisl because 
# the data comes from a distribution that is not fiathfull to the graph meaning that
# it does not show the same independences as the graph. The graph has an edge from A->C but looking at the data we see that C is
# marginally independent from A. The same goes for A,B and C. Consequently the hc algorithm will not be able 
# to find any dependent variables and link them with the edge. The problem lies in that the hc alogrithm
# only addes 1 adge at the time. 

########################################Assignment 3############################################
###HMM

library(HMM)
set.seed(567)

startProbs = rep(1/10, 10)
states = 1:10
symbols = 1:11


transProbs<-matrix(rep(0,length(states)*length(states)), nrow=length(states), ncol=length(states), byrow = TRUE)
for(i in 1:10){
  if(i %% 10 == 0){
    transProbs[i,1] = 0.5
    transProbs[i,i] = 0.5
  }
  else{
    transProbs[i,i]<-.5
    transProbs[i,i+1]<-.5
  }
}

emissionProbs = matrix(rep(0,length(states)*length(symbols)), nrow=length(states), ncol =length(symbols), byrow= TRUE)
for(i in states){
  emissionProbs[i,i:min(i+2,max(states))] = 0.1
  emissionProbs[i,max(i-2, min(states)):i] = 0.1
  emissionProbs[i,11] = 0.5
  if(i<2){
    emissionProbs[i,9:10] = 0.1
  }
  if(i<3){
    emissionProbs[i,10] = 0.1
  }
  
  if(i > 8){
    emissionProbs[i,1] = 0.1
  }
  if(i > 9){
    emissionProbs[i,1:2] = 0.1
  }
}


hmm = initHMM(States = states, 
              Symbols = symbols, 
              startProbs = startProbs,
              transProbs = transProbs, 
              emissionProbs = emissionProbs)

obs_sequence = c(1,11,11,11)

viterbi = viterbi(hmm, obs_sequence)
smoothed = posterior(hmm, obs_sequence)
smoothed = apply(smoothed, 2, which.max)
smoothed



#################################Assignment 4#########################################################
### Gaussian processes

tempData = read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/
Code/TempTullinge.csv", header=TRUE, sep=";")


y = tempData$temp
x = 1:length(y)

subset = seq(1, length(time), 5)
x = x[subset]
y = y[subset]



SEKernel2 <- function(par=c(20,0.2),x1,x2){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- (par[1]^2)*exp(-0.5*( (x1-x2[i])/par[2])^2 )
  }
  return(K)
}

LM <- function(par=c(20,0.2),X,y,k,sigmaNoise){
  n <- length(y)
  L <- t(chol(k(par,X,X)+((sigmaNoise^2)*diag(n))))
  a <- solve(t(L),solve(L,y))
  logmar <- -0.5*(t(y)%*%a)-sum(diag(L))-(n/2)*log(2*pi)
  return(logmar)
}


regression_fit = lm(scale(y) ~ scale(x) + I(x)^2)
sigmaNoise = sd(regression_fit$residuals)

foo<-optim(par = c(1,0.1), fn = LM, X=scale(x),y=scale(y),k=SEKernel2,sigmaNoise=sigmaNoise, method="L-BFGS-B",
           lower = c(.Machine$double.eps, .Machine$double.eps),control=list(fnscale=-1))

foo$par


####################################### Via grid search ####################################################
bestLM = LM(par = c(20,0.2), X = scale(x), y = scale(y), k = SEKernel2, sigmaNoise)
bestLM
besti = 20
bestj = 0.2
for(i in seq(1,50,1)){
  for(j in seq(0.1,10,0.1)){
    accuracy = LM(par = c(i,j), X = scale(x), y = scale(y), k = SEKernel2, sigmaNoise = sigmaNoise)
    if(bestLM < accuracy){
      bestLM = accuracy
      besti = i
      bestj = j
    }
  }
}

bestLM
besti
bestJ



################################################################################
library(kernlab)
data <-
  read.csv(
    "https://github.com/STIMALiU/AdvMLCourse/raw/master/
GaussianProcess/Code/banknoteFraud.csv",
    header = FALSE,
    sep = ","
  )



set.seed(111)
selectTraining = sample(1:dim(data)[1], size = 1000, replace = FALSE)
y = as.factor(data[,5]) #predictor
x = as.matrix(data[,1:4]) # covariates
yTrain = y[selectTraining]
yTest = y[-selectTraining]
xTrain = x[selectTraining,]
xTest = x[-selectTraining,]
#From the training data use 200 points for validation
selectVal = sample(1:1000, size = 200, replace = FALSE)
yVal = yTrain[selectVal]
xVal = xTrain[selectVal,]
yTrain = yTrain[-selectVal]
xTrain = xTrain[-selectVal,]




acVal = function(par = c(0.1)){
  GP.fit_fraud = gausspr(x = xTrain, y = yTrain, kernel = "rbfdot", kpar = list(sigma=par[1]))
  predVal = predict(GP.fit_fraud, xVal)
  table(predVal, yVal)
  accuracyVal = sum(predVal == yVal)/length(yVal)
  return(accuracyVal)
}

selVars = c(1,2,3,4)
GP.fit_fraud = gausspr(x = xTrain, y = yTrain, kernel = "rbfdot", kpar = 'automatic')
GP.fit_fraud
predVal = predict(GP.fit_fraud, xVal)
table(predVal, yVal)
accuracyVal = mean(predVal == yVal)


#GridSearch for classification problem
bestVal = accuracyVal
for ( j in seq(0.1, 10, 0.1)){
  accuracy = acVal(j)
  if(bestVal < accuracy){
    bestVal = accuracy
    bestJ = j
  } 
}

GP.fit_fraud = gausspr(x = xTrain, y = yTrain, kernel = "rbfdot", kpar = list(sigma=bestJ))
predTest = predict(GP.fit_fraud, xTest[,selVars])
table(predTest, yTest)
mean(predTest == yTest)



########################## Alternative solution using optim#############################
foo = optim(par = c(0.1), fn = acVal, method = "L-BFGS-B", lower = c(.Machine$double.eps), control = list(fnscale=-1))
acVal(foo$par)

GP.fit_fraud = gausspr(x = xTrain[,selVars], y = yTrain, kernel = "rbfdot", kpar = list(sigma=foo$par))
predTest = predict(GP.fit_fraud, xTest[,selVars])
table(predTest, yTest)
mean(predTest == yTest)
