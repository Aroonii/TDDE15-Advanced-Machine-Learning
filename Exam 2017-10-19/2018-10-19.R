#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#install.packages("bnlearn")
#BiocManager::install("RBGL")
#BiocManager::install("Rgraphviz")
#BiocManager::install("gRain")


#################################### Assignment 1 Part 1###########################################################
### Graphical models
library(bnlearn)
library(gRain)

data(asia)
View(asia)

#Learn a BN from the Asia dataset. Learn both the structure and the parameters. Use any algorithm setting that
#you consider appropriate. Identify a d-separation in the BN learned and show that it indeed corresponds
# to an independence in the probability distribution represented by the BN. To do so you may want to use exact or
#approximate inference with the help of the bnlearn and gRain packages


### Learn the structure
set.seed(123)
BN_structure = hc(asia, score = "bde", iss = 10)

plot(BN_structure)
cpdag(BN_structure)
score(BN_structure, asia)
#plotting the structure we can see that b is independent from E is independent


### Learn the parameters
# Using bn.fit we get the conditional probability tables for the different nodes. All these conditional probabilities 
# define the joint probability. In conclusion we know have the joint probability but now need the marginal probability
# in order to be able to compute the marginal probability distribution for teh searched node
BN_fit = bn.fit(BN_structure, data = asia, method = "bayes")
BN_fit

# in order to make inference we have to modify the network to a grain object. This makes it possible for setEvidence and 
# querygrain to perform posterior inference via belief propagation. This only works for discrete networks. Also important
# that as.grain foes not allow for conditional probabilities to be Nan. This happens when estimating with maximum
# likelihood and some parent configurations are not observed in the data. In this case the as.grain will give a warning 
# and replace the nan with a uniform distribution, much like the baysian posterior would

fit_grain = as.grain(BN_fit)

#compile the network = create the junction tree and establishing the clique potentials. It is now an undirected MN
compiled_network = compile(fit_grain)
plot(compiled_network)

# Using querygrain we get the marginal distribution of a set of variable given evidence of set of variable
# Evidence is added by using setEvidence which extracts the probability tables for node(S) which can be in states x
# setEvidence allows for specification of hard evidence and likelihood evidence for variables. In this case we want to
# prove that there is a d-separation present on node s and that there is an independence in the probability distribution
# We therfor want to extract get the evidence (the probability table) for node s which can be in state x.

### vad ger querygrain?? marginal distribution?? Hur gör man approximate inference??
prop.table(table(querygrain(compiled_network)))
querygrain(compiled_network)

# In order to prove independence between two nodes and therefor prove the presens of d-separation we want to get the 
# conditional probabilities for Independent node | parentnode1, parentnode2, independent_node. For the same pair of 
# states in parentnodes the conditional probability should therfore be the same.


# nodes = parent nodes + comparing node, states = all combinations of states
hc7<-setFinding(compiled_network,nodes=c("S","T","E"),states=c("yes","yes","yes"))
querygrain(hc7,c("B"))
hc7<-setFinding(compiled_network,nodes=c("S","T","E"),states=c("yes","yes","no"))
querygrain(hc7,c("B"))
hc7<-setFinding(compiled_network,nodes=c("S","T","E"),states=c("yes","no","yes"))
querygrain(hc7,c("B"))
hc7<-setFinding(compiled_network,nodes=c("S","T","E"),states=c("yes","no","no"))
querygrain(hc7,c("B"))
hc7<-setFinding(compiled_network,nodes=c("S","T","E"),states=c("no","yes","yes"))
querygrain(hc7,c("B"))
hc7<-setFinding(compiled_network,nodes=c("S","T","E"),states=c("no","yes","no"))
querygrain(hc7,c("B"))
hc7<-setFinding(compiled_network,nodes=c("S","T","E"),states=c("no","no","yes"))
querygrain(hc7,c("B"))
hc7<-setFinding(compiled_network,nodes=c("S","T","E"),states=c("no","no","no"))
querygrain(hc7,c("B"))

#In the example above we can see that for the same combination of yes/no for the nodes S, T we get the same conditional
# probability for B. This implies that B is independent from E and that d-separation is present.


##################################Assignemt 1 Part 2###############################################################

#Compute approximately the fraction of the DAGs that are essential. An essential DAG is a DAG which is not Markov
#equivalent  to any other DAG. 

# to solve the problem of finding the fraction of graphs which coincide with the true graph we simply sample new random
#graphs for the set of noded and see how many of these coincide with the essential graph

#essential_graph = cpdag(BN_structure)
#plot(essential_graph)
#
#library(bnlearn)
#set.seed(123)
#ss<-50000
#x<-random.graph(c("A","B","C","D","E"),num=ss,method="melancon",every=50,burn.in=30000)
#
#y = unique(x)
#z = lapply(y, cpdag)
#
#r=0
#
#for(i in 1:length(y)) {
#  if(all.equal(y[[i]],z[[i]])==TRUE){
#    r<-r+1
#  }
#}
#length(y)/r

# By creating random samples of the graph and comparing if they correspont to the essential graph we can then compute
# the fraction
##################??? varför length(y)/r?


################################ Assignment 2 Part 1###################################################
###  Hidden markov models

#install.packages("HMM")
library(HMM)
startprobs = rep(1/100, 100)
States = 1:100
transProbs<-matrix(rep(0,length(States)*length(States)), nrow=length(States), ncol=length(States), byrow = TRUE)
for(i in 1:99){
  transProbs[i,i]<-.1
  transProbs[i,i+1]<-.9
}


# Symbols accoding to the probabilities we have which differ door/not door
Symbols<-1:2 
emissionProbs<-matrix(rep(0,length(States)*length(Symbols)), nrow=length(States), ncol=length(Symbols), byrow = TRUE)
for(i in States){
  if(i %in% c(10,11,12,20,21,22,30,31,32)){
    emissionProbs[i,] = c(0.9, 0.1)
  }
  else{
    emissionProbs[i,] = c(0.1, 0.9)
  }
}

# Initialize the hmm 
hmm<-initHMM(States,Symbols,startProbs,transProbs,emissionProbs)

# Add a unimodal sequence of observations. Observatios are referring to observations of what symbol is present
# in our case, door/not door. When theo model in this case see a long sequence of only door coming it
# therefore becomes sure that we must have passed the third doorand therfore continue to guess on segments
# after the third door

obs<-c(1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
alpha = exp(forward(hmm, obs))
pt = prop.table(alpha, margin = 2)

which.maxima<-function(x){ # This function is needed since which.max only returns the first maximum.
  return(which(x==max(x)))
}

apply(pt,2,which.maxima)



############################################Assignment 3##################################################
#library(kernlab)
#library(mvtnorm)
#
## From KernelCode.R: 
## Squared exponential, k
#k <- function(sigmaf = 1, ell = 1)  
#{   
#  rval <- function(x, y = NULL) 
#  {       
#    r = sqrt(crossprod(x-y))       
#    return(sigmaf^2*exp(-r^2/(2*ell^2)))     
#  }   
#  class(rval) <- "kernel"   
#  return(rval) 
#}  



k <- function(sigmaf = 1, ell = 1)  
{   
  rval <- function(x, y = NULL) 
  {       
    r = sqrt(crossprod(x-y))       
    return(sigmaf^2*exp(-r^2/(2*ell^2)))     
  }   
  class(rval) <- "kernel"   
  return(rval) 
}  




### Task a
# Let f~GP(0, k(x, x')) a priori and simulate 5 realizations from the prior distribution
# of f over the grid. We can solve this problem since we know the distribution
# of the resulting f* which will be multivariate nromally distributed with
# mean =  0 and sigma = k(x, x')

xGrid = seq(-1,1,0.1)

#Create the kernel for computing the covariance
# Simulating from the prior for ell = 0.2
kernel02 <- k(sigmaf = 1, ell = 0.2) # This constructs the covariance function
xGrid = seq(-1,1,by=0.1) 
# Compute the covariance matrix. Is necessary to do it this way since the 
# squared exponetial kernel in this case does not return a matrix
K = kernelMatrix(kernel = kernel02, xGrid, xGrid)

colors = list("black","red","blue","green","purple")
f = rmvnorm(n = 1, mean = rep(0,length(xGrid)), sigma = K)
plot(xGrid,f, type = "l", ylim = c(-3,3), col = colors[[1]])
for (i in 1:4){
  f = rmvnorm(n = 1, mean = rep(0,length(xGrid)), sigma = K)
  lines(xGrid,f, col = colors[[i+1]])
}

# Simulating from the prior for ell = 1
kernel1 <- k(sigmaf = 1, ell = 1) # This constructs the covariance function
xGrid = seq(-1,1,by=0.1) 
K = kernelMatrix(kernel = kernel1, xGrid, xGrid)

colors = list("black","red","blue","green","purple")
f = rmvnorm(n = 1, mean = rep(0,length(xGrid)), sigma = K)
plot(xGrid,f, type = "l", ylim = c(-3,3), col = colors[[1]])
for (i in 1:4){
  f = rmvnorm(n = 1, mean = rep(0,length(xGrid)), sigma = K)
  lines(xGrid,f, col = colors[[i+1]])
}


# compute the corr(f(0), f(0.1)) and corr(f(0), f(0.5)). This is the same 
# as computing the correlation between the input values for where f is evaluated
# l = 0.2
# kernel0.2
kernel02(0, 0.1)
kernel02(0, 0.5)


# l = 1
kernel1(0, 0.1)
kernel1(0, 0.5)

# Since we are using a sigmaf = 1 the correlation will be the same as the
# the covariance. Comparing the two scenarios the 2 points which are further
# away will have less correlation. When l is larges the decrease in correlation
# due to distance will not decac as rapidly. 


### Question 3b

# Compute the posterior distribution of f in the model 
# y ~ f(x) + epsilon, epsilon ~ N(0, 0.2^2)
load("GPdata.RData")

### ell = 0.2

sigmaNoise = 0.2

# Set up the kernel function
kernelFunc <- k(sigmaf = 1, ell = 0.2)

# Plot the data and the true 
plot(x, y, main = "", cex = 0.5)

GPfit <- gausspr(x, y, kernel = kernelFunc, var = sigmaNoise^2)
# Alternative: GPfit <- gausspr(y ~ x, kernel = k, kpar = list(sigmaf = 1, ell = 0.2), var = sigmaNoise^2)
xs = seq(min(x),max(x), length.out = 100)
meanPred <- predict(GPfit, data.frame(x = xs)) # Predicting the training data. To plot the fit.
lines(xs, meanPred, col="blue", lwd = 2)

# Compute the covariance matrix Cov(f)
n <- length(x)
Kss <- kernelMatrix(kernel = kernelFunc, x = xs, y = xs)
Kxx <- kernelMatrix(kernel = kernelFunc, x = x, y = x)
Kxs <- kernelMatrix(kernel = kernelFunc, x = x, y = xs)
Covf = Kss-t(Kxs)%*%solve(Kxx + sigmaNoise^2*diag(n), Kxs)

# Probability intervals for f
lines(xs, meanPred - 1.96*sqrt(diag(Covf)), col = "red")
lines(xs, meanPred + 1.96*sqrt(diag(Covf)), col = "red")

# Prediction intervals for y 
lines(xs, meanPred - 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "purple")
lines(xs, meanPred + 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "purple")


legend("topright", inset = 0.02, legend = c("data","post mean","95% intervals for f", "95% predictive intervals for y"), 
       col = c("black", "blue", "red", "purple"), 
       pch = c('o',NA,NA,NA), lty = c(NA,1,1,1), lwd = 2, cex = 0.55)



