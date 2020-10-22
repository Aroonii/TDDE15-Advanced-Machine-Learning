##Install the necessary 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#install.packages("bnlearn")
#BiocManager::install("RBGL")
#BiocManager::install("Rgraphviz")
#BiocManager::install("gRain")
library(bnlearn)
library(gRain)







set.seed(567)
data('asia')
ind <- sample(1:5000, 4000)
tr <- asia[ind,]
te <- asia[-ind,]

train = c(10,20,50,100,2000)
accuracy = c()
for (i in train){
  BN_bayes =  model2network("[S][A|S][T|S][L|S][B|S][E|S][X|S][D|S]")
  BN_fit_bayes = bn.fit(BN_bayes, data = tr[1:i,], method = "bayes")
  BN_fit_bayes = as.grain(BN_fit_bayes)
  BN_fit_bayes = compile(BN_fit_bayes)
  
  prediction_bayes = predict_test(BN_fit_bayes, 
                                  testData = te,
                                  observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                                  prediction_variable = c("S"))
  
  accuracy = append(accuracy, mean(prediction_bayes == te$S))
}



### Use an alternative model which is the opposite of the NB classifier

BN_bayes =  model2network("[A][X][T][B][D][E][L][S|A:X:T:B:D:E:L]")

accuracy_alt = c()
for (i in train){
  BN_bayes =  model2network("[A][X][T][B][D][E][L][S|A:X:T:B:D:E:L]")
  BN_fit_bayes = bn.fit(BN_bayes, data = tr[1:i,], method = "bayes")
  BN_fit_bayes = as.grain(BN_fit_bayes)
  BN_fit_bayes = compile(BN_fit_bayes)
  
  prediction_bayes = predict_test(BN_fit_bayes, 
                                  testData = te,
                                  observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                                  prediction_variable = c("S"))
  
  accuracy_alt = append(accuracy_alt, mean(prediction_bayes == te$S))
}



###################################Assignement 2 ###############################################
library(HMM)
set.seed(567)

startProbs = rep(1/10, 10)
states = 1:10
symbols = 1:10


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
  emissionProbs[i,i:min(i+2,max(states))] = 0.2
  emissionProbs[i,max(i-2, min(states)):i] = 0.2
  if(i<2){
    emissionProbs[i,9:10] = 0.2
  }
  if(i<3){
    emissionProbs[i,10] = 0.2
  }
  
  if(i > 8){
    emissionProbs[i,1] = 0.2
  }
  if(i > 9){
    emissionProbs[i,1:2] = 0.2
  }
}


hmm = initHMM(States = states, 
              Symbols = symbols, 
              startProbs = startProbs,
              transProbs = transProbs, 
              emissionProbs = emissionProbs)

set.seed(123)
HmmSimulation = simHMM(hmm = hmm, length = 100)
simObservations = HmmSimulation$observation
simStates = HmmSimulation$states



#compute alpha(z0) using x0 and z0 which are the emission probabilities and initial state probabilities
#p(x|z) = The probabilties ob observing a certain symbol. Is consequently computed by finding the emissionprob
#for a certain row which is set by the first observed value from the sequence and taking that whole row.
#p(z0) = for z0 this will simply be the inital starting probabilities since we don't know where we start
#alpha(z0) = p(x0|z0)*p(z0) This is then simply done by multiplying the emission probability row of the observed
#value and multiplying with the inital transition matrix and consequently gives us the probabilities for the 
#first time step. For computing the next columns of alpha we apply a similar technique but the transition
# is instead computed by for each for going through each state transition [,j] and multiplying by the previous time 
#steps alpha values and summing. each sum of each transition will then be the input for a row in the new alpha value
# for the next step. In this way we can forward compute the probability of being in a certain state while also taking
#into account the previous states.

a = matrix(data = 0, nrow = length(states), ncol = 100)
xz_0 = emissionProbs[simObservations[1],1:10]
z_0 = startProbs 
temp_xz*temp_z
a[1:10,1] = xz_0*z_0

for (j in 2:length(simObservations)){
  for (i in 1:10){
    p_xz = emissionProbs[simObservations[j], i]
    alpha = a[1:10, j-1]
    p_zt_1 = transProbs[1:10, i]
    a[i,j] = p_xz * sum(alpha*p_zt_1)
  }
}

## Alternative solution
#for(i in states){
#  p_z0 = startProbs[i]
#  p_xz = emissionProbs[simObservations[1], i]
#  a[i,1] = p_xz*p_z0
#}
#
#for (j in 2:length(simObservations)){
#  for(i in states){
#    p_xz = emissionProbs[simObservations[j],i]
#    summa = sum(t(transProbs[,i])*a[,j-1] )
#    a[i,j] = p_xz*summa
#  }
#}

alpha = prop.table(a, 2)
maxa = apply(alpha, 2, which.max)
table(maxa == simStates)
mean(maxa == simStates)

##### Backward algorithm
b = matrix(0,nrow = length(states), ncol = 100)

#computeB(zT)
#according to the algorithm we should set the last value to zero. This makes sense since all the last states should
#have equal probability of getting reached. 
b[,100] = 1

###compute b(zt)
#Simliar compuations as in previous step but since we are going backward we are instead looking by column for
#the emission matrix and by row for transition matrix. Since the sum is computed over all values we can take
# the whole vector immedeately on the emission matrix instead of doing value per value as in alpha. Since we
# are computing backwards we also have to make the multiplications backwards
for (j in 99:1){
  for (i in states){
    p_zt1 = transProbs[i,1:10]
    p_xt1 = emissionProbs[1:10,simObservations[j+1]]
    B_t1 = b[,j+1]
    b[i,j] = sum(p_zt1 * p_xt1 * B_t1)
  }
}

beta = prop.table(b,2)
maxa = apply(b, 2, which.max)
table(maxa == simStates)

#verify
beta = exp(backward(hmm, simObservations))
beta = prop.table(beta, 2)
maxa = apply(beta, 2, which.max)
table(maxa == simStates)


### Forward backward algorithm
forwardBackward = a*b
forwardBackward = prop.table(forwardBackward, 2)
maxa = apply(forwardBackward, 2, which.max)
table(maxa == simStates)

#Control forward backba
fb = posterior(hmm, simObservations)
maxa = apply(fb, 2, which.max)
table(maxa == simStates)



################################## Question 4#########################################
### Gaussian Processes

Matern32 <- function(sigmaf = 1, ell = 1) 
{
  rval <- function(x, y = NULL) {
    r = sqrt(crossprod(x-y));
    return(sigmaf^2*(1+sqrt(3)*r/ell)*exp(-sqrt(3)*r/ell))
  }
  class(rval) <- "kernel"
  return(rval)
} 

sigmaf = 1
ell = 0.5
zGrid = seq(0.01, 1, 0.01)

k = Matern32(sigmaf = sigmaf, ell = ell)
k_vals = c()
for (i in zGrid){
  k_vals = append(k_vals, k(0,i) )
}
plot(zGrid, k_vals, ylab = "f")

# Update sigma value
sigmaf = sqrt(0.5)
k = Matern32(sigmaf = sigmaf, ell = ell)
k_vals = c()
for (i in zGrid){
  k_vals = append(k_vals, k(0,i) )
}
plot(zGrid, k_vals, ylab = "f")

#As we lower the value of sigmaf from 1 to 0.5 we limit how large we allow the variance to be. Since the value is halfed the 
# resulting kernel value will then be limited to half of the previous value. sigma f has an impact of the amplitude of f where a lower value
# will result in a decreased amplitude meaning that the spread of the f values is decreased. As we can see in the plot, as the distance grows
# larger the corresponding covariane will be decreased and go towards zero. Sigmaf will determine the variance of f

lidarData <- read.table('https://raw.githubusercontent.com/STIMALiU/AdvMLCourse/master/GaussianProcess/Code/LidarData', 
                        header = T)
LogRatio <- lidarData$LogRatio
Distance <- lidarData$Distance
library(kernlab)
#Distance = scale(Distance)
#LogRatio = scale(LogRatio)

y = LogRatio
x<-Distance

xs<-seq(min(x), max(x), length.out = length(x))
xs = as.integer(xs)

GP_fit = gausspr(x = Distance,
                 y = LogRatio,
                 kernel = Matern32(sigmaf = 1, ell = 1),
                 var = 0.05^2)

predMean = predict(GP_fit, xs)

plot(Distance, LogRatio)
lines(xs, predMean, col = "red", lwd = 2)

matern = Matern32(sigmaf = 1, ell = 1)


n <- length(x)
Kss <- kernelMatrix(kernel = matern, x = xs, y = xs)
Kxx <- kernelMatrix(kernel = matern, x = x, y = x)
Kxs <- kernelMatrix(kernel = matern, x = x, y = xs)
Covf = Kss-t(Kxs)%*%solve(Kxx + 0.05^2*diag(n), Kxs) # Covariance matrix of fStar.


# Probability intervals for fStar.
lines(xs, predMean - 1.96*sqrt(diag(Covf)), col = "blue", lwd = 2)
lines(xs, predMean + 1.96*sqrt(diag(Covf)), col = "blue", lwd = 2)


lines(xs, predMean - 1.96*sqrt(diag(Covf) + 0.05^2), col = "orange", lwd = 2)
lines(xs, predMean + 1.96*sqrt(diag(Covf) + 0.05^2), col = "orange", lwd = 2)

###########################

GPPosterior = function(x, y, xs,kernel,...){
  k = kernel(...)
  n <- length(x)
  Kss <- kernelMatrix(kernel = k, x = xs, y = xs)
  Kxx <- kernelMatrix(kernel = k, x = x, y = x)
  Kxs <- kernelMatrix(kernel = k, x = x, y = xs)
  Covf = Kss-t(Kxs)%*%solve(Kxx + 0.05^2*diag(n), Kxs)
  return(diag(Covf))
}

# Do we have to scale or not ????
var = GPPosterior(scale(x), scale(y), scale(xs), Matern32, 1, 1)
lines(xs, predMean + 1.96*sqrt(var*sd(LogRatio^2)))
lines(xs, predMean - 1.96*sqrt(var*sd(LogRatio^2)))

lines(xs, predMean + 1.96*sqrt(var*sd(LogRatio^2) + 0.05^2))
lines(xs, predMean - 1.96*sqrt(var*sd(LogRatio^2) + 0.05^2))

#fast iof så borde väl y värdet täcka alla data punkter 