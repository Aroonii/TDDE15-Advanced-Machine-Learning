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

HmmSimulation = simHMM(hmm = hmm, length = 100)
simObservations = HmmSimulation$observation
simStates = HmmSimulation$states

a = matrix(data = 0, nrow = length(states), ncol = 100)

#compute alpha(z0) using x0 and z0 which are the 
#compute p(x0|z0)*p(z0)
for(i in states){
  p_z0 = startProbs[i]
  p_xz = emissionProbs[simObservations[1],i]
  a[i,1] = p_xz*p_z0
}

for (j in 2:length(simObservations)){
  
  for(i in states){
    p_xz = emissionProbs[simObservations[j],i]
    summa = sum(t(transProbs[,i])*a[,j-1] )
    a[i,j] = p_xz*summa
  }
  
}

alpha = prop.table(a, 2)
maxa = apply(alpha, 2, which.max)
table(maxa == simStates)
mean(maxa == simStates)











