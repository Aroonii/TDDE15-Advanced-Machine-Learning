###########################Graphical models##################################
### Find fraction of BN whose independence model can be represented with a MN
# Separators in a BN will only be able to be represented if there are no 
# unshielded colliders. If there is an unshielded collider then the moral graph
# will add another 

library(bnlearn)
library(gRain)
set.seed(567)

ss<-5000
x<-random.graph(c("A","B","C","D","E"),num=ss,method="melancon")

x = unique(x)
skeleton = lapply(x, skeleton)
moral = lapply(x, moral)

r = 0
for(i in 1:length(x)) {
  if(all.equal(skeleton[[i]],moral[[i]])==TRUE){
    r<-r+1
  }
}

fraction = r/length(x)
fraction




###############################Hidden Markov Model################################

states = c("SS", "SR", "RS", "RR")
symbols = c("S", "R")
startProbs = rep(1/4, 4)


transitionProbs = matrix(0, nrow = length(states), ncol = length(states), byrow = TRUE)
transitionProbs[1,1] = 0.75
transitionProbs[1,2] = 0.25
transitionProbs[2,3:4] = 0.5
transitionProbs[3,1:2] = 0.5
transitionProbs[4,3] = 0.25
transitionProbs[4,4] = 0.75

emissionProbs = matrix(0, nrow = length(states), ncol = length(symbols), byrow = TRUE)
emissionProbs[1:4,1] = 0.9
emissionProbs[1:4,2] = 0.1

hmm<-initHMM(states,symbols,startProbs,transitionProbs,emissionProbs)
simulation = simHMM(hmm, 10)
simulation$states
simulation$observation

