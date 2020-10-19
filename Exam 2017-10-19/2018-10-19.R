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

essential_graph = cpdag(BN_structure)
plot(essential_graph)

library(bnlearn)
set.seed(123)
ss<-50000
x<-random.graph(c("A","B","C","D","E"),num=ss,method="melancon",every=50,burn.in=30000)

y = unique(x)
z = lapply(y, cpdag)

r=0

for(i in 1:length(y)) {
  if(all.equal(y[[i]],z[[i]])==TRUE){
    r<-r+1
  }
}
length(y)/r

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

# Varför blir emmission probabilities på detta sättet? Varför blev de som de blev i labben
Symbols<-1:2 # 1=door
emissionProbs<-matrix(rep(0,length(States)*length(Symbols)), nrow=length(States), ncol=length(Symbols), byrow = TRUE)
