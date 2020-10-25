##Install the necessary 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#install.packages("bnlearn")
#BiocManager::install("RBGL")
#BiocManager::install("Rgraphviz")
#BiocManager::install("gRain")
library(bnlearn)
library(gRain)

#Read the data to be used
data(asia)
cat("Data structure: ", colnames(asia))
# Variables = A S T L B E X D
variables = colnames(asia)
set.seed(12345)
#Task 1: Run the hill climbing algorithm on the data. This will return a Bayesian
#network structure of class bn. By modifying the parameters score and restart
# we will examine how the bayesian network is affected.

###General network
network = hc(asia)
#Ways of examining the network
plot(network)
arcs(network)
vstructs(network)
cpdag(network)

#Network 1
network1 = hc(asia, restart = 1, score = "bic")
plot(network1)
cpdag(network1)
cat("Score for network 1, restart = 1, score = bci")
score(network1, asia)

#Network 2
network2 = hc(asia, restart = 100, score = "bic")
plot(network2)
cpdag(network2)
score(network2, asia)

print(all.equal(network1, network2))
# By changing the number of restarts we in general avoid getting stuck into a 
# local optimum. By increasing the amount of restart the amount of tests is 
# greatly increased from 91 to 1832. Since the score is however the same 
# then the two classes are score equivalent. This means that the scoring criteria assign the same score
# to equivalent structures. In this example, meerly changing the restart does therfore not create non 
# equivalent BN in this case

#Network 3 
network3 = hc(asia, score = "aic", restart = 100)
plot(network3)
cpdag(network3)
cat("BN loglike, restart = 100 => Score: ", 
score(network3, asia))

#Network 4 
network4 = hc(asia, score = "bic", restart = 100)
plot(network4)
score(network4, asia)
cpdag(network4)

#Comparing network 3 and 4
all.equal(network3, network4)

#Both in the plot and through the above all.equal we see that there are
# 2 different models which have been created where they have a different
#number of directed and undirected arcs. We can moreover verify that
#There are non-equivalent classes created since their scores differ.
# looking at the score we we prefer the bic since this gives a less negative score
# looking at the graph we can also see that A is independent in the bic network and since the 
# hill climbing algorithm does not produce fake independences we know that we can trust this and is therfore prefereed
# since is gives us less parameters to estimate



#Task 2
# Divide data into train/test and use it to learn the structure and parameters
# of a network. Use the created BN to classify the remaining test dataset into two classes
# S = yes and S = no. Using a scoring algorithm this results in computing the posterior probability dsitribution for S
# for each class and classify the most likely class. 

### Functions to be used
#bn.fit
#as.grain

###Functions for approximate inference
#prop.table
#table
#cpdist

### gRain package
# compile
# setFinding
# querygrain

set.seed(12345)
n = nrow(asia)
id = sample(n, floor(n*0.8))
train = asia[id,]
test = asia[-id,]

#Learn BN structure using parameters of your choice
set.seed(12345)
bn_network = hc(train, score = "bic", restart = 100)
plot(bn_network, main = "learned network")

#Create the true BN for the Asia dataset
true_network = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]", ordering = variables)
plot(true_network, main = "True network")

#Fit the parameters of a bayesian notwork conditional on its structure
bn_fit = bn.fit(bn_network, 
                data = train)
true_bn_fit = bn.fit(true_network, 
                     data = train)

library(gRain)
# In order to make inference we have to convert the fit to a BN and compile it
#As.grain conversts the bn.fit object to a grain. The grain class stores a fitted bayesian network as a list of conditional
# probability tables. This makes it possible for setEvidence and querygrain to perform posterior inference via belief
# propagation. The grain does not allow conditional probabilities to be NaN. This happens when estimationg them via a 
# maximum likelihood and some parent ocnfigurations are not observed in the data. Grain will then replace NaN with uniform 
# distribution 
bn_fit_gRain = as.grain(bn_fit)
bn_fit_gRain
#For the true network
bn_fit_true = as.grain(true_bn_fit)
bn_fit_true


#Compile the BN which means creating teh junction tree and establishing clique potentials
junction_tree = compile(bn_fit_gRain)
junction_tree
#For the true network
junction_tree_true = compile(bn_fit_true)

   
   #Using setEvidende we extract the conditional probability tables for node(s) which can be in states x
   # setEvidence allows to specifivction of hard evidence and likelihood evidence for variables.
   # It attache a nodeto all the 
   #(bn_fit, nodes = asia[-2], states = "S")
   
   #query an independence network i.e obtain the conditional distribution of a set of variables given evidnce on other
   # variables (in most cases)
querygrain(junction_tree)

#prop.table(table(querygrain(junction_tree)))

# Make predictin on the test data 

predict_test = function(BN_model, testData, observation_nodes, prediction_variable){
  prediction = rep(0, length(testData))
  featureData = testData[,c(observation_nodes)]
   for(i in 1:nrow(testData)){
   node_state = t(featureData[i,])
   #for (k in observation_nodes){
   #  node_state[k] = if(testData[i,k] == "yes") "yes" else "no"
   #}
   
   
   #Assign the new data points for each observation to a node
   # observation nodes should be all the nodes except the conditional node we are doing inference on
   # In this case this is S
   evidence = setEvidence(object = BN_model, 
                          nodes = observation_nodes,
                          states = node_state)
   
   #obtain the conditional distribution for a node given the evidence obtained above. In this case we get the conditional
   #distribution for S which is the prediction variable
   conditional_distribution = querygrain(object = evidence, 
                                         nodes = prediction_variable)$S
   
   #For each observation (row in the data set) we classify the prediction variable as yes/no depending on what
   #probability is larger
    
   prediction[i] = if(conditional_distribution["yes"] > 0.5) "yes" else "no"
 } 
  return(prediction)
}


var = c("A", "T", "L", "B", "E", "X", "D")

test_prediction = predict_test(BN_model = junction_tree, 
             testData = test, 
             observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
             prediction_variable = c("S"))

table(test_prediction, test$S)

#For the true graph
test_prediction_true = predict_test(BN_model = junction_tree_true, 
                               testData = test, 
                               observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                               prediction_variable = c("S"))

table(test_prediction_true, test$S)




### Task 3####
# The same excercise as before but classify the S only for the MArakov Blanket of S. i.e its parents plus its children plus the
# parents of of its children minus S itseld. Reports again the confusion matrix

#Useful function
# mb => miscelalaneous utilities: assign or extract various quantities from an object of class bn of bn.fit

# Extract Marakov blaket from the fitted network structure generated by the hill climb algorithm as well as
# the fitted network structure created from the true model

mb_hc = mb(x = bn_fit, 
           node = c("S"))

mb_true = mb(x = true_bn_fit,
             node = c("S"))

#Make the predictions and create the new confusion matrices
test_prediction_mb = predict_test(BN_model = junction_tree, 
                               testData = test, 
                               observation_nodes = mb_hc,
                               prediction_variable = c("S"))

table(test_prediction_mb, test$S)
mean(test_prediction_mb != test$S)

#For the true graph
test_prediction_true_mb = predict_test(BN_model = junction_tree_true, 
                                    testData = test, 
                                    observation_nodes = mb_true,
                                    prediction_variable = c("S"))

table(test_prediction_true, test$S)
mean(test_prediction_true_mb != test$S)

#Using the markov blanket we only make computations using the data belonging to the nodes B and L since they are
# the parents/children. We see that the resulting matrix is the same and this makes sense. By looking at the graph we can
# see that the node S is independent from all the nodes given B, L. Adding thhe other nodes would therfore only create
# more computation but would not gain any moreinformation. 



### TASK 4###
#Repeat the previous excercise but insted using a
#bayes classifier. We are assuming that the predictive
#variables are independent given the class variable
#Model the NB classifier as a BN

# Create the naive bayes network structur under the asumption
# That all predictor variables are independent given
# S. 
BN_bayes =  model2network("[S][A|S][T|S][L|S][B|S][E|S][X|S][D|S]")

# Fit the parameters to the network
BN_fit_bayes = bn.fit(BN_bayes, data = train)
plot(BN_bayes)

#Convert the BN to a gRain object and compile
# to a junction tree. 
BN_fit_bayes_grain = as.grain(BN_fit_bayes)
junction_tree_bayes = compile(BN_fit_bayes_grain)

#Classify the test set using the naive baye classifier
prediction_bayes = predict_test(junction_tree_bayes, 
                                testData = test,
                                observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                                prediction_variable = c("S"))

table(Observartion = test$S, Classification = prediction_bayes)
mean(test$S != prediction_bayes)
# Results in a missclassificartio rate of 0.334 

prediction_true = predict_test(junction_tree_true, 
                                testData = test,
                                observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                                prediction_variable = c("S"))

table(Observartion = test$S, Classification = prediction_true)
mean(test$S != prediction_true)
# Looking at the true structure we know that the predictor varibles are actually
# not independent given the class variables. Since this asumption is wrong it will
# therfore not classify as well as the true network structure