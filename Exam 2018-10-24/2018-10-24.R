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

tr1 = tr[1:10,] 
tr2 = tr[1:20,]
tr3 = tr[1:50,]
tr4 = tr[1:100,]
tr5 = tr[1:1000,]
tr6 = tr[1:2000,]

####1
BN_bayes =  model2network("[S][A|S][T|S][L|S][B|S][E|S][X|S][D|S]")
BN_fit_bayes = bn.fit(BN_bayes, data = tr1, method = "bayes")
BN_fit_bayes = as.grain(BN_fit_bayes)
BN_fit_bayes = compile(BN_fit_bayes)

prediction_bayes = predict_test(BN_fit_bayes, 
                                testData = te,
                                observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                                prediction_variable = c("S"))
mean(prediction_bayes==te$S)

### 2
BN_bayes =  model2network("[S][A|S][T|S][L|S][B|S][E|S][X|S][D|S]")
BN_fit_bayes = bn.fit(BN_bayes, data = tr2, method = "bayes")
BN_fit_bayes = as.grain(BN_fit_bayes)
BN_fit_bayes = compile(BN_fit_bayes)


prediction_bayes = predict_test(BN_fit_bayes, 
                                testData = te,
                                observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                                prediction_variable = c("S"))
mean(prediction_bayes==te$S)

## 3
BN_bayes =  model2network("[S][A|S][T|S][L|S][B|S][E|S][X|S][D|S]")
BN_fit_bayes = bn.fit(BN_bayes, data = tr3, method = "bayes")
BN_fit_bayes = as.grain(BN_fit_bayes)
BN_fit_bayes = compile(BN_fit_bayes)


prediction_bayes = predict_test(BN_fit_bayes, 
                                testData = te,
                                observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                                prediction_variable = c("S"))
mean(prediction_bayes==te$S)

## 4
BN_bayes =  model2network("[S][A|S][T|S][L|S][B|S][E|S][X|S][D|S]")
BN_fit_bayes = bn.fit(BN_bayes, data = tr4, method = "bayes")
BN_fit_bayes = as.grain(BN_fit_bayes)
BN_fit_bayes = compile(BN_fit_bayes)


prediction_bayes = predict_test(BN_fit_bayes, 
                                testData = te,
                                observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                                prediction_variable = c("S"))
mean(prediction_bayes==te$S)

##5
BN_bayes =  model2network("[S][A|S][T|S][L|S][B|S][E|S][X|S][D|S]")
BN_fit_bayes = bn.fit(BN_bayes, data = tr5, method = "bayes")
BN_fit_bayes = as.grain(BN_fit_bayes)
BN_fit_bayes = compile(BN_fit_bayes)


prediction_bayes = predict_test(BN_fit_bayes, 
                                testData = te,
                                observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                                prediction_variable = c("S"))
mean(prediction_bayes==te$S)

##6
BN_bayes =  model2network("[S][A|S][T|S][L|S][B|S][E|S][X|S][D|S]")
BN_fit_bayes = bn.fit(BN_bayes, data = tr6, method = "bayes")
BN_fit_bayes = as.grain(BN_fit_bayes)
BN_fit_bayes = compile(BN_fit_bayes)


prediction_bayes = predict_test(BN_fit_bayes, 
                                testData = te,
                                observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                                prediction_variable = c("S"))
mean(prediction_bayes==te$S)


predict_test = function(BN_model, testData, observation_nodes, prediction_variable){
  prediction = rep(0, length(testData)) 
  feature_data = testData[c(observation_nodes)]
  for(i in 1:nrow(testData)){
    node_state = t(feature_data[i,])

    evidence = setEvidence(object = BN_model, 
                           nodes = observation_nodes,
                           states = node_state)
    
    conditional_distribution = querygrain(object = evidence, 
                                          nodes = prediction_variable)$S
    
    prediction[i] = if(conditional_distribution["yes"] >= 0.5) "yes" else "no"
  } 
  return(prediction)
}

#####################Reverse and change the relationship########################


BN_bayes =  model2network("[A][X][T][B][D][E][L][S|A:X:T:B:D:E:L]")
BN_fit_bayes = bn.fit(BN_bayes, data = tr6, method = "bayes")
BN_fit_bayes = as.grain(BN_fit_bayes)
BN_fit_bayes = compile(BN_fit_bayes)


prediction_bayes = predict_test(BN_fit_bayes, 
                                testData = te,
                                observation_nodes = c("A", "T", "L", "B", "E", "X", "D"),
                                prediction_variable = c("S"))
mean(prediction_bayes==te$S)



