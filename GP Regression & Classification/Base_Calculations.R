############################################Compute mean and covariance######################################################
GPPosterior = function(x, y, xs, sigmaNoise, kernel,...){
  k = kernel(...)
  n <- length(x)
  Kss <- kernelMatrix(kernel = k, x = xs, y = xs)
  Kxx <- kernelMatrix(kernel = k, x = x, y = x)
  Kxs <- kernelMatrix(kernel = k, x = x, y = xs)
  Covf = Kss-t(Kxs)%*%solve(Kxx + sigmaNoise^2*diag(n), Kxs)
  return(diag(Covf))
}

#Alternative way usign algorithm 2.1
PosteriorGP = function(X, y, XStar, sigmaNoise, k,...){
  k = k(...)
  n = length(X)
  L = t(chol(k(X,X) + ((sigmaNoise^2)*diag(n))))
  a = solve(t(L), solve(L,y))
  kStar = k(X, XStar)
  mean = t(kStar)%*%a
  v = solve(L, kStar)
  var = k(XStar, XStar) - (t(v)%*%v)
  return(list("mean" = mean, "var" = var))
}

########## Compute covariance immedeately########################################
# If not scaling any data this can then be used to immedeately compute compute the prediction bands with the mean from gauuspr
xs = seq(min(x),max(x), length.out = length(x))
n <- length(x)
Kss <- kernelMatrix(kernel = kernelFunc, x = xs, y = xs)
Kxx <- kernelMatrix(kernel = kernelFunc, x = x, y = x)
Kxs <- kernelMatrix(kernel = kernelFunc, x = x, y = xs)
Covf = Kss-t(Kxs)%*%solve(Kxx + sigmaNoise^2*diag(n), Kxs)



















###########################Hyperparameter Learning por Classification Problem#######################################################
#Split data into train, validation and Test

#800 points for training, 200 points for validation and the rest for test
set.seed(111)
selectTraining = sample(1:dim(data)[1], size = 1000, replace = FALSE)
y = data[,5] #predictor
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


#function for computing the accuracy, the correctly classified points
# funciton has 1 input parameter whihc is the input for the rbfdot kernel
acVal = function(par = c(0.1)){
  GP.fit_fraud = gausspr(x = xTrain, y = yTrain, kernel = "rbfdot", kpar = list(sigma=par[1]))
  predVal = predict(GP.fit_fraud, xVal)
  table(predVal, yVal)
  accuracyVal = sum(predVal == yVal)/length(yVal)
  return(accuracyVal)
}

selVars = c(1,2,3,4)
GP.fit_fraud = gausspr(x = xTrain, y = yTrain, kernel = "rbfdot", 'automatic', scaled = FALSE, type = "classification")
GP.fit_fraud
predVal = predict(GP.fit_fraud, xVal[,selVars])
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

GP.fit_fraud = gausspr(x = xTrain[,selVars], y = yTrain, kernel = "rbfdot", kpar = list(sigma=bestJ))
predTest = predict(GP.fit_fraud, xTest[,selVars])
table(predTest, yTest)
mean(predTest == yTest)

## Alternative solution using optim
foo = optim(par = c(0.1), fn = acVal, method = "L-BFGS-B", lower = c(.Machine$double.eps), control = list(fnscale=-1))
acVal(foo$par)

GP.fit_fraud = gausspr(x = xTrain[,selVars], y = yTrain, kernel = "rbfdot", kpar = list(sigma=foo$par))
predTest = predict(GP.fit_fraud, xTest[,selVars])
table(predTest, yTest)
mean(predTest == yTest)


#################################################Learning HyperParameters for regression problem########################################
# Similar as before but we instead compute for example the residual sum and minimize this. Another way could be to implement the logmar
# which gives us a weighted value

#Implemented in PosteriorGP algorithm 2.1
PosteriorGP = function(X, y, XStar, sigmaNoise, k){
  n = length(X)
  L = t(chol(k(X,X) + ((sigmaNoise^2)*diag(n))))
  a = solve(t(L), solve(L,y))
  kStar = k(X, XStar)
  mean = t(kStar)%*%a
  v = solve(L, kStar)
  var = k(XStar, XStar) - (t(v)%*%v)
  logmar -0.5*(t(y)%*%a)-sum(diag(L)) - (n/2)*log(2*pi)
  return(list("mean" = mean, "var" = var, "logmar" = logmar))
}

#Only the logmarginal part of algorithm 2.1 and finding the best hyperparameter
#Could alternatively also include sigma as a search parameter and same for choice of kernel. Want to choose the one which maximizes the
#marginal likelihood
LM = function(X, y, sigmaNoise, k,...){
  k = k(...)
  n = length(X)
  L = t(chol(k(X,X) + ((sigmaNoise^2)*diag(n))))
  a = solve(t(L), solve(L,y))
  logmar -0.5*(t(y)%*%a)-sum(diag(L)) - (n/2)*log(2*pi)
}

bestLM = LM(X = scale(X), y = scale(y), sigmaNoise = sigmaNoise, k = SEkernel, par = c(20,0.2))
bestLMs
besti = 20
bestj = 0.2
for(i in seq(1,50,1)){
  for(j in seq(0.1,10,0.1)){
    accuracy = LM(X = scale(X), y = scale(y), sigmaNoise = sigmaNoise, k = SEkernel, par = c(i,j))
    if(bestLM < accuracy){
      bestLM = accuracy
      besti = i
      bestj = j
    }
  }
}

#Alternatively using optim
foo = optim(par = c(1,0.1), fn = LM, X = scale(x), y = scale(y), sigmaNoise = sigmaNoise, k = SEkernel,
            method = "L-BFGS-B",lower = c(.Machine$double.eps, .Machine$double.eps), control = list(fnscale = -1))
