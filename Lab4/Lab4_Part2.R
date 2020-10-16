tempData = read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/
Code/TempTullinge.csv", header=TRUE, sep=";")

library(kernlab)

time = seq(1, nrow(tempData))

day = c()
counter = 1
for (i in 1:nrow(tempData)){
  
 if(counter > 365){
   counter = 1
 }
  
 day = append(day, counter)
 counter = counter +1
  
}


five_sequence = seq(1, nrow(tempData), 5)

time_selection = time[five_sequence]
day_selection = day[five_sequence]

temperature = tempData$temp

temperature_selection = temperature[five_sequence]


### 2.1
# Define your own square exponential kernel function (with parameters ` (ell) and 
# σf (sigmaf)), evaluate it in the point x = 1, x′ = 2, and use the kernelMatrix function
# to compute the covariance matrix K(X, X∗) for the input vectors X = (1, 3, 4)
# T and X∗ = (2, 3, 4)T

x = 1
x_prime = 2
X = c(1,3,4)
X_prime = c(2,3,4)

SE_Kernel <- function(sigmaF=1,l=1){
  rval = function(x, y = NULL){
    n1 <- length(x)
    n2 <- length(y)
    #K <- matrix(NA,n1,n2)
    #for (i in 1:n2){
     # for(j in 1:n1){
     #   K[j,i] <- sigmaF^2*exp(-0.5*( (x[i]-y[j])/l)^2 )
     # }
    #}
    #return(K)
    res = sigmaF^2*exp(-0.5*( (x-y)/l)^2 )
    return(res)
    
  }
  class(rval) = "kernel"
  return(rval)
  
}

## Compute the covariance for (x, x')
# Initialize the kernel, values ell and sigmaF = 1
kernel = SE_Kernel()

# Evalute the kernel in x = 1 and x' = 2
covariance = kernel(x = 1, y = 2)
covariance

# Compute the covariance matrix for the input vectors X, X_prime K[X, X_prime]
# X = c(1,3,4)
# X_prime = c(2,3,4)

cov_matrix = kernelMatrix(kernel = SE_Kernel(), x = X, y = X_prime)
cov_matrix



### 2.2
#Consider the following model:
#temp = f(time) + epsilon with epsilon ∼ N (0, σ2n) and f ∼ GP(0, k(time, time′))

#Let σ2n be the residual variance from a simple quadratic regression fit (using the lm function in R). 
#Estimate the above Gaussian process regression model using the squaredexponential function from (1) with σf = 20 and ` = 0.2. 
#Use the predict function in R to compute the posterior mean at every data point in the training dataset. Make
#a scatterplot of the data and superimpose the posterior mean of f as a curve (use
#type="l" in the plot function). Play around with different values on σf and ` (no needto write this in the report though).

sigmaf = 20
ell = 0.2

# Fit quadratic regression with the scaled data
regression_fit = lm(temperature_selection ~ time_selection + I(time_selection)^2)

#Compute the residual variance
sigmaNoise = sd(regression_fit$residuals)

hyperparam = c(sigmaf, ell)

#Compute GP regression

GP_fit = gausspr(x = time_selection, 
                 y = temperature_selection,
                 kernel = SE_Kernel(sigmaF = sigmaf, l = ell),
                 var = sigmaNoise^2)

#Compute the posterior mean at every data point in the training
# dataset
meanPred = predict(GP_fit, time_selection)

plot(time_selection, temperature_selection, main = "Posterior mean",
     ylab = "Temperature", xlab = "Time")
lines(time_selection, meanPred, col = "red", lwd = 3)
legend("bottomright", legend = c("posterior mean", "observations"),col=c("red", "black"), 
       lty = c(1, 1), lwd = c(2,2))



### 2.3
# Make own computations to obtain the posterior variance of f and plot the
# 95 % probability bands for f. To do this we can use the prviously computed
# function PosteriorGP

sigmaf = 20
ell = 0.2
hyperparam = c(sigmaf, ell)

SquaredExpKernel <- function(x1,x2,sigmaF=1,l=3){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}

PosteriorGP = function(X, y, XStar, sigmaNoise, hyperParameter){

  n = length(X)
  #Compute the covarance matrix [K = K(X,X)]
  K = SquaredExpKernel(X,X, hyperParameter[1], hyperParameter[2])

  L_trans = chol(K + sigmaNoise^2*diag(n))
  #Need to take the transpose since L in the algorithm is a lower triangular matrix where as
  # the R function returns an uper triangular matrix
  L = t(L_trans)
  alpha = solve(L_trans, solve(L,y))
  
  ##Compute f_bar*
  K_X_Xstar = SquaredExpKernel(X, XStar, hyperParameter[1], hyperParameter[2])
  
  f_bar_star = t(K_X_Xstar) %*% alpha
  
  v = solve(L, K_X_Xstar)
  
  K_XStar_XStar = SquaredExpKernel(XStar, XStar, hyperParameter[1], hyperParameter[2])
  #Compute V_f*
  v_f_star = K_XStar_XStar - t(v) %*% v
  #To draw from the posterior we only need the variance of f
  v_f_star = diag(v_f_star)
  
  result = list("Predictive mean" = f_bar_star,
                "Predicitive variance" = v_f_star)
}


posterior  = PosteriorGP(X = scale(time_selection),
                                 y = scale(temperature_selection),
                                 XStar = scale(time_selection),
                                 sigmaNoise = sigmaNoise,
                                 hyperParameter = hyperparam)


## Compute the variance for f
posterior_variance = posterior$`Predicitive variance`
posterior_mean = posterior$`Predictive mean`
posterior_mean = posterior_mean*sd(temperature_selection) + mean(temperature_selection)

L = posterior_mean - 1.96*sqrt(posterior_variance)
U = posterior_mean + 1.96*sqrt(posterior_variance)

# Plot the meanPred, and the prediction bands for the posterior variance
plot(time_selection, temperature_selection, main = "Posterior mean",
     ylab = "Temperature", xlab = "Time")
lines(time_selection, posterior_mean, col = "red", lwd = 2)
legend("bottomright", legend = c("posterior mean", "observations"),col=c("red", "black"), 
       lty = c(1, 1), lwd = c(2,2))
lines(time_selection, meanPred, col = "green")
lines(time_selection, meanPred + 1.96*sqrt(posterior_variance), col = "blue")
lines(time_selection, meanPred - 1.96*sqrt(posterior_variance), col = "green")





### Part 2.4

# temp = f(day) + epsilon with epsilon ∼ N (0, σ2n) and f ∼ GP(0, k(day, day′))
#Estimate the model using the squared exponential funciton with simgaf= 20
# and l = 0.2. Superimpose the posterior mean from this model on the posterior
# mean from the model in (2). Compare the results

# start by doing same as in exercise 2 and fit the quadratic regression

regression_fit_4 = lm(temperature_selection ~ day_selection + I(day_selection)^2)
sigmaNoise_day = sd(regression_fit_4$residuals)

ell = 0.2
sigmaf = 20

GP_fit_day = gausspr(x = day_selection,
                     y = temperature_selection,
                     kernel = SE_Kernel(sigmaF = sigmaf, l = ell),
                     var = sigmaNoise_day)

meanPred_day = predict(GP_fit_day, day_selection)

plot(time_selection, temperature_selection, main = "posterior mean",
     ylab = "temperature", xlab = "time")
polygon(c(time_selection, rev(time_selection)),
        c(L, rev(U)), col = "darkgray")
points(time_selection, temperature_selection)
lines(time_selection, meanPred, col = "red")
lines(time_selection, meanPred_day, col = "blue")



### Part 2.5
#Implement a generalization of the periodic kernel given in the lectures
# Note that we have two different l which controls the correlation between
# the same day in different years. Estimate the GP mmodel using the time variable
# with this kernel and the hyperparamerers:

sigmaf = 20 
ell1 = 1
ell2 = 10
d = 365/sd(time_selection)

#Create the periodic kernel
periodic_kernel = function(sigmaf, l1, l2, d){
  Periodic_K = function(x,y){
    result = (sigmaf^2)*exp(-2*((sin(pi*abs(x-y)/d)^2)/(l1^2)))*
                              exp(-0.5*((x-y)^2)/(l2^2))
    return(result)
  }
  class(Periodic_K) = "kernel"
  return(Periodic_K)
}

GP_periodic = gausspr(x = time_selection,
                      y = temperature_selection,
                      kernel = periodic_kernel(sigmaf, ell1, ell2, d),
                      var = sigmaNoise)

meanPred_periodic = predict(GP_periodic, time_selection)

plot(time_selection, temperature_selection, main = "Posterior mean")
polygon(c(time_selection, rev(time_selection)),
        c(L, U), col = "darkgray")
lines(time_selection, meanPred, col = "red")
lines(time_selection, meanPred_day, col = "green")
lines(time_selection, meanPred_periodic, col = "blue")
