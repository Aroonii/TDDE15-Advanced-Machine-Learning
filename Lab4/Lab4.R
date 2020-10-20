library(kernlab)



### Assignment 1### 

# Covariance function
# The function takes in input values x1, x2 and computes the exponential kernel which gives
# the resulting covariance between the two input values
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
  #n is how many functions we will need, as many as the nr of inputs
  n = length(X)
  
  #Compute the covarance matrix [K = K(X,X)]
  K = SquaredExpKernel(X,X, hyperParameter[1], hyperParameter[2])
  sigmaNoise = sigmaNoise^2
  
  #Compute L by using the cholesky decomposition
  L_trans = chol(K + sigmaNoise*diag(n))
  #Need to take the transpose since L in the algorithm is a lower triangular matrix where as
  # the R function returns an uper triangular matrix
  L = t(L_trans)
  
  ### Predictive mean f_bar*
  ##Comute the predictive mean by solving the equations
  # L\y means the vector x that solves the equation Lx = y. On paper it can be solved by
  # multiplying by the inverse but is better solved by using the function solve
  # [alpha = t(L)\(L\y)]
  alpha = solve(L_trans, solve(L,y))
  
  ##Compute f_bar*
  #f_bar* = alpha * t(K*)
  # [K* = K(X, X*) => t(K*) = K(X*, X)]
  K_X_Xstar = SquaredExpKernel(X, XStar, hyperParameter[1], hyperParameter[2])
  f_bar_star = t(K_X_Xstar) %*% alpha
  
  ### Predictive variance f_star
  #Compute v, [v = L\K*]
  v = solve(L, K_X_Xstar)
  
  ## Compute the variance of f*
  #V[f_star] = K(X*, X*) - t(v)*v
  # First need to compute K[X*, X*]
  K_XStar_XStar = SquaredExpKernel(XStar, XStar, hyperParameter[1], hyperParameter[2])
  #Compute V_f*
  v_f_star = K_XStar_XStar - t(v) %*% v
  #To draw from the posterior we only need the variance of f
  v_f_star = diag(v_f_star)
  
  result = list("Predictive mean" = f_bar_star,
                "Predicitive variance" = v_f_star)
}


### Part 1.2
# Let the hyperparameters be the following: sigmaf = 1, el = 0,3, using a singel observation (x,y) = (0.4, 0.719), sigmanoise = 0.1
# Plot the posterior


#Function for plotting the result
plotGP = function(mean,variance,grid,x,y){
  plot(grid,mean,ylim = c(min(mean-1.96*sqrt(variance))
                          ,max(mean+1.96*sqrt(variance))),
       type = "l")
  lines(grid,
        mean+1.96*sqrt(variance), 
        col = "red")
  lines(grid,
        mean-1.96*sqrt(variance), 
        col = "red")
  points(x,y)
}

sigmaF = 1
ell = 0.3
obs = data.frame(0.4, 0.719)
sigmanoise = 0.1
xGrid = seq(-1,1,0.01)

GP = PosteriorGP(obs[,1], obs[,2], xGrid, sigmanoise,c(sigmaF, ell))
plotGP(GP$`Predictive mean`, GP$`Predicitive variance`, xGrid, obs[1,1], obs[1,2])



### Part 1.3
#Update the posterior with another observation (x,y) = (-0.6, -0.044) and plot the posterior mean of f and the probability bands



newobs = c(-0.6, -0.044)
obs_1.3 = rbind(obs, newobs)
GP = PosteriorGP(obs_1.3[,1], obs_1.3[,2], xGrid, sigmanoise,c(sigmaF, ell))
plotGP(GP$`Predictive mean`, GP$`Predicitive variance`, xGrid, obs_1.3[,1], obs_1.3[,2])


### Part 1.4
#Compute now the posterior distirbution of f using all available data points

obs_1.4 = data.frame(x=c(-1,-0.6, -0.2, 0.4, 0.8), y=c(0.768, -0.044, -0.904, 0.719,-0.664))
GP = PosteriorGP(obs_1.4[,1], obs_1.4[,2], xGrid, sigmanoise,c(sigmaF, ell))
plotGP(GP$`Predictive mean`, GP$`Predicitive variance`, xGrid, obs_1.4[,1], obs_1.4[,2])


#Part 1.5
#Repeat the exercise, this time wiht hyperparameters sigmaf = 1 and ell = 1
# Compare the results
sigmaF = 1
ell = 1
GP = PosteriorGP(obs_1.4[,1], obs_1.4[,2], xGrid, sigmanoise,c(sigmaF, ell))
plotGP(GP$`Predictive mean`, GP$`Predicitive variance`, xGrid, obs_1.4[,1], obs_1.4[,2])

