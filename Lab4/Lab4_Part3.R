#Assignment 3 - GP Classification with kernlab

library(kernlab)
library(AtmRay)

data <-
  read.csv(
    "https://github.com/STIMALiU/AdvMLCourse/raw/master/
GaussianProcess/Code/banknoteFraud.csv",
    header = FALSE,
    sep = ","
  )
names(data) <-
  c("varWave", "skewWave", "kurtWave", "entropyWave", "fraud")

data[, 5] <- as.factor(data[, 5])

#Select 1000 datapoints
set.seed(111)
SelectTraining <-
  sample(1:dim(data)[1], size = 1000, replace = FALSE)
train = data[SelectTraining, ]
test = data[-SelectTraining, ]

# Part 3.1
GP.fit = gausspr(fraud ~ varWave + skewWave, data = train)

x1 <- seq(min(train$varWave), max(train$varWave), length = 100)
x2 <- seq(min(train$skewWave), max(train$skewWave), length = 100)
gridPoints <- meshgrid(x1, x2)
gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))
gridPoints <- data.frame(gridPoints)
names(gridPoints) <- names(data)[1:2]
probPreds <- predict(GP.fit, gridPoints, type = "probabilities")

contour(x1, x2, matrix(probPreds[, 1], 100, byrow = TRUE), 20)
points(train[train$fraud == 1, "varWave"], train[train$fraud == 1, "skewWave"], col = "blue")
points(train[train$fraud == 0, "varWave"], train[train$fraud == 0, "skewWave"], col = "red")

### Part 3.2
# Compute the accuracy on predicting the test data
test_prediction = predict(GP.fit, test, type = "response")
accuracy = mean(test_prediction == test$fraud)
accuracy

### Part 3.3
GP.fit_all = gausspr(fraud ~ ., data = train)
test_prediction_all = predict(GP.fit_all, test, type = "response")
accuracy_all = mean(test_prediction_all == test$fraud)
accuracy_all
