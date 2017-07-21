
# CLASSIFICATION USING NN

###############################################
##########       neuralnet        #############
###############################################
# https://www.packtpub.com/mapt/book/big_data_and_business_intelligence/9781783982042/6/ch06lvl1sec72/training-a-neural-network-with-neuralnet

library(neuralnet)
library(nnet) # for its
library(caret)
data(iris)

## make a 70%/30% train/test partition
trIdx <- createDataPartition(iris$Species, p=.7, list=F)
trainset <- iris[trIdx, ]
testset <- iris[-trIdx, ]
table(trainset$Species)
table(testset$Species)


## Add dummy variables for discrete target (i.e. CLASSIFICATION)
train <- cbind(irisTrain[, 1:4], nnet::class.ind(as.factor(irisTrain[,5])))
test <- cbind(irisTest[, 1:4], nnet::class.ind(as.factor(irisTest[,5])))

# Scale data
scl <- function(x){ (x - min(x))/(max(x) - min(x)) }
train[, 1:4] <- data.frame(lapply(train[, 1:4], scl))

## fit model
f <- as.formula("setosa + versicolor + virginica ~ Sepal.Length + Sepal.Width + Petal.Width + Petal.Length")
network <- neuralnet(f, train, hidden = 3, stepmax=1e6)
network$result.matrix

# GENERALIZED WEIGHTS
# If all the generalized weights are close to zero on the plot, it means the covariate has little effect.
# However, if the overall variance is greater than one, it means the covariate has a nonlinear effect.
par(mfrow=c(2,2))
gwplot(network, selected.covariate="Petal.Width")
gwplot(network, selected.covariate="Petal.Length")
gwplot(network, selected.covariate="Sepal.Width")
gwplot(network, selected.covariate="Petal.Length")

## prediction accuracy
nnout <- compute(network, test[,1:4])
results <- nnout$net.result
head(nnout$neurons)
head(nnout$net.result)
head(test[,5:7])

# Accuracy
original_values <- max.col(test[, 5:7])
preds <- max.col(results)
mean(preds == original_values)
t <- table(original_values, preds)
confusionMatrix(t)


###############################################
##########       nnet             #############
###############################################
## nnet supports classification which makes our life easy!!!
library(nnet)
data(iris)

trIdx <- createDataPartition(iris$Species, p=.7, list=F)
trainset <- iris[trIdx, ]
testset <- iris[-trIdx, ]
table(trainset$Species)
table(testset$Species)

# fit model
fit <- nnet(Species~., data=trainset, size=4, decay=0.0001, maxit=500)
summary(fit)

# make predictions
predictions <- predict(fit, testset[,1:4], type="class")
# summarize accuracy
table(predictions, testset$Species)

## NOTE:
## in neural net, we used 'train/test', but in nnet, we used 'trainset/testset'

###############################################
######## Tuning Parameters in caret ###########
###############################################

# NOTE: Neural net from caret only deals with regression and takes 3 params i.e. layers 1-3.

#Grid of tuning parameters to try:
grid <- expand.grid(size=c(2,3,4,5,6,7,8,9,10,15), decay=c(0.1,0.01,0.001))
model <- train(Species~., trainset, method='nnet', trace = FALSE, tuneGrid = grid) 
plot(model)
model

best_preds <- predict(model, testset[,1:4])
table(testset[,5], best_preds)
mean(testset[,5] == best_preds)



f <- as.formula("setosa + versicolor + virginica ~ Sepal.Length + Sepal.Width + Petal.Width + Petal.Length")
network <- neuralnet(f, train, hidden = 3, stepmax=1e6)

model <- train(f, train, method='neuralnet', trace = FALSE) 
plot(model)
model












