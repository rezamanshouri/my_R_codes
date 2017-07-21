
### read data
boshoff  <- read.table("boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)

library(randomForest)
library(caret)

##### partition data to train and test #####
trainIndex <- createDataPartition(X$Gene, p=.7, list=F)
trainData <- X[trainIndex, ]
testData <- X[-trainIndex, ]
table(testData$Gene)
table(trainData$Gene)


###############################
## Howevr, if you want to tune both 'mtry' and 'ntree', you have to create a customRF first: 
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes


# train model
control <- trainControl(method="repeatedcv", number=10, repeats=3)
tunegrid <- expand.grid(.mtry=c(3:22), .ntree=c(100,500,1000,1500))
set.seed(8)
custom <- train(Gene~., data=trainData, method=customRF, metric="Accuracy", tuneGrid=tunegrid, trControl=control)
summary(custom)

preds <- predict(custom, testData[,-1])
t <- table(testData$Gene, preds)

print(custom)
cat("\n------------------------------------------\n\n")
print(confusionMatrix(t))
cat("\n------------------------------------------\n\n")
cat("Use plot(custom) to plot grid-search of hyper-paramerts!\n\n")













