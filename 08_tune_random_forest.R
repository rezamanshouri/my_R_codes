
### Tnseq
X  <- read.csv("Documents/data/tnseq_data/2_resampling_ready_dups_removed", sep="\t", header= TRUE)
X <- X[,-2] # remove col 2
X[1:10, 1:15]
dim(X)

## boshoff
X  <- read.table("Documents/data/boshoff/new/boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)


library(randomForest)
library(caret)


#### remove all classes with <120 examples
ii = X[,1]=="C" | X[,1]=="E" | X[,1]=="G" | X[,1]=="H" | X[,1]=="I" | X[,1]=="J" | X[,1]=="K" | X[,1]=="L" | X[,1]=="Q"
X <- X[ii,]
X$Gene = factor(X$Gene)
######################


##### partition data to train and test #####
set.seed(123)
trainIndex <- createDataPartition(X$Gene, p=.8, list=F)
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
tunegrid <- expand.grid(.mtry=c(1,3,5,7,9,11,13,15,17,19), .ntree=c(500,1000))
set.seed(8)
custom2 <- train(Gene~., data=trainData, method=customRF, metric="Accuracy", tuneGrid=tunegrid, trControl=control)
summary(custom2)

preds2 <- predict(custom2, testData[,-1])
t2 <- table(testData$Gene, preds2)

print(custom2)
cat("\n------------------------------------------\n\n")
print(confusionMatrix(t2))
cat("\n------------------------------------------\n\n")
cat("Use plot(custom) to plot grid-search of hyper-paramerts!\n\n")













