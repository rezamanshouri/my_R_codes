
### read data
boshoff  <- read.table("boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)
table(boshoff$Gene)
str(boshoff)
##shuffle rows
if(FALSE){
  set.seed(1987)
  g <- runif(nrow(boshoff))
  boshoff <- boshoff[order(g),]
}
###############################################################################################
library(randomForest)
#help(randomForest)
?randomForest

############
#the final classification is the majority "vote" of trees in the forest
#Each tree is grown as follows:
#If the number of cases in the training set is N, sample N cases at random - but with replacement, from the original data. This sample will be the training set for growing the tree.
#If there are M input variables, a number m<<M is specified such that at each node, m variables are selected at random out of the M and the best split on these m is used to split the node. The value of m is held constant during the forest growing.
#Each tree is grown to the largest extent possible. There is no pruning.

#see below link for more info
#https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm#overview
###########



### split data to test and train
k <- nrow(X)
s <- sample(k,k/5) #set 20% of data aside for test data
testData <- boshoff[s,]
table(testData$Gene)
trainData <- boshoff[-s,]
table(trainData$Gene)

#using 'gini gain' (default) OR 'information gain' as splitting criteria in decision trees
#actually using __parms = list(split = 'information')__ in rpart make it use 'information index', I assume(!) it does same thing here too
rf <- randomForest(Gene ~ ., data = trainData, ntree = 5001 , mtry = 8, importance = T, parms = list(split = 'information') )
#rf <- randomForest(Gene ~ ., data = trainData, ntrees = 1000, mtry = 8 ,importance = T)
plot(rf)
print(rf)
#Dotchart of variable importance as measured by a Random Forest
#varImpPlot(rf)


### prediction accuracy
pred <- predict(rf, testData, type = 'class') #NOTE: for rpar you need to specify "type = 'class'"
table(testData[,1], pred)
#(i,i) is the correct predictions of class i, the (i,j) is the number of missclassification where class "i" is (miss)classified as "j"
z <- pred == testData[,1]
table(z)
mean(z)



###### better plot #####
#see 06_pca.R for more details on "mycol"
mycol = c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")
layout(matrix(c(1,2),nrow=1), width=c(4,1)) 
par(mar=c(5,4,4,0)) #No margin on the right side
plot(rf, log="y", col = mycol)  ## OOB is overall error rate
par(mar=c(5,0,4,2)) #No margin on the left side
plot(c(0,1),type="n", axes=F, xlab="", ylab="")
legend("top", colnames(rf$err.rate), col=mycol, cex=0.8, fill=mycol)



###############################################
######## Tuning Parameters in caret ###########
###############################################
# http://machinelearningmastery.com/tune-machine-learning-algorithms-in-r/

##### partition data to train and test #####
trainIndex <- createDataPartition(X$Gene, p=.7, list=F)
trainData <- X[trainIndex, ]
testData <- X[-trainIndex, ]
table(testData$Gene)
table(trainData$Gene)

## if you just want to tune 'mtry'
set.seed(8)
control <- trainControl(method="repeatedcv", number=10, repeats=3)
tunegrid <- expand.grid(.mtry=c(5:10))
rf_gridsearch <- train(Gene~., data=trainData, method="rf", metric="Accuracy", tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch)

preds <- predict(rf_gridsearch, testData[,-1])
table(testData$Gene, preds)
mean(testData$Gene == preds)

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
tunegrid <- expand.grid(.mtry=c(1:4), .ntree=c(10, 100, 500, 1000))
set.seed(8)
custom <- train(Gene~., data=trainData, method=customRF, metric="Accuracy", tuneGrid=tunegrid, trControl=control)
summary(custom)
plot(custom)

preds <- predict(custom, testData[,-1])
table(testData$Gene, preds)

