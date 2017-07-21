X  <- read.table("boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)
table(X$Gene)

library(e1071)
library(naivebayes)
library(caret)

########## 10-fold CV in caret ############
train_control<- trainControl(method="cv", number=10, savePredictions = TRUE)
grid <- expand.grid(.fL=c(0), .usekernel=c(FALSE))  # fix the parameters of the algorithm

# model_fit<- train(Gene ~., data=X, trControl=train_control, method="naive_bayes", preProcess=c("pca"))
model_fit<- train(Gene ~., data=X, 
                  trControl=train_control, 
                  method="naive_bayes",
                  trControl = trainControl(method="none"),
                  tuneGrid = data.frame(fL=0, usekernel=FALSE))
model_fit
model_fit$pred
model_fit$resample
table(model_fit$pred$obs, model_fit$pred$pred)
mean(model_fit$pred$obs == model_fit$pred$pred)
mean(model_fit$resample$Accuracy)
model_fit



################ e1071 naiveBayes #################
set.seed(321)
k <- nrow(X)
s <- sample(k,k/5) #set 20% of data aside for test data
testData <- X[s,]
table(testData$Gene)
trainData <- X[-s,]
table(trainData$Gene)

fit <- naiveBayes(Gene ~., data=trainData)
pred <- predict(fit, testData[,2:64])
mean(pred == testData$Gene)


############ OR do k-fold yourself #################
# k: k in k-fold cv
k_folds <- function(k, myData) {
  n <- nrow(X)
  # shuffle rows
  X <- myData[order(runif(n)),]
  accuracies <- c()
  m <- floor(n/k) # size of test set
  for (i in 0:(k-1)) {
    idx <- i*m  # start index of test set in this fold
    testData <- X[idx:(idx+m),]
    trainData <- X[-(idx:(idx+m)),]
    
    model <- naiveBayes(Gene ~., data=trainData)
    predictions <- predict(object = model, newdata = testData[,2:64])
    accuracies <- c(accuracies, mean(predictions == testData$Gene))
  }
  accuracies
}
#########################################

set.seed(123)
accuracies <- c()
accuracies <- k_folds(5, X)
accuracies
mean.accuracies <- mean(accuracies)
mean.accuracies














