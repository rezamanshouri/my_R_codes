X  <- read.table("boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)
table(X$Gene)

###############################################################################
################### Cross Validadtion in carret ###############################
#http://www.milanor.net/blog/cross-validation-for-predictive-analytics-using-r/
#machinelearningmastery.com/how-to-estimate-model-accuracy-in-r-using-the-caret-package/
###############################################################################
library(caret)
library(klaR)

# So here is how it works:
# 1- Split 'totality' of data into 80%/20% train1/test1 partitions.
# 2- you are going to use ONLY train1 for CV. for example in 10-fold, you split train1 into
#    90%/10% train2/test2 (DON'T CONFUSE WITH ABOVE), and train model on train2, and test on tets2.
#    Further, test2 is used to adjust parameters fo the algorithm (if there is any).
# 3- At the end, you use test1 to evaluate actual prediction accuracy of your TUNED model in step 2.


# define an 80%/20% train/test split of the dataset
# createDataPartition() samples from within factor levels, i.e. portion 80%/20% is respected for each class! (unlike sample())
split=0.8
trainIndex <- createDataPartition(X$Gene, p=split, list=FALSE)
data_train <- X[ trainIndex,]
data_test <- X[-trainIndex,]



############## some tests ################
lda <- lda(Gene~., data=data_train)

# make predictions
x_test <- data_test[,2:64]
y_test <- data_test[,1]
predictions <- predict(lda, x_test)

# summarize results
confusionMatrix(predictions$class, y_test)
table(predictions$class, y_test)
mean(predictions$class == y_test)



############## k-fold with caret ################
train_control<- trainControl(method="cv", number=10, savePredictions = TRUE)
grid <- expand.grid(.fL=c(0), .usekernel=c(FALSE))  # fix the parameters of the algorithm

model_fit<- train(Gene ~., data=X, trControl=train_control, method="lda", preProcess=c("pca"))
model_fit
model_fit$pred
model_fit$resample
table(model_fit$pred$obs, model_fit$pred$pred)
mean(model_fit$pred$obs == model_fit$pred$pred)
plot(model_fit) #if you have tuning parameters

#prediction accuracy on new (test) data
mprobs <- predict(model_fit, data_test, type = "prob")
confusionMatrix(predictions$class, y_test)
head(mprobs)
head(data_test)


#########################################################
##################  DO it yourself  #####################
# So in my_repeated_cv(), I partition myData into 75%/25% train/test sets, and 
# train model on train, and test it on test. This is done k times. 
# This function is some sort of simplified repeated k-fold
my_repeated_cv <- function(k, myData) {
  n <- nrow(X)
  # shuffle rows
  X <- myData[order(runif(n)),]
  accuracies <- c()
  for (i in 1:k) {
    s <- sample(n,n/4)  # partition myData into 75%/25% train/test
    testData <- X[s,]
    trainData <- X[-s,]
    
    #model <- rpart(Gene ~., data = X[folds[[i]],], method = "class")
    model <- lda (Gene ~., data = trainData )
    predictions <- predict(object = model, newdata = testData, type = "class")
    # for some classifiers only one of the followings will work
    accuracies <- c(accuracies, mean(predictions$class == testData$Gene)) # works for lda
    #accuracies <- c(accuracies, mean(predictions == testData$Gene)) # works for rpart
  }
  accuracies
}

#########################################
# k: k in k-fold cv
k_folds <- function(k, myData) {
  n <- nrow(myData)
  # shuffle rows
  X <- myData[order(runif(n)),]
  accuracies <- c()
  m <- floor(n/k) # size of test set
  for (i in 0:(k-1)) {
    idx <- i*m  # start index of test set in this fold
    testData <- X[idx:(idx+m),]
    trainData <- X[-(idx:(idx+m)),]
    
    #model <- rpart(Gene ~., data = X[folds[[i]],], method = "class")
    model <- lda (Gene ~., data = trainData )
    predictions <- predict(object = model, newdata = testData, type = "class")
    # for some classifiers only one of the followings will work
    accuracies <- c(accuracies, mean(predictions$class == testData$Gene)) # works for lda
    #accuracies <- c(accuracies, mean(predictions == testData$Gene)) # works for rpart
  }
  accuracies
}

#########################################


#### let's do some testing #######
set.seed(123)
accuracies <- c()
accuracies <- k_folds(5, X)
#accuracies <- my_repeated_cv(5, X)
accuracies
mean.accuracies <- mean(accuracies)
mean.accuracies
# CI can be obtained as it's written below


##### repeated k fold #####
# it is a bit dodgy taking a mean of 5 samples,
# splitting our sample into >5 folds would greatly reduce the stability of the estimates
# one solution is repeated k-folds cross-validation
set.seed(567)
accuracies <- replicate(10, k_folds(5, X))
accuracies
mean.accuracies <- mean(accuracies)
mean.accuracies
lci <- mean(accuracies) - sd(accuracies) * 1.96
uci <- mean(accuracies) + sd(accuracies) * 1.96
cat("95% CI: (", lci, ",", uci, ")")

hist(accuracies, col = grey(0.8))
abline(v=mean(accuracies),col="red", lty = 2)
abline(v=lci,col="black")
abline(v=uci,col="black")
legend(x="topright", lty=c(1,1), c("95% CI", "mean"), col=c("black","red") )

# http://t-redactyl.io/blog/2015/10/using-k-fold-cross-validation-to-estimate-out-of-sample-accuracy.html
#########################################################



