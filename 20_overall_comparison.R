
# Here I will perform 5-fold cross validation on all algorithms
# I will shuffle (appropriately) the data firs, and the same test/train sets will be used for each algorithm
# This ensures fair comparison!


boshoff  <- read.table("Documents/my_R/boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)
set.seed(1)
randidx <- sample(nrow(boshoff))
X <- boshoff[randidx,] # shuffle rows
table(X$Gene)


## make sure in all folds the portion of test and train set is reasonable
## supposing that we are going to do 5-fold CV:
k <- 5
n <- nrow(X)
m <- floor(n/k) # size of test set
for(i in 0:4) {
  idx <- i*m  # start index of test set in this fold
  cat("\n---------- Fold #", i+1, "\n")
  testData <- X[idx:(idx+m),]
  print(table(testData$Gene))
  trainData <- X[-(idx:(idx+m)),]
  print(table(trainData$Gene))
}


### 'abs' and 'square' values
X_abs <- cbind(X[,1], abs(X[,2:64]))
colnames(X_abs) <- colnames(X)
X_sqr <- as.data.frame(X[,2:64]^2)
X_sqr <- cbind(X[,1],X_sqr)
colnames(X_sqr) <- colnames(X)


############ Naive Bayes ###############
library(e1071)
#####
k_folds_nb <- function(k, X) {
  n <- nrow(X)
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
#####

acc_nb <- k_folds_nb(5, X_sqr)
acc_nb
mean.acc_nb <- mean(acc_nb)
mean.acc_nb


############ Decision Tree ###############
library(rpart)
#####
k_folds_dt <- function(k, X) {
  n <- nrow(X)
  accuracies <- c()
  m <- floor(n/k) # size of test set
  for (i in 0:(k-1)) {
    idx <- i*m  # start index of test set in this fold
    testData <- X[idx:(idx+m),]
    trainData <- X[-(idx:(idx+m)),]
    
    dt <- rpart(Gene ~ ., data = trainData, method="class", parms = list(split = 'information'), cp = -1 )
    # prune tree with 1 SE rule
    errs <- printcp(dt)
    i1 <- which.min(errs[,4]) # inedx of min xerr
    one_SE <- errs[i1,4] + errs[i1,5]
    i2 <- min(which(errs[,4] < one_SE)) # index of smallest tree within 1SE
    best_cp <- errs[i2,1]
    # now prune
    pt <- prune(dt, cp = best_cp)

    model <- pt 
    predictions <- predict(object = model, newdata = testData[,2:64], type = 'class')
    accuracies <- c(accuracies, mean(predictions == testData$Gene))
  }
  accuracies
}
#####

acc_dt <- k_folds_dt(5, X)
acc_dt
mean.acc_dt <- mean(acc_dt)
mean.acc_dt



############ Random Forest ###############
library(randomForest)
#####
k_folds_rf <- function(k, X) {
  n <- nrow(X)
  accuracies <- c()
  m <- floor(n/k) # size of test set
  for (i in 0:(k-1)) {
    idx <- i*m  # start index of test set in this fold
    testData <- X[idx:(idx+m),]
    trainData <- X[-(idx:(idx+m)),]
    
    model <- randomForest(Gene ~ ., data = trainData, ntree = 1001 , mtry = 15, importance = T, parms = list(split = 'information') )
    predictions <- predict(object = model, newdata = testData[,2:64])
    accuracies <- c(accuracies, mean(predictions == testData$Gene))
  }
  accuracies
}
#####

acc_rf <- k_folds_rf(5, X)
acc_rf
mean.acc_rf <- mean(acc_rf)
mean.acc_rf



########### LDA and QDA ################
library(MASS)
#####
k_folds_lda <- function(k, X) {
  n <- nrow(X)
  accuracies <- c()
  m <- floor(n/k) # size of test set
  for (i in 0:(k-1)) {
    idx <- i*m  # start index of test set in this fold
    testData <- X[idx:(idx+m),]
    trainData <- X[-(idx:(idx+m)),]
    
    model <- lda(Gene ~., trainData)
    #model <- qda(Gene ~., trainData)
    predictions <- predict(object = model, newdata = testData[,2:64], type = 'class')
    accuracies <- c(accuracies, mean(predictions$class == testData$Gene))
  }
  accuracies
}
#####

acc_lda <- k_folds_lda(5, X)
acc_lda
mean.acc_lda <- mean(acc_lda)
mean.acc_lda



########## MY RF #################
library(rpart)
##########
##########
random_forest <- function(X, n, m) {
  trees <- c()
  
  for(i in 1:n) {
    
    # select m random colums (plus firs column)
    k <- ncol(X)
    s <- sample(2:k,m)
    Y <- cbind(X[,1], X[,s])
    colnames(Y)[1] <- "Gene"
    
    use_default <- FALSE
    
    if(use_default) {  ## use default
      T <- rpart(Gene ~ ., data = Y, method="class", parms = list(split = 'information'))
      trees[[paste0("Tree_",i)]] <- T
    }
    else {  ## overfit and then prune 
      T <- rpart(Gene ~ ., data = Y, method="class", parms = list(split = 'information'), cp = -1 )  # "cp = -1" makes it overfit, but we next prune
      
      # find cp whithin 1SE of min xerr (error on CV for validation set)
      errs <- printcp(T)
      i1 <- which.min(errs[,4]) # inedx of min xerr
      one_SE <- errs[i1,4] + errs[i1,5]
      i2 <- min(which(errs[,4] < one_SE)) # index of smallest tree within 1SE
      best_cp <- errs[i2,1]
      
      # now prune
      PT <- prune(T, cp = best_cp)
      trees[[paste0("Tree_",i)]] <- PT
    }
  }
  
  return(trees)
}
##########
##########
make_prediction <- function(RF, testData) {
  
  predictions <- c()
  n_correct_predictions <- 0
  for(e in 1:nrow(testData)) {
    preds <- c()
    for(T in RF) {
      p <- as.character( predict(T, testData[e,2:64], type = 'class') )
      preds <- append(preds, p)
    }
    
    final_prediction <- names(which.max(table(preds)))
    predictions[e] <- final_prediction
    
    if(final_prediction == testData[e,1]) {
      n_correct_predictions <- n_correct_predictions+1
    }
  }
  
  #cat("\n The number of correct predictions is ", n_correct_predictions, " out of ", nrow(testData), " test examples")
  #cat("\n Accuracy is: ", (n_correct_predictions/ nrow(testData))  )
  
  return(predictions)
}
##########
##########
k_folds_my_rf <- function(k, X) {
  n <- nrow(X)
  accuracies <- c()
  m <- floor(n/k) # size of test set
  for (i in 0:(k-1)) {
    idx <- i*m  # start index of test set in this fold
    testData <- X[idx:(idx+m),]
    trainData <- X[-(idx:(idx+m)),]
    
    model <- random_forest(trainData, 500, 15)
    predictions <- make_prediction(RF, testData)
    accuracies <- c(accuracies, mean(predictions == testData$Gene))
  }
  accuracies
}
#####

acc_my_rf <- k_folds_my_rf(5, X)
acc_my_rf
mean.acc_my_rf <- mean(acc_my_rf)
mean.acc_my_rf



###########################












