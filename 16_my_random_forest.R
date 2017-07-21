library(party)
library(rpart)

X <- read.table("boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)

####################################################################


############## random_fores() ###############
# X : data frame (testSet)
# n : number of trees to be built
# m : number of feature to be sampled randomly
random_forest <- function(X, n, m) {
  trees <- c()
  
  for(i in 1:n) {
    
    # select m random colums (plus firs column)
    k <- ncol(X)
    s <- sample(2:k,m)
    Y <- cbind(X[,1], X[,s])
    colnames(Y)[1] <- "Gene"
    
    
    use_default <- TRUE
    
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


###############################################


############## make_prediction() ###############
# RF : Random Forest
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
  
  cat("\n The number of correct predictions is ", n_correct_predictions, " out of ", nrow(testData), " test examples")
  cat("\n Accuracy is: ", (n_correct_predictions/ nrow(testData))  )
  
  return(predictions)
}

#######




############## testing my random forest ###################

k <- nrow(X)
s <- sample(k,k/5) #set 20% of data aside for test data
testData <- X[s,]
table(testData$Gene)
trainData <- X[-s,]
table(trainData$Gene)


rf <- random_forest(trainData, 500, 8)
rf1 <- random_forest(trainData, 100, 15)
rf2 <- random_forest(trainData, 100, 15) # no prune, default, stopping criteria: minsplit, minBucket
rf3 <- random_forest(trainData, 100, 15) # gini index, prune

rf2 <- random_forest(trainData, 50, 20)
rf3 <- random_forest(trainData, 50, 8)
RF <- rf3

preds <- make_prediction(RF, testData)
table(preds, testData$Gene)
mean(preds == testData$Gene)


# number of trees that are just a root node after pruning
num_1node_trees <- 0
for(i in 1:length(RF)) {
  if( nrow(RF[[i]]$cptable) == 1 ) { # probably not the best way tho :)
    num_1node_trees <- num_1node_trees+1
  }
}
num_1node_trees


library(rattle)
library(rpart.plot)
fancyRpartPlot(RF[[2]])
for(i in 1:50) {
  if( nrow(RF[[i]]$cptable) == 1 ) {
    png(paste("z_tree_in_RF_", i, ".png", sep=""))
    fancyRpartPlot(RF[[i]])
    dev.off()
  }
}


######## perform repeated k fold CV  ###########
k_folds <- function(k, X) {
  folds <- createFolds(X$Gene, k = k, list = TRUE, returnTrain = TRUE)
  accuracies <- c()
  for (i in 1:k) {
    rf <- random_forest(X[folds[[i]],], 500, 8)
    preds <- make_prediction(rf, X[-folds[[i]],])
    accuracies <- c(accuracies, mean(preds == X[-folds[[i]], ]$Gene)) 
  }
  accuracies
}

set.seed(567)
v <- c()
v <- replicate(20, k_folds(5, X))
accuracies <- c()
for (i in 1 : 20) { 
  accuracies <- c(accuracies, v[,i])
}

v
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




