library(lattice)
library(ggplot2)
library(caret)
library(MASS)


X <- read.table("boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)

#### remove classes with less than 100 examples #####
table(X$Gene)
levels(X$Gene)
X <- X[ X[,1]!="D" , ]
X <- X[ X[,1]!="F" , ]
X <- X[ X[,1]!="U" , ]
X <- X[ X[,1]!="V" , ]
table(X$Gene)
X$Gene <- factor(X$Gene)
levels(X$Gene)



########### just for testing ##############
testData <- X[1:5,]
cl <- c("C", "E", "G", "H", "I")
pairwise_models <- sapply(cl, function(x) sapply(cl, function(y) binary_classifier(X, x,y) ) )
P <- pairwise_models
###########################################


########################################################################
########## multiclass classification using a binary classifier #########


######## binary_classifier() #########
binary_classifier <- function(X, c1, c2) {
  if(c1 == c2) return(NA)
  
  ii = X[,1]==c1  | X[,1]==c2
  Y <- X[ii,]
  Y$Gene <- factor(Y$Gene)  
  
  model <- qda(Gene ~., Y)
  return(model)
}
#######################################

######## make_prediction() ###########
# P: pairwise binary classifiers matrix
predict_multinomial_from_binary_classifiers <- function(P, testData) {
  
  predictions <- c()
  cl <- c("C", "E", "G", "H", "I", "J", "K", "L", "M", "O", "P", "Q", "T")
  
  n_correct_predictions <- 0
  for(e in 1:nrow(testData)) {
    
    max_wins <- 0
    final_prediction <- ""
    
    for(i in 1:nrow(P)) { #iterate rows
      
      n_row_wins <- 0
      for(j in 1:ncol(P)) { #iterate cols
        if(i==j) next
        
        pred <- predict(P[i,j][[1]], testData[e,2:64])$class
        pred <- factor(pred)
        
        #print(paste0("e_i_j = ", e, "_", i, "_", j))
        #print(paste0("pred == cl[[i]] ___", pred ,"==", cl[[i]]))
        
        if(pred == cl[[i]]) {
          #print("--------corect!!")
          n_row_wins = n_row_wins+1
        }
        
      }
      
      #print(paste0("-----------------------n_row_wins= ", n_row_wins))
      if(n_row_wins > max_wins) {
        max_wins = n_row_wins
        final_prediction <- cl[[i]]
        #print(paste0("------final_prediction: ", final_prediction))
      }
    }
    
    predictions[e] <- final_prediction
    if(final_prediction == testData[e,1]) {
      n_correct_predictions <- n_correct_predictions+1
    }
  }
  
  cat("The number of correct predictions is ", n_correct_predictions, " out of ", nrow(testData), " test examples")
  cat("Accuracy is: ", (n_correct_predictions/ nrow(testData))  )

  return(predictions)
}
#######################################



##### partition data to train and test #####
set.seed(123)
trainIndex <- createDataPartition(X$Gene, p=.7, list=F)
trainData <- X[trainIndex, ]
testData <- X[-trainIndex, ]
table(testData$Gene)
table(trainData$Gene)


cl <- c("C", "E", "G", "H", "I", "J", "K", "L", "M", "O", "P", "Q", "T")
pairwise_models <- sapply(cl, function(x) sapply(cl, function(y) binary_classifier(trainData, x,y) ) )
# Note, for example, P[1,3] and P[3,1] are exactly same models

preds <- predict_multinomial_from_binary_classifiers(pairwise_models, testData)
table(preds, testData$Gene)
mean(preds == testData$Gene)





############# K-Fold #################

k_folds <- function(k, X) {
  folds <- createFolds(X$Gene, k = k, list = TRUE, returnTrain = TRUE)
  accuracies <- c()
  for (i in 1:k) {
    pairwise_models <- sapply(cl, function(x) sapply(cl, function(y) binary_classifier(X[folds[[i]],], x,y) ) )
    preds <- predict_multinomial_from_binary_classifiers(pairwise_models, X[-folds[[i]],])
    accuracies <- c(accuracies, mean(preds == X[-folds[[i]], ]$Gene)) 
  }
  accuracies
}
#############

set.seed(123)
accuracies <- c()
accuracies <- k_folds(10, X)
accuracies
mean.accuracies <- mean(accuracies)
mean.accuracies
#######################################
















