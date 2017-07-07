library(caret)
library(MASS)


### binary_classifier() ###
binary_classifier <- function(X, c1, c2) {
  if(c1 == c2) return(NA)
  
  ii = X[,1]==c1  | X[,1]==c2
  Y <- X[ii,]
  Y$Gene <- factor(Y$Gene)  
  
  model <- lda(Gene ~., Y)
  return(model)
}
#######################


X <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE)
k <- nrow(X)
s <- sample(k,k/3) #pick 1/3 of data for test data
testData <- X[s,]
table(testData$Gene)
trainData <- X[-s,]
table(trainData$Gene)

X <- trainData


cl <- c("C", "E", "G", "H", "I", "J", "K", "L", "M", "O", "P", "Q", "T")
pairwise_models <- sapply(cl, function(x) sapply(cl, function(y) binary_classifier(X, x,y) ) )
P <- pairwise_models





######testing
testData <- X[1:5,]
cl <- c("C", "E", "G", "H", "I")
pairwise_models <- sapply(cl, function(x) sapply(cl, function(y) binary_classifier(X, x,y) ) )
P <- pairwise_models
######
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
  
  if(final_prediction == testData[e,1]) {
    n_correct_predictions <- n_correct_predictions+1
  }
}

accuracy <- (n_correct_predictions/ nrow(testData))





a <- mean(pred == testData[,1])
round(a,digits=2)














## print lowest and highest values in P
P
which(P == max(P, na.rm = T) , arr.ind = TRUE)
which(P == min(P, na.rm = T) , arr.ind = TRUE)
mean(P, na.rm = T)





















