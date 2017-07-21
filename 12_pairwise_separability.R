library(caret)
library(MASS)


############ consider all pairwise classes  ###########
## build pairwise models and return accuracy using k-fold CV

### pairwise_model_accuracy() ###
pairwise_model_accuracy <- function(X, c1, c2) {
  if(c1 == c2) return(NA)

  ii = X[,1]==c1  | X[,1]==c2
  X <- X[ii,]
  X$Gene <- factor(X$Gene)  

  k <- nrow(X)
  s <- sample(k,k/10) #pick 1/3 of data for test data
  testData <- X[s,]
  table(testData$Gene)
  trainData <- X[-s,]
  table(trainData$Gene)
  
  lda <- lda(Gene ~., trainData)
  
  if(FALSE) {
    png(paste("pairwise_model_accuracy_", c1, "_", c2, ".png", sep=""))
    plot(lda)
    dev.off()
  }
  
  pred <- predict(lda, testData[,2:64])$class
  a <- mean(pred == testData[,1])
  
  round(a,digits=2)


}
##################



X <- read.table("Documents/my_R/boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)
X <- X[ X[,1]!="D" , ]
X <- X[ X[,1]!="F" , ]
X <- X[ X[,1]!="U" , ]
X <- X[ X[,1]!="V" , ]
table(X$Gene)
X$Gene <- factor(X$Gene)
levels(X$Gene)

cl <- c("C", "E", "G", "H", "I", "J", "K", "L", "M", "O", "P", "Q", "T")
P <- sapply(cl, function(x) sapply(cl, function(y) pairwise_model_accuracy(X, x,y) ) )
# Note, for example, P[1,3] and P[3,1] are exactly same models


## print lowest and highest values in P
P
which(P == max(P, na.rm = T) , arr.ind = TRUE)
max(P, na.rm = T)
which(P == min(P, na.rm = T) , arr.ind = TRUE)
min(P, na.rm = T)

mean(P, na.rm = T)











