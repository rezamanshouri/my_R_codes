library(caret)
library(MASS)


############ consider all pairwise classes  ###########
## I am wondering if there are some pair of classes that can be weel classified!

### lda_pairwise() ###
lda_pairwise <- function(X, c1, c2) {
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
    png(paste("lda_pairwise_", c1, "_", c2, ".png", sep=""))
    plot(lda)
    dev.off()
  }
  
  pred <- predict(lda, testData[,2:64])$class
  a <- mean(pred == testData[,1])
  
  round(a,digits=2)


}
##################



X <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE)
cl <- c("C", "E", "G", "H", "I", "J", "K", "L", "M", "O", "P", "Q", "T")
P <- sapply(cl, function(x) sapply(cl, function(y) lda_pairwise(X, x,y) ) )

## print lowest and highest values in P
P
which(P == max(P, na.rm = T) , arr.ind = TRUE)
which(P == min(P, na.rm = T) , arr.ind = TRUE)
mean(P, na.rm = T)











