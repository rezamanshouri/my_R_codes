
library(neuralnet)
library(nnet) 
library(caret)

X  <- read.table("boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)


#### remove all classes with <120 examples
ii = X[,1]=="C" | X[,1]=="E" | X[,1]=="G" | X[,1]=="H" | X[,1]=="I" | X[,1]=="J" | X[,1]=="K" | X[,1]=="L" | X[,1]=="Q"
X <- X[ii,]
X$Gene = factor(X$Gene)
######################


## make a 70%/30% train/test partition
set.seed(123)
trainIndex <- createDataPartition(X$Gene, p=.7, list=F)
trainset <- X[trainIndex, ]
testset <- X[-trainIndex, ]
table(testset$Gene)
table(trainset$Gene)


# Scale data
scl <- function(x){ (x - min(x))/(max(x) - min(x)) }
trainset[, 2:64] <- data.frame(lapply(trainset[, 2:64], scl))




###############################################
######## Tuning Parameters in caret ###########
###############################################

# NOTE: Neural net from caret only deals with regression and takes 3 params i.e. layers 1-3.

control <- trainControl(method="repeatedcv", number=10, repeats=3)
#Grid of tuning parameters to try:
grid <- expand.grid(size=c(5,7,9,11,13,15), decay=c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
model3 <- train(Gene~., trainset, method='nnet', trace=FALSE, tuneGrid=grid, trControl=control, maxit=1000)


preds3 <- predict(model3, testset[,2:64])
t3 <- confusionMatrix(table(testset[,1], preds3))

print(model3)
cat("\n------------------------------------------\n\n")
print(t3)
cat("\n------------------------------------------\n\n")
cat("Use 'plot(model)' to plot grid-search of hyper-paramerts!\n\n")

plot(model3)















