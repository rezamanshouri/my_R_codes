
library(neuralnet)
library(nnet) 
library(caret)

X  <- read.table("boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)


## make a 70%/30% train/test partition
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
grid <- expand.grid(size=c(5:20), decay=c(0,0.1,0.01,0.001,0.0005))
model <- train(Gene~., trainset, method='nnet', trace=FALSE, tuneGrid=grid, trControl=control)


preds <- predict(model, testset[,2:64])
t <- confusionMatrix(table(testset[,1], preds))

print(model)
cat("\n------------------------------------------\n\n")
print(t)
cat("\n------------------------------------------\n\n")
cat("Use 'plot(model)' to plot grid-search of hyper-paramerts!\n\n")

plot(model)















