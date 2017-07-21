X  <- read.table("Documents/my_R/boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)
table(X$Gene)
factor(X$Gene)

library(neuralnet)
library(caret)
################## usin neuralnet package #########################

###### convert class labels to neumeric values ######
#replace each letter in column 1 (Gene) with a number between 1 and 17
#X <- within(X, Gene <- as.numeric(factor(Gene, labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) )))
#factor(X$Gene)

##### min-max normalization ######
##this is a very important step!
#scale the data in the interval [0,1]
for(i in 2:64) {
  min <- min(X[,i])
  max <- max(X[,i])
  X[,i] = (X[,i]-min)/(max-min)
}


##### partition data to train and test #####
trainIndex <- createDataPartition(X$Gene, p=.7, list=F)
trainData <- X[trainIndex, ]
testData <- X[-trainIndex, ]
table(testData$Gene)
table(trainData$Gene)


# Binarize the categorical output
trainData <- cbind(trainData, trainData$Gene == 'C')
trainData <- cbind(trainData, trainData$Gene == 'D')
trainData <- cbind(trainData, trainData$Gene == 'E')
trainData <- cbind(trainData, trainData$Gene == 'F')
trainData <- cbind(trainData, trainData$Gene == 'G')
trainData <- cbind(trainData, trainData$Gene == 'H')
trainData <- cbind(trainData, trainData$Gene == 'I')
trainData <- cbind(trainData, trainData$Gene == 'J')
trainData <- cbind(trainData, trainData$Gene == 'K')
trainData <- cbind(trainData, trainData$Gene == 'L')
trainData <- cbind(trainData, trainData$Gene == 'M')
trainData <- cbind(trainData, trainData$Gene == 'O')
trainData <- cbind(trainData, trainData$Gene == 'P')
trainData <- cbind(trainData, trainData$Gene == 'Q')
trainData <- cbind(trainData, trainData$Gene == 'T')
trainData <- cbind(trainData, trainData$Gene == 'U')
trainData <- cbind(trainData, trainData$Gene == 'V')

names(trainData)[65:81] <- c('C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'O', 'P', 'Q', 'T', 'U', 'V')


#####################################################
################# using neuralnet package ###########
#http://horicky.blogspot.de/2012/06/predictive-analytics-neuralnet-bayesian.html
#####################################################


n <- names(X)
#f <- as.formula(paste("Gene ~", paste(n[!n %in% "Gene"], collapse = " + ")))
f <- " C + D + E + F + G + H + I + J + K + L + M + O + P + Q + T + U + V ~
    ARP4 + Amikacin + Amp + Antimycin + Ascididemin_Nat_prod + 
    BZA + CCCP + CPZ + Cap + Cephalexin + Cerulenin + Clofazimine + 
    Clotrimazole + DCCD + DIPED + DNP + DTNB + DTT + Deferoxamine + 
    Dipyridyl + EMB + Econazole + Eta + GSNO + GSNO_CFZ + GSNO_CPZ + 
    GSNO_KCN + GSNO_menadione + H2O2 + INH + KCN + Levo + Min_med__succinate_ + 
    Mtm + NRP_1 + NaN3 + Nam + Nigericin + Novobiocin + Oflox + 
    PA_1 + PA_21 + PA824 + PZA + Palmitate + Procept_6776 + Procept_6778 + 
    Rif + Rifp + Rox + SM + Succinate + TLM + TRC + TRZ + Tet + 
    UV + Valinomycin + Verapamil + ZnSO4 + menadione + mercaptoethanol + 
    methoxatin"

#for example, 'hidden=c(5,3)' means 2 hidden layers with 5 and 3 neurons in each, respectively.
# "linear.output=FALSE" says to use activation function
---DONT FORGET PUTTING linear.output = FALSE, OTHERWISE IT DOES REGRESSION
nn1 <- neuralnet(f, data = trainData, hidden=20, stepmax=1e6, linear.output=FALSE, learningrate=0.0001, algorithm = "backprop")
nn2 <- neuralnet(f, data = trainData, hidden=20, stepmax=1e6, linear.output=FALSE, algorithm = "rprop+")
nn3 <- neuralnet(f, data = trainData, hidden=c(15,15), stepmax=1e6, linear.output=FALSE, algorithm = "rprop+")
nn4 <- neuralnet(f, data = trainData, hidden=20, stepmax=1e6, linear.output=FALSE, algorithm = "sag")
nn5 <- neuralnet(f, data = trainData, hidden=20, stepmax=1e6, linear.output=FALSE, algorithm = "slr")

#Note that I set the argument linear.output to FALSE in order to tell the model that
#I want to apply the activation function act.fct and that I am not doing a regression task.




##### prediction accuracy #####
nn <- nn2

print(nn)
plot(nn)

nnpred <- compute(nn, testData[,2:64])
#cbind(testData[,1], nnpred)
nnpred.weights <- nnpred$net.result

# Put multiple binary output to categorical output
idx <- apply(nnpred.weights, c(1), which.max)
prediction <- c('C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'O', 'P', 'Q', 'T', 'U', 'V')[idx]
table(prediction, testData$Gene)
mean(prediction == testData$Gene)




################################################################
############ Doing same thing simpler #########################

## There are 2 much easier ways to tranform categorical variable to numerical:
## 1- function model.matrix() creates dummy varibales for factor variables (Gene)
## 2- nnet::class.ind() ! very easy! Introduce a new varibale, func, and class.ind() will do the job!

train <- cbind(trainData[, 2:64], nnet::class.ind(as.factor(trainData[,1])))
test <- cbind(testData[, 2:64], nnet::class.ind(as.factor(testData[,1])))
f <- " C + D + E + F + G + H + I + J + K + L + M + O + P + Q + T + U + V ~
    ARP4 + Amikacin + Amp + Antimycin + Ascididemin_Nat_prod + 
    BZA + CCCP + CPZ + Cap + Cephalexin + Cerulenin + Clofazimine + 
    Clotrimazole + DCCD + DIPED + DNP + DTNB + DTT + Deferoxamine + 
    Dipyridyl + EMB + Econazole + Eta + GSNO + GSNO_CFZ + GSNO_CPZ + 
    GSNO_KCN + GSNO_menadione + H2O2 + INH + KCN + Levo + Min_med__succinate_ + 
    Mtm + NRP_1 + NaN3 + Nam + Nigericin + Novobiocin + Oflox + 
    PA_1 + PA_21 + PA824 + PZA + Palmitate + Procept_6776 + Procept_6778 + 
    Rif + Rifp + Rox + SM + Succinate + TLM + TRC + TRZ + Tet + 
    UV + Valinomycin + Verapamil + ZnSO4 + menadione + mercaptoethanol + 
    methoxatin"

## fit model
nn2 <- nnet(f, train, size = 5, maxit = 1000)
nn2 <- neuralnet(f, data = trainData, hidden=c(15,15), stepmax=1e6, linear.output=FALSE, algorithm = "rprop+")

## prediction accuracy
nnout <- compute(nn2, test[,1:4])
results <- nnout$net.result

# Accuracy
original_values <- max.col(test[, 5:7])
preds <- max.col(results)
mean(preds == original_values)



#################################################################
#################  Parameter Tuning  ############################
################# using caret with 'nnet' #######################
#################################################################
# https://stackoverflow.com/questions/17457028/working-with-neuralnet-in-r-for-the-first-time-get-requires-numeric-complex-ma
# https://stackoverflow.com/questions/7743768/using-nnet-for-prediction-am-i-doing-it-right
# https://gist.github.com/mick001/5973654a443a79b5d1b911a22c00e487

library(nnet)
library(caret)


trainIndex <- createDataPartition(X$Gene, p=.7, list=F)
trainData <- X[trainIndex, ]
testData <- X[-trainIndex, ]

train <- cbind(trainData[, 2:64], nnet::class.ind(as.factor(trainData[,1])))
test <- cbind(testData[, 2:64], nnet::class.ind(as.factor(testData[,1])))

# Scale data
scl <- function(x){ (x - min(x))/(max(x) - min(x)) }
train[, 2:64] <- data.frame(lapply(train[, 2:64], scl))

f <- as.formula(" C + D + E + F + G + H + I + J + K + L + M + O + P + Q + T + U + V ~
    ARP4 + Amikacin + Amp + Antimycin + Ascididemin_Nat_prod + 
    BZA + CCCP + CPZ + Cap + Cephalexin + Cerulenin + Clofazimine + 
    Clotrimazole + DCCD + DIPED + DNP + DTNB + DTT + Deferoxamine + 
    Dipyridyl + EMB + Econazole + Eta + GSNO + GSNO_CFZ + GSNO_CPZ + 
    GSNO_KCN + GSNO_menadione + H2O2 + INH + KCN + Levo + Min_med__succinate_ + 
    Mtm + NRP_1 + NaN3 + Nam + Nigericin + Novobiocin + Oflox + 
    PA_1 + PA_21 + PA824 + PZA + Palmitate + Procept_6776 + Procept_6778 + 
    Rif + Rifp + Rox + SM + Succinate + TLM + TRC + TRZ + Tet + 
    UV + Valinomycin + Verapamil + ZnSO4 + menadione + mercaptoethanol + 
    methoxatin")

## Fit model
#Grid of tuning parameters to try:
grid <- expand.grid(.size=seq(5,21,by=2),.decay=c(0,0.001,0.1))
model <- train(f, train, method='', trace = FALSE, tuneGrid= grid) 
ps <- predict(model, te)

## prediction accuracy
nnout <- compute(nn2, test[,1:4])
results <- nnout$net.result

# Accuracy
original_values <- max.col(test[, 5:7])
preds <- max.col(results)
mean(preds == original_values)








#Examine results
model
plot(y)
lines(ps, col=2)



















