X  <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE)
table(X$Gene)
factor(X$Gene)

#####################################################


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
s <- sample(2551,850) #pick 850 random numbers (i.e. 1/3 of data for test data)
testData <- X[s,]
table(testData$Gene)
trainData <- X[-s,]
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
library(neuralnet)

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
nn <- neuralnet(f, data = trainData, hidden=20, stepmax=1e6, linear.output=FALSE)
#for example, 'hidden=c(5,3)' means 2 hidden layers with 5 and 3 neurons in each, respectively.
# "linear.output=FALSE" says to use activation function
print(nn)
plot(nn)




##### prediction accuracy #####
nnpred <- compute(nn, testData[,2:64])
#cbind(testData[,1], nnpred)


# Put multiple binary output to categorical output
maxidx <- function(arr) {
  return(which(arr == max(arr)))
}
idx <- apply(nnpred, c(1), maxidx)
prediction <- c('C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'O', 'P', 'Q', 'T', 'U', 'V')[idx]
table(prediction, testData$Gene)
mean(prediction == testData$Gene)


########
nn2 <- neuralnet(f, data = trainData, hidden=c(15,15), stepmax=1e6, linear.output=FALSE)
plot(nn2)

nnpred2 <- compute(nn2, testData[,2:64])$net.result

idx2 <- apply(nnpred2, c(1), maxidx)
prediction2 <- c('C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'O', 'P', 'Q', 'T', 'U', 'V')[idx2]
table(prediction2, testData$Gene)
mean(prediction2 == testData$Gene)













