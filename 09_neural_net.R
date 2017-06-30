boshoff  <- read.table("Desktop/boshoff_ready_4_NB", sep=",", header= TRUE)
table(boshoff$Gene)
str(boshoff)

#####################################################


###### convert class labels to neumeric values ######
#replace each letter in column 1 (Gene) with a number between 1 and 17
boshoff_n <- within(z, Gene <- as.numeric(factor(Gene, labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) )))


##### min-max normalization ######
##this is a very important step!
#scale the data in the interval [0,1]
for(i in 1:64) {
  min <- min(boshoff_n[,i])
  max <- max(boshoff_n[,i])
  boshoff_n[,i] = (boshoff_n[,i]-min)/(max-min)
}




##### partition data to train and test #####
s <- sample(2551,850) #pick 850 random numbers (i.e. 1/3 of data for test data)
testData <- boshoff_n[s,]
table(testData$Gene)
trainData <- boshoff_n[-s,]
table(trainData$Gene)


##### build NN ######
library(neuralnet)
n <- names(boshoff_n)
f <- as.formula(paste("Gene ~", paste(n[!n %in% "Gene"], collapse = " + ")))
nn <- neuralnet(f, data = trainData, hidden=20, stepmax=1e6, linear.output=FALSE)
#for example, 'hidden=c(5,3)' means 2 hidden layers with 5 and 3 neurons in each, respectively.
plot(nn)

##### prediction accuracy #####
pr <- compute(nn, testData[1:20,2:64])








