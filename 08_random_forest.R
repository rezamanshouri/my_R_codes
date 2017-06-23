
### read data
boshoff  <- read.table("Desktop/boshoff_ready_4_NB", sep=",", header= TRUE)
table(boshoff$Gene)
str(boshoff)
##shuffle rows
if(FALSE){
  set.seed(1987)
  g <- runif(nrow(boshoff))
  boshoff <- boshoff[order(g),]
}
###############################################################################################
library(randomForest)
#help(randomForest)
?randomForest

### split data to test and train
s <- sample(2551,850) #pick 850 random numbers (i.e. 1/3 of data for test data)
testData <- boshoff[s,]
table(testData$Gene)
trainData <- boshoff[-s,]
table(trainData$Gene)

rf <- randomForest(Gene ~ ., data = trainData, ntrees = 1000, mtry = 8 ,importance = T )
plot(rf)
print(rf)

### prediction accuracy
pred <- predict(rf, testData, type = 'class') #NOTE: for rpar you need to specify "type = 'class'"
table(testData[,1], pred)
#(i,i) is the correct predictions of class i, the (i,j) is the number of missclassification where class "i" is (miss)classified as "j"
z <- pred == testData[,1]
table(z)
mean(z)

rf$err.rate

###### better plot #####
mycol <- colorRampPalette(c('red','blue','green', "black", "yellow", "brown"))(17)
## OOB is overall error rate
layout(matrix(c(1,2),nrow=1), width=c(4,1)) 
par(mar=c(5,4,4,0)) #No margin on the right side
plot(rf, log="y", col = mycol)
#par(mar=c(5,0,4,2)) #No margin on the left side
plot(c(0,1),type="n", axes=F, xlab="", ylab="")
legend("top", colnames(rf$err.rate),col=mycol,cex=0.8,fill=mycol)










