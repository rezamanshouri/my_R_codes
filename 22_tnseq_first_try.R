

library(class)
library(caret)
# knn method from class package:
# 1- Euclidean distance
# 2- Classification by majority vote, with ties broken at random.


## essentiality hits
X <- read.csv("Documents/my_R/tnseq_data/2_gumbel_hit_no_dups", sep=",", header= TRUE)
X <- read.csv("Documents/my_R/tnseq_data/2_resampling_hit_no_dups", sep=",", header= TRUE)

X<- read.csv("Documents/data/tnseq_data/working_dir/2_gumbel_hit_zbar_ready_no_dups", sep = "\t", header = T)
X<- read.csv("Documents/data/tnseq_data/working_dir/2_gumbel_hit_3bit_ready_no_dups", sep = "\t", header = T)
X <- X[,-2] # remove col 2
dim(X)


################# if you want:
## Remove rows with all 0 or all 1 in HIT TABLE (Ess) from Entr
remove_idx <- c()
dim(X)
for(i in 1:nrow(X)) {
  msum <- sum(X[i,-1])
  if(msum==0 || msum==84) {
    remove_idx <- c(remove_idx, i)
  }
}
X <- X[-remove_idx,]
dim(X)
###############




################################
#######      KNN      ##########
################################


###############################
## below steps are done to prevent getting errors like "too many ties in knn"
# add a small noise to X
vec <- rnorm(2262,0,0.00001)
Y <- X[,-1]+vec
Y <- cbind(X[,1],Y)
names(Y) <- names(X)
##############################

##### partition data to train and test #####
set.seed(123)
Y <- X
trainIndex <- createDataPartition(Y$Gene, p=.8, list=F)
trainData <- Y[trainIndex, ]
testData <- Y[-trainIndex, ]
table(testData$Gene)
table(trainData$Gene)



### using "class" package
accs <- numeric()
for( i in 1:100) {
  knn <- knn(trainData[,-1], testData[,-1], cl=trainData$Gene, k=i)
  accs <- c(accs, mean(knn == testData$Gene))
}

plot(1-accs, type="l", xlab="K", ylab="Error Rate", main="KNN Error Rate")






################################
####### Decision tree ##########
################################
library(rpart)
library(rattle)
library(rpart.plot)

dt <- rpart(Gene ~ ., data = X, method="class", parms = list(split = 'information'), cp = -1 )
plot(dt)
printcp(dt)
plotcp(dt)

pt <- prune(dt, cp = 0.00293830)

png(filename = "ptree", units = "in", width = 10, height = 10, res = 300)
fancyRpartPlot(pt)
dev.off()

######### Predict Accuracy #############
k <- nrow(X)
set.seed(123)
s <- sample(k,k/5)
testData <- X[s,]
trainData <- X[-s,]

dt<-rpart(Gene ~ ., data = trainData, method="class", parms = list(split = 'information'), cp = -1 )
printcp(dt)
plotcp(dt)
pt <- prune(dt, cp = 0.00370599)
fancyRpartPlot(pt)

pred <- predict(pt, testData, type = 'class') #NOTE: for rpar you need to specify "type = 'class'"
table(testData[,1], pred)
z <- pred == testData[,1]
table(z)
mean(z)




####################################
######### Lets try NB  #############
####################################

# for NB you don't need noise, so use X itself in test/train sets
set.seed(123)
trainIndex1 <- createDataPartition(X$Gene, p=.8, list=F)
trainData1 <- X[trainIndex1, ]
testData1 <- X[-trainIndex1, ]
table(testData1$Gene)
table(trainData1$Gene)

library(e1071)
fit <- naiveBayes(Gene ~., data=trainData1)
pred <- predict(fit, testData1[,-1])
mean(pred == testData1$Gene)




####################################
######### Random Forest  ###########
####################################
library(randomForest)
library(caret)

set.seed(123)
Y <- XX[,-c(1,3)]
trainIndex <- createDataPartition(Y$Function, p=.8, list=F)
trainData <- Y[trainIndex, ]
testData <- Y[-trainIndex, ]
table(testData$Function)
table(trainData$Function)

rf <- randomForest(Function ~ ., data = trainData, ntree = 501 , mtry = 5, importance = T, parms = list(split = 'information') )
print(rf)

### prediction accuracy
pred <- predict(rf, testData, type = 'class') #NOTE: for rpar you need to specify "type = 'class'"
table(testData[,1], pred)
z <- pred == testData[,1]
table(z)
mean(z)













####################################
########## visualization ###########
####################################

library(gplots)

Z <- read.table("../data/tnseq_data/resampling/2_pval_sanger_ready", sep=" ", skip=0, header= TRUE, row.names=NULL)
X <- Z[,-c(1,3)]
dim(X)


## remove all rows with same value on all conditions
remove_idx <- c()
for( c in CLASS ) {
  if( length(which(X$Function==c)) < 10 ) {
    remove_idx <- c(remove_idx, which(X$Function==c))
  }
}
X <- X[-remove_idx, ]
dim(X)
X$Function <- factor(X$Function)
table(X$Function)



## remove all rows with same value on all conditions
remove_idx <- c()
for(i in 1:nrow(X) ) {
  if( sd(X[i,-c(1:3)])==0 ) {
    remove_idx <- c(remove_idx, i)
  }
}
length(remove_idx)
dim(X)
X <- X[-remove_idx, ]
dim(X)



CLASS = c("J", "K", "L", "D", "V", "T", "M", "U", "O", "C", "G", "E", "F", "H", "P", "I", "Q")

CLASS = c("I.A", "I.B", "I.C","I.D","I.E",
          "I.F", "I.G","I.H","I.I","II.A","II.B","II.C","III.A", "III.B", "III.C", "III.D",
           "III.E", "III.F", "I.J", "IV.A", "IV.B", "IV.C","IV.D","IV.E","IV.F", "IV.G","IV.H","IV.I","IV.J","IV.K")



for (c in CLASS){
  
  ind = which(X$Function==c)
  if( length(ind) < 5 ){
    next
  }
  Y <- X[ind,-1]
  Y <- data.matrix(Y)
  #Y <- data.matrix(Y[1:50,])
  
  png(paste("z/heatmap_of_calss_", c, ".png", sep=""), units="in", width=13, height=7, res=300)
  my_heatmap <- heatmap.2(Y, scale = "none", Colv=FALSE, dendrogram = "row",labRow = Z[ind,1] ,trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
  # remove 'Colv=FALSE' to group on columns as well
  #my_heatmap <- heatmap(Y, scale = "none", col=my_col)
  #my_heatmap <- heatmap(Y, scale = "row", col=my_col)
  dev.off()
  
}





#####################################
### remove highly correlated cols ###
#####################################


tmp <- cor(X[,-c(1:3)])
tmp[upper.tri(tmp)] <- 0 #upper triagnle
diag(tmp) <- 0 # diagonal

XX <- X[,-c(1:3)]
XX <- XX[,!apply(tmp,2,function(x) any(x > 0.9))]
XX <- cbind(X[,c(1:3)], XX)
dim(XX)
XX[1:10,1:14]












