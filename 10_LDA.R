rm(list = setdiff(ls(), lsf.str()))

boshoff  <- read.table("boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)
table(boshoff$Gene)
#str(boshoff)

#################maybe you want to restrict data####################
ii = boshoff[,1]=="I"  | boshoff[,1]=="Q"  | boshoff[,1]=="G" 
X <- boshoff[ii,]
levels(X$Gene)
table(X$Gene)

#NOTE "levels" in X stil remains that of boshoff
#how to actually make "levels" be restricted to the present labels in X$Gene
X$Gene = factor(X$Gene)
levels(X$Gene)
table(X$Gene)

plot(X$BZA, X$Tet , col=X$Gene)

##############################################################
###################  L.D.A  ##################################
##############################################################
#https://medium.com/towards-data-science/is-lda-a-dimensionality-reduction-technique-or-a-classifier-algorithm-eeed4de9953a
#http://sebastianraschka.com/Articles/2014_python_lda.html
#https://tgmstat.wordpress.com/2014/01/15/computing-and-visualizing-lda-in-r/
#https://codesachin.wordpress.com/2015/08/25/linear-and-quadratic-discriminant-analysis-for-ml-statistics-newbies/


################ LDA usin ade4 vs MASS ##################

library(ade4)
pca1 <- dudi.pca(B[, 2:64], scannf = FALSE)
dis1 <- discrimin(pca1, B$Gene, scannf = FALSE)
names(dis1)
dis1
plot(dis1)



library(MASS)
dis2 <- lda(as.matrix(B[, 2:64]), B$Gene)
names(dis2)
dis2
plot(dis2)

#########################################################


##################  LDA in MASS  #######################
library(MASS)

table(X$Gene)
lda <- lda(Gene ~., X, prior = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/17)
lda <- lda(Gene ~., X, prior = c(1,1,1)/3)
lda <- lda(Gene ~., X)
plot(lda)


##determine how well the model fits.
##elements on main diagonal are true predictions, and for e.g. element (1,2) is how many times model classifies class1 as class2 mistakenly
pred <- predict(lda, X[,2:64])$class
table(pred, X[,1])
mean(pred == X[,1])


## or QDA
qda <- qda(Gene ~., X)
pred <- predict(qda, X[,2:64])$class
table(pred, X[,1])
mean(pred == X[,1])


###### Dimensinality Reduction ########
projected_data <- as.matrix(X[,2:64]) %*% lda$scaling  #"coef(lda)" is the same as "lda$scaling"
plot( projected_data, col = X[,1], pch = 19 )



########## lets dig into LDA a little bit more ##########
lda$counts
lda$means  #mean of eachvariable in each class
lda$scaling
lda$svd #the ratio of the between- and within-group standard deviations on the linear discriminant variables.

#We can use the singular values to compute the amount of the between-group variance that is explained by each linear discriminant
prop = lda$svd^2/sum(lda$svd^2)
prop
#the first linear discriminant explains more than {81%} of the between-group variance 

lda.values <- predict(lda)
ldahist(data = lda.values$x[,1], g=X$Gene)
ldahist(data = lda.values$x[,2], g=X$Gene)
#the following line plots the 'predictions', i.e. what lda "thinks" our data look like
#as we might expect, it will probably have nice separation between data points 
plot(lda.values$x[,1], lda.values$x[,2], col=lda.values$class)




##############################################################
###################  Q.D.A  ##################################
##############################################################

## for QDA, the number of examples in each class should be higher than a min value
## thus, I just remove the 4 classes with <100 examples in them
X <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE)
table(X$Gene)
levels(X$Gene)
X <- X[ X[,1]!="D" , ]
X <- X[ X[,1]!="F" , ]
X <- X[ X[,1]!="U" , ]
X <- X[ X[,1]!="V" , ]
table(X$Gene)
X$Gene <- factor(X$Gene)
levels(X$Gene)

qda <- qda(Gene ~., X)

qpred <- predict(qda, X[,2:64])$class
table(qpred, X[,1])
mean(qpred == X[,1])

#QDA is non-linear, idk how to use it for dimensionality reduction
#Idk how to plot(?) it either




#########################################
############# LDA vs QDA ################
#########################################
#removing all variables except functions
rm(list = setdiff(ls(), lsf.str()))

X <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE)
X <- X[ X[,1]!="D" , ]
X <- X[ X[,1]!="F" , ]
X <- X[ X[,1]!="U" , ]
X <- X[ X[,1]!="V" , ]
table(X$Gene)
X$Gene <- factor(X$Gene)
levels(X$Gene)


k <- nrow(X)
s <- sample(k,k/10) #set 20% of data aside for test data
testData <- X[s,]
table(testData$Gene)
trainData <- X[-s,]
table(trainData$Gene)


lda <- lda(Gene ~., trainData)
lpred <- predict(lda, testData[,2:64])$class
table(lpred, testData[,1])
mean(lpred == testData[,1])

#if you get 'some group is too small for qda', then do sampling for test/train again
qda <- qda(Gene ~., trainData)
qpred <- predict(qda, testData[,2:64])$class
table(qpred, testData[,1])
mean(qpred == testData[,1])





