
X  <- read.table("Documents/my_R/boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)

########## PCA as pre-processing #############
pca <- prcomp(X[,2:64], center=TRUE, scale=TRUE)
plot(pca, type="l")
#summary(pca)
trunc2 <- as.matrix(X[,2:64]) %*% pca$rotation[, 1:2]  #dimension reduction to 3
dim(trunc2)
# pca2 <- cbind(X[,1], as.data.frame(trunc2))
# colnames(pca2)[1] <- "Gene"
# head (pca2)
#########################################

## FIRST look at data yourself, how would you cluster it into 17 classes?
plot(trunc2)
## well, now you might expect already poor results :|

km <- kmeans(trunc2, 10)
table(X$Gene, km$cluster)
plot(trunc2, col=km$cluster)  # this plot shows how km THINKS our data looks like, so there should acuatlly be a nice clustering!
plot(trunc2, col=X$Gene)  # how the data actually looks like


#########################################
## or lets restrict data into only 2 or 3 classes
X <- read.csv("Documents/my_R/boshoff_ready_duplicates_removed.csv", sep=",",skip=0, row.names=NULL)
ii = boshoff[,1]=="K" | boshoff[,1]=="C" | boshoff[,1]=="Q"
X <- boshoff[ii,]
X$Gene = factor(X$Gene)
dim(X)


## use pca and then kmeans on PCs
pca <- prcomp(X[,2:64], center=TRUE, scale=TRUE)
trunc2 <- as.matrix(X[,2:64]) %*% pca$rotation[, 1:20]  #dimension reduction
## FIRST look at data yourself, how would you cluster it into 3 classes?
plot(trunc2)
km <- kmeans(trunc2, 3)
table(X$Gene, km$cluster)
plot(trunc2, col=km$cluster) # what clustering looks like
plot(trunc2, col=X$Gene)  # what the data actually looks like

## use raw data
km <- kmeans(X[,2:64], 3)
table(X$Gene, km$cluster)
plot(X[,2:4], col=km$cluster) # what clustering looks like
plot(X[,2:4], col=X$Gene)  # what the data actually looks like
## plot along "nice" features
plot(X[,c(51,57,59)], col=km$cluster) # what clustering looks like
plot(X[,c(51,57,59)], col=X$Gene)  # what the data actually looks like


km <- kmeans(X[,c(51,57,59)], 3)
table(X$Gene, km$cluster)
plot(X[,c(51,57,59)], col=km$cluster) # what clustering looks like
plot(X[,c(51,57,59)], col=X$Gene)  # what the data actually looks like



#################################
########## on iris ##############
#################################
pca <- prcomp(iris[,1:4], center=TRUE, scale=TRUE)
plot(pca, type="l")
trunc2 <- as.matrix(iris[,1:4]) %*% pca$rotation[, 1:2]  #dimension reduction to 3

## FIRST look at data yourself, how would you cluster it into 3 classes?
plot(trunc2)

km <- kmeans(trunc2, 3)
table(km$cluster, iris$Species)
plot(trunc2, col=km$cluster)  # this plot shows how km THINKS our data looks like, so there should acuatlly be a nice clustering!
plot(trunc2, col=iris$Species) # how the dataset ACTUALLY looks like


## using actual data
km <- kmeans(iris[,1:4], centers=3)
table(iris$Species, km$cluster)
plot(iris[,-5], col=km$cluster)  # this plot shows how km THINKS our data looks like, so there should acuatlly be a nice clustering!
plot(iris[,3:4], col=iris$Species) # how the dataset ACTUALLY looks like

## nicer plots
with(iris, pairs(iris[,-5], col=c(1:3)[km$cluster])) 

library(fpc)
clusplot(iris[,-5], km$cluster) #, color=TRUE, shade=TRUE, labels=2, lines=0)





