boshoff  <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE)
table(boshoff$Gene)
#str(boshoff)


#################maybe you want to restrict data####################
ii = boshoff[,1]=="E" | boshoff[,1]=="G"| boshoff[,1]=="J"
X <- boshoff[ii,]
table(X$Gene)

#NOTE "levels" in X stil remains that of boshoff
#how to actually make "levels" be restricted to the present labels in X$Gene
X$Gene
X$Gene = factor(X$Gene)
X$Gene

plot(X$BZA, X$Tet , col=X$Gene)

##############################################################
#http://sebastianraschka.com/Articles/2014_python_lda.html
#https://tgmstat.wordpress.com/2014/01/15/computing-and-visualizing-lda-in-r/

library(MASS)
X <- boshoff #comment this if you want restricted data to be used

table(X$Gene)
lda <- lda(Gene ~., X, prior = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/17)
lda <- lda(Gene ~., X, prior = c(1,1,1)/3)
lda <- lda(Gene ~., X)
plot(lda)


lda$counts
lda$means
lda$scaling
lda$svd #the ratio of the between- and within-group standard deviations on the linear discriminant variables.

#We can use the singular values to compute the amount of the between-group variance that is explained by each linear discriminant
prop = lda$svd^2/sum(lda$svd^2)
prop
#the first linear discriminant explains more than {81%} of the between-group variance 

X_mean <- X[,2:64]-lda$means
X_mean <- as.matrix(X_mean)
lda_matrix <- X_mean %*% coef(lda)
dim(lda_matrix)

#add labels to be used in plot as color
lda_matrix <- cbind(X[,1], lda_matrix)
plot(lda_matrix[,2:3] , col=lda_matrix[,1])


lda.values <- predict(lda)
ldahist(data = lda.values$x[,1], g=lda_matrix[,1])
ldahist(data = lda.values$x[,2], g=lda_matrix[,1])
plot(lda.values$x[,1], lda.values$x[,2], col=lda.values$class)



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

########################################################




