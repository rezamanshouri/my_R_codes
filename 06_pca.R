library("scatterplot3d")

boshoff  <- read.table("Desktop/boshoff_ready_4_NB", sep=",", header= TRUE)


####for testing purposes, you might want to restrict data in to a few categories only###
ii = boshoff[,1]=="I"  | boshoff[,1]=="K" | boshoff[,1]=="J" | boshoff[,1]=="P"
boshoff <- boshoff[ii,]


X <- boshoff[,2:64]
functions <- boshoff[,1]
###########PCA using prcomp()####################
pca <- prcomp(X, center=TRUE, scale=TRUE)
names(pca)
#pca$center
plot(pca, type="l")
summary(pca)
trunc1 <- as.matrix(X) %*% pca$rotation[, 1:3]  #dimension reduction to 3

######plot
DF <- as.data.frame(trunc1)
DF[,"Functions"] <- functions
head (DF)
plot(DF$PC1, DF$PC2, col = DF$Functions)

####plot on 4 dimensions: first 3 PCs#####
col20 <- palette(rainbow(17)) #generate 20 colors
col20 <- col20[as.numeric(DF$Functions)]
scatterplot3d(DF[,1:3], pch = 16, color = col20)

####biplot
pca$rotation=-pca$rotation
pca$x=-pca$x
biplot (pca , scale =0)

#########################################
############OR do it yourself!###########
X <- scale(X, center=TRUE, scale=TRUE)
R <- cov(X)
ee <- eigen(R)
eval <- ee$values
eval
evec <- ee$vectors

X <- as.matrix(X)
dim(X)
dim(evec)
W <- X %*% evec
trunc <- W[,1:8]
dim(trunc)


###### plot instances on best PSc ######
DF <- as.data.frame(W)
DF[,"Functions"] <- functions
head (DF)
plot(DF$V1, DF$V2, col = DF$Functions)

####plot on 4 dimensions: first 3 PCs#####
col20 <- palette(rainbow(17)) #generate 20 colors
col20 <- col20[as.numeric(DF$Functions)]
scatterplot3d(DF[,1:3], pch = 16, color = col20)


#############################################

#### write truncated data into file ########
write.table(trunc,file="boshoff_pca_33",sep=",", col.names = F, row.names = F)
