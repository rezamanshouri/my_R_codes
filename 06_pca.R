
boshoff  <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE)


####for testing purposes, you might want to restrict data in to a few categories only###
ii = boshoff[,1]=="I"  | boshoff[,1]=="K" | boshoff[,1]=="J"
bb <- boshoff[ii,]
table(bb$Gene)

#NOTE "levels" in X stil remains that of boshoff
#how to actually make "levels" be restricted to the present labels in X$Gene
bb$Gene
bb$Gene = factor(bb$Gene)
bb$Gene

plot(bb$ARP4, bb$DIPED, col=bb$Gene)

###########PCA using prcomp()####################
X <- bb[,2:64]
functions <- bb[,1]

pca <- prcomp(X, center=TRUE, scale=TRUE)
names(pca)
#pca$center
plot(pca, type="l")
summary(pca)
trunc1 <- as.matrix(X) %*% pca$rotation[, 1:3]  #dimension reduction to 3
dim(trunc1)
View(trunc1)

######### plot on PC1&PC2 #########
DF <- as.data.frame(trunc1)
DF[,"Functions"] <- functions
head (DF)
#I found it easier to plot each category separetly on the same plot using points() and defining my own colors
#this way I am able o plot only a desired subset of functions,
#and also I can use better color scheme:))
#Here are 20 distinct colors
#they are:  red       Green        Yellow     Blue     Orange      Purple     Cyan      Magenta     Lime        Pink      Teal      Lavender    Brown      Beige      Maroon      Mint      Olive     Coral       Navy        Grey      White       Black)
mycol = c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")
funcs = c(unique(array(bb[,1])))

first_func <- DF[DF$Functions ==funcs[1],]
plot(first_func$PC1, first_func$PC2, col=mycol[1], pch=19) #plot firs function

n <- length(funcs)
for(i in 2:n) { #add next functions in order to the plot
  f <- DF[DF$Functions ==funcs[i],]
  points(f$PC1, f$PC2, col=mycol[i], pch=19) #Note I am plotting one class at a time, that's why I set col to a color in my list
}
legend(x="topright", legend = funcs, col=mycol[1:n], pch=19)
#########################################


####### plot on PC1&PC2$PC3 ###########
#same way as 2D
library("scatterplot3d")

funcs
first_func <- DF[DF$Functions == funcs[1],]
my3d <- scatterplot3d(first_func[,1:3], color = mycol[1], pch=19, angle = 50) #plot firs function

n <- length(funcs)
for(i in 2:n) { #add next functions in order to the plot
  f <- DF[DF$Functions == funcs[i],]
  #scatterplot3d(f[,1:3], color = mycol[2], pch=19)
  my3d$points3d(f[,1:3], col = mycol[i], pch=19)
}
legend(x="topright", legend = funcs, col=mycol[1:n], pch=19)
#########################################


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
