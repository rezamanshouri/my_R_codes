
#most of this code come from:
#http://www.sthda.com/english/wiki/static-and-interactive-heatmap-in-r-unsupervised-machine-learning


library(gplots)


X <- read.csv("Documents/data/new_RNASeg_data/nrp/", sep=",",skip=0, row.names=NULL)
#X <- X[order(X$Gene),]

CLASS = c("J", "K", "L", "D", "V", "T", "M", "U", "O", "C", "G", "E", "F", "H", "P", "I", "Q")


my_col<- colorRampPalette(c("red", "white", "blue"))(256)


#simple heatmap
#################################################################
if(TRUE){
  for (c in CLASS)
  {
    
    ind = X[,1]==c
    Y <- X[ind,2:64]
    Y <- data.matrix(Y)
    #Y <- data.matrix(Y[1:50,])
    
    png(paste("heatmap_of_calss_", c, ".png", sep=""), units="in", width=13, height=7, res=300)
    my_heatmap <- heatmap.2(Y, scale = "none", Colv=FALSE, dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
    # remove 'Colv=FALSE' to group on columns as well
    #my_heatmap <- heatmap(Y, scale = "none", col=my_col)
    #my_heatmap <- heatmap(Y, scale = "row", col=my_col)
    dev.off()
    
  }
}
#################################################################














head(iris)
ind <- iris$Species == "virginica"
Y <- iris[ind,1:5]
Y$Species <- factor(Y$Species)
Y$Species
Y <- data.matrix(Y)
heatmap.2(Y[,1:4], scale = "col", Colv=FALSE, dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))

## whole data
heatmap.2(as.matrix(iris[,1:2]), scale = "col", Colv=FALSE, dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
heatmap.2(as.matrix(iris[,1:2]), scale = "col", labRow=iris[,5], Colv=FALSE, dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))


heatmap.2(as.matrix(iris[,3:4]), scale = "col", Colv=FALSE, dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
heatmap.2(as.matrix(iris[,3:4]), scale = "col", Colv=FALSE, labRow =iris[,5], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))

heatmap.2(as.matrix(iris[,1:4]), scale = "col", Colv=FALSE, labRow =iris[,5], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))


############################
Y <- iris
plot(Y$Sepal.Length,Y$Sepal.Width, col=Y$Species) # bad features
plot(Y$Petal.Length,Y$Petal.Width, col=Y$Species) # good features



############################
############################
############################
X <- read.csv("Documents/my_R/boshoff_ready_duplicates_removed.csv", sep=",",skip=0, row.names=NULL)
ind = X[,1]=="C"

Y <- X[,c(1,2,3,10,51,57,59)]
Y <- X[,1:64]
Y <- Y[ind,]
Y$Gene = factor(Y$Gene)

my_heatmap <- heatmap.2(as.matrix(Y[,-1]), scale = "none", Colv=FALSE, labRow =Y[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))


## whole data
X <- read.csv("Documents/my_R/boshoff_ready_duplicates_removed.csv", sep=",",skip=0, row.names=NULL)
my_heatmap <- heatmap.2(as.matrix(X[1:100,2:44]), scale = "none", Colv=FALSE, labRow =X[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
my_heatmap <- heatmap.2(as.matrix(X[100:200,2:64]), scale = "none", Colv=FALSE, labRow =X[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
my_heatmap <- heatmap.2(as.matrix(X[150:250,2:64]), scale = "none", Colv=FALSE, labRow =X[,1], trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))




###########################################################################
# let's remove dendogram, and force same labels to show next to each other
###########################################################################
X <- read.csv("Documents/my_R/boshoff_ready_duplicates_removed.csv", sep=",",skip=0, row.names=NULL)
idx <- X[,1]=='C' | X[,1]=='E' | X[,1]=='Q' | X[,1]=='I' | X[,1]=='J'
Y <- X[idx,]
Y$Gene = factor(Y$Gene)
## take a sample cuz its huge to be showd on heatmap
idx <- createDataPartition(Y$Gene, p=.2, list=F)
Z <- Y[idx, ]
table(Z$Gene)

Z <- Z[order(Z$Gene),]
my_heatmap <- heatmap.2(as.matrix(Z[,2:64]), scale = "none", Rowv="FALSE", labRow =testData[,1], dendrogram = "none", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))






############################
############################
############################
## or lets restrict data into only 2 or 3 classes
X <- read.csv("Documents/my_R/boshoff_ready_duplicates_removed.csv", sep=",",skip=0, row.names=NULL)
ii = boshoff[,1]=="K" | boshoff[,1]=="C"
X <- boshoff[ii,]
X$Gene = factor(X$Gene)
dim(X)
my_heatmap <- heatmap.2(as.matrix(X[1:100,2:64]), scale = "none", Colv=FALSE, labRow =X[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
my_heatmap <- heatmap.2(as.matrix(X[100:200,2:64]), scale = "none", Colv=FALSE, labRow =X[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
my_heatmap <- heatmap.2(as.matrix(X[200:300,2:64]), scale = "none", Colv=FALSE, labRow =X[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))

## wow!! scaling...
my_heatmap <- heatmap.2(as.matrix(X[200:300,c(15,51,57,59)]), scale = "none", Colv=FALSE, labRow =X[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
my_heatmap <- heatmap.2(as.matrix(X[200:300,c(2,51,57,59)]), scale = "none", Colv=FALSE, labRow =X[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
my_heatmap <- heatmap.2(as.matrix(X[200:300,c(2,51,57,59)]), scale = "col", Colv=FALSE, labRow =X[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))


## 2 things:
## 1- if a group of genes (which have same functionality)) don't response to a condition/drug, we can say:
##    a- ignore that condition
##    b- OR, it could actually be an indicator of that class

## I'm saying this because now in the heatmaps we are looking for red/blue areas,
## but maybe "whiter" areas should be carefully examined as well!


## whole data with scaling
X <- read.csv("Documents/my_R/boshoff_ready_duplicates_removed.csv", sep=",",skip=0, row.names=NULL)
my_heatmap <- heatmap.2(as.matrix(X[1:2000,2:64]), scale = "col", Colv=FALSE, labRow =X[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
my_heatmap <- heatmap.2(as.matrix(X[200:400,2:64]), scale = "col", Colv=FALSE, labRow =X[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
my_heatmap <- heatmap.2(as.matrix(X[1000:1500,2:64]), scale = "col", Colv=FALSE, labRow =X[,1], dendrogram = "row", trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))














#discretize values: if -1<x<1, then x=1; if x<=-1, then x=-1; if x>=1, then x=1
#################################################################
if(FALSE){
  for (c in CLASS)
  {
    
    ind = X[,1]==c
    Y1 <- X[ind,3:64]  #, row1 is row.names(?), row2 is Gene
    Y1 <- data.matrix(Y1)
    Y1[Y1 >= 1] <- 1
    Y1[Y1 <= -1] <- -1
    Y1[Y1<1 & Y1>-1] <- 0
    
    png(paste("discretized_heatmap_of_calss_", c, ".png", sep=""), units="in", width=11, height=8.5, res=300)
    my_heatmap <- heatmap.2(Y1, scale = "none", dendrogram = "row", col=bluered(100), trace = "none", density.info = "none")
    dev.off()
    
    
  }
}
#################################################################



################## create heatmap for all files ########################
files <- list.files(path="Documents/data/TIGR_and_Sanger_functionalcategories/func_categories/TIGR/out_files/", pattern="*_boshoff", full.names=T, recursive=FALSE)
lapply(files, function(x) {
  X <- read.csv(x, sep=",", header=T) # load file
  
  if(nrow(X) < 3) {
    return()
  }
  cat(x, "\n")
  Y <- data.matrix(X[,2:64])
  fname <- sub(".*role", "", x)  #extract name of file
  fname <- sub("__boshoff", "", fname)
  png(paste("Documents/my_R/z/heatmap_of_calss_", fname, ".png", sep=""), units="in", width=13, height=7, res=300)
  my_heatmap <- heatmap.2(Y, scale = "none", Colv=FALSE, dendrogram = "row", labRow =X[,1], trace = "none", col=bluered(100), density.info = "none", margin=c(10, 10))
  dev.off()
})

#####################################################################






