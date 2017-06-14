#most of this code come from:
#http://www.sthda.com/english/wiki/static-and-interactive-heatmap-in-r-unsupervised-machine-learning


library(gplots)


X <- read.table("boshoff_ready_4_NB", sep=" ",skip=0, row.names=NULL)
#X <- X[order(X$Gene),]

CLASS = c("J", "K", "L", "D", "V", "T", "M", "U", "O", "C", "G", "E", "F", "H", "P", "I", "Q")


my_col<- colorRampPalette(c("red", "white", "blue"))(256)


#simple heatmap
#################################################################
if(TRUE){
  for (c in CLASS)
  {

    ind = X[,1]==c
    Y <- X[ind,3:64]  #, row1 is row.names(?), row2 is Gene
    Y <- data.matrix(Y)
    #Y <- data.matrix(Y[1:50,])

    png(paste("heatmap_of_calss_", c, ".png", sep=""), units="in", width=11, height=8.5, res=300)
    my_heatmap <- heatmap.2(Y, scale = "none", dendrogram = "row", col=bluered(100), trace = "none", density.info = "none")
    #my_heatmap <- heatmap(Y, scale = "none", col=my_col)
    #my_heatmap <- heatmap(Y, scale = "row", col=my_col)
    dev.off()

  }
}
#################################################################


#discretize values: if -1<x<1, then x=1; if x<=-1, then x=-1; if x>=1, then x=1
#################################################################
if(TRUE){
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
