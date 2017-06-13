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
ie = X[,1]=="E" ; Y = X[ie,2:64]; Y = data.matrix(Y[1:100,]); my_hm <- heatmap(Y, scale = "row")
    ind = X[,1]==c
    Y <- X[ind,2:64]
    Y <- data.matrix(Y)
    #Y <- data.matrix(Y[1:50,])

    png(paste("heatmap_of_calss_", c, ".png", sep=""))
    my_heatmap <- heatmap.2(Y, scale = "none", col=bluered(100), trace = "none", density.info = "none")
    #my_heatmap <- heatmap(Y, scale = "none", col=my_col)
    #my_heatmap <- heatmap(Y, scale = "row", col=my_col)
    dev.off()

  }
}
#################################################################
