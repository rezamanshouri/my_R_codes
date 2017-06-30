
boshoff  <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE)
table(boshoff$Gene)
#str(boshoff)


#################maybe you want to restrict data####################
ii = boshoff[,1]=="E" | boshoff[,1]=="G"| boshoff[,1]=="J" | boshoff[,1]=="I"
X <- boshoff[ii,]
table(X$Gene)

#NOTE "levels" in X stil remains that of boshoff
#how to actually make "levels" be restricted to the present labels in X$Gene
X$Gene
X$Gene = factor(X$Gene)
X$Gene


##############################################################
#https://pbil.univ-lyon1.fr/R/pdf/course5.pdf
library(ade4)

B <- X
B <- boshoff


###### bivariate scatter plot ######
par(mfcol = c(3, 4))
for (k in 2:5) {
  j0 <- names(B)[k]
  br0 <- seq(min(B[, k]), max(B[, k]), le = 11)
  x0 <- seq(min(B[, k]), max(B[, k]), le = 50)
  for (i in 1:3) {
    i0 <- levels(B$Gene)[i]
    x <- B[B$Gene == i0, j0]
    hist(x, br = br0, proba = T, col = grey(0.8), main = i0,
         xlab = j0)
    lines(x0, dnorm(x0, mean(x), sd(x)), col = "red", lwd = 2)
  }
}  


###### bivariate scatter plot #######
par(mar = c(0, 0, 0, 0))
pan1 <- function(x, y, ...) {
  xy <- cbind.data.frame(x, y)
  s.class(xy, B$Gene, include.ori = F, add.p = T, clab = 1.5,
          col = mycol, cpoi = 2, csta = 0.5)
}
pairs(B[, 2:3], panel = pan1)
pairs(B[,7:10], panel = pan1)



###### 3-dimensional scatter plot #######
library(scatterplot3d)
par(mfrow = c(2, 3)) #number of rows, and number of plots in each row 
mar0 = c(2, 3, 2, 3)

#pick 3 features foe which you want to see 3D scatterplot
scatterplot3d(B[, 2], B[, 3], B[, 4], mar = mar0, color = c("blue",
                                                            "black", "red")[B$Gene], pch = 19)
scatterplot3d(B[, 3], B[, 4], B[, 6], mar = mar0, color = c("blue",
                                                            "black", "red")[B$Gene], pch = 19)
scatterplot3d(B[, 3], B[, 5], B[, 25], mar = mar0, color = c("blue",
                                                            "black", "red")[B$Gene], pch = 19)

#you can use names if you want
scatterplot3d(B$ARP4, B$Amikacin, B$Amp, mar = mar0, color = c("blue",
                                                            "black", "red")[B$Gene], pch = 19)


#LOL! note using class label as a dimension will make 3 distcinct gropus, obviously!
scatterplot3d(B$Gene, B$ARP4, B$Amikacin, mar = mar0, color = c("blue",
                                                                "black", "red")[B$Gene], pch = 19)






########### plot multiple lines with same X axis ##########
##I used this to compare MSE on training vs validation set in NeuralNet 

nn  <- read.table("zz2nnnnn", sep=",", header= TRUE)
plot(range(nn[,1]), range(c(nn[,3],nn[,2])), type='n')
lines(nn[,1], nn[,2], type='l', col='green')
lines(nn[,1], nn[,3], type='l', col='red')
legend(x="topright", lty=c(1,1), c("MSE on training set", "MSE on validation se"), col=c("green","red") )









