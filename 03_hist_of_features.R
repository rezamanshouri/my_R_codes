

X = read.table("boshoff_ready_4_NB", sep=" ")

for (i in 2:64)
{

    png(paste("hist_of_feature_", i, ".png", sep=""))
    hist(X[,i], breaks=100)
    dev.off()
}
