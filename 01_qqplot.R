

X = read.table("boshoff_ready_4_NB", sep=" ")
iq = X[,1] == "Q"
ii = X[,1] == "I"
iNq = X[,1] != "Q"
iNqi = X[iNq,1] != "I"

ie = X[,1] == "E"
ip = X[,1] == "P"
ih = X[,1] == "H"
ig = X[,1] == "C"
io = X[,1] == "O"
ix = X[,1] == "X"


for (i in 2:64)
{

    png(paste("qqplot_QvsE_feature_", i, ".png", sep=""))
    qqplot(X[iq,i], X[ie,i])
    dev.off()
}
