
X <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE, row.names=NULL)

#discretize all columns except 1st ("row.names") and 2nd ("Functions" or "Gene") columns according to the following
X[-1] [ X[-1] >= 1] <- 1
X[-1] [ X[-1] <= -1] <- -1
X[-1] [ X[-1]<1 & X[-1]>-1] <- 0


########## entropy function ##########
entropy <- function(X) {
#f is going to be a list of label frequencies
f <- table(X[,1])
f = f[ f>0 ] #exclude those with count 0 to avoid undefined behavior in "log"
p <- f/sum(f) # now p is the list of priors
sum(-p*log2(p))
}
#######################################


########## printing ##########
e <- entropy(X)
cat("\n------------------------------------------
      \nEntropy: ", e,
      "\n------------------------------------------\n")
cat("IG\t\t","Feature\n")
cat("-------\t\t","-------\n")
#######################################


########## Calculate IGs ##########
IG_list <- list()

for(i in 2:64) {
  e <- entropy(X)

  sum = 0
  vals = unique(X[,i])  #taking all unique values of this column
  for(v in vals) {
    ind <- X[,i]==v
    S_v <- X[ind,]

    m <- nrow(S_v)
    n <- nrow(X)
    r = m/n

    entropy_S_v <- entropy(S_v)
    sum = sum + r*entropy_S_v
  }

  IG = e - sum
  IG_list[[paste0(colnames(X)[i], i)]] <- IG
  cat(IG, "\t", i-1, "_", colnames(X)[i], "\n")
}


######### print sorted IGs (ascending) ###########

#IG_list[order(sapply(IG_list,'[[',1))]

#######################################
