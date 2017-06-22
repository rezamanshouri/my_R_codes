

X <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE, row.names=NULL)


######################## Discretization ########################
#discretize all columns except 1st ("row.names") and 2nd ("Functions" or "Gene") columns according to the following
if(FALSE) {
  X[-1] [ X[-1] >= 1] <- 1
  X[-1] [ X[-1] <= -1] <- -1
  X[-1] [ X[-1]<1 & X[-1]>-1] <- 0
}


######################## discretize() ########################
discretize <- function(X) {

}



######################## entropy function ########################
entropy <- function(X) {
  #f is going to be a list of label frequencies
  f <- table(X[,1])
  f = f[ f>0 ] #exclude those with count 0 to avoid undefined behavior in "log"
  p <- f/sum(f) # now p is the list of priors
  sum(-p*log2(p))
}



######################## printing ########################
e <- entropy(X)
cat("\n------------------------------------------
      \nEntropy: ", e,
      "\n------------------------------------------\n")
cat("IG\t\tbest_threshold\t", "n_lower\tmean_lowers\t", "n_higher\tmean_highers\t","Feature\n")
cat("------------\t------------\t------------\t------------\t------------\t------------\t","------------\n")


######################### Calculate IGs ########################
IG_list <- list()
e <- entropy(X)

for(i in 2:64) {

  vals = unique(X[,i]) #taking all unique values of this column
  best_gain <- -1 #keep track of best IG for each j
  best_val <- 0
  for(j in vals) {
    Z <- X[order(X[,i]),] # Z is X when sorted based on column "i"
    Z[,i] [ Z[,i] > j] <- 1
    Z[,i] [ Z[,i] <= j] <- 0

    sum = 0
    for(v in 0:1){
      ind <- Z[,i]==v
      S_v <- Z[ind,]

      m <- nrow(S_v)
      n <- nrow(Z)
      r = m/n

      entropy_S_v <- entropy(S_v)
      sum = sum + r*entropy_S_v
    }

    IG = e - sum
    if(IG > best_gain) {
      best_gain = IG
      best_val <- j
    }

  }


  lower <- X[ X[,i] <= best_val , ] #all rows whose i's column value is less <= best_val
  higher <- X[ X[,i] > best_val , ] #all rows whose i's column value is less > best_val

  IG_list[[paste0(colnames(X)[i], "_f",i)]] <- best_gain
  cat(best_gain, "\t", best_val, "\t", nrow(lower), "\t\t", mean(lower[,i]), "\t", nrow(higher), "\t", mean(higher[,i]), "\t", i-1, "_", colnames(X)[i],"\n")
}


######################## print sorted IGs (ascending) ########################

#IG_list[order(sapply(IG_list,'[[',1))]
