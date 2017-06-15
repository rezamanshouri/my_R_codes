
X <- read.table("z_play_tennis_data", sep=",", header= TRUE, row.names=NULL)

entropy <- function(X) {
#f is going to be a list of label frequencies
f <- table(X[,1])
f = f[ f>0 ] #exclude those with count 0 to avoid undefined behavior in "log"
p <- f/sum(f) # now p is the list of priors
sum(-p*log2(p))
}

e <- entropy(X)
cat("\n---------------------\nEntropy: ", e, "\n---------------------\n")
cat("Feature\t","IG\n")

for(i in 2:5) {

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
  cat(i, "\t", IG, "\n")

}
