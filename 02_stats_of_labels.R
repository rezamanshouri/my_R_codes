
X = read.table("boshoff_ready_4_NB", sep=" ")

ie = X[,2] == "E"





L = c("J", "K", "L", "D", "V", "T", "M", "U", "O", "C", "G", "E", "F", "H", "P", "I", "Q")

cat("func \t min \t\t max \t\t mean \t\t sd\n")


for(i in 10:30) {
  cat("----------------------------------feature ", i , "----------------------------------------\n")

  for(l in L) {
    #mean = colMeans()

    rows = X[,1] == l

    min = min(X[rows,i])
    max = max(X[rows,i])
    mean = mean(X[rows,i])
    sd = sd(X[rows,i])

    cat(l, "\t" ,min, "\t" ,max , "\t" , sd, "\t" ,mean, "\n")
  }

  cat("\n")
}
