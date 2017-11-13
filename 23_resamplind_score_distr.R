

## The goal here was:
## Finding the distributioin of 00, 01, and 11 across all experiments for each possibe pair of genes

# "2_binary_pvalue_sanger_ready" has gene, func, gene_name, <values>. 'V' and 'VI' are already excluded!
Y <- read.table("../data/tnseq_data/resampling/result_resampling_hit_table_merged_pval_09_04_2017.dat", sep = "\t", header = T)
Y[1:10,1:10]
X <- Y[,-c(1,2)]
X[X < .05] <- 0
X[X >= .05] <- 1
X[1:10,1:10]
X <- (X+1)%%2 # essentials are 1 now
X[1:10,1:10]
dim(X)

zz <- 0 # 00
zo <- 0 # 01
oo <- 0 # 11

n <- nrow(X)
for(i in 1:n) {
  for(j in i:n) {
    if( i!= j) {
      m = X[i,] + X[j,]
      zz = zz + length(which(m==0))
      zo = zo + length(which(m==1))
      oo = oo + length(which(m==2)) 
    } 
  }
}

zz
zo
oo
X








