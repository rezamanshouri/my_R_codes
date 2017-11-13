
## Shannon Entropy
############################################
## v is of type "numeric" of occurrances: e.g. v <- c(1,2,1,3,1,2,3,1,2)
my.entropy <- function (v) {
  f <- table(v)  # #table of label frequencies
  f = f[ f>0 ] #exclude those with count 0 to avoid undefined behavior in "log"
  entropy <- 0
  if ( length(f)==1 ) {
    entropy
  }else {
    p <- f/sum(f) # now p is the list of priors
    entropy <- sum(-p*log(p,length(f)))  
  }
  
  entropy  
}
############################################




## Input : 2 vectors/lists 
#########################################################
###########          correlation        #################
my.cor <- function(a, b) {
  mcor = 0
  # if a or b are all 0 or 1, i.e. stdev=0, we add 2 columns to both: one 1, and one 0. 
  if( sd(a)==0 | sd(b)==0 ) {
    a <- c(a,1,0)
    b <- c(b,1,0)
    mcor <- cor(a,b)
  }else {
    mcor <- cor(a,b)
  }
  
  mcor*mcor # to capture negative cor as well
}


# test
a <- as.numeric(X[X$Gene=="Rv3001c",-1])
b <- as.numeric(X[X$Gene=="Rv1078",-1])
my.cor(a, b)
################################################################




## Input : 2 vectors/lists of bit vectors 
####################################################################
###########      entropy weighted correlation      #################
entropy_weighted_corr_2bit <- function(a, b) {
  
  correl <- my.cor(a,b)
  entrp.a <- my.entropy(as.numeric(a))
  entrp.b <- my.entropy(as.numeric(b))
  ewcor <- correl*(entrp.a + entrp.b)/2
  
  ewcor*ewcor
}


# test
a <- Ess[Ess$Gene=="Rv3001c",-1]
b <- Ess[Ess$Gene=="Rv1078",-1]
entropy_weighted_corr_2bit(a, b)
################################################################

#########################################
## Remove rows with entropy < 0.2
## 1st row is gene name, and the rest are values
remove_entr_leq <- function(X, threshold) {
  remove_idx <- c()
  for(i in 1:nrow(X)) {
    entr <- my.entropy(as.numeric(X[i,-c(1,2,3)]))
    if( entr < threshold) {
      remove_idx <- c(remove_idx, i)
    }
  }
  X[-remove_idx,]
}
#########################################




################################################################
################        KNN Leave 1 Out       #################
## Each row of Train and Test is like 'Rv0001 1 1 0 1 ...', i.e. gene name followed by its values for each feature
## gfm: 'gene_function_map' which is a 2 column DF and is used to look up gene's functions
## returns list of predictions for Test
## Train: all the examples
## It does K-fold CV where K=N, where N is the size of Train
library(infotheo) # for mutinformation()

my.knn.leave1out <- function(Train, k, gfm) {
  
  all.preds <- c()
  num_correct_pred = 0
  
  for( i in 1:nrow(Train)) { # i iterates on 'Test Example' 
    cat("Query: ",as.character(Train[i,1]), as.character(gene_names[gene_names[,1]==as.character(Train[i,1]),2]), ", entropy:", my.entropy(as.numeric(Train[i,-1])) ,"\n")
    distances <- c() # this will contain all the distances/similarities of row i of Train to all rows in Train
    for( j in 1:nrow(Train) ) {
      if( i==j ) { #ignore itself
        next
      }
      #d <- length(which(Train[i,-1] == Train[j,-1])) #bits in common
      #d <- my.similarity(as.numeric(Train[i,-1]), as.numeric(Train[j,-1]))  
      d <- entropy_weighted_corr_2bit(Train[i,-1], Train[j,-1])
      distances[as.character(Train[j,1])] = d
    }
    
    distances <- sort(distances, decreasing = TRUE) # it's decreasing cuz we have correlation (not distance)
    top_k <- names(distances[1:k])
    
    
    cat(k, " NN: \n")
    for( m in 1:length(top_k)) {
      cat(top_k[m], as.character(gene_names[gene_names[,1]==top_k[m],2]),"\t")
    }
    cat("\n")
    
    
    cat("\t")
    # find all functions corresponding to top_k genes
    all_funcs <- c()
    for( g in top_k ) {
      idx <- which(gfm$Gene == g)
      for( l in idx ) {
        cat( as.character(gfm$Function[[l]]),",")
        all_funcs <- c(all_funcs, as.character(gfm$Function[[l]]) )
      }
      cat("\t\t")
    }
    cat("\n")
    
    
    # find majority function (NOTE it could be NON-unique)
    tbl <- table(all_funcs)
    all_maxs <- names(which(tbl == max(tbl)))
    
    
    # determine prediction: if there are ties, pick one at random
    if(FALSE) {
      pred <- ""
      cat("Prediction: ")
      if( length(all_maxs)==1 ) {
        pred <- all_maxs
      }else {
        pred <- "uncertain"
      }
      cat(pred, "\n")
      all.preds <- c(all.preds, pred)
    }
    
    
    # correct pred?
    ii <- which(gfm$Gene == as.character(Train[i,1]))
    actual_funcs <- as.character(gfm[ii,2])
    
    
    cat("Actual_funcs: ")
    for(f in actual_funcs){
      cat(f, " ")
    }
    cat("\n--------------------\n")
    
    
    cmmn <- intersect(all_maxs,actual_funcs)
    if(length(cmmn) > 0){   # correct pred if 1- appears among majority funcs
      for(ff in cmmn) {
        if(tbl[ff]>2) {  # 2- the majority count is at least 2 (Note this condition works only for sanger where func of each gene is unique)
          num_correct_pred = num_correct_pred + 1;
          cat("#of matches: ", tbl[actual_funcs], "\n")
          cat("+++++++++++++++++++++++\n")
          break
        }
      }  
    }
    
  }
  
  cat("Accuracy is: ", (num_correct_pred/nrow(Train))*100, " %\n")
  all.preds
}

################################################################


























X<- read.csv("Documents/data/tnseq_data/original_files/new/2_gumble_hit_sanger_ready", sep = ",", header = T)
dim(X)
X[1:10,1:10]

## remove highly correlated columns (cor > 95 %)
tmp <- cor(X[,-c(1,2,3)])
tmp[upper.tri(tmp)] <- 0 #upper triagnle
diag(tmp) <- 0 # diagonal

XX <- X[,-c(1,2,3)]
XX <- XX[,!apply(tmp,2,function(x) any(x > 0.95))]
XX <- cbind(X[,c(1,2,3)], XX)
dim(XX)
XX[1:10,1:14]
X[1:10,1:14]
setdiff(colnames(X), colnames(XX))


## remove genes with function frquency less than 25
tbl = table(XX$Function)
nn = names(tbl[tbl<25])
XX <- XX[! XX$Function %in% nn,]
XX$Function <- factor(XX$Function)
dim(XX)
XX[1:10,1:14]
table(XX$Function)


## remove the ones with very low entropy
entropies <- c()
for(i in 1:nrow(XX)) {
  ent <- my.entropy(as.numeric(XX[i,-c(1,2,3)]))
  entropies <- c(entropies, ent)
}
hist(entropies, breaks = 20)

dim(XX)
XX <- remove_entr_leq(XX, 0.2)
dim(XX)
XX[1:10, 1:10]




###################
#### Testing  #####
###################
gene_func_map <- read.table("Documents/data/tnseq_data/original_files/new/H37Rv_sanger_level_2.txt", sep = "\t", header = T)

E <- XX[,-c(2,3)]
preds <- my.knn.leave1out(E, 7, gene_func_map)






