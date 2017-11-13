
## read gene function mapping file
gene_func_map <- read.table("../data/tnseq_data/working_dir/1_gene_function_mapping", sep = ",", header = T)
gene_func_map <- read.table("../data/tnseq_data/original_files/new/H37Rv_sanger_level_2.txt", sep = "\t", header = T)


# "2_binary_pvalue_sanger_ready" has gene, func, gene_name, <values>. 'V' and 'VI' are already excluded!
X <- read.table("../data/tnseq_data/resampling/2_binary_pvalue_sanger_ready", sep = " ", header = T)
gene_names <- X[,c(1,3)]
dim(X)

## remove genes with function frquency less than 25
tbl = table(X$Function)
nn = names(tbl[tbl<25])
X <- X[! X$Function %in% nn,]
X$Function <- factor(X$Function)
dim(X)





###################################################################
## Remove genes whose finction is only EITHER of {"S","R","A","N"}
## i.e. Do not remove gene if its functions are, say, "S" and "K"
## Step 1: exclude genes from hit table which has only on function, and its function is one of {"S","R","A","N"}
## X is essentiality hit table
## gfm : gene_func_map
exclude_RSNA <- function(X, gfm) {
  excld <-c("S","R","A","N")
  remove_idx <- c()
  for(i in 1:nrow(X)) {
    # find function(s) of this gene and if it is ONLY {S or R or A or N} add idx to remove_idx
    idx <- which(gfm$Gene == as.character(X[i,1]))
    if(length(idx) == 1) {
      if(as.character(gfm$Function[[idx]]) %in% excld) {
        remove_idx <- c(remove_idx, i)
      }
    }
  }
  length(remove_idx)
  
  X_excld_RSAN <- X[-remove_idx,]
  dim(X_excld_RSAN)
  dim(X)
  
  X_excld_RSAN
}

## Step 2: exclude all the rows in 'gene_func_map' which contain "S" or"R"
## X: gene_func_map
exclude_S_R_from_gene_func_map <- function(X) {
  ii <- X$Function=="S" | X$Function=="R" 
  new_gene_func_map <- X[!ii,]
  dim(new_gene_func_map)
  
  new_gene_func_map
}

#######################################################
## for Sanger categories
## Step 1: exclude genes from hit table which has only one function, and its function is "V" or "VI"
## X: essentiality hit table
## gfm : gene_func_map
exclude_V_VI <- function(X, gfm) {
  excld <-c("VI","V")
  remove_idx <- c()
  for(i in 1:nrow(X)) {
    # find function(s) of this gene and if it is ONLY {S or R or A or N} add idx to remove_idx
    idx <- which(gfm$Gene == as.character(X[i,1]))
    if(length(idx) == 1) {
      if(as.character(gfm$Function[[idx]]) %in% excld) {
        remove_idx <- c(remove_idx, i)
      }
    }
  }
  length(remove_idx)
  
  X_excld_VVI <- X[-remove_idx,]
  dim(X_excld_VVI)
  dim(X)
  
  X_excld_VVI
}


## Step 2: exclude all the rows in 'gene_func_map' which contain "V" or"VI"
## X: gene_func_map
exclude_V_VI_from_gene_func_map <- function(X) {
  ii <- X$Function=="V" | X$Function=="VI" 
  new_gene_func_map <- X[!ii,]
  dim(new_gene_func_map)
  
  new_gene_func_map
}

#######################################################






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
  
  correl <- my.cor(as.numeric(a),as.numeric(b))
  entrp.a <- my.entropy(as.numeric(a))
  entrp.b <- my.entropy(as.numeric(b))
  ewcor <- correl*(entrp.a + entrp.b)/2
  
  ewcor*ewcor
}


# test
a <- X[X$Gene=="Rv3001c",-1]
b <- X[X$Gene=="Rv1078",-1]
entropy_weighted_corr_2bit(a, b)
################################################################



############################################
###########   my.similarity   ##############
## input two vectors containing p-values
my.similarity <- function(a,b) {
  m = 0
  for(i in 1:length(a) ) {
    if( max(a[i],b[i]) < 0.5 ) {
      s <- (1/max(a[i],b[i])) * abs(log2(max(a[i],b[i]) + 0.0001)) # smallness
      d <- abs(log10(abs(a[i]-b[i]) + 0.0001))  # difference
      cat(a[i], b[i], " *** ", s, "+", d, "= ", s+d,  "\n")
      m = m + ( s + d )
    }
  }
  m
}
a <- as.numeric(X[1,-c(1,2)])
b <- as.numeric(X[18,-c(1,2)])
my.similarity(a,b)
############################################






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




#########################################
## Remove rows with entropy < 0.2
## 1st row is gene name, and the rest are values
remove_entr_leq <- function(E, threshold) {
  remove_idx <- c()
  for(i in 1:nrow(E)) {
    entr <- my.entropy(as.numeric(E[i,-1]))
    if( entr < threshold) {
      remove_idx <- c(remove_idx, i)
    }
  }
  X[-remove_idx,]
}
#########################################

#########################################
## X: Train/Whatever passed to KNN
## gfm : gene_func_map
find_actual_functions <- function(X, gfm) {
  actual_funcs <- c()
  for(i in 1:nrow(X)) {
    p <- as.character(gfm[gfm$Gene==as.character(EE[i,1]),2])
    actual_funcs <- c(actual_funcs, p)
  }
  
  actual_funcs
}
#########################################

#########################################
## Remove rows with entropy < 0.2
## 1st row is gene name, and the rest are values
remove_genes_with_funct_freq_leq <- function(E, threshold) {
  remove_idx <- c()
  for(i in 1:nrow(E)) {
    entr <- my.entropy(as.numeric(E[i,-1]))
    if( entr < threshold) {
      remove_idx <- c(remove_idx, i)
    }
  }
  X[-remove_idx,]
}
#########################################



#################################
##### Prediction Accuracy  ######
#################################

X[1:10,1:10]
E <- X 
new_gene_func_map <- exclude_V_VI_from_gene_func_map(gene_func_map)
dim(E)

## my.knn.leave1out(Train)
EE <- remove_entr_leq(E, 0.2)
table(E$Function)
table(EE$Function)
dim(EE)
preds <- my.knn.leave1out(EE[,-c(2,3)], 7, new_gene_func_map)




act_funcs <- find_actual_functions(EE, new_gene_func_map)
table(act_funcs)


ee <- c()
for(i in 1:nrow(E)) {
  entr <- my.entropy(as.numeric(E[i,-1]))
  ee <- c(ee, entr)
}
hist(ee, breaks = 50)
ee <- sort(ee)
length(ee)
length(ee[ee >0.1])



