

## In this file I am implementing KNN when correlation/distance matrix of genes is ready, i.e.
## no calculation of distance is done, and I just look up that matrix to find K NN.


################ Reading Files  ##################
##################################################
## read gene function mapping file
gene_func_map <- read.table("../data/tnseq_data/working_dir/1_gene_function_mapping", sep = ",", header = T)
gene_func_map <- read.table("../data/tnseq_data/original_files/new/H37Rv_sanger_level_2.txt", sep = "\t", header = T)

## read gene*gene distance matrix
X<- read.csv("../data/tnseq_data/original_files/new/result_gumbel_entropy_weighted_correlation_08_25_2017.dat", sep = "\t", header = T)

gene_names <- X[,c(1,2)]

dim(X)
X <- X[,-c(2:5)] # remove unwanted columns
temp <- abs(X[,-1]) ## we are interested in a absolute values
X <- cbind(X[,1], temp)
names(X)[1] <- "Gene"
Dist.matrix <- X

##################################################



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



##################################################################
# Note in "correl.matrix" all the weighted distance are calculated, thus
# I do not do any calculation for distance here. I just take a 
# set of genes, and look up the corresponding row, and find top_k
# nearest neighbor, and use voting among the functions of those
# top_k to provide a prediction. 
##################################################################
## predicts solely based on entropy
## returns a 'list' of predictions.
## correl.matrix is the entropy (correl.matrix) matrix.
## Test is just a list of gene names (included in correl.matrix)
## for printing information, we use 'gene_names' which is global varible defined above
## To look up gene functions, we use 'gene_func_map' which is obtained from 'gene_func_map' by removing {"S","R","A","N"}
knn_prediction_given_dist_matrix <- function(Dist.matrix, Test.idx, k, exclude_Test) {
  
  Test.genes <- Dist.matrix[Test.idx,1]
  if(exclude_Test) {
    # To make Train data, just remove 'Test' genes' column from E! Note E is similarity/distacce matrix
    ii <- idx +1
    Dist.matrix <- Dist.matrix[,-ii]
  }
  
  
  all.preds <- c()
  num_correct_pred = 0
  for(gene in Test.genes) {
    iii <- which(h.table[,1]==gene)
    cat("Query: ",gene, as.character(gene_names[gene_names[,1]==gene,2]), ", entropy: ", my.entropy(as.numeric(h.table[iii,-c(1,2)])),"\n")
    # simply find the row corresponding to gene in E
    idx <- which(Dist.matrix$Gene == gene)
    row <- Dist.matrix[idx,]
    row <- sort(row, decreasing = TRUE) # it's decreasing cuz we have correlation (not distance)
    if(!exclude_Test) {
      row <- row[-which(colnames(row) == gene)] #exclude gene itself
    }
    top_k <- names(row)[2:(k+1)] #start from 2 since 1st col is "Gene"
    
    
    cat("KNN: ")
    for( i in 1:length(top_k)) {
      cat(top_k[i], "_", as.character(gene_names[gene_names[,1]==top_k[i],2]),"\t")
    }
    cat("\n")
    
    
    # cat("\t")
    # find all functions corresponding to top_k genes
    all_funcs <- c()
    for( g in top_k ) {
      idx <- which(gene_func_map$Gene == g)
      for( i in idx ) { # add function(s) of this gene g to "all_funcs"
        cat( as.character(gene_func_map$Function[[i]]),",")
        all_funcs <- c(all_funcs, as.character(gene_func_map$Function[[i]]) )
      }
      cat("\t\t")
    }
    cat("\n")

    
    # find majority function (NOTE it could be NON-unique)
    tbl <- table(all_funcs)
    all_maxs <- names(which(tbl == max(tbl)))
    
    # determine prediction: if there are ties, pick one at random
    #pred <- sample(all_maxs, 1)
    #all.preds <- c(all.preds, pred)
    
    # correct pred?
    ii <- which(gene_func_map$Gene == gene)
    actual_funcs <- as.character(gene_func_map[ii,2])
    
    cat("Actual_funcs: ")
    for(f in actual_funcs){
      cat(f, " ")
    }
    cat("\n--------------------\n")
    
    
    cmmn <- intersect(all_maxs,actual_funcs)
    if(length(cmmn) >2){   # correct pred if 1- appears among majority funcs
      for(ff in cmmn) {
        if(tbl[ff]>2) {  # 2- the majority count is at least 2 (Note this condition works only for sanger where func of each gene is unique)
          num_correct_pred = num_correct_pred + 1;
          cat("#of matches: ", tbl[actual_funcs], "\n")
          cat("+++++++++++++++++++++++\n")
          all.preds <- c(all.preds, ff)
          break
        }
      }  
    }
  }
  
   cat("Accuracy is: ", num_correct_pred, "/", length(Test.genes), "= ",(num_correct_pred/length(Test.genes))*100, " %\n")
  all.preds
}
###################################


###################################################################
## Remove genes with finction "S","R","A","N" (from BOTH row&column)
excld <-c("S","R","A","N")
remove_idx <- c()
for(i in 1:nrow(Dist.matrix)) {
  # find function(s) of this gene and if it is one of R,R,A,N add idx to remove_idx (I remove gene if its funcs are, say, R and L)
  idx <- which(gene_func_map$Gene == as.character(Dist.matrix[i,1]))
  for( j in idx ) {
    if(as.character(gene_func_map$Function[[j]]) %in% excld) {
      remove_idx <- c(remove_idx, i)
      next
    }
  }
}
length(remove_idx)

Dist.matrix_excld_RSAN <- Dist.matrix[-remove_idx,]
# exclude the same genes from columns (Note gene at row i appears at col i+1)
remove_idx_col <- remove_idx +1
Dist.matrix_excld_RSAN <- Dist.matrix_excld_RSAN[,-remove_idx_col]
dim(Dist.matrix_excld_RSAN)
dim(Dist.matrix)
###################################################################

#######################################################
## for Sanger ## categories
excld <-c("VI","V")
remove_idx <- c()
for(i in 1:nrow(Dist.matrix)) {
  # find function(s) of this gene and if it is ONLY {S or R or A or N} add idx to remove_idx
  idx <- which(gene_func_map$Gene == as.character(X[i,1]))
  if(length(idx) == 1) {
    if(as.character(gene_func_map$Function[[idx]]) %in% excld) {
      remove_idx <- c(remove_idx, i)
    }
  }
}
length(remove_idx)

Dist.matrix_excld_VVI <- Dist.matrix[-remove_idx,]
# exclude the same genes from columns (Note gene at row i appears at col i+1)
remove_idx_col <- remove_idx +1
Dist.matrix_excld_VVI <- Dist.matrix_excld_VVI[,-remove_idx_col]
dim(Dist.matrix_excld_VVI)
dim(Dist.matrix)
#######################################################






##########################
##### Lets Test Them #####
##########################

exclude_Test <- FALSE
idx <- c(1:nrow(Dist.matrix_excld_RSAN))

preds <- knn_prediction_given_dist_matrix(Dist.matrix_excld_RSAN, idx, 7, exclude_Test)




#####################################
#####################################

###################################
find_idx_of_non_zero_sd <- function(hit_table, cor_table) {
  ii <- c()
  for(i in 1:nrow(cor_table)) {
    idx <- which(hit_table[,1] == as.character(cor_table[i,1]) )
    if( sd(hit_table[idx,-c(1,2)]) != 0 ) {
      ii <- c(ii, i)
    }
  }
  
  ii
}
###################################

h.table<- read.csv("../data/tnseq_data/original_files/new/result_gumbel_hit_table_binary_08_25_2017.dat", sep = "\t", header = T)
idx <- find_idx_of_non_zero_sd(h.table, Dist.matrix_excld_VVI)
exclude_Test <- FALSE


preds <- knn_prediction_given_dist_matrix(Dist.matrix_excld_VVI, idx, 11, exclude_Test)
table(preds)
actual_classes <- c()
for(i in idx) {
  ii <- which(gene_func_map[,1] == as.character(Dist.matrix_excld_VVI[i,1]))
  actual_classes <- c(actual_classes, as.character(gene_func_map[ii,2]))
}
table(actual_classes)



#################################
## exclude from Test and Train
################################
idx_col <- idx + 1
EE <- Dist.matrix_excld_VVI[idx,c(1,idx_col)]
dim(EE)
iix <- c(1:nrow(EE))

preds <- knn_prediction_given_dist_matrix(EE, iix, 11, exclude_Test)
table(preds)
actual_classes <- c()
for(i in iix) {
  ii <- which(gene_func_map[,1] == as.character(EE[i,1]))
  actual_classes <- c(actual_classes, as.character(gene_func_map[ii,2]))
}
table(actual_classes)
length(preds)
length(actual_classes)
