

## In this file I am implementing KNN from the scratch, i.e. I calculate
## the correlation/distance from essentiality hit tables, and then I find
## K NN, and ...


################ Reading Files  ##################
##################################################
## read gene function mapping file
gene_func_map <- read.table("Documents/data/tnseq_data/working_dir/1_gene_function_mapping", sep = ",", header = T)
gene_func_map <- read.table("Documents/data/tnseq_data/original_files/new/H37Rv_sanger_level_2.txt", sep = "\t", header = T)

## read hit table
X<- read.csv("Documents/data/tnseq_data/original_files/new/result_gumbel_hit_table_3bit_08_25_2017.dat", sep = "\t", header = T)
X<- read.csv("Documents/data/tnseq_data/original_files/new/result_gumbel_hit_table_binary_08_25_2017.dat", sep = "\t", header = T)


colnames(X)[1] <- "Gene"
colnames(gene_func_map)[1] <- "Gene"

gene_names <- X[,c(1,2)]
Ess <- X[,-2] # remove col 2
dim(Ess)
##################################################


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






## Input: given two rows of bit-vectors, returns tanimoto score
###############################################################
##################    Tanomoto Score        ###################
tanimoto_score <- function(a, b) {
  # if a or b has majority 1's, then negate them
  if(sum(a) > length(a)/2) {
    a <- (a +1)%%2 # flip bits
  }
  if(sum(b) > length(b)/2) {
    b <- (b +1)%%2 # flip bits
  }  
  
  intrsctn <- sum(a&b)
  unn <- sum(a|b)
  intrsctn/unn
}

# test
a <- Ess[Ess$Gene=="Rv3001c",-1]
b <- Ess[Ess$Gene=="Rv1078",-1]
tanimoto_score(a,b)
################################################################




## Input: given two rows of 3bit-vectors
## I assume: 0:un-essentional, 1:uncertain, 2:essential
## idea: intersection: 2's match, union: non-match 1's AND
###############################################################
###########    Modified Tanomoto Score        #################
modified_tanimoto_score <- function(a, b) {
  # if a or b has majority 2's, then flip 2's to 0, and 0's to 2
  if( length(which(a==2)) > length(a)/2 ) {
    a[which(a==0 | a==2)] <- (a[which(a==0 | a==2)] +2)%%4   # flip 2's to 0's and vice versa
  }
  if( length(which(b==2)) > length(b)/2 ) {
    b[which(b==0 | b==2)] <- (b[which(b==0 | b==2)] +2)%%4   # flip 2's to 0's and vice versa
  }
  
  
  intrsctn <- length(which( (a==2 & b==2) | (a==1 & b==1) ))
  unn <- length( which(a>0 | b>0) )
  intrsctn/unn
}

# test
a <- Ess[Ess$Gene=="Rv0411c",-1]
b <- Ess[Ess$Gene=="Rv1078",-1]
modified_tanimoto_score(a,b)
################################################################



## Input : 2 vectors/lists of 3-bit vectors 
####################################################################
#########   entropy weighted modified Tanimoto score   #############
entropy_weighted_modified_tanomoto_score <- function(a, b) {

  modified_t_score <- modified_tanimoto_score(a,b)
  entrp.a <- my.entropy(as.numeric(a))
  entrp.b <- my.entropy(as.numeric(b))
  score <- modified_t_score*(entrp.a + entrp.b)/2
  
  score
}


# test
a <- Ess[Ess$Gene=="Rv0411c",-1]
b <- Ess[Ess$Gene=="Rv1078",-1]
entropy_weighted_modified_tanomoto_score(a, b)
################################################################






## Shannon Entropy
############################################
## v is of type "numeric" of occurrances: e.g. v <- c(1,2,1,3,1,2,3,1,2)
my.entropy <- function (v) {
  f <- table(v)  # #table of label frequencies
  f = f[ f>0 ] #exclude those with count 0 to avoid undefined behavior in "log"
  p <- f/sum(f) # now p is the list of priors
  sum(-p*log(p,length(v)))
}
############################################




## Input : 2 vectors/lists 
#########################################################
###########          correlation        #################
my.cor <- function(a, b) {
  mcor = 0
  # if a or b are all 0 or 1, i.e. stdev=0, we asseume correlation is zero
  n <- length(a)
  if( sd(a)==0 | sd(b)==0 ) {
    mcor = 0
  }else {
    mcor = cor(t(a), t(b))
  }
  
  mcor*mcor
}


# test
a <- Ess[Ess$Gene=="Rv3001c",-1]
b <- Ess[Ess$Gene=="Rv1078",-1]
my.cor(a, b)
################################################################




## Input : 2 vectors/lists of bit vectors 
####################################################################
###########      entropy weighted correlation      #################
entropy_weighted_corr_2bit <- function(a, b) {
  ewcor = 0
  # if a or b are all 0 or 1, i.e. stdev=0, we asseume correlation is zero
  n <- length(a)
  if( sd(a)==0 | sd(b)==0 ) {
    ewcor <- 0
  }else {
    correl <- cor(t(a), t(b))
    entrp.a <- my.entropy(as.numeric(a))
    entrp.b <- my.entropy(as.numeric(b))
    ewcor <- correl*(entrp.a + entrp.b)/2
  }
  
  ewcor
}


# test
a <- Ess[Ess$Gene=="Rv3001c",-1]
b <- Ess[Ess$Gene=="Rv1078",-1]
entropy_weighted_corr_2bit(a, b)
################################################################



## Input : 2 vectors/lists of 3 bit numbers {-1, 0, 1}
## 0: non-essential, 1: uncertain, 2:essential
## idea is to ignore uncertain as correlation, i.e. if at column i,
## one of them has 0, then just skip it.
####################################################################
###########      entropy weighted correlation      #################
entropy_weighted_corr_3bit_skip_uncertain <- function(a, b) {
  ewcor = 0
  # if stdev=0 for either 'a' or 'b', we asseume correlation is zero
  n <- length(a)
  if( sd(a)==0 | sd(b)==0 ) {
    ewcor <- 0
  }else {
    # remove columns where one of a or b has a 1 (uncertain)
    remove_idx <- which(a==1 | b==1)
    ar <- as.vector(a[-remove_idx])
    br <- as.vector(b[-remove_idx])
    
    correl <- cor(t(ar), t(br))
    entrp.a <- my.entropy(as.numeric(ar))
    entrp.b <- my.entropy(as.numeric(br))
    ewcor <- correl*(entrp.a + entrp.b)/2
  }
  
  ewcor
}


# test
a <- Ess[Ess$Gene=="Rv3001c",-1]
b <- Ess[Ess$Gene=="Rv1078",-1]
entropy_weighted_corr_3bit_skip_uncertain(a, b)
################################################################



## Input : 2 vectors/lists of 3 bit numbers {-1, 0, 1} or {0,1,2}
## idea is to treat uncertain as a state in correlation
####################################################################
###########      entropy weighted correlation      #################
entropy_weighted_corr_3bit_count_uncertain <- function(a, b) {
  ewcor = 0
  # if a or b are all 0 or 1, i.e. stdev=0, we asseume correlation is zero
  n <- length(a)
  if( sd(a)==0 | sd(b)==0 ) {
    ewcor <- 0
  }else {

    
    correl <- cor(t(a), t(b))
    entrp.a <- my.entropy(as.numeric(a))
    entrp.b <- my.entropy(as.numeric(b))
    ewcor <- correl*(entrp.a + entrp.b)/2
  }
  
  ewcor
}


# test
a <- Ess[Ess$Gene=="Rv3001c",-1]
b <- Ess[Ess$Gene=="Rv1078",-1]
entropy_weighted_corr_3bit_count_uncertain(a,b)
################################################################



################################################################
###########        MY K Nearest Neighbor       #################
## Each row of Train and Test is like 'Rv0001 1 1 0 1 ...', i.e. gene name followed by its values for each feature
## gfm: 'gene_function_map' which is a 2 column DF and is used to look up gene's functions
## returns list of predictions for Test

library(infotheo) # for mutinformation()

my.knn <- function(Train, Test, k, gfm) {
  
  all.preds <- c()
  num_correct_pred = 0
  
  for( i in 1:nrow(Test)) {
    cat("Query: ",as.character(Test[i,1]), as.character(gene_names[gene_names[,1]==as.character(Test[i,1]),2]), ", entropy:", my.entropy(as.numeric(Test[i,-1])) ,"\n")
    distances <- c() # this will contain all the distances/similarities of row i of Test to all rows in Train
    for( j in 1:nrow(Train) ) {
      #d <- tanimoto_score(Test[i,-1], Train[j,-1])
      #d <- mutinformation(Test[i,-1], Train[j,-1])
      #d <- entropy_weighted_corr(Test[i,-1], Train[j,-1])
      #d <- entropy_weighted_corr_2bit(Test[i,-1], Train[j,-1])
      #d <- entropy_weighted_corr_3bit_skip_uncertain(Test[i,-1], Train[j,-1])
      #d <- entropy_weighted_corr_3bit_count_uncertain(Test[i,-1], Train[j,-1])
      #d <- my.cor(Test[i,-1], Train[j,-1])
      d <- entropy_weighted_modified_tanomoto_score(Test[i,-1], Train[j,-1])
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
    pred <- sample(all_maxs, 1)
    all.preds <- c(all.preds, pred)
    
    
    # correct pred?
    ii <- which(gfm$Gene == as.character(Test[i,1]))
    actual_funcs <- as.character(gfm[ii,2])
    
    
    cat("Actual_funcs: ")
    for(f in actual_funcs){
      cat(f, " ")
    }
    cat("\n--------------------\n")
    
    
    cmmn <- intersect(all_maxs,actual_funcs)
    if(length(cmmn) > 0){
      num_correct_pred = num_correct_pred + 1;
      cat("--------------------\n")
    }

  }
  
  cat("Accuracy is: ", (num_correct_pred/nrow(Test))*100, " %\n")
  all.preds
}

################################################################

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
      #d <- tanimoto_score(Train[i,-1], Train[j,-1])
      #d <- mutinformation(Train[i,-1], Train[j,-1])
      #d <- entropy_weighted_corr(Train[i,-1], Train[j,-1])
      #d <- entropy_weighted_corr_2bit(Train[i,-1], Train[j,-1])
      #d <- entropy_weighted_corr_3bit_skip_uncertain(Train[i,-1], Train[j,-1])
      #d <- entropy_weighted_corr_3bit_count_uncertain(Train[i,-1], Train[j,-1])
      #d <- my.cor(Train[i,-1], Train[j,-1])
      #d <- entropy_weighted_modified_tanomoto_score(Train[i,-1], Train[j,-1])
      d <- length(which(Train[i,-1] == Train[j,-1])) #bits in common
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
    pred <- ""
    if( length(all_maxs)==1 ) {
      pred <- all_maxs
    }else {
      pred <- "uncertain"
    }
    all.preds <- c(all.preds, pred)
    
    
    # correct pred?
    ii <- which(gfm$Gene == as.character(Train[i,1]))
    actual_funcs <- as.character(gfm[ii,2])
    
    
    cat("Actual_funcs: ")
    for(f in actual_funcs){
      cat(f, " ")
    }
    cat("\n--------------------\n")
    
    
    cmmn <- intersect(all_maxs,actual_funcs)
    if(length(cmmn) > 0){
      num_correct_pred = num_correct_pred + 1;
      cat("--------------------\n")
    }
    
  }
  
  cat("Accuracy is: ", (num_correct_pred/nrow(Train))*100, " %\n")
  all.preds
}

################################################################


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







##########################
##### Lets Test Them #####
##########################
E <- exclude_RSNA(Ess, gene_func_map)  # if COG categories used
new_gene_func_map <- exclude_S_R_from_gene_func_map(gene_func_map)

E <- exclude_V_VI(Ess, gene_func_map)  # if Sanger categories used
new_gene_func_map <- exclude_V_VI_from_gene_func_map(gene_func_map)

dim(E)
dim(Ess)
dim(gene_func_map)

n <- nrow(E)
idx <- sample(n, n/5)
Test <- E[idx,]
Train <- E[-idx,]
dim(Train)
dim(Test)


preds <- my.knn(Train, Test, 7, new_gene_func_map)






dim(E)
EE <- remove_entr_leq_.2(E)

dim(EE)


preds <- my.knn.leave1out(EE, 10, new_gene_func_map)
preds <- my.knn.leave1out(E, 7, new_gene_func_map)

actual_funcs <- find_actual_functions(EE,new_gene_func_map)
tbl <- table(actual_funcs, preds)
mean(actual_funcs == preds)

aplly(tbl,2,sum)
###############################################################














gene_func_map <- read.table("Documents/data/tnseq_data/working_dir/1_gene_function_mapping", sep = ",", header = T)
gene_func_map <- read.table("Documents/data/tnseq_data/original_files/new/H37Rv_sanger_level_2.txt", sep = "\t", header = T)


X<- read.csv("Documents/data/tnseq_data/original_files/new/result_gumbel_hit_table_3bit_08_25_2017.dat", sep = "\t", header = T)
X<- read.csv("Documents/data/tnseq_data/original_files/new/result_gumbel_hit_table_binary_08_25_2017.dat", sep = "\t", header = T)

colnames(X)[1] <- "Gene"
colnames(gene_func_map)[1] <- "Gene"


E <- X[,-2]
E <- exclude_RSNA(Ess, gene_func_map)  # if COG categories used
new_gene_func_map <- exclude_S_R_from_gene_func_map(gene_func_map)

E <- exclude_V_VI(Ess, gene_func_map)  # if Sanger categories used
new_gene_func_map <- exclude_V_VI_from_gene_func_map(gene_func_map)


dim(E)
E <- remove_stdev_0_rows(E)
dim(E)
E <- remove_entr_leq_.2(E)
dim(E)


entropies <- c()
for(i in 1:nrow(E)) {
  ent <- my.entropy(as.numeric(E[i,-1]))
  entropies <- c(entropies, ent)
}
hist(entropies, breaks = 20)

length(entropies)
sorted_entr <- sort(entropies)
length(sorted_entr[sorted_entr>0.02])


#########################################
## Remove rows with all same values
## 1st row is gene name, and the rest are values
remove_stdev_0_rows <- function(X) {
  remove_idx <- c()
  for(i in 1:nrow(X)) {
    if( sd(X[i,-1])==0 ) {
      remove_idx <- c(remove_idx, i)
    }
  }
  X[-remove_idx,]
}
#########################################


#########################################
## Remove rows with entropy < 0.2
## 1st row is gene name, and the rest are values
remove_entr_leq_.2 <- function(X) {
  remove_idx <- c()
  for(i in 1:nrow(X)) {
    entr <- my.entropy(as.numeric(X[i,-1]))
    if( entr < 0.02) {
      remove_idx <- c(remove_idx, i)
    }
  }
  X[-remove_idx,]
}
#########################################








