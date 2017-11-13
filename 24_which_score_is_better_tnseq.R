



## Goal: find a similarity measure that maximizes whithin group similarity,
##       while minimizes netween groups similarity


## "V" and "VI" has already been excluded in below file
X <- read.table("Documents/data/tnseq_data/working_dir/sanger_categories/2_gumbel_hit_3bit_sanger_ready_no_dups", sep=",", skip=0, header= TRUE, row.names=NULL)
X <- X[,-2]
dim(X)

## remove I.E and IV.K  with <3 examples
X <- X[-which(X$Gene=="I.E"),]
X <- X[-which(X$Gene=="IV.K"),]
X$Gene <- factor(X$Gene)

table(X$Gene)
levels(X$Gene)

ii <- which(X$Gene == "I.A")

scores <- sapply(ii, function(a) sapply(ii, function(b) my.cor(X[a,-1], X[b,-1])))

summary(scores)
describ


library(Hmisc)




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











