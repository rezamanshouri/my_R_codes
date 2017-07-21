
X  <- read.table("boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)
##shuffle rows
if(T){
  set.seed(1987)
  g <- runif(nrow(boshoff))
  boshoff <- boshoff[order(g),]
}

###############################################################################################
########################################## rpart ##############################################
### https://gormanalysis.com/decision-trees-in-r-using-rpart/

library(rpart)
###the folowing 2 lines are the same! "." in formula means "all variables not yet used"
dt<-rpart(Gene ~ ARP4 + Amikacin + Amp + Antimycin + Ascididemin_Nat_prod + BZA + CCCP + CPZ + Cap + Cephalexin + Cerulenin + Clofazimine + Clotrimazole + DCCD + DIPED + DNP + DTNB + DTT + Deferoxamine + Dipyridyl + EMB + Econazole + Eta + GSNO + GSNO_CFZ + GSNO_CPZ + GSNO_KCN + GSNO_menadione + H2O2 + INH + KCN + Levo + Min_med__succinate_ + Mtm + NRP_1 + NaN3 + Nam + Nigericin + Novobiocin + Oflox + PA_1 + PA_21 + PA824 + PZA + Palmitate + Procept_6776 + Procept_6778 + Rif + Rifp + Rox + SM + Succinate + TLM + TRC + TRZ + Tet + UV + Valinomycin + Verapamil + ZnSO4 + menadione + mercaptoethanol + methoxatin, data = boshoff, method="class", control = rpart.control(minbucket=100, cp=-1) )
dt<-rpart(Gene ~ ., data = boshoff, method="class" , control = rpart.control(minbucket=100, cp=-1))
dt<-rpart(Gene ~ ., data = boshoff, method="class") #default values
plot(dt)
text(dt, use.n=TRUE, cex=0.8)  #annotate the tree
table(boshoff$Gene)
#the numbers at each node are the number of instances of each class at that node.
#Thus, misclassification rate at a node equals (max_number / sum_of_numbers)


########################
#By default, rpart uses "gini impurty" to select splits when performing classification.
#You can use "information gain" (called information index in their manual) instead by specifying it in the parms parameter by __parms = list(split = 'information')__
#It uses exhaustive search over all possible splits in all possible variables
#There are 2 stopping criteria:
#1- the number of instances in a node must be higher than "minsplit" so thata a split is attempted
#2- If the next best split in growing a tree does not reduce the tree’s overall complexity by a certain amount, rpart will terminate the growing process.
#Complexity measure: a combination of the size of a tree and the ability of the tree to separate the classes
#choosing "cp" to a negative number allows the tree to be full grown, but we will need to prune it later
########################


dt<-rpart(Gene ~ ., data = X, method="class", parms = list(split = 'information'), cp = -1 )
#dt
plot(dt)
text(dt, use.n=TRUE, cex=0.8)  #annotate the tree
#text(dt, pretty=0)

library(rattle)
library(rpart.plot)
fancyRpartPlot(dt)


######## Finding cp for Pruning ########
#When rpart grows a tree it performs 10-fold cross validation on the data. To see the cross validation results use printcp().
#The "rel error" of each iter is the fraction of mislabeled elements in the iter relative to the fraction of mislabeled elements in root.
printcp(dt)
plotcp(dt)
#As a rule of thumb, it’s best to prune a decision tree using the cp of smallest tree that is within one standard deviation of the tree with the smallest xerror.
#in this case, the lowest xerror is:
#8   0.00219106     19   0.88563 0.95355 0.0078374
# ...cp of smnallest tree within one standard deviation... => (0.95355 + 0.0078374) = 0.9613874 
#so we pick cp slightly greater than 0.00292141, that is cp=0.003
#(NOTE: you might get a little different numbes, I guess it is because make a random shuffle of data, but at the end, the cp you get should not be that different)
#NOTE2: for automated cp selection see "16_my_random_forest.R"

############## over-fitting data ###############
### note max_depth is 30, and you might not get "the perfect" tree
t<-rpart(Gene ~ ., data = boshoff, method="class" , control = rpart.control(minbucket = 0, minsplit = 0, cp = -1))


######## prune ########
#As a rule of thumb, it’s best to prune a decision tree using the cp of smallest tree that is within one standard deviation of the tree with the smallest xerror.
pt <- prune(t, cp = 0.003)
fancyRpartPlot(pt)

################ Predict Accuracy ##################
k <- nrow(X)
s <- sample(k,k/5) #set 20% of data aside for test data
testData <- boshoff[s,]
trainData <- boshoff[-s,]

##rebuild the tree on trainData (don't forget!)
dt<-rpart(Gene ~ ., data = trainData, method="class", parms = list(split = 'information'), cp = -1 )
fancyRpartPlot(dt)
printcp(dt)
plotcp(dt)
pt <- prune(dt, cp = 0.0037)
fancyRpartPlot(pt)
text(pt, use.n=TRUE, all=TRUE, cex=0.8)

pred <- predict(pt, testData, type = 'class') #NOTE: for rpar you need to specify "type = 'class'"
table(testData[,1], pred)
#(i,i) is the correct predictions of class i, the (i,j) is the number of missclassification where class "i" is (miss)classified as "j"
z <- pred == testData[,1]
table(z)
mean(z)



############### restrict data to a few classes ################
boshoff  <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE)
ii = boshoff[,1]=="E"  |  boshoff[,1]=="J" |  boshoff[,1]=="K" |  boshoff[,1]=="I"
boshoff <- boshoff[ii,]
#restricting "levels" to only those present in 1st column
boshoff$Gene = factor(boshoff$Gene)

#shuffle rows
set.seed(1987)
g <- runif(nrow(boshoff))
boshoff <- boshoff[order(g),]

nrow(boshoff)
s <- sample(873,290)
testData <- boshoff[s,]
trainData <- boshoff[-s,]

##rebuild the tree on trainData (don't forget!)
dt<-rpart(Gene ~ ., data = trainData, method="class", parms = list(split = 'information'), control = rpart.control(minbucket=5, cp=-1) )
#dt<-rpart(Gene ~ ., data = trainData, method="class") #using default values
fancyRpartPlot(dt)
text(dt, use.n=TRUE, all=TRUE, cex=0.8)
printcp(dt)
plotcp(dt)
pt <- prune(dt, cp = 0.03)
plot(pt)
text(pt, use.n=TRUE, all=TRUE, cex=0.8)

pred <- predict(pt, testData, type = 'class') #NOTE: for rpar you need to specify "type = 'class'"
table(testData[,1], pred)
#(i,i) is the correct predictions of class i, the (i,j) is the number of missclassification where class "i" is (miss)classified as "j"
z <- pred == testData[,1]
table(z)
mean(z)



####################################################################################################################


#############################################################################################################
################################################# tree:tree #################################################
library(tree)
model2<-tree(Gene ~ ., data = boshoff, split = "gini", method="class", control = tree.control(nobs = 2551, mincut = 50))
plot(model2)
text(model2,all=TRUE,cex=0.6) #annotate the tree
model2
############################################################################################################


#############################################################################################################
################################################# party:ctree ##############################################
library(party)
model3 <- ctree(Gene ~ ., data = boshoff, controls=ctree_control(maxdepth=4, mincriterion = 0.80) )
plot(model3)
model3

#no pruning available for ctree


### prediction accuracy
s <- sample(2551,850) #pick 850 random numbers (i.e. 1/3 of data for test data)
testData <- boshoff[s,]
trainData <- boshoff[-s,]
t <- ctree(Gene ~. , data = trainData )
plot(t)
pred <- predict(t, testData)
z <- pred == testData[,1]
table(z)
mean(z) #accuracy (i.e. TRUE/ALL)
############################################################################################################


########## How to Evaluate Accuracy ############
data("iris")
head(iris)

s <- sample(150,50) #pick 50 random numbers (i.e. 1/3 of data for test data)
testData <- iris[s,]
trainData <- iris[-s,]

DT <- ctree(Species ~. , data = trainData )
plot(DT)
pred <- predict(DT, testData)
z <- pred == testData[,5] #compare pred with the "class" column of test data
table(z)
mean(z) #accuracy (i.e. TRUE/ALL)

#no pruning available for ctree

################################################


##### Some nice tutorials #####
if(FALSE) {
http://trevorstephens.com/kaggle-titanic-tutorial/r-part-3-decision-trees/
https://en.wikibooks.org/wiki/Data_Mining_Algorithms_In_R/Classification/Decision_Trees
http://rstatistics.net/decision-trees-with-r/
http://scg.sdsu.edu/ctrees_r/
http://rstatistics.net/decision-trees-with-r/

}


















