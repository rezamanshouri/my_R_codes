
# Use a support vector machine with caret
# We will train 3 different SVM kernels with train() method in caret

#######################
# Note:
# caret uses SVM implementation of 'kernelab' which implements an all-vs-all startegy for multi-class classification.
# All-vs-All strategy: build all N(N-1) pairwise classifiers, where for each class i there are (N-1) classifiers that
#   classify class i as positive and one of other classses as negative. The final prediction of a new example is the 
#   classs with maximum number of wins among its (N-1) classifiers. (same as what you implemented in 15_multiclass_classification_using_binary_classifier).
#######################

# see below link for the tuning paramters of each kernel:
# https://topepo.github.io/caret/available-models.html


library(lattice)
library(ggplot2)
library(kernlab)
library(caret)
library(caTools)

X <- read.table("boshoff_ready_duplicates_removed.csv", sep=",", header= TRUE)


### create train and test sets
set.seed(123)
trainIndex <- createDataPartition(X$Gene, p=.8, list=F)
trainData <- X[trainIndex, ]
testData <- X[-trainIndex, ]
table(testData$Gene)
table(trainData$Gene)

# train data
x <- trainData[,2:64]
y <- trainData[,1]

#Pre-Compute CV folds so we can use the same ones for all models
CV_Folds <- createMultiFolds(y, k = 10, times = 5)

#Fit a Linear SVM
L_model <- train(x,y,method="svmLinear",tuneLength=5,
                 trControl=trainControl(method='repeatedCV',index=CV_Folds, classProbs =  TRUE))

#Fit a Poly SVM
P_model <- train(x,y,method="svmPoly",tuneLength=5,
                 trControl=trainControl(method='repeatedCV',index=CV_Folds, classProbs =  TRUE))

#Fit a Radial SVM
# svmRadial has 2 tuning paramters: 1- tuneLength C, the "cost" of the radial kernel. This parameter controls
# the complexity of the boundary between support vectors. 2- The radial kernel also requires setting a smoothing
# parameter, sigma.
set.seed(1492)
# Use the expand.grid to specify the search space	
grid <- expand.grid(sigma = c(.01, .015, 0.2), tuneGrid = grid,
                    C = c(0.75, 0.9, 1, 1.1, 1.25)
)
R_model <- train(x,y,method="svmRadial",tuneLength=5,
                 trControl=trainControl(method='repeatedCV',index=CV_Folds, classProbs =  TRUE))



########################
#Compare 3 models:
results <- resamples(list(Linear = L_model, Poly = P_model, Radial = R_model))
summary(results)
bwplot(results, metric = "Accuracy")
dotplot(results)
densityplot(results, metric = "Accuracy")


#######################
#Test each model's predictive accuracy Using Area under the ROC curve
x_test <- testData[,2:64]
y_test <- testData[,1]

Linear_preds <- predict(L_model, x_test, type='prob')
colAUC(Linear_preds, y_test, plot=TRUE)

preds <- predict(L_model, testData[,2:64])
table(testData[,1], preds)
mean(testData[,1] == preds)










####################################
##########  test on iris  ##########
####################################

if(FALSE) {
  
  #Choose x and y
  x <- iris[,1:4]
  y <- iris[,5]
  
  #Pre-Compute CV folds so we can use the same ones for all models
  CV_Folds <- createMultiFolds(y, k = 10, times = 5)
  
  #Fit a Linear SVM
  L_model <- train(x,y,method="svmLinear",tuneLength=5,
                   trControl=trainControl(method='repeatedCV',index=CV_Folds, classProbs =  TRUE))
  
  #Fit a Poly SVM
  P_model <- train(x,y,method="svmPoly",tuneLength=5,
                   trControl=trainControl(method='repeatedCV',index=CV_Folds, classProbs =  TRUE))
  
  #Fit a Radial SVM
  R_model <- train(x,y,method="svmRadial",tuneLength=5,
                   trControl=trainControl(method='repeatedCV',index=CV_Folds, classProbs =  TRUE))
  
  #Compare 3 models:
  resamps <- resamples(list(Linear = L_model, Poly = P_model, Radial = R_model))
  summary(resamps)
  bwplot(resamps, metric = "Accuracy")
  densityplot(resamps, metric = "Accuracy")
  
  #Test a model's predictive accuracy Using Area under the ROC curve
  #Ideally, this should be done with a SEPERATE test set
  pSpecies <- predict(L_model,x,type='prob')
  colAUC(pSpecies,y,plot=TRUE)
  
  preds <- predict(L_model, testData[,2:64])
  table(testData[,1], preds)
  mean(testData[,1] == preds)


}














