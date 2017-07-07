X  <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE)

###################################################
############## SVM from e1071 package #############
###################################################
# http://horicky.blogspot.de/2012/06/predictive-analytics-neuralnet-bayesian.html

library(e1071)

tune <- tune.svm(Gene~., data=X, gamma=10^(-6:-1), cost=10^(1:4))
summary(tune)

model <- svm(Gene~., data=X, method="C-classification", kernel="radial", probability=T, gamma=0.001, cost=10000)
prediction <- predict(model, X)
table(X$Gene, prediction)
















