boshoff  <- read.table("boshoff_ready_4_NB", sep=",", header= TRUE)
table(boshoff$Gene)

###############################################################################
################### Cross Validadtion in carret ###############################
#http://www.milanor.net/blog/cross-validation-for-predictive-analytics-using-r/
#machinelearningmastery.com/how-to-estimate-model-accuracy-in-r-using-the-caret-package/
###############################################################################
library(caret)
library(klaR)

X <- boshoff

# define an 80%/20% train/test split of the dataset
split=0.80
trainIndex <- createDataPartition(X$Gene, p=split, list=FALSE)
data_train <- X[ trainIndex,]
data_test <- X[-trainIndex,]

# train a lda model
lda <- lda(Gene~., data=data_train)

# make predictions
x_test <- data_test[,2:64]
y_test <- data_test[,1]
predictions <- predict(lda, x_test)

# summarize results
confusionMatrix(predictions$class, y_test)
table(predictions$class, y_test)
mean(predictions$class == y_test)


## perfform CV
train_control<- trainControl(method="cv", number=10, savePredictions = TRUE)
grid <- expand.grid(.fL=c(0), .usekernel=c(FALSE))  # fix the parameters of the algorithm

model_fit<- train(Gene ~., data=X, trControl=train_control, method="lda", preProcess=c("pca"))
model_fit
print(model_fit)
plot(model_fit) #if you have tuning parameters

#prediction accuracy on new (test) data
mpred <- predict(model_fit, data_test)
str(mpred)

mprobs <- predict(model_fit, data_test, type = "prob")
head(mprobs)




