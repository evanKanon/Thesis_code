

# load glm library
library(glmnet)

# set random seed for reproducibility
set.seed(6)

#-------------------
# Set 1
#-------------------

# read in TRAINING dataset
df<- read.csv('D:\\PhD\\input\\KFold\\GSE19429_train_1.expression',row.names=1,sep="\t")
annot<-read.csv('D:\\PhD\\input\\KFold\\GSE19429_train_1_disease.classes',sep='\t')

# establish x and y var data
x_var<-t(df)
y_var=annot$Class


# fit lasso model onto data
fit.lasso <- glmnet(x_var, y_var,family = "binomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
#! savefig!!!


# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='binomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')

# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
#! this is useful only if saved
# find which attributes are maintained
#which(fit.lasso$beta[,bt]!=0)
# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]


#- predict with all training input values
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
# write output to file for python evaluation
write.csv(common,'D:\\PhD\\Results\\LASSO\\KFold_1_training_disease_prediction.csv')

#==============================

#- read in TEST data
test_df<- read.csv('D:\\PhD\\input\\KFold\\GSE19429_test_1.expression',row.names=1,sep="\t")
test_annot<-read.csv('D:\\PhD\\input\\KFold\\GSE19429_test_1_disease.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr_t<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
common_t<-merge(pr_t,test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
write.csv(common_t,'D:\\PhD\\Results\\LASSO\\KFold_1_testing_disease_prediction.csv')

#-------------------
# Set 2
#-------------------

# read in TRAINING dataset
df<- read.csv('D:\\PhD\\input\\KFold\\GSE19429_train_2.expression',row.names=1,sep="\t")
annot<-read.csv('D:\\PhD\\input\\KFold\\GSE19429_train_2_disease.classes',sep='\t')

# establish x and y var data
x_var<-t(df)
y_var=annot$Class


# fit lasso model onto data
fit.lasso <- glmnet(x_var, y_var,family = "binomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
#! savefig!!!


# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='binomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')

# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
#! this is useful only is saved
# find which attributes are maintained
#which(fit.lasso$beta[,bt]!=0)
# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]


#- predict with all training input values
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
# write output to file for python evaluation
write.csv(common,'D:\\PhD\\Results\\LASSO\\KFold_2_training_disease_prediction.csv')

#==============================

#- read in TEST data
test_df<- read.csv('D:\\input\\KFold\\GSE19429_test_2.expression',row.names=1,sep="\t")
test_annot<-read.csv('D:\\input\\KFold\\GSE19429_test_2_disease.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr_t<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
common_t<-merge(pr_t,test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
write.csv(common_t,'D:\\PhD\\Results\\LASSO\\KFold_2_testing_disease_prediction.csv')

#-------------------
# Set 3
#-------------------

# read in TRAINING dataset
df<- read.csv('D:\\PhD\\input\\KFold\\GSE19429_train_3.expression',row.names=1,sep="\t")
annot<-read.csv('D:\\PhD\\input\\KFold\\GSE19429_train_3_disease.classes',sep='\t')

# establish x and y var data
x_var<-t(df)
y_var=annot$Class


# fit lasso model onto data
fit.lasso <- glmnet(x_var, y_var,family = "binomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
#! savefig!!!


# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='binomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')

# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
#! this is useful only is saved
# find which attributes are maintained
#which(fit.lasso$beta[,bt]!=0)
# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]


#- predict with all training input values
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
# write output to file for python evaluation
write.csv(common,'D:\\PhD\\Results\\LASSO\\KFold_3_training_disease_prediction.csv')

#==============================

#- read in TEST data
test_df<- read.csv('D:\\input\\KFold\\GSE19429_test_3.expression',row.names=1,sep="\t")
test_annot<-read.csv('D:\\input\\KFold\\GSE19429_test_3_disease.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr_t<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
common_t<-merge(pr_t,test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
write.csv(common_t,'D:\\PhD\\Results\\LASSO\\KFold_3_testing_disease_prediction.csv')

#-------------------
# Set 4
#-------------------

# read in TRAINING dataset
df<- read.csv('D:\\PhD\\input\\KFold\\GSE19429_train_4.expression',row.names=1,sep="\t")
annot<-read.csv('D:\\PhD\\input\\KFold\\GSE19429_train_4_disease.classes',sep='\t')

# establish x and y var data
x_var<-t(df)
y_var=annot$Class


# fit lasso model onto data
fit.lasso <- glmnet(x_var, y_var,family = "binomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
#! savefig!!!


# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='binomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')

# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
#! this is useful only is saved
# find which attributes are maintained
#which(fit.lasso$beta[,bt]!=0)
# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]


#- predict with all training input values
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
# write output to file for python evaluation
write.csv(common,'D:\\PhD\\Results\\LASSO\\KFold_4_training_disease_prediction.csv')

#==============================

#- read in TEST data
test_df<- read.csv('D:\\input\\KFold\\GSE19429_test_4.expression',row.names=1,sep="\t")
test_annot<-read.csv('D:\\input\\KFold\\GSE19429_test_4_disease.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr_t<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
common_t<-merge(pr_t,test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
write.csv(common_t,'D:\\PhD\\Results\\LASSO\\KFold_4_testing_disease_prediction.csv')

#-------------------
# Set 5
#-------------------

# read in TRAINING dataset
df<- read.csv('D:\\PhD\\input\\KFold\\GSE19429_train_5.expression',row.names=1,sep="\t")
annot<-read.csv('D:\\PhD\\input\\KFold\\GSE19429_train_5_disease.classes',sep='\t')

# establish x and y var data
x_var<-t(df)
y_var=annot$Class


# fit lasso model onto data
fit.lasso <- glmnet(x_var, y_var,family = "binomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
#! savefig!!!


# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='binomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')

# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
#! this is useful only is saved
# find which attributes are maintained
#which(fit.lasso$beta[,bt]!=0)
# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]


#- predict with all training input values
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
# write output to file for python evaluation
write.csv(common,'D:\\PhD\\Results\\LASSO\\KFold_5_training_disease_prediction.csv')

#==============================

#- read in TEST data
test_df<- read.csv('D:\\PhD\\input\\KFold\\GSE19429_test_5.expression',row.names=1,sep="\t")
test_annot<-read.csv('D:\\input\\KFold\\GSE19429_test_5_disease.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr_t<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
common_t<-merge(pr_t,test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
write.csv(common_t,'D:\\PhD\\Results\\LASSO\\KFold_5_testing_disease_prediction.csv')

