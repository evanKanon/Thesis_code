#! create Rscript file which can run with data, outcome, type=binomial as input...

# load glm library
library(glmnet)

# list of results
res<-vector()
rest<-vector()

# read in prepped training dataset
df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_1.expression',row.names=1,sep="\t")
annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_1_disease.classes',sep='\t')
# establish x and y var data
x_var<-t(df)
y_var=annot$Class
# fit lasso and ridge model onto data (ridge will wait a bit)
fit.lasso <- glmnet(x_var, y_var,family = "binomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='binomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')
# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
# find which attributes are maintained
which(fit.lasso$beta[,bt]!=0)
# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]
#- predict with all training input values
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
op<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(df))

res <- c(res, op)
#========================
# TESTING DATA 
#========================

#- read in test data
test_df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_1.expression',row.names=1,sep="\t")
test_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_1_disease.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
common<-merge(pr,test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
#common$Class <- factor(common$Class, levels=c("myelodysplastic syndrome", "healthy"))
opt<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(test_df))

rest <- c(rest, opt)


# read in prepped training dataset
df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_2.expression',row.names=1,sep="\t")
annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_2_disease.classes',sep='\t')
# establish x and y var data
x_var<-t(df)
y_var=annot$Class
# fit lasso and ridge model onto data (ridge will wait a bit)
fit.lasso <- glmnet(x_var, y_var,family = "binomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='binomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')
# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
# find which attributes are maintained
which(fit.lasso$beta[,bt]!=0)
# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]
#- predict with all training input values
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
op<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(df))
# add value to list
res <- c(res, op)

#========================
# TESTING DATA 
#========================

#- read in test data
test_df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_2.expression',row.names=1,sep="\t")
test_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_2_disease.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
common<-merge(pr,test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
#common$Class <- factor(common$Class, levels=c("myelodysplastic syndrome", "healthy"))
#common$`1` <- factor(common$`1`, levels=c("myelodysplastic syndrome", "healthy"))
opt<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(test_df))

rest <- c(rest, opt)

# read in prepped training dataset
df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_3.expression',row.names=1,sep="\t")
annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_3_disease.classes',sep='\t')
# establish x and y var data
x_var<-t(df)
y_var=annot$Class
# fit lasso and ridge model onto data (ridge will wait a bit)
fit.lasso <- glmnet(x_var, y_var,family = "binomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='binomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')
# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
# find which attributes are maintained
which(fit.lasso$beta[,bt]!=0)
# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]
#- predict with all training input values
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
op<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(df))
# add value to list
res <- c(res, op)

#========================
# TESTING DATA 
#========================

#- read in test data
test_df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_3.expression',row.names=1,sep="\t")
test_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_3_disease.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
common<-merge(pr,test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
#common$Class <- factor(common$Class, levels=c("myelodysplastic syndrome", "healthy"))
#common$`1`<-factor(common$`1`, levels=c("myelodysplastic syndrome", "healthy"))
opt<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(test_df))

rest <- c(rest, opt)

# read in prepped training dataset
df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_4.expression',row.names=1,sep="\t")
annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_4_disease.classes',sep='\t')
# establish x and y var data
x_var<-t(df)
y_var=annot$Class
# fit lasso and ridge model onto data (ridge will wait a bit)
fit.lasso <- glmnet(x_var, y_var,family = "binomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='binomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')
# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
# find which attributes are maintained
which(fit.lasso$beta[,bt]!=0)
# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]
#- predict with all training input values
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
op<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(df))
# add value to list
res <- c(res, op)

#========================
# TESTING DATA 
#========================

#- read in test data
test_df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_4.expression',row.names=1,sep="\t")
test_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_4_disease.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
common<-merge(pr,test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
#common$Class <- factor(common$Class, levels=c("myelodysplastic syndrome", "healthy"))
opt<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(test_df))

rest <- c(rest, opt)










