# load glm library
library(glmnet)

# list of results
res<-vector()
rest<-vector()

#- read in data
df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_1.expression',row.names=1,sep="\t")
oth_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_1_specimen.classes',sep='\t')

# establish y var data
x_var<-t(df)
y_var=oth_annot$Class

# fit lasso and ridge model onto data (ridge will wait a bit)
fit.lasso <- glmnet(x_var, y_var,family = "multinomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='multinomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')
# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
# find which attributes are maintained
#! following line does not work for multinomial!!
#which(fit.lasso$beta[,32]!=0)

# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]
#lambdamin <- fit.lasso$lambda[41]

#- predict with all training input values
pr_init<-predict(fit.lasso, newx = x_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
op<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(df))

res <- c(res, op)

#========================
# TESTING DATA 
#========================

#- read in test data
test_df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_1.expression',row.names=1,sep="\t")
oth_test_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_1_specimen.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr_init<-predict(fit.lasso, newx = tx_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
#common$Class <- factor(common$Class, levels=levels(common$`1`))
opt<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(test_df))

rest<-c(rest,opt)

#- read in data
df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_2.expression',row.names=1,sep="\t")
oth_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_2_specimen.classes',sep='\t')

# establish y var data
x_var<-t(df)
y_var=oth_annot$Class

# fit lasso and ridge model onto data (ridge will wait a bit)
fit.lasso <- glmnet(x_var, y_var,family = "multinomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='multinomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')
# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
# find which attributes are maintained
#! following line does not work for multinomial!!
#which(fit.lasso$beta[,32]!=0)

# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]
#lambdamin <- fit.lasso$lambda[41]

#- predict with all training input values
pr_init<-predict(fit.lasso, newx = x_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
op<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(df))

res <- c(res, op)

#========================
# TESTING DATA 
#========================

#- read in test data
test_df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_2.expression',row.names=1,sep="\t")
oth_test_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_2_specimen.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr_init<-predict(fit.lasso, newx = tx_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
#common$Class <- factor(common$Class, levels=levels(common$`1`))
opt<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(test_df))

rest<-c(rest,opt)

#- read in data
df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_3.expression',row.names=1,sep="\t")
oth_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_3_specimen.classes',sep='\t')

# establish y var data
x_var<-t(df)
y_var=oth_annot$Class

# fit lasso and ridge model onto data (ridge will wait a bit)
fit.lasso <- glmnet(x_var, y_var,family = "multinomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='multinomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')
# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
# find which attributes are maintained
#! following line does not work for multinomial!!
#which(fit.lasso$beta[,32]!=0)

# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]
#lambdamin <- fit.lasso$lambda[41]

#- predict with all training input values
pr_init<-predict(fit.lasso, newx = x_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
op<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(df))

res <- c(res, op)

#========================
# TESTING DATA 
#========================

#- read in test data
test_df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_3.expression',row.names=1,sep="\t")
oth_test_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_3_specimen.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr_init<-predict(fit.lasso, newx = tx_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
#common$Class <- factor(common$Class, levels=levels(common$`1`))
opt<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(test_df))

rest<-c(rest,opt)

#- read in data
df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_4.expression',row.names=1,sep="\t")
oth_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_4_specimen.classes',sep='\t')

# establish y var data
x_var<-t(df)
y_var=oth_annot$Class

# fit lasso and ridge model onto data (ridge will wait a bit)
fit.lasso <- glmnet(x_var, y_var,family = "multinomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='multinomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')
# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
# find which attributes are maintained
#! following line does not work for multinomial!!
#which(fit.lasso$beta[,32]!=0)

# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]
#lambdamin <- fit.lasso$lambda[41]

#- predict with all training input values
pr_init<-predict(fit.lasso, newx = x_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
op<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(df))

res <- c(res, op)

#========================
# TESTING DATA 
#========================

#- read in test data
test_df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_4.expression',row.names=1,sep="\t")
oth_test_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_4_specimen.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr_init<-predict(fit.lasso, newx = tx_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
#common$Class <- factor(common$Class, levels=levels(common$`1`))
opt<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(test_df))

rest<-c(rest,opt)

#- read in data
df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_5.expression',row.names=1,sep="\t")
oth_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_train_5_specimen.classes',sep='\t')

# establish y var data
x_var<-t(df)
y_var=oth_annot$Class

# fit lasso and ridge model onto data (ridge will wait a bit)
fit.lasso <- glmnet(x_var, y_var,family = "multinomial") # same as setting alpha=1
# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")
# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='multinomial')
# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')
# find which is the minimum lambda
bt<-which(fit.lasso$lambda==fit.cv.lasso$lambda.min)
# find which attributes are maintained
#! following line does not work for multinomial!!
#which(fit.lasso$beta[,32]!=0)

# save minimal lambda
lambdamin <- fit.lasso$lambda[bt]
#lambdamin <- fit.lasso$lambda[41]

#- predict with all training input values
pr_init<-predict(fit.lasso, newx = x_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
op<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(df))

res <- c(res, op)

#========================
# TESTING DATA 
#========================

#- read in test data
test_df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_5.expression',row.names=1,sep="\t")
oth_test_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/KFold/GSE19429_test_5_specimen.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr_init<-predict(fit.lasso, newx = tx_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
#common$`1` <- factor(common$`1`, levels=levels(common$Class))
opt<-1-(sum(ifelse(common$`1`==common$Class,0,1))/ncol(test_df))

rest<-c(rest,opt)
