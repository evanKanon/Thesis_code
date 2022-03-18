
# load glm library
library(glmnet)

# set random seed for reproducibility
set.seed(6)

# read in prepped training dataset
df<- read.csv('D:\\input\\GSE19429_train_prepped.expression',row.names=1,sep="\t")
annot<-read.csv('D:\\input\\GSE19429_train_disease.classes',sep='\t')

# establish x and y var data
x_var<-t(df)
y_var=annot$Class

# fit lasso and ridge model onto data (ridge will wait a bit)
fit.lasso <- glmnet(x_var, y_var,family = "binomial") # same as setting alpha=1

# plot lasso beta values as trajectories
plot(fit.lasso, main="Lasso trajectories")

# predict outcome with models
#predict(fit.lasso, newx = x_var[c(140:150),c(1:2500)], type = "response")

# run model with cv in order to identify minimal lambda
fit.cv.lasso <- cv.glmnet(x_var, y_var, family='binomial')

# plot crossvalidated lambda values
plot(fit.cv.lasso,main='Cross-val error curve')

# find which is the minimum lambda
which(fit.lasso$lambda==fit.cv.lasso$lambda.min)

# find which attributes are maintained
which(fit.lasso$beta[,41]!=0)

# save minimal lambda
lambdamin <- fit.lasso$lambda[41]
lambdamin <- fit.cv.lasso$lambda.min

#- predict with all training input values
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
sum(ifelse(common$`1`==common$Class,0,1))
#100% accurate model
# write output to file for python evaluation
write.csv(common,'D:\\PhD\\Results\\LASSO\\GSE19429_training_disease_prediction.csv')


#========================
# TESTING DATA 
#========================

#- read in test data
test_df<- read.csv('D:\\input\\GSE19429_test.expression',row.names=1,sep="\t")
test_annot<-read.csv('D:\\input\\GSE19429_test_disease.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
common<-merge(pr,test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
sum(ifelse(common$`1`==common$Class,0,1))
#85% correct identification
# write output to file for python evaluation
write.csv(common,'D:\\PhD\\Results\\LASSO\\GSE19429_testing_disease_prediction.csv')


#===========================
#===========================
# Different annotation
#===========================
#===========================

#- read in data
oth_annot<-read.csv('D:\\input\\GSE19429_train_specimen.classes',sep='\t')


# establish y var data
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
which(fit.lasso$lambda==fit.cv.lasso$lambda.min)

# find which attributes are maintained
#! following line does not work for multinomial!!
#which(fit.lasso$beta[,32]!=0)

# save minimal lambda
lambdamin <- fit.lasso$lambda[32]
#lambdamin <- fit.lasso$lambda[41]
lambdamin <-fit.cv.lasso$lambda.min

#- predict with all training input values
pr_init<-predict(fit.lasso, newx = x_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
sum(ifelse(common$`1`==common$Class,0,1))
#95% accurate model
#90.6%
# write output to file for python evaluation
write.csv(common,'D:\\PhD\\Results\\LASSO\\GSE19429_training_specimen_prediction.csv')

#========================
# TESTING DATA 
#========================

#- read in test data
oth_test_annot<-read.csv('D:\\input\\GSE19429_test_specimen.classes',sep='\t')

#- predict class participation with optimal model
pr_init<-predict(fit.lasso, newx = tx_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
sum(ifelse(common$`1`==common$Class,0,1))
#60% correct identification
#52.5%
# write output to file for python evaluation
write.csv(common,'D:\\PhD\\Results\\LASSO\\GSE19429_testing_specimen_prediction.csv')

