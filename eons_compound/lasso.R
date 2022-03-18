#! create Rscript file which can run with data, outcome, type=binomial as input...

# read in prepped training dataset
df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/GSE19429_train_prepped.expression',row.names=1,sep="\t")
annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/GSE19429_train_disease.classes',sep='\t')

#df<- read.csv('/Users/s1242130/Desktop/EONS/lasso/OLD/GSE19429_train_prepped.expression',row.names=1,sep="\t")
#annot<-read.csv('/Users/s1242130/Desktop/EONS/lasso/OLD/GSE19429_train.classes',sep='\t')

# load glm library
library(glmnet)

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

#- predict class participation for model with minimal lambda (with subset)
# this will give attribute name
#predict(fit.lasso, newx = x_var[c(140:150),], type = "class",s=lambdamin)
# this will give log-odds instead of attribute name
#predict(fit.lasso, newx = x_var[c(140:150),], type = "class",s=lambdamin)

#- compare output with actual class participation (with subset)
# save prediction outcome
#pr<-predict(fit.lasso, newx = x_var[c(140:150),], type = "class",s=lambdamin)
# merge with known outcome for easier comparison
#common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
# find if the two columns contain discrepancies
#sum(ifelse(common$`1`==common$Class,0,1))

#- predict with all training input values
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
common<-merge(pr,annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
sum(ifelse(common$`1`==common$Class,0,1))
#100% accurate model

#========================
# TESTING DATA 
#========================

#- read in test data
test_df<- read.csv('/Users/s1242130/Desktop/EONS/input/cancer/GSE19429_test_prepped.expression',row.names=1,sep="\t")
test_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/GSE19429_test_disease.classes',sep='\t')

#- establish x var data
tx_var<-t(test_df)

#- predict class participation with optimal model
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
common<-merge(pr,test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
sum(ifelse(common$`1`==common$Class,0,1))
#85% correct identification


#===========================
#===========================
# Different annotation
#===========================
#===========================

#- read in data
oth_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/GSE19429_train_specimen.classes',sep='\t')


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

#- predict with all training input values
pr_init<-predict(fit.lasso, newx = x_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = x_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
sum(ifelse(common$`1`==common$Class,0,1))
#95% accurate model

#========================
# TESTING DATA 
#========================

#- read in test data
oth_test_annot<-read.csv('/Users/s1242130/Desktop/EONS/input/cancer/GSE19429_test_specimen.classes',sep='\t')

#- predict class participation with optimal model
pr_init<-predict(fit.lasso, newx = tx_var, type = "response",s=lambdamin)
pr<-predict(fit.lasso, newx = tx_var, type = "class",s=lambdamin)
rownames(pr)<-rownames(pr_init)
common<-merge(pr,oth_test_annot,by.x="row.names",by.y='X.samples',all.x=TRUE,suffixes=c(".predicted",".actual"))
sum(ifelse(common$`1`==common$Class,0,1))
#60% correct identification

