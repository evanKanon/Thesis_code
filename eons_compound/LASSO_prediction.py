#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------
import pandas as pd
#--------------------------

#----------------
# Simple LASSO
#----------------

#- disease annotation

f_trainin=[]
f_testin=[]

train=pd.read_csv('/mnt/d/PhD/Results/LASSO/GSE19429_training_disease_prediction.csv')
test=pd.read_csv('/mnt/d/PhD/Results/LASSO/GSE19429_testing_disease_prediction.csv')
#- find correct prediction by LASSO in training
correct_tr=0
for j in range(len(train)):
	if train[train.keys()[2]][j]==train[train.keys()[3]][j]:
		correct_tr+=1
	else:
		pass
f_trainin.append(correct_tr/len(train))
#- find correct prediction by LASSO in testing
correct_ts=0
for y in range(len(test)):
	if test[test.keys()[2]][y]==test[test.keys()[3]][y]:
		correct_ts+=1
	else:
		pass
f_testin.append(correct_ts/len(test))

#- specimen annotation

#f_trainin=[]
#f_testin=[]

train=pd.read_csv('/mnt/d/PhD/Results/LASSO/GSE19429_training_specimen_prediction.csv')
test=pd.read_csv('/mnt/d/PhD/Results/LASSO/GSE19429_testing_specimen_prediction.csv')
#- find correct prediction by LASSO in training
correct_tr=0
for j in range(len(train)):
	if train[train.keys()[2]][j]==train[train.keys()[3]][j]:
		correct_tr+=1
	else:
		pass
f_trainin.append(correct_tr/len(train))
#- find correct prediction by LASSO in testing
correct_ts=0
for y in range(len(test)):
	if test[test.keys()[2]][y]==test[test.keys()[3]][y]:
		correct_ts+=1
	else:
		pass
f_testin.append(correct_ts/len(test))


#-----------------
# KFold LASSO
#-----------------

#- run for disease annotation
trainin=[]
testin=[]
for i in range(1,6):
	train=pd.read_csv('/mnt/d/PhD/Results/LASSO/KFold_'+str(i)+'_training_disease_prediction.csv')
	test=pd.read_csv('/mnt/d/PhD/Results/LASSO/KFold_'+str(i)+'_testing_disease_prediction.csv')
	#- find correct prediction by LASSO in training
	correct_tr=0
	for j in range(len(train)):
		if train[train.keys()[2]][j]==train[train.keys()[3]][j]:
			correct_tr+=1
		else:
			pass
	trainin.append(correct_tr/len(train))
	#- find correct prediction by LASSO in testing
	correct_ts=0
	for y in range(len(test)):
		if test[test.keys()[2]][y]==test[test.keys()[3]][y]:
			correct_ts+=1
		else:
			pass
	testin.append(correct_ts/len(test))

#- save result
myres=pd.DataFrame(list(zip(trainin,testin)),columns=['LASSO_training','LASSO_testing'])
myres.to_csv('/mnt/d/PhD/Results/LASSO_disease.csv')


#- run for other annotation
o_trainin=[]
o_testin=[]
for i in range(1,6):
	train=pd.read_csv('/mnt/d/PhD/Results/LASSO/KFold_'+str(i)+'_training_specimen_prediction.csv')
	test=pd.read_csv('/mnt/d/PhD/Results/LASSO/KFold_'+str(i)+'_testing_specimen_prediction.csv')
	#- find correct prediction by LASSO in training
	correct_tr=0
	for j in range(len(train)):
		if train[train.keys()[2]][j]==train[train.keys()[3]][j]:
			correct_tr+=1
		else:
			pass
	o_trainin.append(correct_tr/len(train))
	#- find correct prediction by LASSO in testing
	correct_ts=0
	for y in range(len(test)):
		if test[test.keys()[2]][y]==test[test.keys()[3]][y]:
			correct_ts+=1
		else:
			pass
	o_testin.append(correct_ts/len(test))

#- other result
o_myres=pd.DataFrame(list(zip(o_trainin,o_testin)),columns=['LASSO_training','LASSO_testing'])
o_myres.to_csv('/mnt/d/PhD/Results/LASSO_specimen.csv')
