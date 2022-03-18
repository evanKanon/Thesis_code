#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------
# random allocation of classes to patients
#--------------------------
import os
import random
import pandas as pd
from scipy import stats
#--------------------------

def check_dir(this_dir):
	if not os.path.isdir(this_dir):
		try:
			os.mkdir(this_dir)
		except:
			updir = os.path.split(this_dir)[0]
			os.mkdir(updir)
			fix_permissions(updir)
			os.mkdir(this_dir)
	fix_permissions(this_dir)

def read_classes(thisfile):
	''' patient_id \t class \n'''
	f=open(thisfile)
	lines=[x.strip().split('\t') for x in f.readlines()]
	f.close()
	classtypes = lines[0][1:]
	cdic = {line[0]:line[1:] for line in lines[1:]}
	clists = {}
	for i,t in enumerate(classtypes):
		theseclassnames = list(set([x[i] for x in list(cdic.values())]))
		for n in theseclassnames:
			print("classname:", n)
			clists[n] = [x for x in cdic if cdic[x][0]==n]
	return cdic, clists

def fix_permissions(this_path):
	os.system("/bin/chmod 755 %s"%(this_path))

#--------------------------

#- select and read in classes 
classfile = '/mnt/d/input/KFold/GSE19429_test_5_disease.classes'
outputdir = '/mnt/d/output/KFold/eons_test_EONSv0.1.4_default_train_5_disease/'
nodeclasses, classlists = read_classes(classfile)
chosenclasses = list(classlists.keys())

#- set dictionary for evaluation
evaldir = os.path.join(outputdir, 'random/')
check_dir(evaldir)

#- identify folder of networks used for evaluation
if os.path.isdir(outputdir+'test_networks/')==True:
	networkdir=outputdir+'test_networks/'
else:
	#! update this to output to slack
	print('ERROR: No networks folder found')

#- Read-in input
networks={}
for filename in os.listdir(networkdir):
	s=os.path.abspath(networkdir)
	if filename.endswith('.csv'):
		networks[filename.split('.csv')[0]]=pd.read_csv(s+'/'+filename,index_col=0)


network_dict={}
#- assign and calculate random class allocation for each model
for netkey in networks.keys():
	network=networks[netkey]
	#- initialise relevant variable to store resulting class for each patient
	pat_dict={}
	#- for each patient in the network (here taking columns)
	for i in network.keys():
		#- assign a random class to each patient
		pat_dict[i]=random.sample(chosenclasses,1)
	#- initialise sum of correctly identified patients
	finsum=0
	#- score 1 for correct identification and 0 for incorrect
	for pat in pat_dict.keys():
		if pat_dict[pat][0]==nodeclasses[pat][0]:
			finsum+=1
		else:
			finsum+=0
	#- return final quality value as a percentage of correctly identified patients
	finvalue=float(finsum)/len(pat_dict.keys())
	#- save value in the dictionary
	network_dict[netkey]=finvalue



#- save dictionary in csv format
df=pd.DataFrame(list(network_dict.items()), columns=['Network', 'Random_Predictability'])
df.to_csv(evaldir+"Random_scores.csv",index=False)

#- 
ca=pd.read_csv(outputdir+'outcome/IGP_accross.csv')

#- 
stats.ttest_ind(ca['igp_avg_score_test'],df['Random_Predictability'])


