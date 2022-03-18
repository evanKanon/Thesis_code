#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------
import os
import shutil
import sys
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import distance
from operator import itemgetter, attrgetter
#from eons_funct_mini import *

#-------------------------------------------------------------

def correlation_function(mydata,myoutputfile,transpose=False,correlation='spearman',threshold=-99,er=False):
	''' Function which correlates different distances between columns of input, scale values to range -1,1 and save matrix output '''
	if transpose==True:
		mydata=mydata.T
	else:
		pass
	#- kenny's suggestion to make the code run quicker
	if correlation=='spearman':
		corr_dat=mydata.rank(ascending=True).corr('pearson')
	else:
		try:
			corr_dat=mydata.corr(correlation)
		except:
			corr_dat=pd.DataFrame(distance.cdist(mydata.T,mydata.T,correlation),columns=mydata.columns,index=mydata.columns)
	#- set all values of the dataframe to zero if they are smaller than or equal to the threshold
	corr_dat[corr_dat<=threshold]=0
	#- find range of values in correlation dataframe
	minimum, maximum = (min(corr_dat.min()), max(corr_dat.max()))
	#- see if values are within expected range and if not scale within range (range is applied this way due to rounding errors)
	if minimum<-1.0001 or maximum>1.0001:
		#! incorporate scalariser functionality into its own function !!
		values = corr_dat.values
		scaler = MinMaxScaler(feature_range=(0,1))
		scaler.fit(values)
		#print(scaler.data_max_)
		#print(scaler.transform(data))
		scaled_data=scaler.transform(values)
		corr_dat=pd.DataFrame(scaled_data,columns=corr_dat.columns,index=corr_dat.index)
		#- this calculates similarity, i.e. 0 would mean no similarity and 1 virtual tautology
		corr_dat=(1-corr_dat)
	#- see if diagonal is 0, hence using distance instead of similarity, and inverse it
	#- the round function is added for the same reason as the 1.0001 above, pandas error inflating values in miniscule levels
	elif round(corr_dat[corr_dat.keys()[0]][corr_dat.index[0]],10)==0 and maximum<1.0001 and minimum>=0:
		corr_dat=(1-corr_dat)
	else:
		corr_dat=corr_dat.abs()
		#- this calculates the distance, i.e. 0 would mean complete closeness and 1 total distance
		#corr_dat=(1-corr_dat).abs()
	#- save dataframe as csv file
	corr_dat.to_csv(myoutputfile+'.csv')
	if er==True:
		max_list=list(corr_dat.max())
		min_list=list(corr_dat.min())
		er=(min(min_list),max(max_list))
		return er

def nda_new(mydataframe,myclasselements,permutations=100):
	#- create random sets of nodes, same leng as subset
	randomclassindex = [random.sample(list(range(len(mydataframe))), len(myclasselements)) for x in range(permutations)]
	randomclasses = []
	#- convert indexes (numbers) to patient names
	for ci_list in randomclassindex:
		group=[]
		for ci in ci_list:
			group.append(mydataframe.index[ci])
		randomclasses.append(group)
	#- add target class elements
	randomclasses.append(myclasselements)
	perm_totals=[]
	#- calculate emp pvalue for each set
	for c,randclass in enumerate(randomclasses):
		nodedist={}
		for node in randclass:
			nodedist[node]=list(mydataframe[node])
		#- sort all values for each key within the dict to facilitate the emp pvalue calculation
		nodedist={x:sorted(nodedist[x]) for x in nodedist}
		#- delete last value of nodedist (here because in adjacency matrices have 1 self-loop)
		for nodekey in list(nodedist.keys()):
			del nodedist[nodekey][-1]
		total_p=0
		#- calculate pvalue
		for i, thisnode in enumerate(randclass):
			for j, othernode in enumerate(randclass): #symmetrical because nsp
				#- this line is added because we want to delete the self loop of 1
				if thisnode==othernode:
					pass
				else:
					w = mydataframe[thisnode][othernode]
					p = empiricalp(w, nodedist[thisnode])
					total_p += -math.log(p)
		if c < permutations: #ie this is a permutation and not the real data
			perm_totals.append(total_p)
	return zscore(total_p, perm_totals)

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

def forgivingmsd(x):
	'''Function which finds mean and standard deviation of list'''
	n, mean, std = len(x), 0, 0
	for a in x:
		try:
			mean = mean + a
		except:
			pass # eg if there's an 'x' in there
	mean = mean / float(n)
	for a in x:
		try:
			std = std + (a - mean)**2
		except:
			pass
	std = math.sqrt(std / float(n-1))
	return mean, std

def zscore(thisitem, thislist):
	m,s = forgivingmsd(thislist)
	if s>0:
		return float(thisitem-m)/s
	else:
		return 0

def empiricalp(thisvalue, empiricalcollection):
	if len(empiricalcollection) == 0:
		return 1
	#empiricalcollection = sorted(empiricalcollection) #collection is alredy sorted, saves time
	return 1-float(bisect_left(empiricalcollection,thisvalue))/len(empiricalcollection)

#-------------------------------------------------------------


#inputdir='/mnt/d/PhD/Results/GSE19429_train_23052020_EONSv0.1.7_default_pca_GSE19429_train_disease/'
#inputdir='/mnt/d/PhD/Results/GSE19429_train_23052020_EONSv0.1.7_default_pca_GSE19429_train_specimen/'
#inputdir='/mnt/d/PhD/Results/KFold/GSE19429_train_1_29052020_EONSv0.1.7_default_pca_GSE19429_train_1_disease/'
#inputdir='/mnt/d/PhD/Results/KFold/GSE19429_train_2_29052020_EONSv0.1.7_default_pca_GSE19429_train_2_disease/'
#inputdir='/mnt/d/PhD/Results/KFold/GSE19429_train_3_29052020_EONSv0.1.7_default_pca_GSE19429_train_3_disease/'
inputdir='/mnt/d/PhD/Results/KFold/GSE19429_train_4_29052020_EONSv0.1.7_default_pca_GSE19429_train_4_disease/'
#inputdir='/mnt/d/PhD/Results/KFold/GSE19429_train_5_29052020_EONSv0.1.7_default_pca_GSE19429_train_5_disease/'
#inputdir='/mnt/d/PhD/Results/KFold/GSE19429_train_5_29052020_EONSv0.1.7_default_pca_GSE19429_train_5_specimen/'
inputdir='/mnt/d/PhD/Results/KFold/GSE19429_train_4_29052020_EONSv0.1.7_default_pca_GSE19429_train_4_specimen/'
#inputdir='/mnt/d/PhD/Results/KFold/GSE19429_train_3_29052020_EONSv0.1.7_default_pca_GSE19429_train_3_specimen/'
#inputdir='/mnt/d/PhD/Results/KFold/GSE19429_train_2_29052020_EONSv0.1.7_default_pca_GSE19429_train_2_specimen/'
inputdir='/mnt/d/PhD/Results/KFold/GSE19429_train_1_29052020_EONSv0.1.7_default_pca_GSE19429_train_1_specimen/'

#testfile='/mnt/d/PhD/input/GSE19429_test.expression'
#testfile='/mnt/d/PhD/input/KFold/GSE19429_test_1.expression'
#testfile='/mnt/d/PhD/input/KFold/GSE19429_test_2.expression'
#testfile='/mnt/d/PhD/input/KFold/GSE19429_test_3.expression'
testfile='/mnt/d/PhD/input/KFold/GSE19429_test_4.expression'
#testfile='/mnt/d/PhD/input/KFold/GSE19429_test_5.expression'
#testfile='/mnt/d/PhD/input/KFold/GSE19429_test_4.expression'
#testfile='/mnt/d/PhD/input/KFold/GSE19429_test_3.expression'
#testfile='/mnt/d/PhD/input/KFold/GSE19429_test_2.expression'
testfile='/mnt/d/PhD/input/KFold/GSE19429_test_1.expression'


#classfile='/mnt/d/PhD/input/'+'_'.join(inputdir.split('/')[-2:-1][0].split('.')[-1:][0].split('_')[3:])+'.classes'
classfile='/mnt/d/PhD/input/KFold/'+'_'.join(inputdir.split('/')[-2:-1][0].split('.')[-1:][0].split('_')[3:])+'.classes'
classfile=classfile.replace('train','test')

#- create output directory
outputdir=inputdir.split('_')[0][:-8]+'eons_test_'+'_'.join([x for i,x in enumerate(inputdir.split('_')[4:]) if i!=3])
if not os.path.isdir(outputdir):
	os.mkdir(outputdir)
#- create folder to store test features results
featuredir=outputdir+'test_features/'
if not os.path.isdir(featuredir):
	os.mkdir(featuredir)
#- create folder to store test network results
networkdir=outputdir+'test_networks/'
if not os.path.isdir(networkdir):
	os.mkdir(networkdir)
#- create folder to store evaluation results
evaldir=outputdir+'test_evaluation/'
if not os.path.isdir(evaldir):
	os.mkdir(evaldir)
#- create folder for comparison
outcomedir=outputdir+'outcome/'
if not os.path.isdir(outcomedir):
	os.mkdir(outcomedir)


#- find name of evaluation from folder name
o_evaldir=inputdir+'evaluation_'+'_'.join(inputdir.split('/')[-2:-1][0].split('.')[-1:][0].split('_')[3:])


#- select models from training results
models=list(pd.read_csv(inputdir+'results/best_networks.txt',sep=' ',header=None,index_col=0)[1])

#- find the best networks from the evaluation
evaluation_best=o_evaldir+'/top_evaluation_compound.txt'
e_best=list(pd.read_csv(evaluation_best,sep='\t',header=None)[1])


#- select which features are pertinent
pfeatures_list=[]
for inp in models:
	try:
		pfeatures_list.append(inp.split('expression')[0])
	except:
		print('Check manually: '+inp)


#- keep only unique elements
features_list=[]
for x in pfeatures_list: 
	if x not in features_list: 
		features_list.append(x)


#- for each feature read in the respective feature file from the original folder to keep its index and subset the test features
for feature in features_list:
	try:
		orig_feat=pd.read_csv(inputdir+'features/'+feature+'expression',sep='\t',index_col=0)
		full_test=pd.read_csv(testfile,sep='\t',index_col=0)
		mysubset=full_test.loc[orig_feat.index]
		mysubset.to_csv(featuredir+feature+'expression',sep='\t')
	except:
		if '_'.join(inputdir.split('/')[-2:-1][0].split('.')[0].split('_')[:-2]) in feature:
			shutil.copy(testfile, featuredir)
		else:
			print('Unknown feature, please investigate: '+feature)




#- for each model, create the appropriate one here
for m in models:
	try:
		correlation=m.split('_')[-2:][1].split('-9900')[0]
		mynet=pd.read_csv(featuredir+m.split('expression')[0]+'expression',sep='\t',index_col=0)
		myoutputfile=networkdir+m.split('expression')[0]+'expression'+'_'+correlation+"-9900"
		correlation_function(mynet,myoutputfile=myoutputfile,correlation=correlation)
	except:
		#if '_'.join(inputdir.split('/')[-2:-1][0].split('.')[0].split('_')[:-2]) in m:		
		if 'train' in m:
			f=m.replace("train", "test")
			correlation=f.split('_')[-2:][1].split('-9900')[0]
			mynet=pd.read_csv(testfile,sep='\t',index_col=0)
			myoutputfile=networkdir+f.split('expression')[0]+'expression'+'_'+correlation+"-9900"
			correlation_function(mynet,myoutputfile=myoutputfile,correlation=correlation)




#=========
#==========
#=========


#- inputing classfile
nodeclasses, classlists = read_classes(classfile)

#- select which classes we will compute NDA scores for
#! make this user optionable too!!
chosenclasses = [] #['gram-positive_sepsis', 'gram-negative_sepsis']
chosenclasses = list(set(chosenclasses) & set(classlists.keys())) # classes only count if they were found in the class file
if len(chosenclasses)==0:
	chosenclasses = list(classlists.keys())
print("chosenclasses", chosenclasses)


#- Read-in input
networks={}
for filename in os.listdir(networkdir):
	s=os.path.abspath(networkdir)
	if filename.endswith('.csv'):
		networks[filename.split('.csv')[0]]=pd.read_csv(s+'/'+filename,index_col=0)


#- initiate best dictionary where the names of the top networks from each evaluation will be saved in 
bestdict={}


####################################################################################
# AVERAGE NETWORK NEIGHBOUR
####################################################################################


#- initialise variables for storage of final comparison result
network_dict={}
#- for each network in the networks produced by the analysis
for netkey in list(networks.keys()):
	network=networks[netkey]
	#- initialise relevant variable to store resulting class for each patient
	pat_dict={}
	#- for each patient in the network (here taking columns)
	for i in network.keys():
		#- initialise variables to store maximum sum, class which exhibits it and 
		maxsum=0
		maxclass=None
		#- for each class in the desired classes
		for myclass in chosenclasses:
			#- consider the relevant patient list without including the current patient
			patlist = [x for x in classlists[myclass] if x != i]
			#- calculate the average sum of edgeweights in this subset
			try:
				mysum = network[network.index.isin(patlist)][i].sum()/len(patlist)
			except:
				mysum = 0
			#- see if new value is better than previous maximum, if so update
			if mysum > maxsum:
				maxsum=mysum
				maxclass=myclass
		pat_dict[i]=maxclass
	#- initialise sum of correctly identified patients
	finsum=0
	#- score 1 for correct identification and 0 for incorrect
	for pat in pat_dict.keys():
		if pat_dict[pat]==nodeclasses[pat][0]:
			finsum+=1
		else:
			finsum+=0
	#- return final quality value as a percentage of correctly identified patients
	finvalue=float(finsum)/len(pat_dict.keys())
	#- save value in the dictionary
	network_dict[netkey]=finvalue

#- save dictionary in csv format
#average network neighbour
df=pd.DataFrame(list(network_dict.items()), columns=['Network', 'ANN_Predictability'])
df.to_csv(evaldir+"ANN_evaluation_scores.csv",index=False)

#- save best result from ANN evaluation
bestdict['ANN']=[max(network_dict,key=network_dict.get),network_dict[max(network_dict,key=network_dict.get)]]
#max(network_dict,key=network_dict.get)             # find name of network with highest percentage of patients correctly identified

#- save dictionary in txt format
#dict = {'Python' : '.py', 'C++' : '.cpp', 'Java' : '.java'}
#f = open("dict.txt","w")
#f.write( str(dict) )
#f.close()

'''

####################################################################################
# NODE DENSITY ANALYSIS
####################################################################################

#- COMPUTE NDA - NEW (CORE FUNCTIONALITY)
#networkscores={}
#for network in networks.keys():
#	allscores=[]
#	for thisclass in chosenclasses:
#		nda_score = nda_new(networks[network], classlists[thisclass])
#		allscores.append(nda_score)
#	av_score=forgivingmsd(allscores)[0]
#	networkscores[network]=av_score


#- COMPUTE NDA - NEW
networkscores={}
allnetscores={}
net=[]
for thclass in chosenclasses:
	net.append(thclass)
for network in networks.keys():
	allscores=[]
	#- class header used as sanity check
	class_header=[]
	for thisclass in chosenclasses:
		class_header.append(thisclass)
		nda_score = nda_new(networks[network], classlists[thisclass])
		allscores.append(nda_score)
	av_score=forgivingmsd(allscores)[0]
	networkscores[network]=av_score
	if class_header==net:
		allnetscores[network]=allscores
	else:
		print('Problem with header! Evaluate Manually!')


#- sort average scores by stregth
sorted_x=sorted(networkscores.items(),key=operator.itemgetter(1),reverse=True)


#- plot results to see the distribution of NDA scores
plt.bar(networkscores.keys(), networkscores.values(), color='g')
plt.savefig(evaldir+'NDA_scores.png')
plt.gcf().clear()

#- save results as CSV file
nda_file=pd.DataFrame.from_dict(networkscores,orient='index',columns=['NDA_avg_scores'])
nda_file.to_csv(evaldir+'NDA_avg_evaluation_scores.csv')
#with open(evaldir+'NDA_avg_evaluation_scores.csv', 'w') as f:
#	writer = csv.writer(f)
#	f.write('Network,avg_NDA_score\n')
#	for row in networkscores.items():
#		writer.writerow(row)

#- save sorted results as text file
with open(evaldir+'NDA_avg_evaluation_scores_sorted.txt', 'w') as nf:
	nf.write('Network\tavg_NDA_score\n')
	nf.write('\n'.join('%s %s' % x for x in sorted_x))

#- create df to save full list of NDA scores
ndf=pd.DataFrame.from_dict(allnetscores,orient='index',columns=net)
ndf.to_csv(evaldir+'NDA_all_evaluation_scores.csv')

#- save best result from NDA evaluation
bestdict['NDA']=sorted_x[0]

'''

####################################################################################
# IN GROUP PROPORTION
####################################################################################


#- In Group Proportion (IGP)
#- i.e., how many patients have as their nearest neighbour a patient of the same class
#- idea taken from Kapp and Tibshirani Biostatistics 2007 8(1):9-31


network_dict={}
#- for each network in the networks produced by the analysis
for netkey in networks.keys():
	network=networks[netkey]
	classdict={}
	#- for each class in the given classes
	for myclass in chosenclasses:
		thislist=classlists[myclass]
		#- for each patient in this class
		acc=0
		for patkey in thislist:
			#- identify its nearest neighbour
			for pj in range(len(network[patkey])):
				#- find nearest neighbour that is not hte same patient
				if network[patkey].nlargest(10,keep='all').index[pj]!=patkey:
					nearest_neighbour=network[patkey].nlargest(10,keep='all').index[pj]
					break
			#- check if the nearest neighbour is of the same class
			if nodeclasses[nearest_neighbour]==nodeclasses[patkey]:
				acc+=1
		ig_p=float(acc)/len(thislist)
		classdict[myclass]=ig_p
	network_dict[netkey]=classdict

#- create a header for the dictionary of dictionaries
myheader=['network']
for i in chosenclasses:
	myheader.append(i)
#- create a list of lists containing the relevant data
mydata=[]
for mi in list(network_dict.keys()):
	myvalues=[mi]
	for j in myheader[1:]:
		myvalues.append(network_dict[mi][j])
	mydata.append(myvalues)
#- create a dataframe from the records and save in the appropriate fashion
mydf=pd.DataFrame.from_records(mydata,columns=myheader)
mydf.to_csv(evaldir+'igp_evaluation_scores.csv',index=False)

#- get average igp values for each patient in a dict
netigpscores={}
for item in mydata:
	#netigpscores[item[0]]=forgivingmsd(item[1:])[0]
	netigpscores[item[0]]=np.mean(item[1:])

#- save average igp values
avg_igp=pd.DataFrame.from_dict(netigpscores,orient='index',columns=['igp_avg_score'])
avg_igp.to_csv(evaldir+'igp_average_scores.csv')
#- sort average igp values 
#sorted_y=sorted(netigpscores.items(),key=operator.itemgetter(1),reverse=True)
#- save best result from IGP evaluation
#bestdict['IGP']=sorted_y[0]
max_value = max(netigpscores.values())  # maximum value
max_keys = [k for k, v in netigpscores.items() if v == max_value] # getting all keys containing the `maximum`
bestdict['IGP']=[max_value,max_keys[0]]


#- create and save final top result from each evaluation technique in a single txt file
best_df=pd.DataFrame.from_dict(bestdict,orient='index')
best_df.to_csv(evaldir+'top_evaluation_compound.txt',sep='\t',header=False)


####################################################################################
# FIRST NEAREST NEIGHBOUR
####################################################################################



onetwork_dict={}
#- for each network in the networks produced by the analysis
for netkey in list(networks.keys()):
	network=networks[netkey]
	classdict={}
	acc=0
	tot=0
	#- for each class in the given classes
	for myclass in chosenclasses:
		thislist=classlists[myclass]
		#- for each patient in this class
		#acc=0
		for patkey in thislist:
			#- identify its nearest neighbour
			for pj in range(len(network[patkey])):
				#- find nearest neighbour that is not the same patient
				if network[patkey].nlargest(10,keep='all').index[pj]!=patkey:
					nearest_neighbour=network[patkey].nlargest(10,keep='all').index[pj]
					break
			#- check if the nearest neighbour is of the same class
			if nodeclasses[nearest_neighbour]==nodeclasses[patkey]:
				acc+=1
		#ig_p=float(acc)/len(thislist)
		tot+=len(thislist)
		#classdict[myclass]=ig_p
	classdict=float(acc)/tot
	onetwork_dict[netkey]=classdict


#- save average igp values
fnn_score=pd.DataFrame.from_dict(onetwork_dict,orient='index',columns=['FNN _score'])
fnn_score.to_csv(evaldir+'fnn_scores.csv')
#- sort average igp values 
sorted_y=sorted(netigpscores.items(),key=itemgetter(1),reverse=True)
#- save best result from IGP evaluation
bestdict['FNN']=sorted_y[0]

#--------------------------------------------------
# Save best snippet
#--------------------------------------------------

#- create and save final top result from each evaluation technique in a single txt file
#best_df=pd.DataFrame.from_dict(bestdict,orient='index')
#! this is implemented to run on lx05 server
best_df=pd.DataFrame.from_dict(bestdict).T
best_df.to_csv(evaldir+'top_evaluation_compound.txt',sep='\t',header=False)


####################################################################################
# COMPARISON
####################################################################################


#- create folder for comparison
outcomedir=outputdir+'outcome/'
if not os.path.isdir(outcomedir):
	os.mkdir(outcomedir)


#- load previous evaluation results
ANN_previous=pd.read_csv(o_evaldir+'/ANN_evaluation_scores.csv',index_col=0)
#NDA_previous=pd.read_csv(o_evaldir+'/NDA_avg_evaluation_scores.csv',index_col=0)
IGP_previous=pd.read_csv(o_evaldir+'/igp_average_scores.csv',index_col=0)
FNN_previous=pd.read_csv(o_evaldir+'/fnn_scores.csv',index_col=0)

#- subset previous scores based on curren preferences
msub_ANN_previous=ANN_previous.loc[models]
#msub_NDA_previous=NDA_previous.loc[models+e_best]
msub_IGP_previous=IGP_previous.loc[models]
msub_FNN_previous=FNN_previous.loc[models]


#- change index names from train to test
to_replace={}
for i in msub_ANN_previous.index:
	if 'train' in i:
		newi=i.replace('train','test')
		to_replace[i]=newi
msub_ANN_previous.rename(index=to_replace,inplace=True)
'''
to_replace={}
for i in msub_NDA_previous.index:
	if 'train' in i:
		newi=i.replace('train','test')
		to_replace[i]=newi
msub_NDA_previous.rename(index=to_replace,inplace=True)
'''
to_replace={}
for i in msub_IGP_previous.index:
	if 'train' in i:
		newi=i.replace('train','test')
		to_replace[i]=newi
msub_IGP_previous.rename(index=to_replace,inplace=True)

to_replace={}
for i in msub_FNN_previous.index:
	if 'train' in i:
		newi=i.replace('train','test')
		to_replace[i]=newi
msub_FNN_previous.rename(index=to_replace,inplace=True)


#- prep test values
ann_dat=df.set_index('Network')

#- merge with current scores to see difference
ANN_accross=pd.merge(msub_ANN_previous,ann_dat,right_index=True,left_index=True,suffixes=('_train','_test'))
#NDA_accross=pd.merge(msub_NDA_previous,nda_file,right_index=True,left_index=True,suffixes=('_train','_test'))
IGP_accross=pd.merge(msub_IGP_previous,avg_igp,right_index=True,left_index=True,suffixes=('_train','_test'))
FNN_accross=pd.merge(msub_FNN_previous,fnn_score,right_index=True,left_index=True,suffixes=('_train','_test'))


#- save the output
ANN_accross.to_csv(outcomedir+'ANN_across.csv')
#NDA_accross.to_csv(outcomedir+'NDA_accross.csv')
IGP_accross.to_csv(outcomedir+'IGP_accross.csv')
FNN_accross.to_csv(outcomedir+'FNN_accross.csv')






