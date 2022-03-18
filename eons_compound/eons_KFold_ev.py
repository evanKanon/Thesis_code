#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------
import scipy
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from itertools import permutations 
#--------------------------

def dendrogrammer(df, leaf_labels,outputname='dendrogram.png'):
	''' Function which plots and saves agglomerative clustering dendrogram plot '''
	#import numpy as np
	#from scipy.cluster.hierarchy import dendrogram, linkage
	#from matplotlib import pyplot as plt
	#- keep the values of the input
	D=df.values
	#- transpose data if we to cluster the other way
	if len(leaf_labels) != len(D):
		D = np.transpose(D)
	#- perform clustering
	#Linkage methods could be ‘single’, ‘average’, ‘complete’, ‘median’, ‘centroid’, ‘weighted’, or ‘ward’
	#There are many possible distance metrics (e.g., ‘cityblockk’, ‘yule’, ‘hamming’, ‘dice’, ‘kulsinski’, ‘correlation’, ‘jaccard’, and many more), or you can create your own. See the scipy documentation for pdist for more info.
	Z = linkage(D, method='ward', metric='euclidean')
	#- plot figure
	plt.figure(figsize=(10, 6))
	ax = plt.subplot()
	plt.subplots_adjust(left=0.07, bottom=0.3, right=0.98, top=0.95, wspace=0, hspace=0)
	plt.xlabel('Samples')
	plt.ylabel('Distance')
	dendrogram(Z, leaf_rotation=90., leaf_font_size=10., labels=leaf_labels)
	plt.savefig(outputname)

'''
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h
'''

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return h

#--------------------------
#findict={'KMeans':[],'Agglomerative':[],'Random':[],'EONS':[],'LASSO':[0.875,0.950,0.900,0.950,0.900]}
#findict={'KMeans':[],'Agglomerative':[],'Random':[],'EONS':[],'LASSO':[0.600,0.425,0.550,0.500,0.475]}
findict={'KMeans':[],'Agglomerative':[],'Random':[]}
eonsdict={}
mytitle='Disease State annotation'
#mytitle='Specimen annotation'
foldnumber=5
#--------------------------
#eons=pd.read_csv('/mnt/d/output_old/KFold/eons_test_EONSv0.1.4_default_train_5_disease/test_evaluation/igp_average_scores.csv')
#imp = pd.read_csv('/mnt/d/input/KFold/GSE19429_test_5.expression',sep='\t',index_col=0)
#classes=pd.read_csv('/mnt/d/input/KFold/GSE19429_test_5_disease.classes',sep='\t')
#eons=pd.read_csv('/mnt/d/output_old/KFold/eons_test_EONSv0.1.4_default_train_1_specimen/test_evaluation/igp_average_scores.csv')
#imp = pd.read_csv('/mnt/d/input/KFold/GSE19429_test_1.expression',sep='\t',index_col=0)
#classes=pd.read_csv('/mnt/d/input/KFold/GSE19429_test_1_specimen.classes',sep='\t')
#--------------------------
#las=pd.read_csv('/mnt/d/PhD/Results/LASSO_specimen.csv',index_col=0)
las=pd.read_csv('/mnt/d/PhD/Results/LASSO_disease.csv',index_col=0)
findict['LASSO']=list(las['LASSO_testing'])
for oi in range(1,foldnumber+1):
	#eons=pd.read_csv('/mnt/d/PhD/Results/KFold/eons_test_EONSv0.1.7_default_pca_train_'+str(oi)+'_specimen/test_evaluation/fnn_scores.csv')
	eons=pd.read_csv('/mnt/d/PhD/Results/KFold/eons_test_EONSv0.1.7_default_pca_train_'+str(oi)+'_disease/test_evaluation/fnn_scores.csv')
	imp = pd.read_csv('/mnt/d/input/KFold/GSE19429_test_'+str(oi)+'.expression',sep='\t',index_col=0)
	#classes=pd.read_csv('/mnt/d/input/KFold/GSE19429_test_'+str(oi)+'_specimen.classes',sep='\t')
	classes=pd.read_csv('/mnt/d/input/KFold/GSE19429_test_'+str(oi)+'_disease.classes',sep='\t')
	#- read input
	#df1=imp.drop(columns=['IDENTIFIER'])
	#df1=df1.T
	df1=imp.T
	cells = list(imp.columns.values) 
	#- Run K-Means clustering for the same number of clusters as there are classes
	kmeans = KMeans(n_clusters=len(np.unique(classes['Class']))).fit(df1)
	centroids = kmeans.cluster_centers_
	#print(centroids)
	#- plot clusters
	#plt.scatter(df1, c= kmeans.labels_.astype(float), s=50, alpha=0.5)
	#plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=50)
	#plt.show()
	#- Create the dictionary that defines the order for sorting (based on what was used for K-Maens/clustering)
	sorterIndex = dict(zip(cells,range(len(cells))))
	# Generate a rank column that will be used to sort
	# the dataframe numerically
	#classes['id_Rank'] = classes['id'].map(sorterIndex)
	classes['id_Rank'] = classes['#samples'].map(sorterIndex)
	# sort values based on sorter column
	classes.sort_values(by=['id_Rank'], inplace=True)
	# reset index
	classes.reset_index(drop=True, inplace=True)
	# drop sorter column
	classes.drop('id_Rank', 1, inplace = True)
	# make categories numerical
	#classes['cluster'] = pd.factorize(classes['soft_classifiers'])[0]
	classes['cluster'] = pd.factorize(classes['Class'])[0]
	# set labels
	labels_pred=kmeans.labels_
	labels_true=pd.factorize(classes['Class'])[0]
	#- Evaluation of clustering
	# Get all permutations of [1, 2, 3] 
	#perm = list(permutations([0, 1, 2, 3, 4]))
	perm = list(permutations(list(np.unique(labels_true))))
	# Save final result
	feval={}
	# Print the obtained permutations 
	for item in perm:
		# initiate empty dict for converter
		mdict={}
		# convert 
		for counter, value in enumerate(item):
			mdict[str(counter)]=value
		# make the substitution
		new=[]
		for sampl in labels_pred:
			new.append(mdict[str(sampl)])
		# initiate the main dictionary
		df=pd.DataFrame()
		# add the different values
		df['true']=labels_true
		df['old_pred']=labels_pred
		df['new_pred']=new
		# evaluate prediction
		a=0
		for li in range(len(df)):
			if df['true'][li]==df['new_pred'][li]:
				a+=1
		# save final result to dictionary
		feval[''.join(str(x) for x in item)]=float(a)/len(df)
	# find key with maximum value
	Keymax = max(feval, key=feval.get)
	feval[Keymax]
	#! improve this point to make list-of lists
	# save final top result into ilst
	fin=['KMeans',feval[Keymax]]
	findict['KMeans'].append(feval[Keymax])
	#--------------------------
	#- Agglomerative Clustering
	#--------------------------
	from sklearn.cluster import AgglomerativeClustering
	#clustering = AgglomerativeClustering(n_clusters=4).fit_predict(df1)
	clustering = AgglomerativeClustering(n_clusters=len(np.unique(classes['Class']))).fit(df1)
	labels_true=pd.factorize(classes['Class'])[0]
	labels_pred=clustering.labels_
	#- Evaluation of clustering
	# Get all permutations of [1, 2, 3] 
	#perm = list(permutations([0, 1, 2, 3]))
	# Save final result
	feval={}
	# Print the obtained permutations 
	for item in perm:
		# initiate empty dict for converter
		mdict={}
		# convert 
		for counter, value in enumerate(item):
			mdict[str(counter)]=value
		# make the substitution
		new=[]
		for sampl in labels_pred:
			new.append(mdict[str(sampl)])
		# initiate the main dictionary
		df=pd.DataFrame()
		# add the different values
		df['true']=labels_true
		df['old_pred']=labels_pred
		df['new_pred']=new
		# evaluate prediction
		a=0
		for li in range(len(df)):
			if df['true'][li]==df['new_pred'][li]:
				a+=1
		# save final result to dictionary
		feval[''.join(str(x) for x in item)]=float(a)/len(df)
	# find key with maximum value
	Keymax = max(feval, key=feval.get)
	feval[Keymax]
	# append final result to final list
	#fin.append(['Agglomerative Clustering',feval[Keymax]])
	fin.append(['Agglomerative',feval[Keymax]])
	findict['Agglomerative'].append(feval[Keymax])
	#--------------------------
	#- RANDOM CLASS ALLOCATION
	#--------------------------
	#- set seed for reproducibility
	np.random.seed(42)
	#- randomly allocate classes
	classes['Random'] = np.random.randint(0, len(np.unique(labels_true)), classes.shape[0])
	#- evaluate which ones of them are correct
	s=0
	for r in range(len(classes)):
		if classes['Random'][r]==classes['cluster'][r]:
			s+=1
	#- aooend final percentage result to final list
	fin.append(['Random',float(s)/len(classes)])
	#- append final percentage result to final dictionary
	findict['Random'].append(float(s)/len(classes))
	#--------------------------
	#- EONS SCORE
	#--------------------------
	#eons=pd.read_csv('/mnt/d/output_old/KFold/eons_test_EONSv0.1.4_default_train_1_disease/test_evaluation/igp_average_scores.csv')
	#fin.append(['EONS',max(eons['igp_avg_score'])])
	#fin.append(['EONS',np.mean(eons['igp_avg_score'])])
	#findict['EONS'].append(max(eons['igp_avg_score']))
	#findict['EONS'].append(np.average(eons['igp_avg_score']))
	#- name eons key based on Fold number
	mystrit='eons_'+str(oi)
	#- get all values of the eons analysis for the fold
	eonsdict[mystrit]=eons['FNN _score'].to_list()



#--------------------------
#- PLOT
#--------------------------
mlot=[]
for key in findict.keys():
	avg=np.mean(findict[key])
	ci=mean_confidence_interval(findict[key])
	lil=[key,avg,ci]
	mlot.append(tuple(lil))

#- find average and 95% CI for each of the different models in each eons fold
for nkey in eonsdict.keys():
	navg=np.mean(eonsdict[nkey])
	nci=mean_confidence_interval(eonsdict[nkey])
	nlil=['EONS',navg,nci]
	mlot.append(tuple(nlil))

#- create the Dataframe which will be used to plot
pdf = pd.DataFrame.from_records(mlot, columns=("Method", "Correct prediction value", 'Confidence'))

'''
sns.catplot(x="Method", y="Value", data=pdf)
#plt.title('Disease State')
#plt.savefig('/mnt/d/res_disease.png')
plt.title('Specimen')
plt.savefig('/mnt/d/res_specimen.png')
'''

'''
def errplot(x, y, yerr, **kwargs):
    ax = plt.gca()
    data = kwargs.pop("data")
    data.plot(x=x, y=y, yerr=yerr, ax=ax, **kwargs)


g = sns.FacetGrid(pdf)
g.map_dataframe(errplot, "Method", "Value",0.5)


data.plot(x=x, y=y, yerr=yerr, kind="bar", ax=ax, **kwargs)
'''

#- set name for output file
name='_'.join(mytitle.split(' '))
#- plot and save image
sns.catplot('Method','Correct prediction value',data=pdf,jitter=False)
plt.errorbar(pdf['Method'], pdf['Correct prediction value'], yerr=pdf['Confidence'],ecolor='lightgray', elinewidth=3, capsize=0,ls='none')
plt.title(mytitle)
plt.savefig('/mnt/d/PhD/Text/images/'+name+'.png',bbox_inches='tight')
