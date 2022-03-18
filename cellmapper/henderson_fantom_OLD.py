#!/opt/local/bin/python
# -*- coding: UTF-8 -*-


# Script comparing the fantom5 dataset with the henderson dataset (Specific data-data comparison)

import pandas as pd 
import scipy as sp
from scipy import stats 

from slacker import Slacker 

slack = Slacker('xoxp-56103715248-56114624145-154905509446-72c754c99fabcd5e9e1772648f3ee7c4')
slack.chat.post_message('@evan','Commencing analysis')




#--------------------------------
# Format the Henderson Dataset
#--------------------------------




#- Import henderson dataset
henderson=pd.read_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/Henderson_full_liver_tissue_cleansed.csv')
slack.chat.post_message('@evan','Henderson dataset imported')

#- Import the file with the cluster annotation
clusters=pd.read_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/Henderson_full_liver_tissue_cleansed_clustering.csv')
clusters_list=clusters.x.unique()

#- subset the dataset by the cluster participation of each cell type
workingdir='/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/output/'
subsetdir='/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/output/subsets/'


#- subset gene expression based on cluster participation
for c in clusters_list:
	clustername='my_cluster%s'%(c)
	#list_of_clust=[]
	mydf=pd.DataFrame({'A':[]})
	for ctype in range(len(clusters['Unnamed: 0'])):
		if clusters.x[ctype]==c:
			mydf[henderson.keys()[clusters['Unnamed: 0'][ctype]]]=henderson[henderson.keys()[clusters['Unnamed: 0'][ctype]]]
			#list_of_clust.append(henderson[henderson.keys()[ctype]])
	mydf.to_csv(subsetdir+clustername)


'''
#- original version of writing the averaged expression for each cluster, does not contain the gene names but has everything in order

clust_dict2={}
for c in clusters_list:
	clust_dict2[c]=pd.DataFrame()
	for ctype in range(len(clusters['Unnamed: 0'])):
		if clusters.x[ctype]==c:
			clust_dict2[c][ctype]=henderson[henderson.keys()[ctype]]

average_dict=pd.DataFrame()
for key in clust_dict2.keys():
	clustername='Cluster%s_average'%(key)
	average_dict[clustername]=clust_dict2[key].mean(axis=1)
average_dict.to_csv(workingdir+'cluster_average_expression.csv')

ncorrHenderson=average_dict.corr(method='spearman')
ncorrHenderson.to_csv(workingdir+'henderson_network.csv')

'''

#- create a file averaging all the expression values of the genes within each clustered sample
#! this is problematic, it produces incorrect results, find out why
sumdf=pd.DataFrame()
sumdf['Gene_names']=henderson[henderson.keys()[0]].values
for c in clusters_list:
	clustername='Cluster%s_average'%(c)
	list_of_clust=[]
	for ctype in range(len(clusters['Unnamed: 0'])):
		if clusters.x[ctype]==c:
			list_of_clust.append(henderson[henderson.keys()[clusters['Unnamed: 0'][ctype]]])
	average=sum(list_of_clust)/len(list_of_clust)
	sumdf[clustername]=average
sumdf.to_csv(workingdir+'cluster_sum_expression.csv')


#- create a correlation map of the henderson average expression clusters
corrHenderson=sumdf.corr(method='spearman')
corrHenderson.to_csv(workingdir+'henderson_cluster_network.csv')


#- read in correlation table from csv
#! this would be the same as using corrHenderson
correlation=pd.read_csv(workingdir+'henderson_network.csv')

#correlation.keys()[0]    #use this for generalised function in order to specify the data of the first column





#- convert matrix files into ncol files to be read by graphia

#! UNDER CONSTRUCTION 

node_list=[]
for fnode in range(len(correlation['Unnamed: 0'])):
	for tnode in correlation.keys()[1:]:
		node_list.append([correlation['Unnamed: 0'][fnode], tnode, correlation[tnode][fnode]])



o=open(workingdir+'henderson_network.ncol','w')
for fnode in range(len(correlation['Unnamed: 0'])):
	for tnode in correlation.keys()[1:]:
		o.write('%s\t%s\t%s\n'%(correlation['Unnamed: 0'][fnode], tnode, correlation[tnode][fnode]))
o.close()



#! make this not same node loops

fcorr=pd.read_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/fantom5_sample_corrmatrix.csv')
o=open(workingdir+'fantom_network.txt','w')
for fnode in range(len(fcorr['Unnamed: 0'])):
	for tnode in fcorr.keys()[1:]:
		o.write('%s\t%s\t%s\n'%(fcorr['Unnamed: 0'][fnode], tnode, fcorr[tnode][fnode]))
o.close()



#--------------------------------
# Format the Fantom5 Dataset
#--------------------------------



#- import the metadata file
metadata=pd.read_table('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/METADATA_SYNC015')

#- import the data file
fantom5=pd.read_table('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/nozeros_phase2_nopooled.expression')

#- find names which are different between the dataset and the metadata file
differentnames=set(metadata['00Annotation']).symmetric_difference(set(fantom5['00Annotation']))



#- finding relative expression of each region




#- selecting gene regions based on peak 1
gene_regions=[]
name_list=[]
gene_reg_association=[]
for i in range(len(metadata)):
	try:
		peak,name=metadata.short_description[i].split('@')
		if peak=='p1':
			gene_regions.append(metadata['00Annotation'][i])
			name_list.append(name)
			gene_reg_association.append((metadata['00Annotation'][i],name))
		else:
			pass
	except:
		hold=metadata.short_description[i].split('@')
		all_items=[]
		for obj in hold:
			try:
				one,two=obj.split(',')
				all_items.append(one)
				all_items.append(two)
			except:
				all_items.append(obj)
		for pos in range(len(all_items)):
			if all_items[pos]=='p1':
				sm_name=all_items[pos+1]
				gene_regions.append(metadata['00Annotation'][i])
				name_list.append(sm_name)
				gene_reg_association.append((metadata['00Annotation'][i],sm_name))
			else:
				pass



#- subset fantom dataset by gene regions identified previously
#! evaluate that this does what it should

#- keep only genes that are have peak1
#fantom_subset=fantom5[fantom5['00Annotation'].isin(gene_regions)]

#- subset henderson dataset based on common genes names with genes which have peak1
#henderson_subset=henderson[henderson['Unnamed: 0'].isin(name_list)]
henderson_subset=henderson_average[henderson_average['Gene_names'].isin(name_list)]

#- find common genes in fantom5 from henderson subset 
#hend_names=henderson_subset['Unnamed: 0']
hend_names=henderson_subset['Gene_names']
common_gene=[]
for i in range(len(gene_reg_association)):
	if gene_reg_association[i][1] in set(hend_names):
		common_gene.append(gene_reg_association[i])


#- keep only gene regions to use for subselection of the fantom5 dataset
common_gene_regions=[]
for i in range(len(common_gene)):
	common_gene_regions.append(common_gene[i][0])



#====================EVAL=====================

#- keep only common genes that are have peak1
fantom_subset=fantom5[fantom5['00Annotation'].isin(common_gene_regions)]


#- keep gene regions which have gene equivalent in henderson 
#- common_gene_association is the dictionary for the gene locations that are common between the two datasets
common_gene_association=[]
for i in range(len(gene_reg_association)):
	for jn in hend_names:
		if gene_reg_association[i][1]==jn:
			common_gene_association.append(gene_reg_association[i])

'''
#- break down list to replace gene names with genomic locations for correlation
position=[]
gname=[]
for item in range(len(common_gene_association)):
	position.append(common_gene_association[item][0])
	gname.append(common_gene_association[item][1])
'''

#- replace henderson gene names for fantom5 genomic regions
for i in range(len(common_gene)):
	#for jl in henderson_subset['Unnamed: 0']:
	for jl in henderson_subset['Gene_names']:
		if jl==common_gene[i][1]:
			#henderson_subset['Unnamed: 0']=henderson_subset['Unnamed: 0'].replace(jl,common_gene[i][0])
			henderson_subset['Gene_names']=henderson_subset['Gene_names'].replace(jl,common_gene[i][0])

'''
#- preffered method of replacing
df.loc[index, "A"] = "I am working! {}".format(row["B"])
.loc[row_indexer,col_indexer] = value
'''


#- create dataframes with proper indexes, sorted the same way for comparison
x=pd.DataFrame()
for i in range(len(common_gene_regions)):
	x=x.append(fantom5[fantom5['00Annotation']==common_gene[i][0]])         # ignore_index=True

#! in final implementation will drop the addition of another variable keeping only x and subset 
testdf=x.sort_values('00Annotation', axis=0)
testdf = testdf.reset_index(drop=True)
#testdf2=henderson_subset.sort_values('Unnamed: 0', axis=0)
testdf3=henderson_subset.sort_values('Gene_names', axis=0)
testdf3 = testdf3.reset_index(drop=True)
testdf3=testdf3.drop('Unnamed: 0',axis=1)


#- correlation matrix of the two datasets
matrix_list=[]
for henderson_key in testdf3.keys():
	slack.chat.post_message('@evan','Commencing correlation computation of Henderson cluster: '+str(henderson_key))
	for fantom_key in testdf.keys():
		#if fantom_key=='00Annotation' or henderson_key=='Unnamed: 0':
		if fantom_key=='00Annotation' or henderson_key=='Gene_names':
			pass
		else:
			node1=fantom_key
			node2=henderson_key
			weight=sp.stats.spearmanr(testdf[fantom_key],testdf3[henderson_key])[0]
			matrix_list.append((node1, node2, weight))



outfile=open(workingdir+'henderson_fantom_spearcorr_matrix_correct.ncol','w')
for i in range(len(matrix_list)):
	outfile.write('%s\t%s\t%s\n'%(matrix_list[i]))
outfile.close()





