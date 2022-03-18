#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------
import os
import pandas as pd
#--------------------------

inputfile='/mnt/d/input/GDS5826_prepped.expression'
inputdir='/mnt/d/output/GDS5826_27042020_EONSv0.1.7_default_pca_GDS5826_cl/'
outputfile='/mnt/d/Text/images/Sepsis_pilot_genes.csv'

#- set features dir
featdir=inputdir+'features/'

#- Read-in input
allfeatures={}
for filename in os.listdir(featdir):
	s=os.path.abspath(featdir)
	if filename.endswith('.expression'):
		allfeatures[filename.split('.expression')[0]]=pd.read_csv(s+'/'+filename,index_col=0,sep='\t')

#- Initiate empty directory
features=pd.DataFrame()

#- add data in table with component as column header and gene name as content
for key in allfeatures.keys():
	features[key]=allfeatures[key].index


#- read in correct gene names by reading in inputfile
if inputfile.endswith('.csv'):
	inp=pd.read_csv(inputfile)
else:
	inp=pd.read_csv(inputfile,sep='\t')

#- substitute the identifier for the name in the features file (THIS WILL TAKE LONG !!!)
for i in range(len(inp)):
	#features.loc[(features[list(features.keys())[0]]==inp[list(inp.keys())[0]][i]), 'Comp_10_t10' ] = inp[list(inp.keys())[1]][i]
	features=features.replace(inp[list(inp.keys())[0]][i], inp[list(inp.keys())[1]][i])

'''
#- Alternative Method
#- create a dictionary of termA:termB
idict={}
for i in range(len(inp)):
	idict[inp[list(inp.keys())[0]][i]]=inp[list(inp.keys())[1]][i]
#- create dictionary with key:idict to use for replace one-liner
repl={}
for j in range(len(features)):
	repl[list(features.keys())[j]]=idict
#- substitute values with gene names
features=features.replace(repl)
'''

#- change column names of features 
ndict={}
for num in range(1,len(features.keys())+1):
	ndict['Comp_'+str(num)+'_t10']='Component '+str(num)


features = features.rename(ndict, axis='columns')


#- save final version into csv
features.to_csv(outputfile)


#- split into 10 subparts
f1,f2,f3,f4,f5,f6,f7,f8,f9,f10=np.array_split(features,10)
f1.to_csv('/mnt/d/Sepsis_pilot_genes_p1.csv',index_label='N')
f2.to_csv('/mnt/d/Sepsis_pilot_genes_p2.csv',index_label='N')
f3.to_csv('/mnt/d/Sepsis_pilot_genes_p3.csv',index_label='N')
f4.to_csv('/mnt/d/Sepsis_pilot_genes_p4.csv',index_label='N')
f5.to_csv('/mnt/d/Sepsis_pilot_genes_p5.csv',index_label='N')
f6.to_csv('/mnt/d/Sepsis_pilot_genes_p6.csv',index_label='N')
f7.to_csv('/mnt/d/Sepsis_pilot_genes_p7.csv',index_label='N')
f8.to_csv('/mnt/d/Sepsis_pilot_genes_p8.csv',index_label='N')
f9.to_csv('/mnt/d/Sepsis_pilot_genes_p9.csv',index_label='N')
f10.to_csv('/mnt/d/Sepsis_pilot_genes_p10.csv',index_label='N')

