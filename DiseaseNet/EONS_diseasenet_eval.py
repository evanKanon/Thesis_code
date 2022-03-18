
# script which measured the EONS disease network





import os
import pandas as pd


def write_ncol(myfile):
	'''Function which converts matrix-type file into ncol'''
	#import pandas as pd
	#import os
	data=pd.read_csv(myfile,index_col=0)
	myoutputfile=os.path.split(myfile)[0]+'/'+os.path.split(myfile)[1].split('.')[0]
	newfile=open(myoutputfile+'.ncol','w')
	for pat1 in data.keys():
		for pat2 in data.index:
			newfile.write('%s\t%s\t%s\n'%(pat1, pat2, data[pat1][pat2]))
	newfile.close()


#write_ncol('/mnt/d/PhD/Results/Diseasenet_EONS/Diseasenetcomulative_GWAS_prepped.expression_spearman-9900.csv')

#! NOTE due to the current for of write_ncol, the output name is impractical so it must be run separately for each file
write_ncol('/mnt/d/PhD/Results/Diseasenet_EONS/Diseasenetcomulative_GWAS_prepped.expression_kendall-9900.csv')


ddf=pd.read_csv('/mnt/d/PhD/Results/Diseasenet_EONS/Diseasenetcomulative_GWAS_prepped.ncol',sep='\t',header=None)


# delete values of items when they are for the same 
for index, row in ddf.iterrows():
	if row[0] == row[1]:
		ddf.drop(index, inplace=True)

ddf['ranges']=pd.cut(ddf[2],bins=[-2,-1,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2])


ddf.groupby(['ranges']).agg(['mean','count'])
