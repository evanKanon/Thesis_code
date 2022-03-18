#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# create k_nearest neighbour network from fantom5 correlation network

import pandas as pd
import os
from slacker import Slacker 


slack = Slacker('xoxp-56103715248-56114624145-154905509446-72c754c99fabcd5e9e1772648f3ee7c4')



def k_nearest_neighbour(mydata,k,outname):
	''' Function which creates a ncol file of the k-nearest neighbours of a dataset '''
	#! improve functionnality to include filename and, optionally, raw data input
	#- initiate a list to keep the edges in
	edgelist=[]
	for fnode in mydata.keys():
		#- largest value present in the column for the current key
		maxval=mydata[fnode].sort_values()[max(range(len(mydata[fnode])))]
		#- see if the largest value is with the same node (i.e. if the largest correlation is 1 due to self correlation)
		#! note, this does not accomodate the possibility of having correlations of 1 with nodes other than the same node at the lowest position of the sort
		if mydata[fnode].loc[mydata[fnode]==maxval].index.values[0]==fnode:
			for i in range(k):
				weight=mydata[fnode].sort_values()[max(range(len(mydata[fnode])))-i-1]
				tnode=mydata[fnode].loc[mydata[fnode]==weight].index.values[0]
				edgelist.append((fnode, tnode, weight))
		else:
			for i in range(k):
				weight=mydata[fnode].sort_values()[max(range(len(mydata[fnode])))-i]
				tnode=mydata[fnode].loc[mydata[fnode]==weight].index.values[0]
				edgelist.append((fnode, tnode, weight))
	#- write edgelist into a file
	o=open(outname,'w')
	for edge in edgelist:
		o.write("%s\t%s\t%s\n"%edge)
	o.close()
	#- if we don't want to write a file, we can just return the edge list
	#return edgelist
		

#! could save data subset, delete self looping instance if any and then proceed to find neighbours


#- import data
myfile='/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/fantom5_sample_corrmatrix.csv'
corr_dat=pd.read_csv(myfile,index_col=0)

#- set output name
outname=os.path.abspath(myfile).split('.')[0]+'_3NN.ncol'

#- find 3 nearest neihbour network
k_nearest_neighbour(corr_dat,3,outname)

#- notify me that the creation has been completed
slack.chat.post_message('@evan','Fantom5 k-nearest neighbour matrix creation complete')

