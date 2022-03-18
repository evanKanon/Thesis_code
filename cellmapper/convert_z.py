#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

#- - - Changing z coord to smth related to similarity

import os
import string
import numpy as np
import pandas as pd

coordfile='/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/output/just_coord_slim.txt'
resfile='/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/output/nodes_to_add.txt'


resdf=pd.read_table(resfile, header=None)
coorddf=pd.read_table(coordfile, header=None)


mykey=resdf.keys()[max(range(len(resdf.keys())))]


#- create data dictionary with all edge weights for each cluster
dat_dict={}
for clustkey in coorddf[0]:
	clustdat=[]
	for i in range(len(resdf[0])):
		if clustkey==resdf[0][i]:
			clustdat.append(float(resdf[mykey][i]))
	dat_dict[clustkey]=clustdat


#- create single list with all values for min computation
flat_list = [x for key in dat_dict.keys() for x in dat_dict[key]]
#- select minimum
mymin=min(flat_list)

#- calculate z coordinate for each node
for key in dat_dict:
	#- calculate the average of the 3 edge weights
	avg=np.mean(dat_dict[key])
	#- calculate the factor of max(-log10)-log10, with max(-log10) being the log10 of the min value
	factor=-np.log10(mymin)-np.log10(avg)
	#- calculate the z coord by multiplying the 'distance factor' with the 50 pixel distance
	zcoord=factor*50
	#- append new z score to the relevant column in the coord file
	for ind in range(len(coorddf)):
		if coorddf[0][ind]==key:
			coorddf[3][ind]=zcoord


#- save altered result
coorddf.to_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/output/new_xyz.txt',header=None,index=None,sep='\t')













