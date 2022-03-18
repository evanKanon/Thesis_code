#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

import os
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split

#------------------------------

#- save output
outpuptdir='/mnt/d/PhD/input/KFold/'


# Load the cancer dataset
df = pd.read_csv('/mnt/d/PhD/input/GSE19429_series_matrix_prepped.expression',sep='\t',index_col=0)
y_df=pd.read_csv('/mnt/d/PhD/input/GSE19429.classes',sep='\t',index_col=0)


# select current outcome variable
my_y1=y_df[y_df.keys()[0]]
my_y2=y_df[y_df.keys()[1]]

# ensure x and y have same index 
df=df.T
xdf=df.sort_index()
ydf1=my_y1.sort_index()
ydf2=my_y2.sort_index()

# create training and testing vars
kf = KFold(n_splits=5,random_state=21,shuffle=True)
kf.get_n_splits(xdf)

#- save different 
a=0
for train_index, test_index in kf.split(xdf):
	a+=1
	xdf.loc[xdf.index[train_index]].T.to_csv(outpuptdir+'GSE19429_train_'+str(a)+'.expression',sep='\t') #X_train=
	xdf.loc[xdf.index[test_index]].T.to_csv(outpuptdir+'GSE19429_test_'+str(a)+'.expression',sep='\t') #X_test=
	ydf1.loc[ydf1.index[train_index]].to_csv(outpuptdir+'GSE19429_train_'+str(a)+'_disease'+'.classes',sep='\t',header=['Class']) #Y1_train=
	ydf1.loc[ydf1.index[test_index]].to_csv(outpuptdir+'GSE19429_test_'+str(a)+'_disease'+'.classes',sep='\t',header=['Class']) #Y1_test=
	ydf2.loc[ydf2.index[train_index]].to_csv(outpuptdir+'GSE19429_train_'+str(a)+'_specimen'+'.classes',sep='\t',header=['Class']) #Y2_train=
	ydf2.loc[ydf2.index[test_index]].to_csv(outpuptdir+'GSE19429_test_'+str(a)+'_specimen'+'.classes',sep='\t',header=['Class']) #Y2_test=





