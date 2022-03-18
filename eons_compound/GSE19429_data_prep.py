#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

import os
import pandas as pd
from sklearn.model_selection import train_test_split

#------------------------------

def make_expression_file(file,path=False):
	'''Function which converts .soft or series_matrix.txt files to .expression'''
	input_file=open(file)
	out=[]
	for line in input_file:
		if not line.startswith('#') and not line.startswith('!') and not line.startswith('^'):
			out.append(line)
	input_file.close()
	new_file_path=os.path.splitext(file)[0]+'_prepped.expression'
	output_file=open(new_file_path, 'w')
	output_file.writelines(out)
	output_file.close()
	if path:
		return new_file_path

#------------------------------


# Load the cancer dataset
make_expression_file('/mnt/d/PhD/input/GSE19429_series_matrix.txt')
df = pd.read_csv('/mnt/d/PhD/input/GSE19429_series_matrix_prepped.expression',sep='\t',index_col=0)
y_df=pd.read_csv('/mnt/d/PhD/input/GSE19429.classes',sep='\t',index_col=0)

# select current outcome variable
my_y=y_df[y_df.keys()[0]]

# ensure x and y have same index 
df=df.T
xdf=df.sort_index()
ydf=my_y.sort_index()

# create training and testing vars
X_train, X_test, y_train, y_test = train_test_split(xdf, ydf, test_size=0.2,random_state=15)
print(X_train.shape, y_train.shape)
print(X_test.shape, y_test.shape)

# save train/test split
X_train.T.to_csv('/mnt/d/PhD/input/GSE19429_train.expression',sep='\t')
X_test.T.to_csv('/mnt/d/PhD/input/GSE19429_test.expression',sep='\t')
y_train.to_csv('/mnt/d/PhD/input/GSE19429_train_disease.classes',sep='\t',header=['Class'])
y_test.to_csv('/mnt/d/PhD/input/GSE19429_test_disease.classes',sep='\t',header=['Class'])
# save second class file by split
y_df[y_df.keys()[1]][y_train.index].to_csv('/mnt/d/PhD/input/GSE19429_train_specimen.classes',sep='\t',header=['Class'])
y_df[y_df.keys()[1]][y_test.index].to_csv('/mnt/d/PhD/input/GSE19429_test_specimen.classes',sep='\t',header=['Class'])

