#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

#EONS zeroth script
#Define all functions required for the analysis

#------------------------------
#v0.0.0 initial version of the new code, adapted and augmented from EONS_funct version 0.2
#v0.0.1 updated to keep same version as other files (will change version descrpt convention in current version)
#v0.0.2 added config file, changed EONS_02 defaults, redefined outputdirs, made single EONS_funct and deleted unused functions (version descr now includes changes in any/all scripts)
#v0.0.3 added extreme network selection functionality and the the relevant option for the distance multiplier selection, removed config hardpath and renamed EONS_main to EONS
#v0.0.4 changed MCL feature selection with PCA based on networktools code
#v0.0.5 changed feature selection to only principal components, deleted unsued code
#v0.0.6 merged MCL feature selection and PCA, created multiple feature subsets through different inflation values, made it possible to notify many users simulataneously, added option to run on normalised data as well, automated static deletion of IDENTIFIER column of data
#v0.0.7 created prep script, changed file names from capital to lowercase, changed networks so that all networks are similarity and not distance networks
#v0.0.8 gave two options for feature selection, either MCL or PCA-based feature selection, added classfile option and config path as well as working and useful eval script
#v0.1.0 BASELINE EONS added failed network functionality to eons_02 making all distances runable
#v0.1.1 changed make_expression function to accomodate file names other than .soft, added NA dealing option to eons_prep, improved config option for eons_main, made code python 3 compatible
#v0.1.2 established perpendicular point plots
#v0.1.3 improved folder names, added .layout visualisation, updated environment file
#v0.1.4 made running the MCL component faster and allowed user to select inflation values
#v0.1.5 changed naming of folders and certain default values
#v0.1.6 corrected PCA random selection in eons_01, corrected network names in eons_02 and deleted unused inputfile option in eons_03 with the appropriate corrections in eons_funct and _main
#v0.1.7 added correct evaluation measure

#-----------------------------
# Definition of version
#-----------------------------

__version__ = '0.1.7'

#-----------------------------
# Import relevant packages
#-----------------------------

import string, sys, os, math, random, copy
import time
import json
import igraph
import datetime
import numpy as np
import pandas as pd
from bisect import *
from igraph import *
import networkx as nx
from textwrap import wrap
from itertools import count
from sklearn import preprocessing
from scipy.spatial import distance
from networkx.readwrite import json_graph
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
#from sklearn.preprocessing import Imputer


from sklearn.impute import SimpleImputer

from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer


import matplotlib
#- this has to be done for creation of the network plots at the end when run on the server
#! learn why
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

#-----------------------------
# Setting of paths for MCL
#-----------------------------

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),'config.json')) as json_data_file:
	config=json.load(json_data_file)

pathtomcl = config['File_paths']['pathtomcl']
pathtomcx = config['File_paths']['pathtomcx']
pathtomcxload = config['File_paths']['pathtomcxload']

#-----------------------------
redo_all = True # ANY PROBLEMS, TRY THIS FIRST!
threads = 6
num_perms = 100
sourcedir = "sourcefiles" #! elaborate on this idea or delete


#-----------------------------
# Definition of functions
#-----------------------------

def ncheck(cval, list1, infl=1):
	return(all(x < cval*infl for x in list1))

def pandas_correlation_function(mydata,myoutputfile,transpose=False,correlation='spearman',threshold=-99,er=False):
	''' 'pearson', 'kendall', 'spearman' '''
	if transpose==True:
		corr_dat=mydata.T.corr(correlation)
	else:
		corr_dat=mydata.corr(correlation)
	keylist=list(corr_dat.keys())
	newfile=open(myoutputfile+'.ncol','w')
	for pat1 in range(len(keylist)):
		for pat2 in range(len(keylist)):
			if corr_dat[keylist[pat1]][keylist[pat2]]>=threshold:
				newfile.write('%s\t%s\t%s\n'%(keylist[pat1], keylist[pat2], corr_dat[keylist[pat1]][keylist[pat2]]))
	newfile.close()
	if er==True:
		max_list=list(corr_dat.max())
		min_list=list(corr_dat.min())
		er=(min(min_list),max(max_list))
		return er

#! evaluate proper working of the function
def q_pandas_correlation_function(mydata,myoutputfile,transpose=False,correlation='spearman',threshold=-99,er=False):
	''' 'pearson', 'kendall', 'spearman' '''
	if transpose==True:
		mydata=mydata.T
	else:
		pass
	#- kenny's suggestion to make the code run quicker
	if correlation=='spearman':
		corr_dat=mydata.rank(ascending=True).corr('pearson')
	else:
		corr_dat=mydata.corr(correlation)
	#- set all values of the dataframe to zero if they are smaller than or equal to the threshold
	corr_dat[corr_dat<=threshold]=0
	#- save dataframe as csv file
	corr_dat.to_csv(myoutputfile+'.csv')
	if er==True:
		max_list=list(corr_dat.max())
		min_list=list(corr_dat.min())
		er=(min(min_list),max(max_list))
		return er

def correlation_function(mydata,myoutputfile,transpose=False,correlation='spearman',threshold=-99,er=False,save=True):
	''' Function which correlates different distances between columns of input, scale values to range -1,1 and save matrix output '''
	if transpose==True:
		mydata=mydata.T
		name_add='_T'
	else:
		name_add=''
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
	if save==True:
		#- save dataframe as csv file
		corr_dat.to_csv(myoutputfile+name_add+'.csv')
	else:
		#- return dataframe
		return corr_dat
	if er==True:
		max_list=list(corr_dat.max())
		min_list=list(corr_dat.min())
		er=(min(min_list),max(max_list))
		return er

def scalariser(mydataframe, minim=-1.0001, maxim=1.0001):
	''' Function which scales all values of a dataframe to given range '''
	minimum, maximum = (min(mydataframe.min()), max(mydataframe.max()))
	#- see if values are within expected range and if not scale within range (range is applied this way due to rounding errors)
	if minimum<minim or maximum>maxim:
		values = mydataframe.values
		scaler = MinMaxScaler(feature_range=(minim,maxim))
		scaler.fit(values)
		#print(scaler.data_max_)
		#print(scaler.transform(data))
		scaled_data=scaler.transform(values)
		newframe=pd.DataFrame(scaled_data,columns=mydataframe.columns,index=mydataframe.index)
		return newframe
	else:
		return mydataframe
	#- save dataframe as csv file
	#mydataframe.to_csv(myoutputfile+'.csv')
	
def ppd_plot(mydf,highlight=['Nothing'],outputname='~/my_ppd.png'):
	''' Function which creates a perpendicular point distribution plot '''
	#- create coordinates dictionary for plotting
	coord_dict={'x':[],'y':[],'c':[]}
	#- fill the dictionary in with the relevant points
	x_values=[]
	x_val=3
	for c in mydf.keys():
		x_values.append(x_val)
		for i in range(len(mydf)):
			coord_dict['x'].append(x_val)
			coord_dict['y'].append(mydf[c][i])
			if mydf.index[i] in highlight:
				coord_dict['c'].append(1)
			else:
				coord_dict['c'].append(0)
		x_val+=1
	#- set proper names tuple for plot x-axis
	mycolname=list(mydf.keys())
	mycols=['']+mycolname+['']
	mycols=tuple(mycols)
	#- turn the coordinate dictionary into pandas dataframe from plotting
	plotdf=pd.DataFrame(coord_dict)
	#- then plot using pandas plotting
	#plotdf.plot(kind='scatter',x='x',y='y',color='c') #style='bx' will give x type marks instead of complete circles, 'r+' will give red plusses; label='point'
	#plotdf.plot(x='x',y='y',style='bo') #style='bx' will give x type marks instead of complete circles, 'r+' will give red plusses; label='point'
	#plotdf.plot(x='x',y='y',style='bo') #style='bx' will give x type marks instead of complete circles, 'r+' will give red plusses; label='point'
	#plotdf.plot(x='x',y='y',style='ro') #style='bx' will give x type marks instead of complete circles, 'r+' will give red plusses; label='point'
	#- make scatter plot with different colours for interesting points
	#categories=np.array(cats)
	#colormap=np.array(['b','r'])
	#plt.scatter(a[0], a[1], s=100, c=colormap[categories])
	#- plot through scatter 
	#plt.scatter(plotdf['x'],plotdf['y'],c=colormap[plotdf['c']])
	#- plot irrelevant first as blue and then relevant ones as red
	plt.scatter(plotdf['x'][plotdf['c']==0],plotdf['y'][plotdf['c']==0],c='b')
	plt.scatter(plotdf['x'][plotdf['c']==1],plotdf['y'][plotdf['c']==1],c='r')
	#- write correct x-axis values
	x_values.insert(0,x_values[0]-1)
	x_values.append(x_values[-1]+1)
	plt.xticks(tuple(x_values), mycols)
	#- write x-axis title
	plt.xlabel('Evaluator type')
	#- save the image
	plt.savefig(outputname)

def hamming_similarity(mydataframe1,mydataframe2):
	''' Function which calculates hamming similarity for dataframe matrix data '''
	#- add the two dataframes together
	sumdf=mydataframe1.add(mydataframe2,fill_value=0)
	#- calculate the number of 0s present in dataframe1
	x=sum((mydataframe1==0).sum(axis=1))
	#- calculate the number of 0s present in dataframe2
	y=sum((mydataframe2==0).sum(axis=1))
	#- calculate the number of 0s present in the sum dataframe
	z=sum((sumdf==0).sum(axis=1))
	#- calulate the number of edges that are different between the two dataframes (i.e. hamming distance)
	hd=(x+y-2*z)/2
	#- calculate maximum number of possible edges
	dim=len(sumdf)
	maxedg=((dim*dim)-dim)/2
	similarity=1-float(hd)/maxedg
	return similarity

def SEWD_similarity(mydataframe1,mydataframe2):
	''' Function which calculates Sum Edge Weight Distance similarity for dataframe matrix data '''
	#- subtract the two dataframes from each other
	diffdf=mydataframe1.sub(mydataframe2,fill_value=0)
	#- transfor each value to its absolute value
	diffdf=diffdf.transform(abs)
	#- sum them to calculate SEWD
	#- consider NAs as equal to zero
	sewd=diffdf.sum().sum()/2
	#sewd=diffdf.values.sum()/2
	#- calculate maximum possible difference
	dim=len(diffdf)
	maxdist=((dim*dim)-dim)/2    #- /2 is ommited because edge weights can take all values between -1 and 1, spanning 2 
	#- calculate similarity
	similarity=1-(sewd/maxdist)
	return similarity

def normalise(mydataframe):
	''' Function which normalises values of a pandas dataframe and returns a new, normalised dataframe '''
	x = mydataframe.values
	normalised_x = preprocessing.normalize(x)
	return pd.DataFrame(data=normalised_x, index=mydataframe.index, columns=mydataframe.columns)

def scale(mydataframe):
	''' Function which applies z-score scaling to dataframe data '''
	scaler = preprocessing.StandardScaler()
	scaled_df = scaler.fit_transform(mydataframe)
	scaled_df = pd.DataFrame(scaled_df, index=mydataframe.index, columns=mydataframe.columns)
	return scaled_df

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
	#return zscore(total_p, perm_totals)
	perm_totals.append(total_p)
	return zscore(total_p, perm_totals)

def run_mcxload(thisfile,redo=False):
	if thisfile.endswith('.ncol'):
		matrixfile=thisfile.replace('.ncol','.mcx_matrix')
		tabfile=thisfile.replace('.ncol','.mcx_tab')
	else:
		matrixfile=thisfile+'.mcx_matrix'
		tabfile=thisfile+'.mcx_tab'
	if not os.path.exists(matrixfile) or not os.path.exists(tabfile) or redo:
		cmd = pathtomcxload+"mcxload  -abc %s -o %s -write-tab %s"%(thisfile, matrixfile, tabfile)
		print(cmd)
		os.system(cmd)
	else:
		print("Matrix already exists. Using",matrixfile)
		print("and",tabfile)
	return matrixfile, tabfile

def read_dict(myfile):
	'''function which reads in two column files as dictionaries with first column as key and second as item'''
	inputfile=open(myfile)
	d={}
	for line in inputfile:
		key,item=line.strip().split('\t') # this substituted the =line.split() because it would also split spaces, while we only want tabs split
		d[key]=item
	return d

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

def pandas_to_json(mydataframe,mygroupdict='Nan',myoutname='test.json'):
	''' Function which converts pandas dataframe network to json format for d3 visualisation '''
	#- create overarching json dictionary
	json_dict={}
	json_dict["nodes"]=[]
	json_dict["links"]=[]
	#- set features dictionary
	if mygroupdict!='Nan':
		grdct=mygroupdict
	else:
		grdct={}
		for item in list(mydataframe.keys()):
			grdct[item]='none'
	#- run through the pandas dataframe to select and assort appropriately the various values
	for i in range(len(mydataframe.keys())):
		nodedict={}
		netw1=list(mydataframe.keys())[i]
		nodedict['id']=netw1
		nodedict['group']=grdct[netw1]
		json_dict["nodes"].append(nodedict)
		#- run over only the half table, not the entire thing
		x=max(range(len(mydataframe.keys())))
		while x!=i: #(i+1 possibly?)
			linkdict={}
			netw2=list(mydataframe.keys())[x]
			linkdict["source"]=list(mydataframe.keys())[i]
			linkdict["target"]=list(mydataframe.keys())[x]
			linkdict["value"]=mydataframe[netw1][netw2]
			json_dict["links"].append(linkdict)
			x+=(-1)
	#- save dictionary as json file
	with open(myoutname,'w') as fp:
		json.dump(json_dict,fp)

def csv_to_json(myfilename,mygroupdict='Nan'):
	''' Function which converts csv file network to json format for d3 visualisation '''
	mydataframe=pd.read_csv(myfilename,index_col=0)
	#- create overarching json dictionary
	json_dict={}
	json_dict["nodes"]=[]
	json_dict["links"]=[]
	#- set features dictionary
	if mygroupdict!='Nan':
		grdct=pd.read_csv(mygroupdict)
	else:
		grdct={}
		for item in list(mydataframe.keys()):
			grdct[item]='none'
	#- run through the pandas dataframe to select and assort appropriately the various values
	for i in range(len(mydataframe.keys())):
		nodedict={}
		netw1=list(mydataframe.keys())[i]
		nodedict['id']=netw1
		nodedict['group']=grdct[netw1]
		json_dict["nodes"].append(nodedict)
		#- run over only the half table, not the entire thing
		x=max(range(len(mydataframe.keys())))
		while x!=i: #(i+1 possibly?)
			linkdict={}
			netw2=list(mydataframe.keys())[x]
			linkdict["source"]=list(mydataframe.keys())[i]
			linkdict["target"]=list(mydataframe.keys())[x]
			linkdict["value"]=mydataframe[netw1][netw2]
			json_dict["links"].append(linkdict)
			x+=(-1)
	#- set output file name
	attrs = myfilename.split('.')[:-1]
	myoutname = '.'.join(attrs)+'.json'
	#- save dictionary as json file
	with open(myoutname,'w') as fp:
		json.dump(json_dict,fp)

#! can be upgraded to use multiple classes
def pandas_to_layout(mydataframe,mygroupdf='Nan',myoutputname='default'):
	'''Function which saves pandas df  to .layout'''
	data=mydataframe
	#- set output name or do inplace
	if myoutputname=='default':
		attrs = myfilename.split('.')[:-1]
		myoutname = '.'.join(attrs)+'.layout'
	else:
		myoutname=myoutputname
	#- set attributes file
	if mygroupdf!='Nan':
		if mygroupdf.endswith('.csv'):
			grdct=pd.read_csv(mygroupdf,index_col=0)                 #! index_col might cause problems, harmonise formats!
		else:
			grdct=pd.read_csv(mygroupdf,sep='\t',index_col=0)
	else:
		grdct={}
		for item in list(data.keys()):
			grdct[item]='none'
	#- write s,t,w format of matrix dataframe
	row_list=[]
	for pat1 in data.keys():
		for pat2 in data.index:
			if pat1!=pat2:
				row_list.append([pat1,pat2,data[pat1][pat2]])
	#for classfile in classfiles:
	#	mygroupdf
	#- write classes of nodes
	for i in grdct.index:
		myrow=['//NODECLASS']
		myrow.append(i)
		myrow.append(grdct[grdct.keys()[0]][i])
		row_list.append(myrow)
	#- write rows into file
	with open(myoutname, 'w') as nf:
		for i in row_list:
			nf.write("%s\n"%('\t'.join(str(x) for x in i)))
	nf.close()

#! can be upgraded to use multiple classes
def csv_to_layout(myfilename,mygroupdf='Nan',myoutputname='default'):
	'''Function which converts .csv file to .layout'''
	data=pd.read_csv(myfilename,index_col=0)
	#- set output name or do inplace
	if myoutputname=='default':
		attrs = myfilename.split('.')[:-1]
		myoutname = '.'.join(attrs)+'.layout'
	else:
		myoutname=myoutputname
	#- set attributes file
	if mygroupdf!='Nan':
		if mygroupdf.endswith('.csv'):
			grdct=pd.read_csv(mygroupdf,index_col=0)                 #! index_col might cause problems, harmonise formats!
		else:
			grdct=pd.read_csv(mygroupdf,sep='\t',index_col=0)
	else:
		grdct={}
		for item in list(data.keys()):
			grdct[item]='none'
	#- write s,t,w format of matrix dataframe
	row_list=[]
	for pat1 in data.keys():
		for pat2 in data.index:
			if pat1!=pat2:
				row_list.append([pat1,pat2,data[pat1][pat2]])
	#for classfile in classfiles:
	#	mygroupdf
	#- write classes of nodes
	for i in grdct.index:
		myrow=['//NODECLASS']
		myrow.append(i)
		myrow.append(grdct[grdct.keys()[0]][i])
		row_list.append(myrow)
	#- write rows into file
	with open(myoutname, 'w') as nf:
		for i in row_list:
			nf.write("%s\n"%('\t'.join(str(x) for x in i)))
	nf.close()

def prog(count, total, interval=100, info=""):
	if count in range(0,total,max(int(float(total)/interval),1)): print("%s: done %s of %s"%(info, count,total))

def fix_permissions(this_path):
	os.system("/bin/chmod 755 %s"%(this_path))

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

def clearquotes(thisstring):
	thisstring = thisstring.replace('"',"")
	thisstring = thisstring.replace("'",'')
	return thisstring

def plothist(histlist, filename, label=""):
	title = "Histogram"
	values_to_remove = [-99]
	bindivisions = [100]
	log = "no" 
	c = "blue"
	v = [float(x) for x in histlist if x!="" and x not in values_to_remove]
	v = [x for x in v if x not in values_to_remove]
	print("plotting histogram:", min(v), max(v))
	for bincount in bindivisions:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		if log == "yes":
			ax.set_yscale('log')
		n, bins, patches = ax.hist(v, bincount,  facecolor=c, alpha=0.75)
		bincenters = 0.5*(bins[1:]+bins[:-1])
		ax.set_xlabel('')
		ax.set_ylabel('Frequency')
		ax.set_title(title)
		#if log == "yes":
		#	ax.set_ylim(math.pow(10,-1), len(v))
		#ax.grid(True)
		fig.savefig(filename)
		plt.close()

def run_mcxquery(matrixfile, er, qualitythresh):
	'''#! will require work to be made python 3 compatible '''
	outputfile = '.'.join(string.split(matrixfile,'.')[:-1]) + '.query'
	cmd = pathtomcx+"mcx query --clcf --output-table -imx %s -vary-threshold %d,%d,%d,100 -o %s"%(matrixfile, round(er[0]*100,0), max(1,int((er[1]-er[0])/50)), round(er[1]*100,0), outputfile)
	print(cmd)
	os.system(cmd)
	f=open(outputfile)
	lines=[string.split(string.strip(x),'\t') for x in f.readlines()]
	if len(lines)==1: # then this might be space-delimited. 
		try:
			lines=[string.split(x,' ') for x in lines]
		except:
			print("error parsing lines of mcxquery outputfile")
			print(lines)
	lines=[[string.strip(x) for x in line] for line in lines]
	lines=[x for x in lines if x!='']
	f.close()
	threshold=er[1]
	for line in lines:
		try:
			int(line[0]) # works if this is a data line
		except:
			continue
		if float(line[-3])<qualitythresh:
			try:
				threshold = float(line[-1])
			except:
				threshold=1 #corrupted file.
			break
	return outputfile, threshold

def run_mcl(pearson_file,
			inflation=2,\
			preinflation=3,\
			best_edges='no',\
			cpus=threads,\
			verbose=True, \
			redo=redo_all
			):
	#-- original line was 
	#if best_edges>0: best_edge_name="_B%s"%best_edges
	if best_edges!='no' and best_edges>0: best_edge_name="_B%s"%best_edges
	else: best_edge_name=""
	name_addition="_%s"%(best_edge_name)
	if '.' in pearson_file:
		label = '.'.join(pearson_file.split(".")[:-1])
	else:
		label = pearson_file
	mcl_out=label+name_addition+"m%s.mcl_out"%(int(inflation*10))
	#run mcl
	if not os.path.exists(mcl_out) or redo:
		if best_edges!="no": tf=" -tf '#knn("+str(best_edges)+")'"
		else: tf=""
		command2=(pathtomcl+"mcl "+pearson_file+" -I "+str(inflation)+" -pi "+str(preinflation) \
				   + tf +" -te "+str(cpus)+" -scheme 6 "+" -o " + mcl_out)
		if verbose: print(command2)
		os.system(command2)
	return mcl_out

def readmcx(mcx_correlation_file, rthreshold=0, correctpearson=0):
	'''#! this function will require work to be python 3 compatible '''
	# CONVERT MCX TO GML
	f=open(mcx_correlation_file)
	lines=f.readlines()
	f.close()
	newedges=[]
	innewedges=[]
	linelist = lines[lines.index('begin\n')+1:]       #- cuts off the begining segment of the mcl file (mclmatrix begin)
	edgevalues=[]
	for i,line in enumerate(linelist):
		#prog(i,len(linelist),10, info="readmcx")
		sline = string.split(string.strip(line,' '))
		if line[:5] != '     ' and line[0]!=")":      #- changed [:8] to [:5] and the double tab to five spaces to fit format
			a = string.split(sline[0],':')[0]         #- a stores the name of the begining node if there is no space before it
			sline=sline[1:]
		for item in sline:
			if item == "$": continue
			item = string.split(item,":")             #- item stores the connected node (the second one and its edgeweight)
			if item[0]!=")":
				if a != item[0]: #no point in writing edge between the same node.
					item[1] = float(item[1]) + correctpearson       #! what is the point of 'correctpearson'?
					#print float(item[1]), rthreshold, float(item[1]) >= rthreshold
					if float(item[1]) >= rthreshold:
						newedges.append([int(a),int(item[0]),item[1]])
						innewedges += [int(a),int(item[0])]
						edgevalues.append(item[1])
	try:
		er=(min(edgevalues), max(edgevalues))
	except:
		er=(1,1)
	return newedges, innewedges, er

def readmcl(clusterfile):
	f=open(clusterfile)
	lines=f.readlines()
	f.close()
	theseclusters = {}
	linelist = lines[lines.index('begin\n')+1:]
	for i,line in enumerate(linelist):
		#prog(i,len(linelist),10, info="readmcl")
		newcluster=False
		if line[0]!=" ":
			newcluster=True
		line=[x for x in line.strip().split(' ') if x!='' and x!="$" and x!=")"]
		if len(line)>0:
			if newcluster:
				cluster=line[0]
				theseclusters[cluster]=[int(x) for x in line[1:]]
			else:
				theseclusters[cluster]+=[int(x) for x in line]
	return theseclusters

def readtab(tabfile):
	f=open(tabfile)
	lines=[string.split(string.strip(x),'\t') for x in f.readlines()]
	f.close()
	tabs={}
	decode={}
	for i,line in enumerate(lines):
		if len(line)>1:
			tabs[int(line[0])]=clearquotes(line[1])
			decode[clearquotes(line[1])]=int(line[0])
	return tabs, decode

def make_clustering_list(clustering_dict, total_length='unspecified', minclustersize = 0):
	''' take a cluster dict and make it ready for conversion to an igraph clustering object '''
	''' #! may not be python 3 compatible as is '''
	c = {}
	for cluster in clustering_dict:
		if len(clustering_dict[cluster]) >= minclustersize:
			for vertex in clustering_dict[cluster]:
				c[vertex]=cluster
	vertices = sorted(c.keys()) # vertices are integers already by definition
	ci=[int(c[vertex]) for vertex in vertices]
	return ci

def write_json(edgelist, outputfile, outputdir, nodenames={}, theseclasses={}):
	G=nx.Graph()
	for i, edge in enumerate(edgelist):
		prog(i,len(edgelist),1,"reading nx edgelist:")
		G.add_edge(edge[0],edge[1],weight=float(edge[2])) # specify edge data
	#put in a blank edge-free node for each one that is missing
	allnodes = list(G.nodes(data=False))
	try:
		for missing in set(range(max(allnodes)+1))-set(allnodes):
			G.add_node(missing)
	except:
		pass
	for n in G:
		G.node[n]['nodesize'] = 5000
		try:
			G.node[n]['name'] = nodenames[n]
		except:
			G.node[n]['name'] = 'none'
		try:
			G.node[n]['group'] = theseclasses[nodenames[n]]      #deleted [0] inorder to keep full name of class
		except:
			G.node[n]['group'] = 'none'                          #changed nme from 'class' to 'group' to accomodate for new visualisation software
	d = json_graph.node_link_data(G)
	check_dir(os.path.join(outputdir,'data'))   #changed from resviewdir which got deleted #! could be reason why following line fails, investigate
	jo = open(os.path.join(outputdir,'data',outputfile),'w')    #changed from resviewdir which got deleted 
	json.dump(d, jo)
	jo.close()

def write_index(jsonlist, linknames, outputfile):
	f=open(os.path.join(sourcedir,'template.html'))
	lines=f.readlines()
	f.close()
	o=open(os.path.join(resviewdir,outputfile),'w')
	for line in lines:
		o.write(line)
		if '<select id="networkselection"' in line:
			for i,j in enumerate(jsonlist):
				o.write('<option value="%s">%s</option>\n'%(j,linknames[i]))
	o.close()

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

def write_mcl_matrix(G, filename):
	mcifile = ""
	mclheader = (
		("\n(mclheader\n"),
		("mcltype matrix\n"),
		( "dimensions %sx%s\n)\n" % (G.vcount(), G.vcount()) ),
		("(mclmatrix\n"),
		("begin\n")
	)
	for mhline in mclheader:
		mcifile += mhline
	for v in G.vs:
		mcifile += str(v.index) + " "
		for n in G.neighbors(v.index):
			mcifile += str(n) + " "
		mcifile += " $\n"
	mcifile += ")\n"
	o=open(filename,'w')
	o.write(mcifile)
	o.close()


#----- Potentially delete in future -----

def nda(graphobject, thisclass):
	''' return global NDA score for this network for these classes'''
	global num_perms
	permclasses = [random.sample(list(range(len(graphobject.vs))), len(thisclass)) for x in range(num_perms)]
	permclasses.append(thisclass)
	perm_totals=[]
	for c, runningclass in enumerate(permclasses):
		#make a distribution for each node (for node-specific p value)
		nodedists={nodename:[] for nodename in runningclass}
		for i, edge in enumerate(graphobject.get_edgelist()):
			try:
				nodedists[edge[0]].append(graphobject.es[i]['weight'])
			except:
				pass
			try:
				nodedists[edge[1]].append(graphobject.es[i]['weight'])
			except:
				pass
		#fill in zeros for unmeasured edges (there should be an edge for every other node) # and SORT so empirical p works
		nodedists = {x:sorted(nodedists[x]+[0 for y in range(len(nodes)-len(nodedists[x])-1)]) for x in nodedists}
		total_p=0 # sum log p for this community
		for i, node in enumerate(runningclass):
			for j, othernode in enumerate(runningclass): #symmetrical because nsp
				try:
					w = graphobject.es[graphobject.get_eid(node, othernode)]['weight']
				except:
					continue #equivalent to saying p=0
				p = empiricalp(w, nodedists[node])
				total_p += -math.log(p)
		if c < num_perms: #ie this is a permutation and not the real data
			perm_totals.append(total_p)
	return zscore(total_p, perm_totals)

def igraph_from_df(mydataframe):
	''' Function which converts a pandas dataframe into an igraph network '''
	val=mydataframe.values
	mygraph=igraph.Graph.Adjacency((val>0).tolist()) # alternatively, use val.astype(bool).tolist() or (val/val).tolist()
	mygraph.es['weight']=val[val.nonzero()]
	mygraph.vs['name']=mydataframe.index # alternatively use mydataframe.columns
	return mygraph


#- network

def scale_and_pca(thisdf, n):
	# scale data
	x = thisdf.values
	x = StandardScaler().fit_transform(x)
	# run pca
	pca = PCA(n_components=n)
	principalComponents = pca.fit_transform(x)
	return pd.DataFrame(data = principalComponents)


def write_json_net(G, outputfile, online=False):
	data = json_graph.node_link_data(G)
	jo = open(outputfile,'w')
	if online:
		jo.write('<? header("Access-Control-Allow-Origin: *") ?>\n')
	json.dump(data, jo)
	jo.close()


def write_a2b(G, outputfile, online=False):
	data = json_graph.node_link_data(G)
	with open(outputfile,'w') as o:
		for edge in data['links']:
			s = data['nodes'][edge['source']]['id']
			t = data['nodes'][edge['target']]['id']
			o.write("{}\t{}\t{}\n".format(s,t,edge['weight']))


def pairwise_difference(thisdict):
	outdict = {}
	thesekeys = list(thisdict.keys())
	for i in range(len(thesekeys)):
		for j in range(i+1, len(thesekeys)):
			try:
				outdict[thesekeys[i]]
			except:
				outdict[thesekeys[i]]={}
			outdict[thesekeys[i]][thesekeys[j]] = abs(thisdict[thesekeys[i]] - thisdict[thesekeys[j]])
	return outdict

def newext(filepath, thisext):
	return filepath[:filepath.rfind('.')] + thisext

def make_graph(edgedict, theseclasses={}):
	''' #! may require corrections for python 3 compatibility (emphasis on list(mydict.keys()[x]) principle !!!) '''
	nodenames = []
	for n in edgedict.keys():
		nodenames += list(edgedict[n].keys())
	nodenames = sorted(list(set(nodenames + list(edgedict.keys()))))
	G=nx.Graph()
	for x in edgedict:
		for y in edgedict[x]:
			if y!=x:
				G.add_edge(nodenames.index(x),nodenames.index(y),weight=edgedict[x][y]) # specify edge data
	#put in a blank edge-free node for each one that is missing
	allnodes = list(G.nodes(data=False))
	if len(allnodes)>0:
		for missing in set(range(max(allnodes)+1))-set(allnodes):
			G.add_node(missing)
		for n in G:
			G.node[n]['nodesize'] = 5000
			G.node[n]['name'] = nodenames[n]
			thesecolors = ['red', 'lightgreen', 'gold']
			for i, thisclass in enumerate(theseclasses):
				if theseclasses[thisclass][nodenames[n]]:
					#print (n, thisclass, thesecolors[i], nodenames[n])
					G.node[n]['color'] = thesecolors[i]
	return G, nodenames

def apply_threshold_and_knn(edgedict, t=0.6, knn=3):
	out = {}
	for x in edgedict:
		ycors = [(k, edgedict[x][k]) for k in sorted(edgedict[x], key=edgedict[x].get, reverse=True)][:knn+1]
		for y, c in ycors:
			if y!=x:
				if c >= t:
					try:
						out[x]
					except:
						out[x]={}
					out[x][y]=c
	return out
