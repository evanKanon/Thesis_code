#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

#CellMapper

import datetime, time, os, sys, string, json, operator
import pandas as pd 
import scipy as sp
from scipy import stats 
import subprocess
import urllib2
from bisect import *

from slacker import Slacker 

#-----------------------------
# Definition of functions
#-----------------------------

#-------------------------------------
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

def convert_f5_samplename(thisname):
	n = string.split(thisname,'.')
	if len(n)>1:
		return urllib2.unquote(n[1])
	else:
		return thisname

def write_edgefile(edgelist,outname):
	''' Function which writes edges (f,t,w) from a list of lists into a file '''
	o=open(outname,'w')
	for edge in edgelist:
		o.write("%s\t%s\t%s\n"%edge)
	o.close()

def list_uniq(inputlist):
	''' Function which creates list with unique elements of input list'''
	nodelist_uniq=[]
	for x in inputlist:
		if x not in nodelist_uniq:
			nodelist_uniq.append(x)
	return nodelist_uniq

def k_nearest_neighbour(myfile,k,outname,write=True):
	''' Function which creates a ncol file of the k-nearest neighbours of a dataset '''
	#! improve functionnality to include raw data input
	#- initiate a list to keep the edges in
	edgelist=[]
	#- read in file input (conditional to inputfile type)
	if myfile.endswith('.csv'):
		#- read in matrix file 
		#! could designate specific file type known as .matrix
		mydata=pd.read_csv(myfile, header=None)
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
		if write==True:
			#- write edgelist into a file
			write_edgefile(edgelist,outname)
		else:
			#- if we don't want to write a file, we can just return the edge list
			return edgelist
	elif myfile.endswith('.txt'):
		#- read in ncol-like file instead of matrix
		mydata=pd.read_csv(myfile,sep=" ",header=None)
		nodelist=mydata[1]
		#- create list with unique nodes
		nodelist_uniq=list_uniq(nodelist)
		#! could turn this into a function called k_neighbours docstring 'function which finds k nearest neighbours'
		for fnode in nodelist_uniq:
			templist=mydata.loc[mydata[1]==fnode]
			maxval=max(templist[2].sort_values())
			#! needs improvement with respect to pandas syntax
			if templist[0][templist.loc[templist[2]==maxval].index[0]]==fnode:
				for i in range(k):
					tnode=templist[0][templist.nlargest(k+1,2).index[i+1]]
					weight=templist.nlargest(k+1,2)[2][templist.nlargest(k+1,2).index[i+1]]
					edgelist.append((fnode,tnode,weight))
			else:
				for i in range(k):
					tnode=templist[0][templist.nlargest(k,2).index[i]]
					weight=templist.nlargest(k,2)[2][templist.nlargest(k,2).index[i]]
					edgelist.append((fnode,tnode,weight))
		if write==True:
			#- write edgelist into a file
			write_edgefile(edgelist,outname)
		else:
			#- if we don't want to write a file, we can just return the edge list
			return edgelist
	elif myfile.endswith('.ncol'):
		#- read in ncol file instead of matrix
		mydata=pd.read_table(myfile,header=None)
		nodelist=mydata[1]
		#- create list with unique nodes
		nodelist_uniq=list_uniq(nodelist)
		for fnode in nodelist_uniq:
			templist=mydata.loc[mydata[1]==fnode]
			maxval=max(templist[2].sort_values())
			#! needs improvement with respect to pandas syntax
			if templist[0][templist.loc[templist[2]==maxval].index[0]]==fnode:
				for i in range(k):
					tnode=templist[0][templist.nlargest(k+1,2).index[i+1]]
					weight=templist.nlargest(k+1,2)[2][templist.nlargest(k+1,2).index[i+1]]
					edgelist.append((fnode,tnode,weight))
			else:
				for i in range(k):
					tnode=templist[0][templist.nlargest(k,2).index[i]]
					weight=templist.nlargest(k,2)[2][templist.nlargest(k,2).index[i]]
					edgelist.append((fnode,tnode,weight))
		if write==True:
			#- write edgelist into a file
			write_edgefile(edgelist,outname)
		else:
			#- if we don't want to write a file, we can just return the edge list
			return edgelist

def empiricalp(thisvalue, empiricalcollection):
	if len(empiricalcollection) == 0:
		return 1
	#empiricalcollection = sorted(empiricalcollection) #collection is alredy sorted, saves time
	return 1-float(bisect_left(empiricalcollection,thisvalue))/len(empiricalcollection)

def empiricalp_fract(thisvalue, empiricalcollection):
	if len(empiricalcollection) == 0:
		return 1
	#empiricalcollection = sorted(empiricalcollection) #collection is alredy sorted, saves time
	num=len(empiricalcollection)-bisect_left(empiricalcollection,thisvalue)
	denom=len(empiricalcollection)
	return str(num)+'/'+str(denom)

def empiricalp_both(thisvalue, empiricalcollection1, empiricalcollection2):
	#empiricalcollection = sorted(empiricalcollection) #collection is alredy sorted, saves time
	num1=len(empiricalcollection1)-bisect_left(empiricalcollection1,thisvalue)
	denom1=len(empiricalcollection1)
	num2=len(empiricalcollection2)-bisect_left(empiricalcollection2,thisvalue)
	denom2=len(empiricalcollection2)
	num=num1+num2
	denom=denom1+denom2
	return float(num)/denom

def find_genes(xranks,yranks,mapdf,querydf3,common_genes,maxrank=10,maxgenes=20):
	""" Function which finds genes based on either their rank or quantity 
		input is lists of ranks for either, their original datasets and the dictionary which contains the link between the names and identifiers
		#! this function is NOT generalisable, requires the linker dictionary common_genes need to make it so
	"""
	myrank=0
	over_gene_list=[]
	rank_gene_dict={}
	#- select genes that drive similarity (either number of genes or genes which form part of particular rank differences, or both)
	while len(over_gene_list)<maxgenes and myrank<maxrank:
		#- find list of ranks that are the same by a value of distance (i.e. that are the same [0 distance], that differ by one [1 distance] etc)
		rank_list=[]
		for ri in range(len(xranks)):
			if abs(xranks[ri]-yranks[ri])==myrank:
			#if xranks[ri]==yranks[ri]+myrank or xranks[ri]==yranks[ri]-myrank:
				rank_list.append((ri,xranks[ri],yranks[ri]))
		#- create a list of only the indexes, to be used for name identification
		index_list=[]
		for j in range(len(rank_list)):
			index_list.append(rank_list[j][0])
		#- start adding values to the corresponding dictionary
		rank_gene_dict[myrank]=[]
		#- find which gene names correspond to the particular indexes based on the common_genes dictionary
		for myindex in index_list:
			#! this could be problematic, cuts when list is up to 20, arbitrarily
			if mapdf[mapdf.keys()[0]][myindex]==querydf3[querydf3.keys()[0]][myindex] and len(over_gene_list)<maxgenes:
				gene_region=mapdf[mapdf.keys()[0]][myindex]
				gene_row=common_genes[common_genes[0]==gene_region]
				gene_name=gene_row.iloc[0][1]
				over_gene_list.append(gene_name)
				rank_gene_dict[myrank].append(gene_name)
				#linelist.append((fnode,tnode,weight,myrank,gene_name))   # this is the line that writes
				#myline.append(gene_name)
			elif mapdf[mapdf.keys()[0]][myindex]!=querydf3[querydf3.keys()[0]][myindex]:
				print('ERROR: Manual check needed')
				print(mapdf[mapdf.keys()[0]][myindex])
				print(querydf3[querydf3.keys()[0]][myindex])
			else:
				pass
				#print('Loop over')
		myrank+=0.5
	return over_gene_list, rank_gene_dict


#-------------------------------------

