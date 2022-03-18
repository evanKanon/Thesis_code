#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# Script used to change identifier name
#------------------------------
#v0.0.0 initial version of the new code, developed using python 3.7.2, converts Ensembl identifiers to MGI/HGNC
#------------------------------
__version__ = '0.0.0'
#------------------------------
import os
import json
import mygene #conda install not yet supported, requires pip install
import argparse
import pandas as pd
from slacker import Slacker
#------------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#------------------------------
parser=argparse.ArgumentParser(description='Gene identifier renaming script')
parser.add_argument('-i', '--input', type=str, help='file path for input file')
parser.add_argument('-o', '--output', type=str, help='file path for results and figures to be saved', default='0')
parser.add_argument('-slack','--slack', default=False, action='store_true', help='send messages to user through Slack platform')
parser.add_argument('-c', '--configfile', default=os.path.join(scriptpath,'config.json'), help='custom config file')
args=parser.parse_args()
#------------------------------

mg = mygene.MyGeneInfo()

#- open config file
with open(args.configfile) as json_data_file:
	config=json.load(json_data_file)

#- SET SLACK PARAMETERS
slack_arg=args.slack
slack=Slacker(config['Slack_env']['address'])
user=config['Slack_env']['user']

#- set input filepath
inputfile=args.input

#- keep input folder and filename as variables
inputfolder=os.path.dirname(inputfile)
inputname=os.path.basename(inputfile)

#- set outputdir filepath
if args.output=='0':
	outputdir=inputfolder+'/'
else:
	outputdir=args.output

#- post a slack message to user
if slack_arg==True:
	slack.chat.post_message(user, 'Commencing identifier conversion with '+str(inputname))
else:
	pass

#- read in data
#othdat=pd.read_csv('/Users/s1242130/Desktop/deng.csv',index_col=0)
dat=pd.read_csv(inputfile, index_col=0, sep='\t')
#- transpose data to be in proper direction
#dat=dat.transpose()
#- turn new index into column #! may be unecessary
dat.reset_index(inplace=True)



#-ensembl to symbol
#mg.getgene(name,species='mouse')['symbol']

#-symbol to ensembl
#test=mg.query(dat.keys()[5], scopes='symbol',fields=['ensembl'])
#test['hits'][0]['ensembl']['gene']


#- create list to save combinations of genes
genesymbols=[]
dat['Gene_names']=""
#- save combinations of genes
for i in range(len(dat)):
	#ens=dat['index'][i]
	ens=dat.index[i]
	try:
		symb=mg.getgene(ens,species='mouse')['symbol']
		dat['Gene_names'][i]=symb
		genesymbols.append([ens,symb])
	except:
		#- if entry is depracated, check manually
		dat['Gene_names'][i]='ERROR_CHECK'
		genesymbols.append([ens,'ERROR_CHECK'])

#- drop all genes which don't have a gene name identified
df=dat[dat.Gene_names != 'ERROR_CHECK']
#- drop previous identifier column
df=df.drop(columns=['index'])
#- change order of columns in new dataframe
cols=df.columns.tolist()
cols=cols[-1:]+cols[:-1]
df=df[cols]
#- reset index without saving it
df=df.reset_index(drop=True)
#- save file to csv format
#! naming scheme requires rethinking
df.to_csv(outputdir+inputname.split('.')[0]+'_stdt.'+inputname.split('.')[1])

if slack_arg==True:
	slack.chat.post_message(user, 'Conversion of '+str(inputname)+' completed.')
else:
	pass

'''
#- evaluate which genes were not identified
si=0
for i in range(len(othdat)):
	if othdat['Gene_names'][i]=='ERROR_CHECK':
		si+=1
'''
