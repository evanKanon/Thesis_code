#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

#EONS preparation script 
#Format input data to format appropriate for EONS run

#------------------------------
#v0.0.7 created prep script, changed file names from capital to lowercase, changed networks so that all networks are similarity and not distance networks
#v0.0.8 gave two options for feature selection, either MCL or PCA-based feature selection, added classfile option and config path as well as working and useful eval script
#v0.1.0 BASELINE EONS added failed network functionality to eons_02 making all distances runable, corrected multiplier functionality in eons_03, made eons_eval save proper output, used similarity and not distance networks in eons_funct i.e. 0 is least similar
#v0.1.1 changed make_expression function to accomodate file names other than .soft, added NA dealing option to eons_prep, improved config option for eons_main, made code python 3 compatible
#v0.1.2 established perpendicular point plots
#v0.1.3 improved folder names, added .layout visualisation, updated environment file
#v0.1.4 made running the MCL component faster and allowed user to select inflation values
#v0.1.5 changed naming of folders and certain default values
#v0.1.6 corrected PCA random selection in eons_01, corrected network names in eons_02 and deleted unused inputfile option in eons_03 with the appropriate corrections in eons_funct and _main
#v0.1.7 added correct evaluation measure

#! prep script does not support file type other than tab delimited, need to improve range!!!!!!!!!!!!!!

#------------------------------
__version__ = '0.1.7'
#------------------------------
import argparse
from eons_funct import *
#------------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#------------------------------
parser=argparse.ArgumentParser(description='Exhaustive Observation of Network Space Analysis')
parser.add_argument('-i', '--input', type=str, help='file path for results and figures to be saved', default='0')
#! not actually currently implemented
parser.add_argument('-o', '--output', type=str, help='file path to save output')
#! can add more numbers for more options, i.e. 1 mean, 2 median, 3 MICE
parser.add_argument('-na','--na_comp', type=int, choices=[0,1], help='set whether missing data should be imputed (1) or dropped (0)', default=0)
args=parser.parse_args()

#- Processing arguments 
outputdir=args.output
inputfile=args.input
imput=args.na_comp

#- if the file is a soft file, run appropriate conversion 
new_path=make_expression_file(inputfile,path=True)

#- read in new file to see delete any lines with missing data
mdf=pd.read_csv(new_path,sep='\t')
if imput==0:
	#- just drop lines with missing values
	ndf=mdf.dropna()
else:
	#- impute missing values using MICE
	#imp = IterativeImputer(max_iter=10, random_state=0)
	#data=mdf.values
	#imp.fit(data)
	#tdata=imp.transform(data)
	#ndf=pd.DataFrame(tdata,index=mdf.index,columns=mdf.keys())
	#- impute missing values to the mean
	imp = SimpleImputer(strategy='mean')
	data=mdf.values
	imp.fit(data)
	tdata=imp.transform(data)
	ndf=pd.DataFrame(tdata,index=mdf.index,columns=mdf.keys())
	#imputer = Imputer()
	#X=mdf.values
	#transformed_X = imputer.fit_transform(X)
	#ndf=pd.DataFrame(transformed_X,index=mdf.index,columns=mdf.keys())

#- save file with complete data only
ndf.to_csv(new_path,sep='\t',index=False)

