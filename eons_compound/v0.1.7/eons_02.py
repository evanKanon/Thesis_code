#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

#EONS second script 
#Create all different networks

#------------------------------
#v0.0.0 initial version of the new code, adapted from EONS_02 version 0.2
#v0.0.1 substituted np network creation function for the pandas one and cosine similarity for kendall taf
#v0.0.2 added config file, EONS thresh defaults, redefined outputdirs, made single EONS_funct and deleted unused functions (version descr now includes changes in any/all scripts)
#v0.0.3 added extreme network selection functionality and the the relevant option for the distance multiplier selection, removed config hardpath and renamed EONS_main to EONS
#v0.0.4 changed MCL feature selection with PCA based on networktools code
#v0.0.5 changed feature selection to only principal components, deleted unused code
#v0.0.6 merged MCL feature selection and PCA, created multiple feature subsets through different inflation values, made it possible to notify many users simulataneously, added option to run on normalised data as well, automated static deletion of IDENTIFIER column of data
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

#------------------------------
__version__ = '0.1.7'
#------------------------------
import argparse
from eons_funct import *
#------------------------------

#- Options
parser=argparse.ArgumentParser(description='EONS-Part_II__Exhaustive_Network_Creation')
parser.add_argument('-s', '--save', type=str, help='file path for results and figures to be saved')
parser.add_argument('-f', '--file', type=str, help='file path for input file')
parser.add_argument('-a', '--analysis', type=int, choices=[1,2,3], default=3, help='choice of form of analysis: 1- full analysis, 2- quick analysis, 3- normal analysis')
parser.add_argument('-thr', '--thresholded', type=int, default=1, choices=[0,1], help='Choose whether threshlods will be applied to networks: 0- apply thresholds, 1- run using -99 threshold, making it threshold-free')
args=parser.parse_args()

#- Processing arguments 
outputdir=args.save
inputfile=args.file 
analysis=args.analysis
thresh=args.thresholded


#- ensure program is using proper features directory
print(outputdir)
if os.path.isdir(outputdir+'features/')==True:
	subsetdir=outputdir+'features/'
else:
	#! update this to output to slack
	print('ERROR: No features folder found')

#- Read-in input
mini_expression_files=[]
mini_expression_files.append(inputfile)
for filename in os.listdir(subsetdir):
	s=os.path.abspath(subsetdir)
	if filename.endswith(".expression"):
		mini_expression_files.append(s+'/'+filename)
	#miniexpressionfiles.append(os.path.abspath(filename))
print(mini_expression_files)

#- define folder to write all output of this stage
networkdir=os.path.join(outputdir, 'networks/')
check_dir(networkdir)


#-----------------------------
# CREATE NETWORKS
#-----------------------------


#old cortypes not included yet in current analysis ['dot', 'slowcosine', 'anglenorm', 'acuteanglenorm', 'sine', 'slowsine']:
#- set default parameters initially
allcortypes = ['pearson','spearman','kendall']
allthresholds = [0.4,0.6,0.8]
#- if default analysis, change nothing
if analysis==3:
	pass
#- if quick analysis change everything
elif analysis==2:
	allcortypes = ['spearman']
	allthresholds = [0.8]
	mini_expression_files = mini_expression_files[:5]
#- if full analysis, add all other distance metrics
elif analysis==1:
	#- 'correlation' distance not used because it is simply 1-pearson correlation, hence will not provide information from a really 'new' distance metric
	#- taxi == citycab == manhattan, max == chebyshev
	#-! yule is exculed as produces non-applicable values throughout
	#-! kulsinski, russellrao are excluded as it shows distance between boolean vectors
	#-! 'dice', 'hamming', 'jaccard', 'matching', 'rogerstanimoto', 'sokalmichener', 'sokalsneath' are also excluded as they are for boolean vectors and produce only complete networks with 1.0 similarity throughout
	#-! mahalanobis is excluded as it fails constantly
	othercorrtypes = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'cosine', 'euclidean', 'minkowski', 'seuclidean', 'sqeuclidean']
	allcortypes=allcortypes+othercorrtypes
#- if no thresholds, change thresholds
if thresh==1:
	allthresholds= [-99]

print(allcortypes)

failed_networks=[]
networklist=[]
# ========== create all the networks ========== 
#- for every file in the mini_expression_files
for i, thisexpression in enumerate(mini_expression_files):
	#- readin data based on their type
	thisname = os.path.split(thisexpression)[1]
	if thisexpression.endswith('.csv'):
		mydata=pd.read_csv(thisexpression,index_col=0)
	else:
		mydata=pd.read_csv(thisexpression,index_col=0,sep='\t')
	#- delete any text columns that might be found in the dataset/ only numeric columns kept
	#! this might slow things down, if that proves to be the case, consider implementing earlier in-code
	for y in mydata.columns:
		if(mydata[y].dtype == np.float64 or mydata[y].dtype == np.int64):
			pass
		else:
			mydata=mydata.drop(columns=[y])
	#- for each correlation type in the list we have specified
	for cortype in allcortypes:
		#- for each of the thresholds we have selected
		for threshold in allthresholds:
			#- this try/except is added because of the potential of some distances to fail because D=0, hence we want to keep a list of such failures but let the analysis continue anyway
			try:
				network_name=networkdir+thisname+"_%s%s"%(cortype,int(threshold*100))
				correlation_function(mydata,network_name,transpose=False,threshold=threshold,correlation=cortype)
				networklist.append(network_name)
			except:
				failure=[thisname,cortype]
				failed_networks.append(failure)


#- write log file of networks that failed to be created
with open(outputdir+'failed_networks.txt','w') as fn:
	for ni in failed_networks:
		fn.write(("%s\n")%"\t".join(ni))

