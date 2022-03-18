#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

#EONS first script 
#Create different subset of genes

#------------------------------
#v0.0.0 initial version of the new code, adapted from EONS_01 version 0.2
#v0.0.1 added ability to run csv files (needs more work, this is note yet general) and generalised the code required to write subset files
#v0.0.2 added config file, changed EONS_02 defaults, redefined outputdirs, made single EONS_funct and deleted unused functions (version descr now includes changes in any/all scripts)
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

parser=argparse.ArgumentParser(description='EONS-Part_I__Selection_of_Features')
parser.add_argument('-s', '--save', type=str, help='file path for results and figures to be saved')
parser.add_argument('-f', '--file', type=str, help='file path for input file')
parser.add_argument('-chs','--choose_threshold',type=int, default=1, choices=[0,1], help='use mcxquery to choose gene network threshold')
parser.add_argument('-mcl','--use_mcl',type=int, default=0, choices=[0,1], help='use MCL in order to create desired features')
parser.add_argument('-infl','--infl_values',nargs='*',type=float, default=[x / 10 for x in range(10,32,2)], help='set list of inflation values to use for MCL feature selection method')
parser.add_argument('-feat','--feature',type=int, default=19131261, help='set number of features for PCA feature selection method')
parser.add_argument('-norm', '--normalise', type=int, choices=[0,1], default=0, help='normalise input and use that as additional data input')
args=parser.parse_args()

#- Processing arguments 
inputfile=args.file
outputdir=args.save
if args.choose_threshold==0:
	choose_threshold=False
else:
	choose_threshold=True
#- set number of components
nc=args.feature
#- set variable to normalise data
norm=args.normalise
#- set option to use MCL instead of PCA
mcl_use=args.use_mcl
#- set inflation values for MCL
inflation_list=args.infl_values

#===================================================================
#===================================================================
# COMMENCING ANALYSIS
#===================================================================
#===================================================================


#- CREATE THE DIRECTORIES
#- create output directory
check_dir(outputdir)
#- create a sub-directory for the final features that will be used
featuredir = os.path.join(outputdir, 'features/')
check_dir(featuredir)

#- Read in input file 
#- if the filename doesn't end in .csv it will be treated as a tab delimited file
if inputfile.endswith('.csv'):
	mydata=pd.read_csv(inputfile,index_col=0)
else:
	mydata=pd.read_csv(inputfile,index_col=0,sep='\t')


#- delete any text columns that might be found in the dataset/ only numeric columns kept
for y in mydata.columns:
	if(mydata[y].dtype == np.float64 or mydata[y].dtype == np.int64):
		pass
	else:
		mydata=mydata.drop(columns=[y])

#- save patient identifiers in list
pat_list=list(mydata.columns.values)
#! if in the final implementation we will need rows (unlikely) this can be done by:
#pat_list=list(mydata.index)


#- either use PCA or MCL to create feature subsets
if mcl_use==0:
	#- create appropriate directory for output of this subsection
	prepdir = os.path.join(featuredir,'PCA_prep/')
	check_dir(prepdir)
	#- transpose data for PCA run on features
	pcadata=mydata.transpose()
	#- check if set number of components for pca. if not, do all
	if nc==19131261:
		nc=min(len(pcadata.keys()), len(pcadata.index))
	#- prepare data for pca
	x = pcadata.values
	x = StandardScaler().fit_transform(x)
	#- perform pca
	#-! random selection of random state to maintain the same result across different applications on the same dataset, comes into play with large datasets
	pca = PCA(n_components=nc,random_state=13)
	principalComponents = pca.fit_transform(x)
	#- turn result into dataframe with proper index
	principalDf=pd.DataFrame(data = principalComponents)
	principalDf.index=pcadata.index
	#- save prinicpal components as csv file
	principalDf.to_csv(prepdir+'PCA_result.csv')
	#- create list of components
	comps=[]
	for x in range(1,nc+1):
		mystr='Comp_{}'.format(x)
		comps.append(mystr)
	#- find which feature forms part of which principal component
	table_format=pd.DataFrame(pca.components_,columns=pcadata.columns,index = comps)
	############################################################################################
	#- select components which have more than 0.5 similarity in any given component
	#for this_comp in table_format.index[:1]:
	#	myvariables=[]
	#	for variable in table_format.columns:
	#		if table_format[variable][this_comp]>=0.5:
	#			myvariables.append(variable)
	#		print(table_format[this_comp])
	#- select features based on the additive variance explained (if more than 0.3 here)
	#ncomps=0
	#sum_var=0
	#for i in pca.explained_variance_ratio_:
	#	sum_var+=i
	#	ncomps+=1
	#	if sum_var>=0.3:
	#		break
	############################################################################################
	#- find roughly how many elements correspond to 10% for the particular sample
	#! this should be made user selected
	thr=0.1
	myval=int(round(len(table_format.columns)*thr))
	for compkey in table_format.index:
		#! one option is to simply chop the components according to the top 10%, opt out for now
		#mysubset=table_format.T.nlargest(myval,compkey)
		mysubset=table_format.T.nlargest(myval,compkey).index
		mysdf=mydata[mydata.index.isin(mysubset)]
		mysdf.to_csv(featuredir+'%s_t%s.expression'%(compkey,int(thr*100)),sep='\t')
	#table_format.T.nlargest(3,'Comp 2')
elif mcl_use==1:
	#- create appropriate directory for output of this subsection
	workingdir = os.path.join(featuredir, 'MCL/')
	check_dir(workingdir)
	#- select to normalise data and run analysis on both types
	alldata={}
	alldata['standard']=mydata
	if norm==1:
		normdata=scale(mydata)
		alldata['scaled']=normdata
	#- for each of the different data inputs (i.e. standard and normalised)
	for dtype in alldata:
		#prefix=str(dtype)
		data=alldata[dtype]
		#-----------------------------------------MCL_component------------------------------------------------
		# CREATE INITIAL CLUSTERING OF GENES TO USE FOR FEATURE SELECTION
		#? what is the metric that we are looking for?
		#! the range of inlfation values should be user directed
		for infl_v in inflation_list:
			# FIRST DECIDE ON CLUSTERING PARAMETERS
			if choose_threshold:
				er=correlation_function(data,workingdir+dtype+'_genenet_spear',transpose=True,correlation='spearman',threshold=0.6,er=True)
				write_ncol(workingdir+dtype+'_genenet_spear.csv')
				'''
				#- create spearman feature similarity network
				cordf=data.T.rank(ascending=True).corr('pearson')
				#- set all values smaller than the threshold as equal to zero
				cordf[cordf<=0.6]=0
				#- identify edge range for feature similarity network
				er=(max(cordf.max()),min(cordf.min()))
				#- create a2b semi matrix
				a2b_list=[]
				for i in range(len(cordf.keys())):
					x=max(range(len(cordf.keys())))
					while x!=i:
						a2b_list.append([cordf.keys()[i],cordf.keys()[x],cordf[cordf.keys()[i]][cordf.keys()[x]]])
						x+=(-1)
				'''
				#- turn output into MCL readable format
				speargene_matrix_file, speargene_tab_file = run_mcxload(workingdir+dtype+'_genenet_spear.ncol')
				# CHOOSE A THRESHOLD 
				queryfile, threshold = run_mcxquery(speargene_matrix_file, er, 0.6)
			else:
				threshold=0.85
			print('Selected threshold is:  '+str(threshold))
			# NOW USE THE CHOSEN SPEARMAN THRESHOLD TO CREATE A MASTER CLUSTERING OF CHARACTERISTICS
			#- create teh feature-feature network
			correlation_function(data,workingdir+dtype+'_genenet_pears',transpose=True,correlation='pearson',threshold=threshold)
			write_ncol(workingdir+dtype+'_genenet_pears.csv')
			gene_matrix_file, gene_tab_file = run_mcxload(workingdir+dtype+'_genenet_pears.ncol')
			#- create clusters
			gene_cluster_file = run_mcl(gene_matrix_file,inflation=infl_v)
			#- input clusters as readable object in python
			gene_clusters = readmcl(gene_cluster_file)
			#-----------------------------
			#CREATE MINI EXPRESSION FILES FOR EACH GENE CLUSTER
			# why use clusters for feature selection?
			# ==> one of the defining characteristics of a syndrome is that it is a constellation of features
			# each feature in the cluster supports the other features. 
			#-----------------------------
			#- read gene tabfile for name extraction based on mcl index
			#gene_names, decode_gene_names = readtab(gene_tab_file)
			gene_id_list=read_dict(gene_tab_file)
			#- enter the subselections of genes (the clusters found above)
			for i, clustername in enumerate(gene_clusters):
				#! another 'secret' threshold
				if len(gene_clusters[clustername]) > 10:
					#- name follows reasoning 'mcl infaltion value x subset y'
					newfile = os.path.join(featuredir,dtype+"_imcl%s_subset%s.expression"%(infl_v,clustername))
					o=open(newfile,'w')
					#o.write("Patient\t%s\n"%('\t'.join(pat_list)))
					o.write("Patient\t%s\n"%('\t'.join([str(x) for x in list(pat_list)])))  #- this is added because patient names are numpy array instead of string on their own, evaluate if this can work with normal input and keep
					for gene_id in gene_clusters[clustername]:
						gene_name=gene_id_list[str(gene_id)]
						o.write("%s\t%s\n"%(gene_name, '\t'.join([str(x) for x in list(data.loc[str(gene_name)])])))
					o.close()
