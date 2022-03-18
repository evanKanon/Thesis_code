#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# script which creates disease correlation network
#------------------------------
#v0.0.0 initial version of diseasenet, using folder input and ldsc code to create genomic disease correlation network
#v0.0.1 version of diseasenet, which uses only z-scores to calculate perason similarity between diseases
#v0.1.0 standard version of DiseaseNet, allows for prep based on ldsc munge_stats, choice between running DiseaseNet with LDSC or spearman correlation

#------------------------------
__version__ = '0.1.0'
#------------------------------
import os
import json
import time
import argparse
import subprocess
import pandas as pd
#from scipy.stats import spearmanr
from scipy.stats import pearsonr
#------------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#------------------------------
#- Options
#- General functionality options
parser=argparse.ArgumentParser(description='DiseaseNet: a disease association network')
parser.add_argument('-i', '--input', type=str, help='file path for input directory', default='0')
parser.add_argument('-o', '--output', type=str, help='file path for results and figures to be saved', default='0')
#currently considered as default action, could be added as option
#parser.add_argument('-ldsc', '--ldsc_analysis', action='store_true', help='create DiseaseNet using LDSC genetic correlation to create the network')
parser.add_argument('-pear', '--pearson_analysis', default=False, help='create DiseaseNet using spearman correlation to create the network', action='store_true')
parser.add_argument('-prep', '--data_preparation', default=False, help='run LDSC munge_stats prep script to standardize input', action='store_true')
parser.add_argument('-hp', '--use_hard', default=False, action='store_true', help='use custom config file instead of hardpath for file paths, NOTE all other paths in config will be used')

#parser.add_argument('-log','--log_file', default=False, action='store_true', help='create log file containing the argunents used, the date run, the version and the input file name')
parser.add_argument('-c', '--configfile', default=os.path.join(scriptpath,'config.json'), help='custom config file')
args=parser.parse_args()
#------------------------------

#- open config file
with open(args.configfile) as json_data_file:
	config=json.load(json_data_file)

#- SET MAIN PATHS
if args.use_hard==False:
	#- define inputfile path
	if config['File_paths']['input']:
		if args.input=='0':
			inputdir=str(config['File_paths']['input'])
		else:
			inputdir=str(config['File_paths']['input'])+args.input
	else:
		inputdir=args.input
	#- define outputdir path
	if config['File_paths']['output']:
		if args.output=='0':
			outputdir=str(config['File_paths']['output'])
		else:
			outputdir=str(config['File_paths']['output'])+args.output
	else:
		outputdir=args.output
else:
	inputdir=args.input
	outputdir=args.output


#- define filepath for ldsc scripts
try:
	filepath=config['File_paths']['ldsc']
except:
	#if slack_arg:
	print('ERROR: No filepath detected for LDSC, please provide full filepath in config file')
#- define filepath for reference data files
try:
	datafiles_path=config['File_paths']['reference']
except:
	#if slack_arg:
	print('ERROR: No filepath detected for the reference datafiles, please provide full filepath in config file')
#------------------------------
#datadir='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/marie/data/imputed_gwas_hg38_1.1/'
#filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
#datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'
#outputdir='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/new_input/'
#------------------------------
#! simplify and integrate
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
#------------------------------

#- Create overarching output folder
daterun='_'+str(time.strftime("%d%m%Y"))
#analysislabel = os.path.split(inputfile)[1].split('.')[0]
outputdir = os.path.join(outputdir,'DiseaseNet_v'+str(__version__)+daterun)#+'/'
if args.pearson_analysis==False:
	outputdir=outputdir+'_ldsc/'
else:
	outputdir=outputdir+'_pearson/'
check_dir(outputdir)

#- create results directory to store resulting network
resultsdir = os.path.join(outputdir, 'network/')
check_dir(resultsdir)

#- create comparison directory to store files with comparison values
compdir = os.path.join(outputdir, 'ldsc-logs/')
check_dir(compdir)



#----------------
# DATA PREP
#----------------

if args.data_preparation==True:
	#- read in data
	mydata=[]
	for filename in os.listdir(datadir):
		s=os.path.abspath(datadir)
		if filename.endswith('.txt'):
			mydata.append(s+'/'+filename)

	#- create folder to save prepped output
	inputdir=os.path.join(outputdir+'/prepped_input/')
	check_dir(inputdir)

	#- run sumstats prep from LDSC
	for dataset in mydata:
		cmd='python '+filepath+'munge_sumstats.py --sumstats '+dataset+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+'--out '+inputdir+dataset.split('/')[-1].split('.')[0]+' --snp variant_id --N-col sample_size --a1 effect_allele --a2 non-effect_allele --ignore  n_cases'
		sumstat=subprocess.call(cmd, shell=True)




#-----------------
# Create LDSC network
#-----------------


if args.pearson_analysis==False:
	#- save names of formatted files to list
	formatted=[]
	for filename in os.listdir(inputdir):
		s=os.path.abspath(inputdir)
		if filename.endswith('.sumstats.gz'):
			s=os.path.abspath(inputdir)
			formatted.append(s+'/'+filename)


	#- create seperate genetic correlation files for each comparison
	for i in range(len(formatted)):
		x=max(range(len(formatted)))
		while x!=i: #(i+1 possibly?)
			dis1 = formatted[i]
			dis2 = formatted[x]
			name1 = dis1.split('/')[-1].split('.')[0]
			name2 = dis2.split('/')[-1].split('.')[0]
			mystring=dis1+','+dis2
			out_id=name1+'_'+name2
			cmd2='python '+filepath+'ldsc.py --rg '+mystring+' --ref-ld-chr '+datafiles_path+'eur_w_ld_chr/ '+'--w-ld-chr '+datafiles_path+'eur_w_ld_chr/ '+'--out '+compdir+out_id
			ldsc_corr=subprocess.call(cmd2, shell=True)
			if ldsc_corr==0:
				pass
			else:
				print("ERROR can't compute with file combination: "+mystring)
			x+=(-1)


	#- create single dataframe result of all associations
	res=pd.DataFrame()
	for filename in os.listdir(compdir):
		s=os.path.abspath(compdir)
		if filename.endswith('.log'):
			s=os.path.abspath(compdir)
			#- keep results part of logfile in list format
			mylines = []                                # Declare an empty list.
			with open(s+'/'+filename, 'rt') as myfile:     # Open lorem.txt for reading text.
				for myline in myfile:                   # For each line in the file,
					mylines.append(myline.rstrip('\n')) # strip newline and add to list.
			for r in range(len(mylines)):
				if mylines[r]=='Summary of Genetic Correlation Results':
					mysubset=mylines[r:]
					break
			for ir in range(len(mysubset)):
				if mysubset[ir]=='':
					mysubset=mysubset[:ir]
					break
			#- create list of lists elements in order to create dataframe
			listoflists=[]
			for line in mysubset[1:]:
				templist=line.split(' ')
				mylist=[]
				for i in templist:
					if i=='':
						pass
					else:
						mylist.append(i)
				listoflists.append(mylist)
			#- create pandas dataframe of current result
			df=pd.DataFrame(listoflists[1:],columns=listoflists[0])
			#- merge with previous results
			res=pd.concat([res,df], axis=0, join='outer', join_axes=None, ignore_index=True, names=None, copy=True)

	#- save resulting dataframe as csv for potential further visualisation
	res.to_csv(resultsdir+'result_df.csv',index=False)
else:
	#! add gunzipper functionality 

	#- save filepaths of all relevant files
	filelist=[]
	for file in os.listdir(inputdir):
		s=os.path.abspath(inputdir)
		if file.endswith('.sumstats'):
			filelist.append(s+'/'+file)


	#- calculate similarity by pearson correlating z-scores of SNPs with the same rs
	lines=[['d1','d2','p']]
	for i in range(len(filelist)):
		x=max(range(len(filelist)))
		data=pd.read_csv(filelist[i],sep='\t')
		df1=data.dropna()
		while x!=i:
			line=[]
			data2=pd.read_csv(filelist[x],sep='\t')
			df2=data2.dropna()
			inter=pd.merge(df1,df2,how='inner',on=['SNP'])
			similarity=pearsonr(inter.Z_x,inter.Z_y)[0]
			line.append(filelist[i].split('/')[-1].split('.')[0])
			line.append(filelist[x].split('/')[-1].split('.')[0])
			line.append(similarity)
			lines.append(line)
			x+=(-1)

	#- create dataframe from list of lines
	result=pd.DataFrame.from_records(lines[1:],columns=lines[0])
	#- save result to csv
	result.to_csv(resultsdir+'network.csv')



