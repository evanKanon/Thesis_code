#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import os
import argparse
import subprocess
import pandas as pd
#------------------------------
#- Options
parser=argparse.ArgumentParser(description='DiseaseNet data concatenation')
parser.add_argument('-i', '--input', type=str, help='file path for input directory', default='/mnt/d/PhD/input/Diseasenet/complete_input/')
parser.add_argument('-un', '--unzip', action='store_true', help='use gunzipper to unzip .gz files', default=False)
args=parser.parse_args()
#------------------------------


inputdir=args.input


if args.unzip==True:
	cmd='python gunzipper.py -i '+inputdir
	sub1=subprocess.call(cmd,shell=True)


GWAS=[]
maindf='empty'
old=0
for filename in os.listdir(inputdir):
	s=os.path.abspath(inputdir)
	if filename.endswith(".sumstats"):
		file=s+'/'+filename
		df1=pd.read_csv(file,sep='\t')
		name=file.split('/')[-1].split('.')[0]
		if isinstance(maindf, str):
			old=name
			maindf=df1.loc[:, ['SNP','Z']]
		else:
			subdf=df1.loc[:, ['SNP','Z']]
			maindf=maindf.merge(subdf,on='SNP',how='inner',suffixes=[old,name])
			old=name


outputdir='/'.join(inputdir.split('/')[:-2])

maindf.to_csv(outputdir+'comulative_GWAS.txt', sep='\t', index=False)