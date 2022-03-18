#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

import os
import json
import time
import argparse
import subprocess
import pandas as pd

#------------------------------
workingdir='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/'
networkdir='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/output/'
filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'
#------------------------------


#- save names of formatted files to list
formatted=[]
for filename in os.listdir(workingdir):
	s=os.path.abspath(workingdir)
	if filename.endswith('.sumstats.gz'):
		s=os.path.abspath(workingdir)
		formatted.append(s+'/'+filename)

#- create seperate correlation files for each comparison
for i in range(len(formatted)):
	x=max(range(len(formatted)))
	while x!=i: #(i+1 possibly?)
		dis1 = formatted[i]
		dis2 = formatted[x]
		name1 = dis1.split('/')[-1].split('.')[0]
		name2 = dis2.split('/')[-1].split('.')[0]
		mystring=dis1+','+dis2
		out_id=name1+'_'+name2
		cmd2='python '+filepath+'ldsc.py --rg '+mystring+' --ref-ld-chr '+datafiles_path+'eur_w_ld_chr/ '+'--w-ld-chr '+datafiles_path+'eur_w_ld_chr/ '+'--out '+networkdir+out_id
		ldsc_corr=subprocess.call(cmd2, shell=True)
		if ldsc_corr==0:
			pass
		else:
			print("ERROR can't compute with file combination: "+mystring)
		x+=(-1)

#! ----- DELETE AFTER MERGING CODE GETS WRITTEN ----- !#

#- create single string with all input file names to run the comamnd only once
#! if this fails in implementation, create loop to do it pairwise with half matrix
mystring=''
for item in formatted:
	mystring+=item
	mystring+=','
#- delete last character from string (it's a superfluous comma)
mystring=mystring[:-1]


#- run command to calculate genetic correlations between the inputted diseases
cmd2='python '+filepath+'ldsc.py --rg '+mystring+' --ref-ld-chr '+datafiles_path+'eur_w_ld_chr/ '+'--w-ld-chr '+datafiles_path+'eur_w_ld_chr/ '+'--out '+networkdir+'net'
ldsc_corr=subprocess.call(cmd2, shell=True)
if ldsc_corr==0:
	pass
else:
	print('ERROR with prep of file '+mystring)

