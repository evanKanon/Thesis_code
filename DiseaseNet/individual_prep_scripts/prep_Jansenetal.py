#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# prep script
#------------------------------
import subprocess
import numpy as np
import pandas as pd
#------------------------------

#- correct formatting errors 

with open('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Alzheimers_Jansen_2018/AD_sumstats_Jansenetal.txt') as f:
	content = f.readlines()

content[0]=content[0].split(' ')[0]+'\t'+content[0].split(' ')[-1]

fo = open('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Alzheimers_Jansen_2018/AD_sumstats_Jansenetal_edited.txt', "w")
fo.writelines(content)
fo.close()

#- create sumstats input file required for ldsc
inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Alzheimers_Jansen_2018/AD_sumstats_Jansenetal_edited.txt'
filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+'--ignore BETA --N-col Nsum '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Alzheimers_Jansen_2018'
sumstat1=subprocess.call(cmd1, shell=True)


