#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import subprocess
import pandas as pd

myfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/GAINS/GAINs_LD_hub_format.txt'

df=pd.read_csv(myfile1,sep='\t')

df=df.fillna(0)

df = df.rename(columns={'snpid': 'SNP'})

df.to_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/GAINS/GAINs.txt',sep='\t',index=False)


inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/GAINS/GAINs.txt'

filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/GAINS'
sumstat1=subprocess.call(cmd1, shell=True)
