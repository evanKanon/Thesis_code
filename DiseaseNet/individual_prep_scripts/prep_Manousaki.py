#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# prep script
#------------------------------
import subprocess
#------------------------------

#- unzip file
filename='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Vitamin-D_Manousaki_2017/vitd.gwama.out.gz'
inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Vitamin-D_Manousaki_2017/vit_D_man.txt'
unzip=subprocess.call('gunzip -c '+filename+' > '+inputfile1,shell=True)

#- hard paths to necessary folders
filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist --ignore beta,beta_95L,beta_95U,eaf,q_statistic,q_p-value,i2,effects --N-col n_samples '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Serum-Urate_Kottgen_2013'
sumstat1=subprocess.call(cmd1, shell=True)
