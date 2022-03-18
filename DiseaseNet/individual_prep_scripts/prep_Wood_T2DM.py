#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# prep script
#------------------------------
import subprocess
#------------------------------

#- create sumstats input file required for ldsc
inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Type-2-Diabetes_Wood_2016/t2d_dom_dev.txt'

#- hard paths to necessary folders
filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist --ignore MarkerName '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Serum-Urate_Kottgen_2013'
sumstat1=subprocess.call(cmd1, shell=True)
