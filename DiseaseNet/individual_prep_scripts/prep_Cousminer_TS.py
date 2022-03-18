#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# prep script
#------------------------------
import subprocess
#------------------------------

#- create sumstats input file required for ldsc
inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Tanner-Stage_Cousminer_2014/EGG_TANNER_males_and_females_combined.v2.txt'

#- hard paths to necessary folders
filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Tanner-Stage_Cousminer_2014'
sumstat1=subprocess.call(cmd1, shell=True)
