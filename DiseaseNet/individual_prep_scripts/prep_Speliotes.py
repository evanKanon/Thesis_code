#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# prep script
#------------------------------
import subprocess
#------------------------------

#- create sumstats input file required for ldsc
inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/BMI_Speliotes_2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt'
filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist --a1-inc '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/BMI_Speliotes_2010'
sumstat1=subprocess.call(cmd1, shell=True)


