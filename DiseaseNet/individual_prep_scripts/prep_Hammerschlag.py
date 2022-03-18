#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# prep script
#------------------------------
import subprocess
#------------------------------

#- create sumstats input file required for ldsc
inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Insomnia_Hammerschlag_2017/Hammerschlag_NatGenet2017_insomnia_sumstats-full_090617.txt'
filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist --ignore BETA '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Insomnia_Hammerschlag_2017'
sumstat1=subprocess.call(cmd1, shell=True)


