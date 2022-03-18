#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import subprocess

inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/multiple_sclerosis_formated.txt'

filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/MS'
sumstat1=subprocess.call(cmd1, shell=True)
