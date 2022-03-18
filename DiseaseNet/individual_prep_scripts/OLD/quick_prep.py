#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import subprocess
import pandas as pd

inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/EUR.CD.gwas_info03_filtered.assoc'
inputfile2='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/EUR.IBD.gwas_info03_filtered.assoc'
inputfile3='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/EUR.UC.gwas_info03_filtered.assoc'

filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+' --N 20883 '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/CD'
sumstat1=subprocess.call(cmd1, shell=True)

cmd2='python '+filepath+'munge_sumstats.py --sumstats '+inputfile2+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+' --N 34652 '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/IBD'
sumstat2=subprocess.call(cmd2, shell=True)

cmd3='python '+filepath+'munge_sumstats.py --sumstats '+inputfile3+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+' --N 27432 '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/UC'
sumstat3=subprocess.call(cmd3, shell=True)


