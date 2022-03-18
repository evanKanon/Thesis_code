#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import subprocess
import pandas as pd
#------------------------------

inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Inflammatory-Bowel-Disease_Liu_2015/EUR.CD.gwas_info03_filtered.assoc'
inputfile2='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Inflammatory-Bowel-Disease_Liu_2015/EUR.IBD.gwas_info03_filtered.assoc'
inputfile3='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Inflammatory-Bowel-Disease_Liu_2015/EUR.UC.gwas_info03_filtered.assoc'

filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+' --N 20883 '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Crohns_Disease_Liu_2015'
sumstat1=subprocess.call(cmd1, shell=True)

cmd2='python '+filepath+'munge_sumstats.py --sumstats '+inputfile2+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+' --N 34652 '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Iinflammatory_Bowel_Disease_Liu_2015'
sumstat2=subprocess.call(cmd2, shell=True)

cmd3='python '+filepath+'munge_sumstats.py --sumstats '+inputfile3+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+' --N 27432 '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Ulcerative_Colitis_Liu_2015'
sumstat3=subprocess.call(cmd3, shell=True)


