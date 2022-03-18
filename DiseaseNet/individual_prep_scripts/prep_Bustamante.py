#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# prep script
#------------------------------
import subprocess
#------------------------------

#- create sumstats input file required for ldsc
inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Diarrhoeal-Disease_Bustamante_2016/2016-09-06_EAGLE_diarrhoea_dd1y.txt'
inputfile2='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Diarrhoeal-Disease_Bustamante_2016/2016-09-06_EAGLE_diarrhoea_dd2y.txt'
inputfile3='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Diarrhoeal-Disease_Bustamante_2016/2016-09-06_EAGLE_diarrhoea_d1y.txt'
inputfile4='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Diarrhoeal-Disease_Bustamante_2016/2016-09-06_EAGLE_diarrhoea_d2y.txt'
filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Diarrhoeal-Disease_Bustamante_2016_dd1yr'
sumstat1=subprocess.call(cmd1, shell=True)

cmd2='python '+filepath+'munge_sumstats.py --sumstats '+inputfile2+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Diarrhoeal-Disease_Bustamante_2016_dd2yr'
sumstat2=subprocess.call(cmd2, shell=True)

cmd3='python '+filepath+'munge_sumstats.py --sumstats '+inputfile3+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Diarrhoeal-Disease_Bustamante_2016_d1yr'
sumstat3=subprocess.call(cmd3, shell=True)

cmd4='python '+filepath+'munge_sumstats.py --sumstats '+inputfile4+' --merge-alleles '+datafiles_path+'w_hm3.snplist '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Diarrhoeal-Disease_Bustamante_2016_d2yr'
sumstat4=subprocess.call(cmd4, shell=True)

