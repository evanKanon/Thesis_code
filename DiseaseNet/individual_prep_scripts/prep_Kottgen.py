#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# prep script
#------------------------------
import subprocess
import pandas as pd
#------------------------------

#- create sumstats input file required for ldsc
orig='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/GWAS/Serum-Urate_Kottgen_2013/GUGC_MetaAnalysis_Results_UA.csv'
data=pd.read_csv(orig)
data.to_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/temp/GUGC_MetaAnalysis_Results_UA.txt',sep='\t',index=False)

inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/temp/GUGC_MetaAnalysis_Results_UA.txt'

#- hard paths to necessary folders
filepath='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/scripts/ldsc-master/'
datafiles_path='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/Data/reference/'

cmd1='python '+filepath+'munge_sumstats.py --sumstats '+inputfile1+' --merge-alleles '+datafiles_path+'w_hm3.snplist --p p_gc --signed-sumstats beta,0 --snp MarkerName  --N-col n_total '+'--out /mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/Serum-Urate_Kottgen_2013'
sumstat1=subprocess.call(cmd1, shell=True)
