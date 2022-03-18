#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import pandas as pd


myfile='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/EUR.UC.gwas_info03_filtered.assoc'

data=pd.read_csv(myfile,comment='#',sep='\t')

#int(round(0.0988*21770+0.0969*12882))

data['N']=None
for i in range(len(data)):
	data['N'][i]=int(round(data['FRQ_A_6968'][i]*6968+data['FRQ_U_20464'][i]*20464))



data.to_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/ulcerative_colitis_formated.txt',sep='\t',index=False)