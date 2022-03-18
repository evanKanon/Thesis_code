#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import pandas as pd


myfile='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/EUR.IBD.gwas_info03_filtered.assoc'

data=pd.read_csv(myfile,comment='#',sep='\t')

#int(round(0.0988*21770+0.0969*12882))

data['N']=None
for i in range(len(data)):
	data['N'][i]=int(round(data['FRQ_A_12882'][i]*12882+data['FRQ_U_21770'][i]*21770))



data.to_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/inflammatory_bowel_disease_formated.txt',sep='\t',index=False)