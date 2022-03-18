#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import pandas as pd


myfile='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/EUR.CD.gwas_info03_filtered.assoc'

data=pd.read_csv(myfile,comment='#',sep='\t')

#int(round(0.0988*21770+0.0969*12882))

data['N']=None
for i in range(len(data)):
	data['N'][i]=int(round(data['FRQ_A_5956'][i]*5956+data['FRQ_U_14927'][i]*14927))



data.to_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/crohns_disease_formated.txt',sep='\t',index=False)