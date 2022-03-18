#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import pandas as pd


myfile='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/phs000139.pha002854.txt'

data=pd.read_csv(myfile,comment='#',sep='\t')

data.columns = data.columns.str.replace(' ', '_')

data.rename(columns={'Odds_ratio': 'OR'}, inplace=True)
data.rename(columns={'SNP_ID': 'SNP'}, inplace=True)


#- droping is done in this way because the ldsc virtual env has pandas version 0.20.3
#- proper code would be the following
#data=data.drop(columns=['Submitted SNP ID','Transmit allele1','Transmit allele2'])
data=data.drop('Submitted_SNP_ID',axis=1)
data=data.drop('Transmit_allele1',axis=1)
data=data.drop('Transmit_allele2',axis=1)

#- sample size column added with 950*3, since they said that 950 parent-affected chold trios were genotyped
data['N']='2850'


data.to_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/multiple_sclerosis_formated.txt',sep='\t',index=False)




