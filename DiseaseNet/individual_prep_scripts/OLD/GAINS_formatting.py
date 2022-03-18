#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import pandas as pd

inputfile1='/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/GAINs_LD_hub_format.txt'

df=pd.read_csv(inputfile1,sep='\t')

df=df.fillna(0)

df = df.rename(columns={'snpid': 'SNP'})

df.to_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/input/GAINs.txt',sep='\t',index=False)
