#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
import pandas as pd
#------------------------------

#- read file in
inputfile='/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/Castle_data/Embryo/deng.csv'
data=pd.read_csv(inputfile, index_col=0)

#- format it appropriately
data_t=data.transpose()
data_t.reset_index(inplace=True)
data_t.to_csv(inputfile.split('.')[0]+'_stdt.'+inputfile.split('.')[1])

