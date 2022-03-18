#!/opt/local/bin/python
# -*- coding: UTF-8 -*-

import pandas as pd 
from slacker import Slacker 



slack = Slacker('xoxp-56103715248-56114624145-154905509446-72c754c99fabcd5e9e1772648f3ee7c4')

#- import dataset
fantom=pd.read_table('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/nozeros_phase2_nopooled.expression')
#- create correlation matrix using spearman correlation
corr_dat=fantom.corr(method='spearman')
#- write result into csv file
corr_dat.to_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/fantom5_sample_corrmatrix.csv')
#- notify me that the creation has been completed
slack.chat.post_message('@evan','Fantom5 corr matrix creation complete')
