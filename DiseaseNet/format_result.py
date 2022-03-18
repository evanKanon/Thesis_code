

import pandas as pd
import networkx as nx



ddf=pd.read_csv('/mnt/d/PhD/Results/DiseaseNet_v0.1.0_05102020_ldsc/network/result_df.csv')  

ddf['ranges']=pd.cut(ddf['rg'],bins=[-2,-1,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2])


ddf.groupby(['ranges']).agg(['mean','count'])


G = nx.from_pandas_edgelist(ddf,'p1', 'p2', 'rg')
G=nx.from_pandas_dataframe(ddf, 'p1', 'p2', 'rg')




ddf['p1']=ddf['p1'].str.replace('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/complete_input/','')
ddf['p2']=ddf['p2'].str.replace('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/complete_input/','')

ddf['p1']=ddf['p1'].str.replace('.sumstats.gz','')
ddf['p2']=ddf['p2'].str.replace('.sumstats.gz','')



df3 = ddf[ddf['rg'] > 0.7]  


relevant=df3[['p1','p2','rg']]


relevant.to_csv('/mnt/d/mynet.txt',sep='\t',index=False)
