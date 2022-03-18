#- Script outlining commands used for data prep of GSE59114 for CaSTLe comparison


#- import required packages
import pandas as pd


#- read in relevant fiiles
dat1=pd.read_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/Castle_data/HSCs/GSE59114_C57BL6_GEO_all.csv',header=1)
dat2=pd.read_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/Castle_data/HSCs/GSE59114_DBA_GEO_all.csv',header=1)


#- drop useless column
dat1=dat1.drop(columns=['UCSC transcripts'])

#- remove quotation marks from gene names
for i in range(len(dat1)):
	dat1[dat1.keys()[0]][i]=dat1[dat1.keys()[0]][i].replace("'",'')



#- write output
dat1.to_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/Castle_data/HSCs/GSE59114_C57BL6_stdt.csv')
dat2.to_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/Castle_data/HSCs/GSE59114_DBA_stdt.csv')











