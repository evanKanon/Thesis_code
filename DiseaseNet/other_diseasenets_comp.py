#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------
# diseasenet comparison
#------------------------------
import os
import re
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
#------------------------------




mydnet=pd.read_csv('/mnt/d/PhD/Results/DiseaseNet_v0.1.0_05102020_ldsc/network/result_df.csv')
barrenas=pd.read_csv('/mnt/d/PhD/input/Diseasenet/barrenas/journal.pone.0008090.s002_barrenas.csv',sep=';')
menche=pd.read_csv('/mnt/d/PhD/input/Diseasenet/menche/DataS4_disease_pairs.tsv',sep='\t',comment='#')




mydnet['p1']=mydnet['p1'].str.replace('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/complete_input/','')
mydnet['p2']=mydnet['p2'].str.replace('/mnt/ris-fas1a/linux_groups2/baillie_grp/diseasenet/evan/complete_input/','')

mydnet['p1']=mydnet['p1'].str.replace('.sumstats.gz','')
mydnet['p2']=mydnet['p2'].str.replace('.sumstats.gz','')

mydnet['p1']=mydnet['p1'].str.replace('imputed_','')
mydnet['p2']=mydnet['p2'].str.replace('imputed_','')

mydnet['p1']=mydnet['p1'].str.replace('UKB_20002_','')
mydnet['p2']=mydnet['p2'].str.replace('UKB_20002_','')

mydnet['p1']=mydnet['p1'].str.replace('self_reported_','')
mydnet['p2']=mydnet['p2'].str.replace('self_reported_','')

mydnet['p1']=mydnet['p1'].str.replace('Astle et al 2016','')
mydnet['p2']=mydnet['p2'].str.replace('Astle et al 2016','')


mydnet['p1']=mydnet['p1'].str.replace('DIAGRAM T2D','Type 2 diabetes')
mydnet['p2']=mydnet['p2'].str.replace('DIAGRAM T2D','Type 2 diabetes')


mydnet['p1']=mydnet['p1'].str.replace('Alzheimers',"Alzheimer's")
mydnet['p2']=mydnet['p2'].str.replace('Alzheimers',"Alzheimer's")

mydnet['p1']=mydnet['p1'].str.replace('Alzheimer',"Alzheimer's")
mydnet['p2']=mydnet['p2'].str.replace('Alzheimer',"Alzheimer's")

mydnet['p1']=mydnet['p1'].str.replace('BreastCancer','Breast Cancer')
mydnet['p2']=mydnet['p2'].str.replace('BreastCancer','Breast Cancer')


mydnet['p1']=mydnet['p1'].str.replace('Multiple-Sclerosis', 'Multiple Sclerosis')
mydnet['p2']=mydnet['p2'].str.replace('Multiple-Sclerosis', 'Multiple Sclerosis')


mydnet['p1']=mydnet['p1'].str.replace('_',' ')
mydnet['p2']=mydnet['p2'].str.replace('_',' ')


clean_dnet=pd.DataFrame.dropna(mydnet)
finlist=list(set().union(clean_dnet['p1'],clean_dnet['p1']))

# use only clean mydnet
mydnet=mydnet.dropna()


#------------------------
# format barrenas
#------------------------


mydnet['p1'].unique()
barrenas['Disease/Trait'].unique()

#- create dictionary with genes for each trait in barrenas network
gene_dict={}
for trait in barrenas['Disease/Trait'].unique():
	gene_dict[trait]=[]

#- fill the dictionary with the gene names associated with them
for i in range(len(barrenas)):
	gene_dict[barrenas.loc[i][0]].append(barrenas.loc[i][1])

#- recreate barrenas network
listoflists=[]
for l in range(len(list(gene_dict.keys()))):
	x=max(range(len(list(gene_dict.keys()))))
	while x!=l:
		list1=gene_dict[list(gene_dict.keys())[l]]
		list2=gene_dict[list(gene_dict.keys())[x]]
		common=list(set(list1).intersection(list2))
		if len(common)!=0:
			listoflists.append([list(gene_dict.keys())[l],list(gene_dict.keys())[x],1])
		x+=(-1)


bdf = pd.DataFrame(listoflists, columns = ['From', 'To', 'Connect'])  




#------------------------
# finding common parts of networks
#------------------------


#- find only common diseases accross different networks
common_dict={'barrenas':[],'diseasenet':[]}
lol=[]
for i in barrenas['Disease/Trait'].unique():
	for j in mydnet['p1'].unique():
		if re.search(i, j, re.IGNORECASE):
			common_dict['barrenas'].append(i)
			common_dict['diseasenet'].append(j)
			lol.append([i,j])

common_dict2={'menche':[],'diseasenet':[]}
lol2=[]
for i in menche['disease_A'].unique():
	for j in mydnet['p1'].unique():
		if re.search(i, j, re.IGNORECASE):
			common_dict2['menche'].append(i)
			common_dict2['diseasenet'].append(j)
			lol2.append([i,j])


#! Diabetes melliltus t2/t2; alzheimer; inflamatory bowel disease; parkinson; body weight is common between meche and us but does not show









#- subset networks based on the common elements they have
cdnet1=mydnet[mydnet['p1'].isin(common_dict['diseasenet'])]
cdnet=cdnet1[cdnet1['p2'].isin(common_dict['diseasenet'])]

cbar1=bdf[bdf['From'].isin(common_dict['barrenas'])]
cbar=cbar1[cbar1['To'].isin(common_dict['barrenas'])]

cmen1=menche[menche['disease_A'].isin(common_dict2['menche'])]
cmen=cmen1[cmen1['disease_B'].isin(common_dict2['menche'])]

cdnet2=mydnet[mydnet['p2'].isin(common_dict2['diseasenet'])]
cdneto=cdnet2[cdnet2['p1'].isin(common_dict2['diseasenet'])]


#- clean up !!!
pd.DataFrame.dropna(cdnet)



#- replace names to match similar names
for h in range(len(lol2)):
	cdneto['p1']=cdneto['p1'].str.replace(lol2[h][1],lol2[h][0])
	cdneto['p2']=cdneto['p2'].str.replace(lol2[h][1],lol2[h][0])




#- find common or different links between them
samedfdis=pd.DataFrame()
samedfmen=pd.DataFrame()
for it1 in cdneto.index:
	for it2 in cmen.index:
		if cdneto.loc[it1]['p1']==cmen.loc[it2]['disease_A'] and cdneto.loc[it1]['p2']==cmen.loc[it2]['disease_B']: 
			samedfdis=samedfdis.append(cdneto.loc[it1])
			samedfmen=samedfmen.append(cmen.loc[it2])
		if cdneto.loc[it1]['p2']==cmen.loc[it2]['disease_A'] and cdneto.loc[it1]['p1']==cmen.loc[it2]['disease_B']: 
			samedfdis=samedfdis.append(cdneto.loc[it1])
			samedfmen=samedfmen.append(cmen.loc[it2])


# remove duplicate indexes
samedfmen = samedfmen[~samedfmen.index.duplicated(keep='first')]


# identify which network has which connection uniquely
cdneto[~cdneto.isin(samedfdis)].dropna()
cmen[~cmen.isin(samedfmen)].dropna()




# find links onlu in pne
cd=cdnet.dropna()
cd[cd['rg']>0.3]




##############################

#- save lists for use as appendix tables

myl=list(set().union(mydnet['p1'],mydnet['p2']))

ddf=pd.DataFrame()
ddf['DiseaseNet Nodes']=myl

ddf.to_csv('/mnt/d/PhD/Text/images/mydiseases.csv',index=False)


# find number of unique nodes in barrenas network
myl2=list(set().union(bdf['From'],bdf['To']))
# find number of unique nodes in menche network
myl3=list(set().union(menche['disease_A'],menche['disease_B']))


#############################



#------------------------
# Read in EONS networks
#------------------------

def write_ncol_up(myfile,outputfile='Default'):
	'''Function which converts matrix-type file into ncol, only reads half-matrix'''
	#import pandas as pd
	#import os
	data=pd.read_csv(myfile,index_col=0)
	if outputfile=='Default':
		myoutputfile=os.path.split(myfile)[0]+'/'+os.path.split(myfile)[1].split('.')[0]
	else:
		myoutputfile=outputfile
	newfile=open(myoutputfile+'.ncol','w')
	#for pat1 in data.keys():
	#	for pat2 in data.index:
	#		newfile.write('%s\t%s\t%s\n'%(pat1, pat2, data[pat1][pat2]))
	for i in range(len(data.keys())):
		x=max(range(len(data.keys())))
		while x!=i:
			pat1=data.keys()[i]
			pat2=data.keys()[x]
			newfile.write('%s\t%s\t%s\n'%(pat1, pat2, data[pat1][pat2]))
			x+=(-1)
	newfile.close()


write_ncol_up('/mnt/d/PhD/Results/Diseasenet_EONS/Diseasenetcomulative_GWAS_prepped.expression_spearman-9900.csv','/mnt/d/PhD/Results/Diseasenet_EONS/Diseasenetcomulative_GWAS_prepped.expression_spearman-9900_prep')

ddf=pd.read_csv('/mnt/d/PhD/Results/Diseasenet_EONS/Diseasenetcomulative_GWAS_prepped.expression_spearman-9900_prep.ncol',sep='\t',header=None)








#------------------------
# Find common elements accross all networks
#------------------------

#- establish common names
namelist1=barrenas['Disease/Trait'].unique()
namelist2=menche['disease_A'].unique()
namelist3=mydnet['p1'].unique()
namelist4=ddf[0].unique()


#- turn lists into columns and save for transofrmation
ldf1 = pd.DataFrame({'barrenas':namelist1})
ldf2 = pd.DataFrame({'menche':namelist2})
ldf3 = pd.DataFrame({'diseasenet':namelist3})
ldf4 = pd.DataFrame({'EONS':namelist4})

#ldf34=pd.concat([ldf3,ldf4], ignore_index=True, axis=1)

ldf1.to_csv('/mnt/d/PhD/Results/Diseasenet_comparison/barrenas_namelist.csv', index=False)
ldf2.to_csv('/mnt/d/PhD/Results/Diseasenet_comparison/menche_namelist.csv', index=False)
ldf3.to_csv('/mnt/d/PhD/Results/Diseasenet_comparison/diseasenet_namelist.csv', index=False)
ldf4.to_csv('/mnt/d/PhD/Results/Diseasenet_comparison/EONS_namelist.csv', index=False)


#- read in manually altered lists

ndf1=pd.read_csv('/mnt/d/PhD/Results/Diseasenet_comparison/barrenas_namelist_edited.csv',sep=';')
ndf3=pd.read_csv('/mnt/d/PhD/Results/Diseasenet_comparison/diseasenet_namelist_edited.csv',sep=';')
ndf4=pd.read_csv('/mnt/d/PhD/Results/Diseasenet_comparison/EONS_namelist_edited.csv',sep=';')


ndict1={}
for i in range(len(ndf1)):
	ndict1[ndf1['barrenas'][i]]=ndf1['Unnamed: 0'][i]


ndict3={}
for i in range(len(ndf3)):
	ndict3[ndf3['diseasenet'][i]]=ndf3['Unnamed: 0'][i]

ndict4={}
for i in range(len(ndf4)):
	ndict4[ndf4['EONS'][i]]=ndf4['Unnamed: 0'][i]


#- change names of nodes in networks

bdf['From']=bdf['From'].replace(ndict1)
bdf['To']=bdf['To'].replace(ndict1)


menche.replace(',','', regex=True, inplace=True)



mydnet['p1']=mydnet['p1'].replace(ndict3)
mydnet['p2']=mydnet['p2'].replace(ndict3)



ddf[0]=ddf[0].replace(ndict4)
ddf[1]=ddf[1].replace(ndict4)



common_dict3={'barrenas':[],'diseasenet':[],'menche':[],'EONS':[]}
lol_running=[]
for i in list(bdf['From'].unique())+list(bdf['To'].unique()):
	for j in list(menche['disease_A'].unique())+list(menche['disease_B'].unique()):
		if re.search(i, j, re.IGNORECASE):
			common_dict3['barrenas'].append(i)
			common_dict3['menche'].append(j)
			lol_running.append([i,j])
			break


lol_2=[]
for x in common_dict3['barrenas']:
	for y in list(mydnet['p1'].unique())+list(mydnet['p2'].unique()):
		if re.search(x,y,re.IGNORECASE):
			common_dict3['diseasenet'].append(y)
			lol_2.append([x,y])
			break
			


lol_3=[]
for z in common_dict3['diseasenet']:
	for w in list(ddf[0].unique())+list(ddf[1].unique()):
		if re.search(z,w,re.IGNORECASE):
			common_dict3['EONS'].append(w)
			lol_3.append([z,w])
			break




myset = set(common_dict3['EONS']+common_dict3['diseasenet']+common_dict3['menche']+common_dict3['barrenas'])


#- create common network
G=nx.Graph()
G.add_nodes_from(myset)



fndf=pd.DataFrame()
for i in range(len(myset)):
	x=len(myset)-1
	while i!=x:
		fndf=pd.concat([fndf,bdf.loc[(bdf['From']==list(myset)[i]) & (bdf['To']==list(myset)[x])]], ignore_index=True)
		fndf=pd.concat([fndf,bdf.loc[(bdf['To']==list(myset)[i]) & (bdf['From']==list(myset)[x])]], ignore_index=True)		
		#bdf.loc[(bdf['To']==list(myset)[i]) & (bdf['From']==list(myset)[x])]
		x-=1


fndf2=pd.DataFrame()
for i in range(len(myset)):
	x=len(myset)-1
	while i!=x:
		fndf2=pd.concat([fndf2,mydnet.loc[(mydnet['p1']==list(myset)[i]) & (mydnet['p2']==list(myset)[x])]], ignore_index=True)
		fndf2=pd.concat([fndf2,mydnet.loc[(mydnet['p2']==list(myset)[i]) & (mydnet['p1']==list(myset)[x])]], ignore_index=True)		
		#bdf.loc[(bdf['To']==list(myset)[i]) & (bdf['From']==list(myset)[x])]
		x-=1


fndf3=pd.DataFrame()
for i in range(len(myset)):
	x=len(myset)-1
	while i!=x:
		fndf3=pd.concat([fndf3,ddf.loc[(ddf[0]==list(myset)[i]) & (ddf[1]==list(myset)[x])]], ignore_index=True)
		fndf3=pd.concat([fndf3,ddf.loc[(ddf[1]==list(myset)[i]) & (ddf[0]==list(myset)[x])]], ignore_index=True)		
		#bdf.loc[(bdf['To']==list(myset)[i]) & (bdf['From']==list(myset)[x])]
		x-=1


#- delete duplicate connections which arise from having multiple studies of the same disease as unique elements
#fndf3=fndf3.drop_duplicates(subset=[0, 1], keep='first')
#fndf3.reset_index(inplace=True)

#- delete duplicate connections which arise from having multiple studies of the same disease as unique elements
#- use list of unique elements present to find duplicate lines
lines_to_delete=[]
theset=list(myset)
for i in range(len(theset)):
	w=len(theset)-1
	while i!=w:
		item1=theset[i]
		item2=theset[w]
		a=0
		myline=[]
		for line in range(len(fndf3)):
			if fndf3[0][line]==item1 and fndf3[1][line]==item2:
				a+=1
				myline.append(line)
			elif fndf3[0][line]==item2 and fndf3[1][line]==item1:
				a+=1
				myline.append(line)
		if a>1:
			lines_to_delete.append(myline[1])
		w-=1

#- actually delete duplicated
fndf3=fndf3.drop(fndf3.index[lines_to_delete])
fndf3.reset_index(inplace=True)


fndf4=pd.DataFrame()
for i in range(len(myset)):
	x=len(myset)-1
	while i!=x:
		fndf4=pd.concat([fndf4,menche.loc[(menche['disease_A']==list(myset)[i]) & (menche['disease_B']==list(myset)[x])]], ignore_index=True)
		fndf4=pd.concat([fndf4,menche.loc[(menche['disease_B']==list(myset)[i]) & (menche['disease_A']==list(myset)[x])]], ignore_index=True)		
		#bdf.loc[(bdf['To']==list(myset)[i]) & (bdf['From']==list(myset)[x])]
		x-=1



for i in range(len(fndf)):
	if G.has_edge(fndf['From'][i], fndf['To'][i]):
		w=G.edges[fndf['From'][i], fndf['To'][i]]['weight']
		G.edges[fndf['From'][i], fndf['To'][i]]['weight']=w+1
	else:
		G.add_edge(fndf['From'][i], fndf['To'][i], weight=1 )


for i in range(len(fndf2)):
	if G.has_edge(fndf2['p1'][i], fndf2['p2'][i]):
		w=G.edges[fndf2['p1'][i], fndf2['p2'][i]]['weight']
		G.edges[fndf2['p1'][i], fndf2['p2'][i]]['weight']=w+1
	else:
		G.add_edge(fndf2['p1'][i], fndf2['p2'][i], weight=1 )


for i in range(len(fndf3)):
	if G.has_edge(fndf3[0][i], fndf3[1][i]):
		w=G.edges[fndf3[0][i], fndf3[1][i]]['weight']
		G.edges[fndf3[0][i], fndf3[1][i]]['weight']=w+1
	else:
		G.add_edge(fndf3[0][i], fndf3[1][i], weight=1 )


for i in range(len(fndf4)):
	if G.has_edge(fndf4['disease_A'][i], fndf4['disease_B'][i]):
		w=G.edges[fndf4['disease_A'][i], fndf4['disease_B'][i]]['weight']
		G.edges[fndf4['disease_A'][i], fndf4['disease_B'][i]]['weight']=w+1
	else:
		G.add_edge(fndf4['disease_A'][i], fndf4['disease_B'][i], weight=1 )



elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] == 4]
esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] == 1]
esmali = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] == 2]
elargi = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] == 3]

pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility
nx.draw_networkx_nodes(G, pos, node_size=700)
nx.draw_networkx_edges(G, pos, edgelist=elarge, width=6, alpha=0.5, edge_color="r")
nx.draw_networkx_edges(G, pos, edgelist=elargi, width=6, alpha=0.5, edge_color="orange")
nx.draw_networkx_edges(G, pos, edgelist=esmali, width=6, alpha=0.5, edge_color="c")
nx.draw_networkx_edges(G, pos, edgelist=esmall, width=6, alpha=0.5, edge_color="b")


nx.draw_networkx_labels(G, pos, font_size=9, font_family="sans-serif")

ax = plt.gca()
ax.margins(0.08)
plt.axis("off")
plt.tight_layout()

l,r = plt.xlim()
plt.xlim(l-0.45,r+0.45)

plt.savefig('/mnt/d/PhD/Text/images/disease_comp3.png')

plt.gcf().clear()




for (u, v, wt) in G.edges.data('weight'):
    print(u,v,wt)