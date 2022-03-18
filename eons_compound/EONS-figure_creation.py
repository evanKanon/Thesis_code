#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#--------------------------

#inputdir='/mnt/d/PhD/Results/GDS5826_22052020_EONSv0.1.7_default_pca_GDS5826_cl/'
#plotttitle='Cancer Data Set - Cell Line Annotation'
#inputdir='/mnt/d/PhD/Results/GDS5826_22052020_EONSv0.1.7_default_pca_GDS5826_ct/'
#plotttitle='Cancer Data Set - Cell Treatment Annotation'
#inputdir='/mnt/d/PhD/Results/GDS5826_22052020_EONSv0.1.7_default_pca_GDS5826_co/'
#plotttitle='Cancer Data Set - Original Cell Line Annotation'
#inputdir='/mnt/d/PhD/Results/GDS6063_22052020_EONSv0.1.7_default_pca_GDS6063_individual/'
#plotttitle='Influenza Data Set - Donor Annotation'
#inputdir='/mnt/d/PhD/Results/GDS6063_22052020_EONSv0.1.7_default_pca_GDS6063_infection/'
#plotttitle='Influenza Data Set - Disease State Annotation'
inputdir='/mnt/d/PhD/Results/GDS3085_full_22052020_EONSv0.1.7_extensive_pca_GDS3085_full/'
plotttitle='Sepsis Patients Data Set - Extensive Run'
#inputdir='/mnt/d/PhD/Results/GDS3085_full_22052020_EONSv0.1.7_default_pca_GDS3085_full/'
#plotttitle='Sepsis Patients Data Set - Default Run'


o_evaldir=inputdir+'evaluation_'+'_'.join(inputdir.split('/')[-2:-1][0].split('.')[-1:][0].split('_')[3:])
of_interest=list(pd.read_csv(inputdir+'results/best_networks.txt',sep=' ',header=None,index_col=0)[1])
ann=pd.read_csv(o_evaldir+'/ANN_evaluation_scores.csv')
nda=pd.read_csv(o_evaldir+'/NDA_avg_evaluation_scores.csv')
igp=pd.read_csv(o_evaldir+'/igp_average_scores.csv')




numbs=ann[ann.keys()[1]].to_list()
attrs=['ANN'] * len(ann[ann.keys()[1]])
typ=['Predictive accuracy'] * len(ann[ann.keys()[1]])

#numbs=numbs+nda[nda.keys()[1]].to_list()
#attrs=attrs+['NDA']*len(nda[nda.keys()[1]])
#typ=typ+['Network structure'] * len(nda[nda.keys()[1]])

numbs=numbs+igp[igp.keys()[1]].to_list()
attrs=attrs+['FNN']*len(igp[igp.keys()[1]])
typ=typ+['Predictive accuracy'] * len(igp[igp.keys()[1]])

impot=[]
for x in ann[ann.keys()[0]]:
	if x in of_interest:
		impot.append('EONS-identified')
	else:
		impot.append('non-identified')

#for x in nda[nda.keys()[0]]:
#	if x in of_interest:
#		impot.append('EONS-identified')
#	else:
#		impot.append('non-identified')

for x in igp[igp.keys()[0]]:
	if x in of_interest:
		impot.append('EONS-identified')
	else:
		impot.append('non-identified')

#! ADD SANITY CHECK

#myfdf = pd.DataFrame({'Correct prediction' : numbs,'Evaluation metrics' : attrs,'Selection':impot, 'Metric type':typ})
myfdf = pd.DataFrame({'Score' : numbs,'Evaluation metrics' : attrs,'Selection':impot, 'Metric type':typ})


#sns.catplot(x="day", y="total_bill", data=tips);
#sns.catplot(x="day", y="total_bill", hue="sex", kind="swarm", data=tips) # for hue add third column with info
#sns.catplot(x="cats", y="data", hue="imp", kind='swarm', data=myfdf) 
'''
sns.catplot(x="Evaluation metrics", y="Correct prediction", hue="Selection", kind='swarm', data=myfdf) 
sns.catplot(x="Evaluation metrics", y="Score", hue="Selection", kind='swarm', col='Metric type', data=myfdf) 


g = sns.FacetGrid(myfdf, hue="Selection", col="Metric type", height=4)
#g.map(catplot, "Evaluation metrics", "Correct prediction")
g.map(sns.catplot(x="Evaluation metrics", y="Correct prediction", hue="Selection", kind='swarm', data=myfdf))
g.add_legend()

fig, axs = plt.subplots(ncols=2)
sns.catplot(x='Evaluation metrics', y='Score', hue="Selection", kind='swarm', data=myfdf, ax=axs[0])
sns.catplot(x='Evaluation metrics', y='Score', hue="Selection", kind='swarm', data=myfdf2, ax=axs[1])
'''
#- seperate NDA score plot (to integrate with subplotting)

numbs2=nda[nda.keys()[1]].to_list()
attrs2=['NDA']*len(nda[nda.keys()[1]])

impot2=[]
for x in nda[nda.keys()[0]]:
	if x in of_interest:
		impot2.append('EONS-identified')
	else:
		impot2.append('non-identified')


#myfdf2 = pd.DataFrame({'data' : numbs2,'cats' : attrs2,'imp':impot2})
myfdf2 = pd.DataFrame({'Score' : numbs2,'Evaluation metrics' : attrs2,'Selection':impot2})

'''
#===============
# plot with jitter
#===============

#f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
f, (a0, a1) = plt.subplots(1, 2)
a0=sns.stripplot(x='Evaluation metrics', y='Score', hue="Selection", jitter=True, data=myfdf, ax=a0)
a1=sns.stripplot(x='Evaluation metrics', y='Score', hue="Selection", jitter=True, data=myfdf2, ax=a1)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

handles, labels = a0.get_legend_handles_labels()
#fig.legend(handles, labels, loc='upper center')
#f.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

#a0.legend(handles[:0], labels[:0])
a0.get_legend().remove()
plt.legend(handles[2:4], labels[2:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('/mnt/d/test.png',bbox_inches='tight')
plt.gcf().clear()

#==============
# two plots
#==============

sns.catplot(x="Evaluation metrics", y="Score", hue="Selection", kind='swarm', data=myfdf) 
plt.savefig('/mnt/d/plot1.png')
plt.gcf().clear()

sns.catplot(x="Evaluation metrics", y="Score", hue="Selection", kind='swarm', data=myfdf2)
plt.savefig('/mnt/d/plot2.png')
plt.gcf().clear()

'''
#================
# swarm plots
#================

#! SET CONSTANT SCALE FOR NDA ACCROSS DIFFERENT FIGURES

f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
a0=sns.swarmplot(x='Evaluation metrics', y='Score', hue="Selection", data=myfdf, ax=a0, linewidth=1)
a1=sns.swarmplot(x='Evaluation metrics', y='Score', hue="Selection", data=myfdf2, ax=a1, linewidth=1)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

handles, labels = a0.get_legend_handles_labels()
#fig.legend(handles, labels, loc='upper center')
#f.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

#a0.legend(handles[:0], labels[:0])
a0.get_legend().remove()
a1.set_ylabel('')
a0.set_ylim([-0.05, 1.05])
a1.set_ylim([-1, 8])
f.suptitle(plotttitle, fontsize=16)
plt.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
name='_'.join(plotttitle.split(' '))
plt.savefig('/mnt/d/PhD/Text/images/'+name+'.png',bbox_inches='tight')
plt.gcf().clear()


#palette=['#91bfdb','#fc8d59']
#edgecolor='black'
