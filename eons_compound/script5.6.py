#---------------------------------------------------------
# version 4.0.0.py is the normal analysis (1000 iter, 10 nodes)
# version 4.0.1.py is the normal analysis with 20 nodes
# version 4.0.py accepts automated input but cannot be nohuped
# version 5.0.py is the most automated version to date
# version 5.1.py has different namings of files and corrected figures
# version 5.2.py saves distribution images and network snapshots
# version 5.3.py first creates 0 , 1 networks and then adds noise
# version 5.4.py substitutes WHD for EDES and saves network edges
# version 5.5.py creates histograms of distributions of HD, WHD and random values
# version 5.6.py substitutes EDES for SEWD, improves options settings and runs on python3
#---------------------------------------------------------
p1 = '/nfs_netapp/s1242130/Results/' #server path
p2 = '/Users/s1242130/Desktop/Test/' #mac path
p3 = 'C:\\Users\\User\\Desktop\\Test' #pc path
p4 = '' #website path
p5 = '/nfs_netapp/s1242130/Results/Fig/'

import numpy as np
import networkx as nx
import random, sys
import matplotlib.pyplot as plt
import scipy.stats
import time
import datetime
import argparse
import math

#print the date for log purposes
print(datetime.datetime.today())
print()

#NEW HAMMING DISTANCE FUNCTION WITH THRESHOLD INSIDE
def newhamming(N1, N2, thr):
	list1=[]
	list2=[]
	for (f,t,d) in N1:
		b=(f,t,d)
		if d<=thr: b=b[:2]+(0,)
		list1.append(b)
	for (f,t,d) in N2:
		b=(f,t,d)
		if d<=thr: b=b[:2]+(0,)
		list2.append(b)
	df=0
	for (f,t,d) in list1:
		for (f2,t2,d2) in list2:
			if f==f2 and t==t2 and d!=0 and d2==0 : df+=1
			if f==f2 and t==t2 and d==0 and d2!=0 : df+=1
	return df

#WEIGHT HAMMING
def weighthamming(N1,N2):
	mylist=[]
	for (f,t,d) in N1:
		for (f2,t2,d2) in N2:
			if f==f2 and t==t2: mylist.append(abs(d-d2))
	return sum(mylist)


#PROGRAM FUNCTIONALITY
parser=argparse.ArgumentParser(description='Creation and comparison of complete networks')
parser.add_argument('-n','--nodes',type=int, help="number of network's nodes")
parser.add_argument('-i','--iterations', default=1000, type=int, help='number of iterations')
parser.add_argument('-s', '--save', type=str, help='file path for results and figures to be saved')
parser.add_argument('-hist', '--histograms', default=False, help='creation of histograms of distributions', action='store_true')
parser.add_argument('-snap','--snapshots', default=False, help='save snapshots of networks and edge weights', action='store_true')
parser.add_argument('-cr','--crazy', default=False, help='save all networks created by the analysis', action='store_true')
parser.add_argument('-d','--data', default=False, help='save data of metrics (HD and WHD) from the analysis', action='store_true')

# add verbose -v paramater which will show snapshots of the networks actualy used not just top ones

args=parser.parse_args()



# MAKE ONE FIGURE WITH NETWORKS, FINAL FIGURES ETC IN MANY PANELS


nodes=args.nodes
hmi=args.iterations
save=args.save

hist=args.histograms
snap=args.snapshots
crazy=args.crazy
data=args.data

#this print put is for log purposes
print('Nodes: '+str(nodes)+' Iterations: '+str(hmi))
print

#Creation of networks
print('COMMENCING ANALYSIS')
print('Creating Networks...')
start_time1=time.time()
#list1=[[] for _ in range(11)]
#list2=[[] for _ in range(11)]

standdev=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

#Creation of dictionaries to store networks
dict1={}
dict2={}
for i in standdev:
	dict1[i]={}
for i in standdev:
	for n in range(hmi):
		dict1[i][n]={}
for i in standdev:
	dict2[i]={}
for i in standdev:
	for n in range(hmi):
		dict2[i][n]={}

#Creation of the original networks (no noise)
N1=nx.complete_graph(nodes)
N0=nx.complete_graph(nodes)
for (f,t,d) in N1.edges(data='weight'):
	N1[f][t]['weight']=1
for (f,t,d) in N0.edges(data='weight'):
	N0[f][t]['weight']=0

for stdv in standdev:
	iteration=0
	while iteration!=hmi:
		x1=np.random.normal(0, stdv, 1000)
		#x2=np.random.normal(1, stdv, 1000)
		norm1=[]
		for v in x1:
			if v>0 and v<1:
				norm1.append(v)
		#norm2=[]
		#for v in x2:
		#	if v<1 and v>0:
		#		norm2.append(v)
		T1=N0
		T2=N1
		for (f,t,d) in T1.edges(data='weight'):
			T1[f][t]['weight']=norm1[random.randrange(len(norm1))]
		for (f,t,d) in T2.edges(data='weight'):
			T2[f][t]['weight']=1-norm1[random.randrange(len(norm1))]
		dict1[stdv][iteration]=T1.edges(data='weight')
		dict2[stdv][iteration]=T2.edges(data='weight')
		if crazy==True:
			nx.write_edgelist(T1, 'Net_edges_T1_noise_'+str(stdv)+'_iter_'+str(iteration))
			nx.write_edgelist(T2, 'Net_edges_T2_noise_'+str(stdv)+'_iter_'+str(iteration))
		iteration+=1
#saving snapshots of networks and edge weight histograms
	if iteration==hmi:
		if snap==True:
			h1=sorted(norm1, reverse=True)
			nor=[]
			for i in norm1:
				nor.append(1-i)
			h2=sorted(nor)
			#plt.plot(h1)
			#plt.plot(h1)
			plt.hist(h1, histtype='step')
			plt.hist(h2, histtype='step')
			plt.savefig('Distributions_'+str(stdv)+'.png')
			plt.gcf().clear()
			pos=nx.spring_layout(T1, k=5/math.sqrt(T1.order()))   #k could be revised
			edge_list_t1=dict([((u,v,),d['weight'])
				for u,v,d in T1.edges(data=True)])
			edge_list2_t1={}
			keys=edge_list_t1.keys()
			for i in keys:
				edge_list2_t1[i]=round(edge_list_t1[i], 3)
			nx.draw(T1, pos)
			nx.draw_networkx_edge_labels(T1, pos, edge_labels=edge_list2_t1)
			plt.savefig('Network1_'+str(stdv)+'_'+str(hmi)+'.png')
			plt.gcf().clear()
			pos2=nx.spring_layout(T2, k=5/math.sqrt(T2.order()))   #k could be revised
			edge_list_t2=dict([((u,v,),d['weight'])
				for u,v,d in T2.edges(data=True)])
			edge_list2_t2={}
			keys=edge_list_t2.keys()
			for i in keys:
				edge_list2_t2[i]=round(edge_list_t2[i], 3)		
			nx.draw(T2, pos2)
			nx.draw_networkx_edge_labels(T2, pos2, edge_labels=edge_list2_t2)
			plt.savefig('Network2_'+str(stdv)+'_'+str(hmi)+'.png')
			plt.gcf().clear()


#Creation of completely random networks
print('...Creating Random Networks...')
dict3={}
dict4={}
for i in range(hmi):
	dict3[i]={}
for i in range(hmi):
	dict4[i]={}
for iteration in range(hmi):
	x1=np.random.normal(0, 20, 1000)
	x2=np.random.normal(1, 20, 1000)
	norm1=[]
	for v in x1:
		if v>0 and v<1:
			norm1.append(v)
	norm2=[]
	for v in x2:
		if v<1 and v>0:
			norm2.append(v)
	T1=nx.complete_graph(nodes)
	T2=nx.complete_graph(nodes)
	for (f,t,d) in T1.edges(data='weight'):
		T1[f][t]['weight']=norm1[random.randrange(len(norm1))]
	for (f,t,d) in T2.edges(data='weight'):
		T2[f][t]['weight']=norm2[random.randrange(len(norm2))]
	dict3[iteration]=T1.edges(data='weight')
	dict4[iteration]=T2.edges(data='weight')

print('......Networks Created......')
print('--- %s seconds ---' % (time.time()-start_time1))


#NETWORK ANALYSIS
print('Initialising Network Analysis...')
start_time2=time.time()

list_hamm=[[[],[],[],[],[],[],[],[],[],[],[]] for _ in range(11)]
for a in standdev:
	for b in range(hmi):
		for c in range(11):
			list_hamm[int(a*10)][c].append(newhamming(dict1[a][b], dict2[a][b], (float(c)/10)))
del list_hamm[0]

list_weighthamm=[[] for _ in range(11)]
for a in standdev:
	for b in range(0,hmi):
		list_weighthamm[int(a*10)].append(weighthamming(dict1[a][b],dict2[a][b]))
del list_weighthamm[0]

#---------------------------------------------------------------------
randomHD=[[] for _ in range(11)]
for c in range(11):
	for b in range(0,hmi):
		randomHD[c].append(newhamming(dict3[b],dict4[b],(float(c)/10)))

randomWHD=[]
for b in range(0,hmi):
	randomWHD.append(weighthamming(dict3[b],dict4[b]))

mean_WHD=[]
for a in range(0,10):
	mean_WHD.append(np.mean(list_weighthamm[a]))
mean_HD=[[] for _ in range(10)]
for a in range(0,10):
	for c in range(11):
		mean_HD[a].append(np.mean(list_hamm[a][c]))

#randmean_WHD=np.mean(randomWHD)
#randmean_HD=np.mean(randomHD)

#calculates number of edges in network, for edge weights in (0,1) equiv to total possible distance for both metrics
totdist=(nodes*(nodes-1))/2

#normalised distances
normal_WHD=[]
for a in range(10):
	normal_WHD.append(mean_WHD[a]/totdist)
normal_HD=[[] for _ in range(10)]
for a in range(10):
	for c in range(11):
		normal_HD[a].append(mean_HD[a][c]/totdist)

if data==True:
	with open(str(save)+"HD_Results_("+str(nodes)+'n)_('+str(hmi)+'iter).txt', "w") as f:
		for s in list_hamm:
			f.write(str(s) +"\n")
	with open(str(save)+"WHD_Results_("+str(nodes)+'n)_('+str(hmi)+'iter).txt', "w") as f:
		for s in list_weighthamm:
			f.write(str(s) +"\n")
	with open(str(save)+"Random_HD_Results_("+str(nodes)+'n)_('+str(hmi)+'iter).txt', "w") as f:
		for s in randomHD:
			f.write(str(s) +"\n")
	with open(str(save)+"Random_WHD_Results_("+str(nodes)+'n)_('+str(hmi)+'iter).txt', "w") as f:
		for s in randomWHD:
			f.write(str(s) +"\n")

'''
#to open a file saved as above use (this will read strings as integers, may need to convert to floats)
with open("file.txt", "r") as f:
  for line in f:
    score.append(int(line.strip()))
'''

#randomHD_norm=randmean_HD/totdist
#randomWHD_norm=randmean_WHD/totdist

print('...Network Analysis Completed...')
print('--- %s seconds ---' % (time.time() - start_time2))


#PLOT
print('Creating Comparison Plot...')
print('...Calculating IQR...')


def iqr_up(input):
	value=np.subtract(*np.percentile(input, [75,50]))
	return value
def iqr_down(input):
	value=np.subtract(*np.percentile(input, [50,25]))
	return value

iqrHD_up=[[] for _ in range(10)]
for a in range(10):
	for b in range(11):
		iqrHD_up[a].append(iqr_up(list_hamm[a][b]))

iqrHD_down=[[] for _ in range(10)]
for a in range(10):
	for b in range(11):
		iqrHD_down[a].append(iqr_down(list_hamm[a][b]))

niqrHD_up=[[] for _ in range(10)]
for a in range(10):
	for b in range(11):
		niqrHD_up[a].append(iqrHD_up[a][b]/totdist)

niqrHD_down=[[] for _ in range(10)]
for a in range(10):
	for b in range(11):
		niqrHD_down[a].append(iqrHD_down[a][b]/totdist)

iqrWHD_up=[]
for a in range(10):
	iqrWHD_up.append(iqr_up(list_weighthamm[a]))

iqrWHD_down=[]
for a in range(10):
	iqrWHD_down.append(iqr_down(list_weighthamm[a]))

niqrWHD_up=[]
for a in range(10):
	niqrWHD_up.append(iqrWHD_up[a]/totdist)

niqrWHD_down=[]
for a in range(10):
	niqrWHD_down.append(iqrWHD_down[a]/totdist)

print('......IQR acquired......')

#FINAL PLOTS
from textwrap import wrap
listx=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

#FIRST PLOT (NORMALISED)
title1='Comparison Between Hamming Distance (HD) and Sum Edge Weight Distance (SEWD) [Thresholds 0.0-0.4] ('+str(nodes)+'Nodes) ('+str(hmi)+'Iterations)'
title2='Comparison Between Hamming Distance (HD) and Sum Edge Weight Distance (SEWD) [Thresholds 0.5-1.0] ('+str(nodes)+'Nodes) ('+str(hmi)+'Iterations)'
trialist=np.asarray(niqrHD_down).T.tolist()
trialist2=np.asarray(niqrHD_up).T.tolist()
trialist1=np.asarray(normal_HD).T.tolist()
for a in range(5):
	plt.errorbar(listx, trialist1[a], yerr=[trialist[a], trialist2[a]], label='HD_'+str(float(a)/10))
plt.errorbar(listx, normal_WHD, yerr=[niqrWHD_down, niqrWHD_up], linewidth=5, label='SEWD', color='k')
#plt.ylim(0,1)
plt.xlabel('Value of SD')
plt.ylabel('% of difference found over maximum possible')
plt.title('\n'.join(wrap(title1,60)))
lgd=plt.legend(bbox_to_anchor=(1,1), loc=2, fontsize=10)
plt.savefig(str(save)+'A!=B_Norm_Comparison(0-4)_('+str(nodes)+'n)_('+str(hmi)+'iter).png',bbox_extra_artists=(lgd,),bbox_inches='tight')
plt.gcf().clear()

#SECOND PLOT (NORMALISED)
for a in range(5,11):
	plt.errorbar(listx, trialist1[a], yerr=[trialist[a], trialist2[a]], label='HD_'+str(float(a)/10))

plt.errorbar(listx, normal_WHD, yerr=[niqrWHD_down, niqrWHD_up], linewidth=5, label='SEWD', color='k')
plt.xlabel('Value of SD')
plt.ylabel('% of difference found over maximum possible')
plt.title('\n'.join(wrap(title2,60)))
lgd=plt.legend(bbox_to_anchor=(1,1), loc=2, fontsize=10)
plt.savefig(str(save)+'A!=B_Norm_Comparison(5-10)_('+str(nodes)+'n)_('+str(hmi)+'iter).png',bbox_extra_artists=(lgd,),bbox_inches='tight')


print('.........Comparison Plot Created.........')

#Mann-Whitney U-test
print('Performing Mann-Whitney U-test...')

resultsWHD=[]
i=1
for a in range(10):
	local=[]
	local2=[]
	local.append(str(scipy.stats.mannwhitneyu(list_weighthamm[a],randomWHD)))
	local2.append(local[0].split())
	resultsWHD.append(float(local2[0][1][7:-1]))
	local=[]
	local2=[]
	i+=1

resultsHD=[[] for _ in range(12)]
resultsHD[11].append('thr')
for thre in range(10):
	resultsHD[11].append('SD_'+str(float(thre+1)/10))
for why in range(11):
	resultsHD[why].append('thr_'+str(float(why)/10))

for a in range(10):
	for c in range(11):
		try:
			local=[]
			local2=[]
			local.append(str(scipy.stats.mannwhitneyu(list_hamm[a][c],randomHD[c])))
			local2.append(local[0].split())
			resultsHD[c].append(float(local2[0][1][7:-1]))
			local=[]
			local2=[]		
		except:
			resultsHD[c].append('NA')

#sresultsHD=np.asarray(resultsHD).T.tolist()

with open(str(save)+'pvalHD_Results_('+str(nodes)+'n)_('+str(hmi)+'iter).txt', 'w') as f:
	for s in resultsHD:
		f.write(str(s)+'\n')

with open(str(save)+'pvalWHD_Results_('+str(nodes)+'n)_('+str(hmi)+'iter).txt', 'w') as f:
	for s in resultsWHD:
		f.write(str(s)+'\n')

#Summary results information
cresultsHD=[[] for _ in range(11)]
for a in range(0,10):
	for c in range(11):
		try:
			cresultsHD[c].append('Comparison_of_HD_between_SD_'+str(float(a+1)/10)+'_thr_'+str(float(c)/10)+'_and_rand_thr_'+str(float(c)/10)+'_: '+str(scipy.stats.mannwhitneyu(list_hamm[a][c],randomHD[c])))
		except:
			cresultsHD[c].append('Comparison_of_HD_between_SD_'+str(float(a+1)/10)+'_thr_'+str(float(c)/10)+'_and_rand_thr_'+str(float(c)/10)+'_: '+'Not_operable')

cresultsWHD=[]
i=1
for a in range(10):
	cresultsWHD.append('Comparison_of_WHD_between_SD_'+str(float(i)/10)+'_and_rand:\n'+str(scipy.stats.mannwhitneyu(list_weighthamm[a],randomWHD)))
	i+=1

with open(str(save)+'compWHD_Results_('+str(nodes)+'n)_('+str(hmi)+'iter).txt', 'w') as f:
	for s in cresultsWHD:
		f.write(str(s)+'\n')

with open(str(save)+'compHD_Results_('+str(nodes)+'n)_('+str(hmi)+'iter).txt', 'w') as f:
	for s in cresultsHD:
		f.write(str(s)+'\n')

print('...Mann-Whitney U-test performed...')


#_____!!!! log plot x-axis should be changed to (0.1-1) instead of (0-9)

#-log(pval) plot
print('......Creating -log Plot......')
title1='-ln(pvalue) for Hamming Distance (HD) and Sum Edge Weight Distance (SEWD) ('+str(nodes)+' Nodes) ('+str(hmi)+' Iterations)'
plt.gcf().clear()
#plt.xticks(x, listx)
ax=[0,1,2,3,4,5,6,7,8,9]
col=['b','g','r','c','m','y','b','orange','brown','lime']
for a in range(11):
	try:
		plt.plot(ax,-np.log(resultsHD[a][1:len(resultsHD[a])]), label='pval_HD_'+str(resultsHD[a][0]), color=col[a])
		plt.xticks(ax,listx)
	except:
		plt.plot([0,0], linewidth='3', label='HD_'+str(resultsHD[a][0]))
plt.plot(-np.log(resultsWHD), linewidth='5', label='pval_SEWD', color='k')
plt.xlabel('Value of SD')
plt.ylabel('-ln(pvalue)')
plt.title('\n'.join(wrap(title1,60)))
#will need to make legend a bit smaller
lgd=plt.legend(bbox_to_anchor=(1,1), loc=2, fontsize=10)
plt.savefig(str(save)+'-log(pvalue)_plot_('+str(nodes)+'n)_('+str(hmi)+'iter).png', bbox_extra_artists=(lgd,),bbox_inches='tight')
plt.gcf().clear()

print('.........-log Plot Created.........')

#Histograms of distributions
if hist==True:
	print('Creating Histograms...')
	for i in range(10):
		for j in range(11):
			title_hist='Distribution form stddev_'+str(float(i)/10)+' of HD_'+str(float(j)/10)+' vs RandomHD_'+str(float(j/10))+' vs WHD_'+str(float(i/10))+' vs RandomWHD'
			plt.hist(list_hamm[i][j],label='HD_'+str(float(j)/10))
			plt.hist(randomHD[j], label='randomHD_'+str(float(j)/10))
			plt.hist(list_weighthamm[i], label='WHD_'+str(float(i)/10))
			plt.hist(randomWHD, label='randomWHD')
			lgd=plt.legend(bbox_to_anchor=(1,1),loc=2,fontsize=10)
			plt.title('\n'.join(wrap(title_hist,60)))
			plt.savefig(str(save)+'hist_plot_stdev_'+str(float(i)/10)+'_(HD_'+str(float(j)/10)+')_(WHD_'+str(float(i)/10)+').png', bbox_extra_artists=(lgd,),bbox_inches='tight')
			plt.gcf().clear()
	print('...Histograms Created')


print('ANALYSIS COMPLETED')
print('--- %s seconds ---' % (time.time()-start_time1))
