#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

#EONS third script 
#Create network-of-networks and cluster it

#------------------------------
#v0.0.0 initial version of the new code, adapted from EONS_03 version 0.2
#v0.0.1 updated in EONS_03 to keep same version as other files (will change version descrpt convention in current version)
#v0.0.2 added config file, simplified EONS_02, redefined outputdirs, made single EONS_funct and deleted unused functions, (version descr now includes changes in any/all scripts)
#v0.0.3 added extreme network selection functionality and the the relevant option for the distance multiplier selection, removed config hardpath and renamed EONS_main to EONS
#v0.0.4 changed MCL feature selection with PCA based on networktools code
#v0.0.5 changed feature selection to only principal components, deleted unused code
#v0.0.6 merged MCL feature selection and PCA, created multiple feature subsets through different inflation values, made it possible to notify many users simulataneously, added option to run on normalised data as well, automated static deletion of IDENTIFIER column of data
#v0.0.7 created prep script, changed file names from capital to lowercase, changed networks so that all networks are similarity and not distance networks
#v0.0.8 gave two options for feature selection, either MCL or PCA-based feature selection, added classfile option and config path as well as working and useful eval script
#v0.1.0 BASELINE EONS added failed network functionality to eons_02 making all distances runable, corrected multiplier functionality in eons_03, made eons_eval save proper output, used similarity and not distance networks in eons_funct i.e. 0 is least similar
#v0.1.1 changed make_expression function to accomodate file names other than .soft, added NA dealing option to eons_prep, improved config option for eons_main, made code python 3 compatible
#v0.1.2 established perpendicular point plots
#v0.1.3 improved folder names, added .layout visualisation, updated environment file
#v0.1.4 made running the MCL component faster and allowed user to select inflation values
#v0.1.5 changed naming of folders and certain default values
#v0.1.6 corrected PCA random selection in eons_01, corrected network names in eons_02 and deleted unused inputfile option in eons_03 with the appropriate corrections in eons_funct and _main
#v0.1.7 added correct evaluation measure

#------------------------------
__version__ = '0.1.7'
#------------------------------
import argparse
from pandas import DataFrame
from eons_funct import *
#------------------------------

#- Options
parser=argparse.ArgumentParser(description='EONS-Part_III__Network-of-Networks_Analysis')
parser.add_argument('-s', '--save', type=str, help='file path for results and figures to be saved')
parser.add_argument('-abs','--absolute', default=0, type=float, help='use absolute value for network retention i.e. if cluster contains more than 10 networks keep it')
parser.add_argument('-rel','--relative', default=0, type=float, help='use relative value for network retention i.e. if cluster contains more than 10 per cent of all networks keep it')
parser.add_argument('-thr', '--thresholded', type=int, default=1, choices=[0,1], help='Choose whether to run with thresholded or unthresholded netowrks: 0- run with thresholded and use HD 1- run with unthresholded and use SEWD')
parser.add_argument('-mult','--multiplier', nargs='*', type=float, help='multiplier value for the distance threshold of the extreme network selection algorithm')
args=parser.parse_args()

#- Processing arguments 
outputdir=args.save
absolute=args.absolute
relative=args.relative
thresh=args.thresholded

#- Set extreme network multiplier/s
m=args.multiplier
if not m:
	m=[1]


#- ensure program is using proper networks directory
if os.path.isdir(outputdir+'networks/')==True:
	networkdir=outputdir+'networks/'
else:
	#! update this to output to slack
	print('ERROR: No networks folder found')


#- Read-in input
networks={}
for filename in os.listdir(networkdir):
	s=os.path.abspath(networkdir)
	if filename.endswith('.csv'):
		networks[filename.split('.csv')[0]]=pd.read_csv(s+'/'+filename,index_col=0)


#- create results directory
resultsdir = os.path.join(outputdir, 'results/')
check_dir(resultsdir)


#----------------------------
# Create network-of-networks
#----------------------------

#! unsure whether histogram and hammingnet json work correctly, evaluate !!!!

hammingnet=[] #write this as a list and plot it in Miru

#- create network nodes
nodes = list(networks.keys())
meta = Graph()
meta.add_vertices(range(len(nodes)))
meta.vs["name"]=[nodes[x] for x in range(len(nodes))]

if thresh==0:
	#- compute hamming distance for all networks (create network edges)
	for i in range(len(networks.keys())):
		x=max(range(len(networks.keys())))
		while x!=i: #(i+1 possibly?)
			netw1=networks[list(networks.keys())[i]]
			netw2=networks[list(networks.keys())[x]]
			h=hamming_similarity(netw1,netw2)
			hammingnet.append([list(networks.keys())[i],list(networks.keys())[x],h])	
			x+=(-1)
	#! this part with the corresponding functions may be deleted in future versions
	#- create hist of Hamming distances in network
	hs = [x[2] for x in hammingnet] # all the h scores
	plothist(hs, os.path.join(resultsdir, "hist.png"))
	#! maybe will delete that too...
	#- select hamming threshold and add edges that pass it #! here addition of running threshold
	#! another 'hidden' threshold
	for hamlist in hammingnet:
		h=hamlist[2]
		if h > 0.5:
			meta.add_edge(hamlist[0], hamlist[1], weight=h) # specify edge data
elif thresh==1:
	#- compute SEWD for all networks (create network edges)
	for i in range(len(networks.keys())):
		x=max(range(len(networks.keys())))
		while x!=i: #(i+1 possibly?)
			netw1=networks[list(networks.keys())[i]]
			netw2=networks[list(networks.keys())[x]]
			h=SEWD_similarity(netw1,netw2)
			hammingnet.append([list(networks.keys())[i],list(networks.keys())[x],h])	
			x+=(-1)
	#- create hist of Hamming distances in network
	hs = [x[2] for x in hammingnet] # all the h scores
	#- changed resviewdir to outputdir here too
	plothist(hs, os.path.join(resultsdir, "hist.png"))
	#- caclulate critical value of 90th percentile for threshold 
	#! make this user selected
	#! another 'hidden' threshold
	cv=np.percentile(hs,90)
	print("Critical value selected : ",cv)
	#- create hamming network
	for hamlist in hammingnet:
		th=hamlist[2]
		if th>=cv:        #! make this running (based on previously defined histogram quantiles)
			meta.add_edge(hamlist[0], hamlist[1], weight=th) # specify edge data
	#! alternatively will add:
	#for hamlist in hammingnet:
	#	meta.add_edge(hamlist[0],hamlist[1],hamlist[2])

#- the next lines have been kept from original script
#- save matrix in mcl readable format
write_mcl_matrix(meta, os.path.join(resultsdir, "meta.matrix"))

#- save hamming net as json file for visualisation
meta_nodes={}
for i in range(len(meta.vs)):
	meta_nodes[i]=meta.vs[i]['name']
write_json(hammingnet, 'hammingnet.json', resultsdir, meta_nodes)

#- save hamming net as network file for alternative visualisation
meta.write_ncol(resultsdir+'/hammingnet')

#- run mcl on hammingnet to find clusters within
hamsterfile = run_mcl(os.path.join(resultsdir, "meta.matrix"), 2, redo=True) # always need to redo this.
hamcl = readmcl(hamsterfile)
print("hamcl:", hamcl)


#----------------------------
# Select best network/s
#----------------------------


connectednesswithincluster={}
clusterstats={}
best_networks={}


#! could improve to tell user the attributes of the cluster the resulting network belonged to
#- for each of the clusters of the hammingnet identified by MCL
for thiscluster in hamcl:
	print(thiscluster, len(hamcl[thiscluster]))
	clusterstats[thiscluster]=[]	
	for node in hamcl[thiscluster]:
		n = meta.vs[node]['name']
		connectednesswithincluster[n] = 0
		for othernode in hamcl[thiscluster]:
			if node==othernode:
				continue
			try:
				#- if edge exists in this direction
				w=meta.es.select(_source=node,_target=othernode)['weight'][0]
			except:
				try:
					#- if edge exists in the other direction
					w=meta.es.select(_source=othernode,_target=node)['weight'][0]
				except:
					#- if edge does not exist
					w=0
			#w = hammingdic[n][meta.vs[othernode]['name']]
			connectednesswithincluster[n]+=w		
		#print(n, connectednesswithincluster[n])
		clusterstats[thiscluster].append((n, connectednesswithincluster[n]))
	if absolute!=0:
		if len(hamcl[thiscluster])>absolute:
			#chosen_clustersizes[n] = len(hamcl[thiscluster])
			our_metric=0
			our_node=''
			for i in range(len(clusterstats[thiscluster])):
				if float(clusterstats[thiscluster][i][1])>float(our_metric):
					our_metric=clusterstats[thiscluster][i][1]
					our_node=clusterstats[thiscluster][i][0]
			best_networks[thiscluster]=(our_node,our_metric)
	#! this could be cumulative, i.e., include all clusters starting from biggest to smallest until a certain percentage
	if relative!=0:
		if len(hamcl[thiscluster])/float(len(networks)) > relative:
			#chosen_clustersizes[n] = len(hamcl[thiscluster])
			our_metric=0
			our_node=''
			for i in range(len(clusterstats[thiscluster])):
				if float(clusterstats[thiscluster][i][1])>float(our_metric):
					our_metric=clusterstats[thiscluster][i][1]
					our_node=clusterstats[thiscluster][i][0]
				#elif float(clusterstats[thiscluster][i][1])==float(our_metric): #! should consider what will happen if two values are equal
			best_networks[thiscluster]=(our_node,our_metric)


#- write best network names (and scores) into file
with open(resultsdir+'best_networks.txt', 'w') as nf:
	nf.write('\n'.join('%s %s' % (x, best_networks[x][0]) for x in list(best_networks.keys())))
nf.close()


#- create complete hammingnet
nodes = list(networks.keys())
mega = Graph()
mega.add_vertices(range(len(nodes)))
mega.vs["name"]=[nodes[x] for x in range(len(nodes))]
for hamlist in hammingnet:
	mega.add_edge(hamlist[0], hamlist[1], weight=hamlist[2]) # specify edge data


#- list of inflation values for d
#- constant list of inflation values
#infl_list=[round(x*0.1,2) for x in range(10,4,-1)]         # the round function is added because of aberrant values for 0.7 and 0.6
infl_list=m

#! if the list above gets deleted, loop below along with two in-loop instances need to be deleted too
#- identification of extreme networks (requires preceding EONS_03 code)
result_list=[]
for infl in infl_list:
	for ken in list(best_networks.keys()):
		#- find central node index in igraph n-o-n object
		cn=0
		cn_name = best_networks[ken][0]
		for indx in range(len(mega.vs)):
			if mega.vs[indx]['name']==cn_name:
				cn=indx
		# create accessible dictionary of edgeweights between central node and all other nodes
		lendict={}
		for n1 in range(len(mega.vs)):
			try:
				lendict[mega.vs[n1]['name']]=mega.es.select(_source=n1, _target=cn)['weight'][0]
			except:
				try:
					lendict[mega.vs[n1]['name']]=mega.es.select(_source=cn, _target=n1)['weight'][0]
				except:
					pass
		#- create sorted list of above
		sorted_w=sorted(lendict.items(),key=operator.itemgetter(1))
		#- assign first pole
		poles={}
		poles[sorted_w[0][0]]=sorted_w[0][1]
		distances=[]
		distances.append(sorted_w[0][1])
		#- set initial length
		nlength=len(poles)
		#- start iterations
		i=-1
		diff=500
		while diff!=0:
			plength=nlength
			for n in range(len(sorted_w)):
				dt=SEWD_similarity(networks[list(poles.keys())[i+1]],networks[sorted_w[n][0]])
				if ncheck(dt,distances,infl):
					poles[sorted_w[n][0]]=dt
					distances.append(dt)
					break
			nlength=len(poles)
			i+=1
			diff=abs(plength-nlength)
		#- format final line which will be added to final output list
		finline=[infl, cn_name]
		for final_key in list(poles.keys()):
			finline.append(final_key)
		result_list.append(finline)


#- write best and pole network names into file
with open(resultsdir+'extreme_networks.txt', 'w') as extreme_file:
	for newline in result_list:
		extreme_file.write("%s\n"%('\t'.join([str(x) for x in result_list])))
extreme_file.close()

