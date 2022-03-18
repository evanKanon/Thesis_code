#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

#EONS evaluation script 
#Evaluate fit of interesting networks with respect to clinical outcome

#------------------------------
#v0.0.2 first implementation of eval script in EONSv0.0.2
#v0.0.8 baseline script for evaluation, incorporation of average neighbour score, NDA and  score
#v0.1.0 BASELINE EONS added failed network functionality to eons_02 making all distances runable, corrected multiplier functionality in eons_03, made eons_eval save proper output, used similarity and not distance networks in eons_funct i.e. 0 is least similar
#v0.1.1 changed make_expression function to accomodate file names other than .soft, added NA dealing option to eons_prep, improved config option for eons_main, made code python 3 compatible
#v0.1.2 established perpendicular point plots
#v0.1.3 improved folder names, added .layout visualisation, updated environment file, deleted unecessary options of eons_eval
#v0.1.4 made running the MCL component faster and allowed user to select inflation values
#v0.1.5 changed naming of folders and certain default values
#v0.1.6 corrected PCA random selection in eons_01, corrected network names in eons_02 and deleted unused inputfile option in eons_03 with the appropriate corrections in eons_funct and _main
#v0.1.7 added correct evaluation measure

#------------------------------
__version__ = '0.1.7'
#------------------------------
import csv
import argparse
from eons_funct import *
#------------------------------

#- Options
parser=argparse.ArgumentParser(description='EONS-Evaluation_based_on_outcome')
parser.add_argument('-s', '--save', type=str, help='file path for results and figures to be saved')
parser.add_argument('-class', '--class_file', type=str, help='file path for class file')
parser.add_argument('-thr', '--thresholded', type=int, default=1, choices=[0,1], help='Choose whether threshlods will be applied to networks: 0- apply thresholds, 1- run using -99 threshold, making it threshold-free')
#! ADD OPTIONS TO SELECT EVALUATION METHOD !! DEFAULT == TRUE TO ALL
#! ADD OPTIONS TO CHOOSE CLASSES
args=parser.parse_args()

#- Processing arguments 
outputdir=args.save
classfile=args.class_file
thrs=args.thresholded

#- create results directory
#! suffix potentially obsolete
evaldir = os.path.join(outputdir, 'evaluation_'+classfile.split('/')[-1].split('.')[0]+'/')
check_dir(evaldir)

#- identify folder of networks used for evaluation
if os.path.isdir(outputdir+'networks/')==True:
	networkdir=outputdir+'networks/'
else:
	#! update this to output to slack
	print('ERROR: No networks folder found')

#- inputing classfile
nodeclasses, classlists = read_classes(classfile)

#- select which classes we will compute NDA scores for
#! make this user optionable too!!
chosenclasses = [] #['gram-positive_sepsis', 'gram-negative_sepsis']
chosenclasses = list(set(chosenclasses) & set(classlists.keys())) # classes only count if they were found in the class file
if len(chosenclasses)==0:
	chosenclasses = list(classlists.keys())
print("chosenclasses", chosenclasses)


#- Read-in input
networks={}
for filename in os.listdir(networkdir):
	s=os.path.abspath(networkdir)
	if filename.endswith('.csv'):
		networks[filename.split('.csv')[0]]=pd.read_csv(s+'/'+filename,index_col=0)


#- initiate best dictionary where the names of the top networks from each evaluation will be saved in 
bestdict={}


#! this should be made subsettable, to run analysis only on networks that are desired

#- create list of attributes to drop
#dropclass
droppat=[x for x in set(classlists.keys()) if x not in chosenclasses]
#- subset networks if there are fewer classes
if len(droppat)!=0:
	#- subset networks droping patients which are of the unchosen classes
	for onekey in networks.keys():
		for notclass in droppat:
			networks[onekey]=networks[onekey].loc[~networks[onekey].index.isin(classlists[notclass])]
			networks[onekey]=networks[onekey].filter(networks[onekey].index)
else:
	pass


####################################################################################
# AVERAGE NETWORK NEIGHBOUR
####################################################################################


#- initialise variables for storage of final comparison result
network_dict={}
if thrs==1:
	#- for each network in the networks produced by the analysis
	for netkey in networks.keys():
		network=networks[netkey]
		#- initialise relevant variable to store resulting class for each patient
		pat_dict={}
		#- for each patient in the network (here taking columns)
		for i in network.keys():
			#- initialise variables to store maximum sum, class which exhibits it and 
			maxsum=0
			maxclass=None
			#- for each class in the desired classes
			for myclass in chosenclasses:
				#- consider the relevant patient list without including the current patient
				patlist = [x for x in classlists[myclass] if x != i]
				#- calculate the average sum of edgeweights in this subset
				mysum = network[network.index.isin(patlist)][i].sum()/len(patlist)
				#- see if new value is better than previous maximum, if so update
				if mysum > maxsum:
					maxsum=mysum
					maxclass=myclass
			pat_dict[i]=maxclass
		#- initialise sum of correctly identified patients
		finsum=0
		#- score 1 for correct identification and 0 for incorrect
		for pat in pat_dict.keys():
			if pat_dict[pat]==nodeclasses[pat][0]:
				finsum+=1
			else:
				finsum+=0
		#- return final quality value as a percentage of correctly identified patients
		finvalue=float(finsum)/len(pat_dict.keys())
		#- save value in the dictionary
		network_dict[netkey]=finvalue
else:
	#- for each network in the networks produced by the analysis
	for netkey in networks.keys():
		network=networks[netkey]
		#- initialise relevant variable to store resulting class for each patient
		pat_dict={}
		#- for each patient in the network (here taking columns)
		for i in network.keys():
			#- initialise variables to store maximum sum, class which exhibits it and 
			maxsum=0
			maxclass=None
			#- for each class in the desired classes
			for myclass in chosenclasses:
				#- consider the relevant patient list without including the current patient
				patlist = [x for x in classlists[myclass] if x != i]
				#- calculate the sum of edgeweights in this subset
				mysum = network[network.index.isin(patlist)][i].sum()
				#- see if new value is better than previous maximum, if so update
				if mysum > maxsum:
					maxsum=mysum
					maxclass=myclass
			pat_dict[i]=maxclass
		#- initialise sum of correctly identified patients
		finsum=0
		#- score 1 for correct identification and 0 for incorrect
		for pat in pat_dict.keys():
			if pat_dict[pat]==nodeclasses[pat][0]:
				finsum+=1
			else:
				finsum+=0
		#- return final quality value as a percentage of correctly identified patients
		finvalue=float(finsum)/len(pat_dict.keys())
		#- save value in the dictionary
		network_dict[netkey]=finvalue


#- save dictionary in csv format
#average network neighbour
df=pd.DataFrame(list(network_dict.items()), columns=['Network', 'ANN_Predictability'])
df.to_csv(evaldir+"ANN_evaluation_scores.csv",index=False)

#- save best result from ANN evaluation
bestdict['ANN']=[max(network_dict,key=network_dict.get),network_dict[max(network_dict,key=network_dict.get)]]
#max(network_dict,key=network_dict.get)             # find name of network with highest percentage of patients correctly identified

#- save dictionary in txt format
#dict = {'Python' : '.py', 'C++' : '.cpp', 'Java' : '.java'}
#f = open("dict.txt","w")
#f.write( str(dict) )
#f.close()


####################################################################################
# NODE DENSITY ANALYSIS
####################################################################################

#- COMPUTE NDA - NEW (CORE FUNCTIONALITY)
#networkscores={}
#for network in networks.keys():
#	allscores=[]
#	for thisclass in chosenclasses:
#		nda_score = nda_new(networks[network], classlists[thisclass])
#		allscores.append(nda_score)
#	av_score=forgivingmsd(allscores)[0]
#	networkscores[network]=av_score


#- COMPUTE NDA - NEW
networkscores={}
allnetscores={}
net=[]
for thclass in chosenclasses:
	net.append(thclass)
for network in networks.keys():
	allscores=[]
	#- class header used as sanity check
	class_header=[]
	for thisclass in chosenclasses:
		class_header.append(thisclass)
		nda_score = nda_new(networks[network], classlists[thisclass])
		allscores.append(nda_score)
	av_score=forgivingmsd(allscores)[0]
	networkscores[network]=av_score
	if class_header==net:
		allnetscores[network]=allscores
	else:
		print('Problem with header! Evaluate Manually!')


#- sort average scores by stregth
sorted_x=sorted(networkscores.items(),key=operator.itemgetter(1),reverse=True)


#- plot results to see the distribution of NDA scores
plt.bar(networkscores.keys(), networkscores.values(), color='g')
plt.savefig(evaldir+'NDA_scores.png')
plt.gcf().clear()

#- save results as CSV file
nda_file=pd.DataFrame.from_dict(networkscores,orient='index',columns=['NDA_avg_scores'])
nda_file.to_csv(evaldir+'NDA_avg_evaluation_scores.csv')
#with open(evaldir+'NDA_avg_evaluation_scores.csv', 'w') as f:
#	writer = csv.writer(f)
#	f.write('Network,avg_NDA_score\n')
#	for row in networkscores.items():
#		writer.writerow(row)

#- save sorted results as text file
with open(evaldir+'NDA_avg_evaluation_scores_sorted.txt', 'w') as nf:
	nf.write('Network\tavg_NDA_score\n')
	nf.write('\n'.join('%s %s' % x for x in sorted_x))

#- create df to save full list of NDA scores
ndf=pd.DataFrame.from_dict(allnetscores,orient='index',columns=net)
ndf.to_csv(evaldir+'NDA_all_evaluation_scores.csv')

#- save best result from NDA evaluation
bestdict['NDA']=sorted_x[0]

####################################################################################
# IN GROUP PROPORTION
####################################################################################


#- In Group Proportion (IGP)
#- i.e., how many patients have as their nearest neighbour a patient of the same class
#- idea taken from Kapp and Tibshirani Biostatistics 2007 8(1):9-31


network_dict={}
#- for each network in the networks produced by the analysis
for netkey in list(networks.keys()):
	network=networks[netkey]
	classdict={}
	#- for each class in the given classes
	for myclass in chosenclasses:
		thislist=classlists[myclass]
		#- for each patient in this class
		acc=0
		for patkey in thislist:
			#- identify its nearest neighbour
			for pj in range(len(network[patkey])):
				#- find nearest neighbour that is not the same patient
				if network[patkey].nlargest(10,keep='all').index[pj]!=patkey:
					nearest_neighbour=network[patkey].nlargest(10,keep='all').index[pj]
					break
			#- check if the nearest neighbour is of the same class
			if nodeclasses[nearest_neighbour]==nodeclasses[patkey]:
				acc+=1
		ig_p=float(acc)/len(thislist)
		classdict[myclass]=ig_p
	network_dict[netkey]=classdict

#- create a header for the dictionary of dictionaries
myheader=['network']
for i in chosenclasses:
	myheader.append(i)
#- create a list of lists containing the relevant data
mydata=[]
for mi in list(network_dict.keys()):
	myvalues=[mi]
	for j in myheader[1:]:
		myvalues.append(network_dict[mi][j])
	mydata.append(myvalues)
#- create a dataframe from the records and save in the appropriate fashion
igpdf=pd.DataFrame.from_records(mydata,columns=myheader)
igpdf.to_csv(evaldir+'igp_evaluation_scores.csv',index=False)

#- get average igp values for each patient in a dict
netigpscores={}
for item in mydata:
	netigpscores[item[0]]=forgivingmsd(item[1:])[0]
#- save average igp values
avg_igp=pd.DataFrame.from_dict(netigpscores,orient='index',columns=['igp_avg_score'])
avg_igp.to_csv(evaldir+'igp_average_scores.csv')
#- sort average igp values 
sorted_y=sorted(netigpscores.items(),key=operator.itemgetter(1),reverse=True)
#- save best result from IGP evaluation
bestdict['IGP']=sorted_y[0]



####################################################################################
# FIRST NEAREST NEIGHBOUR
####################################################################################



onetwork_dict={}
#- for each network in the networks produced by the analysis
for netkey in list(networks.keys()):
	network=networks[netkey]
	classdict={}
	acc=0
	tot=0
	#- for each class in the given classes
	for myclass in chosenclasses:
		thislist=classlists[myclass]
		#- for each patient in this class
		#acc=0
		for patkey in thislist:
			#- identify its nearest neighbour
			for pj in range(len(network[patkey])):
				#- find nearest neighbour that is not the same patient
				if network[patkey].nlargest(10,keep='all').index[pj]!=patkey:
					nearest_neighbour=network[patkey].nlargest(10,keep='all').index[pj]
					break
			#- check if the nearest neighbour is of the same class
			if nodeclasses[nearest_neighbour]==nodeclasses[patkey]:
				acc+=1
		#ig_p=float(acc)/len(thislist)
		tot+=len(thislist)
		#classdict[myclass]=ig_p
	classdict=float(acc)/tot
	onetwork_dict[netkey]=classdict


#- save average igp values
fnn_score=pd.DataFrame.from_dict(onetwork_dict,orient='index',columns=['FNN _score'])
fnn_score.to_csv(evaldir+'fnn_scores.csv')
#- sort average igp values 
sorted_y=sorted(netigpscores.items(),key=operator.itemgetter(1),reverse=True)
#- save best result from IGP evaluation
bestdict['FNN']=sorted_y[0]

#--------------------------------------------------
# Save best snippet
#--------------------------------------------------

#- create and save final top result from each evaluation technique in a single txt file
#best_df=pd.DataFrame.from_dict(bestdict,orient='index')
#! this is implemented to run on lx05 server
best_df=pd.DataFrame.from_dict(bestdict).T
best_df.to_csv(evaldir+'top_evaluation_compound.txt',sep='\t',header=False)



####################################################################################
# VISUALISATION
####################################################################################

#! IMPROVE VISUALISATION INTO A SINGLE (OR MULTIPLE INTERLOCKING) COHERENT FUNCTION

#- create common df of evaluation metrics
mydf=pd.merge(nda_file,avg_igp,right_index=True,left_index=True)
ann_file=df.set_index('Network')
mydf=pd.merge(ann_file,mydf,right_index=True,left_index=True)
mydf=pd.merge(fnn_score,mydf,right_index=True,left_index=True)

#- read in interesting (best) networks resulting from EONS analysis
#res=pd.read_csv(outputdir+'results/best_networks.txt',sep=' ',header=None,index_col=0)
#of_interest=list(res[1])
of_interest=list(pd.read_csv(outputdir+'results/best_networks.txt',sep=' ',header=None,index_col=0)[1])

#- plot outcome
#ppd_plot(mydf,highlight=of_interest,outputname=evaldir+'Evaluator_value_distributions.png')


coord_dict={'x':[],'y':[],'c':[],'s':[]}
#- fill the dictionary in with the relevant points
x_values=[]
x_val=3
for c in mydf.keys():
	x_values.append(x_val)
	for i in range(len(mydf)):
		coord_dict['x'].append(x_val)
		coord_dict['y'].append(mydf[c][i])
		if mydf.index[i] in of_interest:
			coord_dict['c'].append(1)
		else:
			coord_dict['c'].append(0)
		if c.startswith('NDA'):
			coord_dict['s'].append(2)
		else:
			coord_dict['s'].append(1)
	x_val+=1
#- set proper names tuple for plot x-axis
mycolname=list(mydf.keys())
mycols=['']+mycolname+['']
mycols=tuple(mycols)
#- turn the coordinate dictionary into pandas dataframe from plotting
plotdf=pd.DataFrame(coord_dict)


#- create two subplots with 2:1 ratio
f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
a0.scatter(plotdf['x'][plotdf['c']==0][plotdf['s']==1],plotdf['y'][plotdf['c']==0][plotdf['s']==1],c='b')
a0.scatter(plotdf['x'][plotdf['c']==1][plotdf['s']==1],plotdf['y'][plotdf['c']==1][plotdf['s']==1],c='r')

a1.scatter(plotdf['x'][plotdf['c']==0][plotdf['s']==2],plotdf['y'][plotdf['c']==0][plotdf['s']==2],c='b')
a1.scatter(plotdf['x'][plotdf['c']==1][plotdf['s']==2],plotdf['y'][plotdf['c']==1][plotdf['s']==2],c='r')

#- write correct x-axis values
a0.set_xticklabels(['','ANN_score','','','','FNN_score'])	
a1.set_xticklabels(['','','NDA_score'])

f.suptitle(classfile.split('/')[-1:][0], fontsize=16)

f.savefig(evaldir+'Evaluator_value_distributions.png')

plt.gcf().clear()
