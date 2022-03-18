#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------

#CellMapper
#A script that maps unknown samples to a reference atlas (fantom5) using spearman correlation

#------------------------------
#v0.0.0 initial version of the new code, using elements from original scripts
#v0.0.1 added pvalue selection method
#v0.0.2 made more generalisable, added corr option and added config file
#v0.0.3 defined seperate functions file, deleted unnecessary code, made subset selection option, made log option working, created visualisation files&folder and set proper (not hard) config file path

#-----------------------------
# Definition of version
#-----------------------------

__version__ = '0.0.3'

#-----------------------------
# Import relevant packages
#-----------------------------
import argparse
from cellmapfunctions import *
#-----------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#-----------------------------

#- Options
parser=argparse.ArgumentParser(description='Mapping cell types to reference map')
parser.add_argument('-f', '--file', type=str, help='file path for input file')
parser.add_argument('-m', '--mapfile', type=str, help='file path for file to be used as a map', default='0')
parser.add_argument('-o', '--output', type=str, help='file path for results and figures to be saved', default='0')
parser.add_argument('-a', '--ann_file', type=str, help='file path for annotation file', default='0')
parser.add_argument('-map_a', '--map_annot', type=str, help='file path for annotation file for the map', default='0')
parser.add_argument('-cor', '--correlation', default=False, action='store_true', help='create correlation matrix for the map file')
parser.add_argument('-sub', '--save_subset', action='store_true', default=False, help='save files of subsets of inputfile')
parser.add_argument('-coord', '--coordinates', type=str, help='file path for map coordiates file', default='0')
parser.add_argument('-c', '--configfile', default=os.path.join(scriptpath,'config.json'), help='custom config file')
#arser.add_argument('-fnn', '-file_near', type=int, default=3, help='number of nearest neighbours for input file')
#parser.add_argument('-mnn', '-map_near', type=int, default=5, help='number of nearest neighbours for map')
parser.add_argument('-log', '--log_file', default=False, action='store_true', help='create log file containing the argunents used, the date run, the version and the input file name')
parser.add_argument('-slack', '--slack', default=False, action='store_true', help='send messages to user through Slack platform')
args=parser.parse_args()


#- open config file
#! need to see how to make config file path
with open(args.configfile) as json_data_file:
	config=json.load(json_data_file)

#- SET SLACK PARAMETERS
slack=Slacker(config['Slack_env']['address'])
#! make allowance for more than one user to be notified
user=config['Slack_env']['user1']
slack_arg=args.slack

#- SET MAIN PATHS
#- define inputfile path
if config['File_paths']['input']:
	inputfile=str(config['File_paths']['input'])+args.file
else:
	inputfile=args.file
#- define outputdir path
if config['File_paths']['output']:
	if args.output=='0':
		outputdir=str(config['File_paths']['output'])
	else:
		outputdir=str(config['File_paths']['output'])+args.save
else:
	outputdir=args.save
#- define mapfile path
if config['File_paths']['map_file']:
	if args.mapfile=='0':
		mapfile=str(config['File_paths']['map_file'])
	else:
		mapfile=str(config['File_paths']['map_file'])+args.mapfile
else:
	mapfile=args.mapfile
#- define annotation for query file path
if config['File_paths']['ann_file']:
	if args.ann_file=='0':
		annofile=str(config['File_paths']['ann_file'])
	else:
		annofile=str(config['File_paths']['ann_file'])+args.ann_file
else:
	annofile=args.ann_file
#- define annotation for map file path
#! add else clause to ensure ability to run with no annotation file, find example files that can be used like that
if config['File_paths']['map_annot']:
	if args.map_annot=='0':
		mapannfile=str(config['File_paths']['map_annot'])
	else:
		mapannfile=str(config['File_paths']['map_annot'])+args.map_annot
else:
	mapannfile=args.map_annot
#- define map coordinates file path
if config['File_paths']['map_coord']:
	if args.coordinates=='0':
		ncoord=str(config['File_paths']['ann_file'])
	else:
		ncoord=str(config['File_paths']['ann_file'])+args.coordinates
else:
	ncoord=args.coordinates

#- Set all other arguments
corr_args=args.correlation

#- Create overarching output folder
daterun='_'+str(time.strftime("%d%m%Y"))
analysislabel = os.path.split(inputfile)[1].split('.')[0]
outputdir = os.path.join(outputdir,analysislabel+daterun+'_CellMv'+str(__version__)+'/')
if not os.path.isdir(outputdir):
	os.mkdir(outputdir)


#- Optionally create log file for run
logfile=args.log_file

if logfile==True:
	if not os.path.isdir(outputdir):
		try:
			os.mkdir(outputdir)
		except:
			updir = os.path.split(outputdir)[0]
			os.mkdir(updir)
			os.system("/bin/chmod 755 %s"%(updir))
			os.mkdir(this_dir)
	if outputdir.endswith('/'):
		log=open(outputdir+'logfile.txt','w')
	else:
		log=open(outputdir+'/logfile.txt','w')
	log.write("Date run: "+time.strftime("%d%m%Y")+'\n')
	log.write("On datafile: "+inputfile+'\n')
	log.write('Using EONS version: '+__version__+'\n')
	log.write('With arguments:\n')
	arglist=args
	#print(sys.argv)
	#print(args)
	log.write(str(arglist)+'\n')
	log.write('\n')
	log.write('Full command used was:\n')
	log.write("%s\n"%(' '.join([str(x) for x in sys.argv])))
	log.close()
else:
	pass

#- Inform user that analysis is commencing
if slack_arg == True:
	slack.chat.post_message(user,'Commencing analysis')


#--------------------------------
# Format the Query Dataset
#--------------------------------


#- Import query dataset (henderson)
query=pd.read_csv(inputfile)
if slack_arg == True:
	slack.chat.post_message(user,'Query dataset imported')

#- Import the file with the cluster annotation
#! this needs improvement, clusters.x.unique will fail if annotation file has different headers, should generalise
clusters=pd.read_csv(annofile)
clusters_list=clusters.x.unique() #! is also created later, check why
#! this provides better generalisation, overcomes header name issu but requires specific format of annofile
#clusters_list=clusters[clusters.keys()[1]].unique()

#! this seems unecessary in it current form, should delete
#! add '_' to cluster names if they have space
for c in clusters_list:
	if " " in c:
		c.replace(' ','_')



#- save copies of subsets of data in designated folder
if args.save_subset:
	#-Create a directory for the subsets
	subsetdir = os.path.join(outputdir, 'subsets/')
	check_dir(subsetdir)
	#- subset gene expression based on cluster participation
	for c in clusters_list:
		clustername='my_cluster%s'%(c)
		#list_of_clust=[]
		mydf=pd.DataFrame({'A':[]})
		for ctype in range(len(clusters['Unnamed: 0'])):
			if clusters.x[ctype]==c:
				mydf[query.keys()[clusters['Unnamed: 0'][ctype]]]=query[query.keys()[clusters['Unnamed: 0'][ctype]]]
				#list_of_clust.append(query[query.keys()[ctype]])
		mydf.to_csv(subsetdir+'/'+clustername)
else:
	pass


#- create a file averaging all the expression values of the genes within each clustered sample
sumdf=pd.DataFrame()
sumdf['Gene_names']=query[query.keys()[0]].values
for c in clusters_list:
	clustername='Cluster%s_average'%(c)
	list_of_clust=[]
	for ctype in range(len(clusters['Unnamed: 0'])):
		if clusters.x[ctype]==c:
			list_of_clust.append(query[query.keys()[clusters['Unnamed: 0'][ctype]]])
	average=sum(list_of_clust)/len(list_of_clust)
	sumdf[clustername]=average
#! improve naming convention
sumdf.to_csv(outputdir+'query_cluster_sum_expression.csv')




#- create a correlation map of the query average expression clusters
corrquery=sumdf.corr(method='spearman')
#! naming could be improved
corrquery.to_csv(outputdir+'query_network.csv')



#- readin query average clusters (would be same as using sumdf but values are slightly different, because when written some decimals get chopped off)
query_average=pd.read_csv(outputdir+'query_cluster_sum_expression.csv',index_col=0)
#- calculate the average expression across the clusters
query_average['mean']=query_average.mean(axis=1)
#- delete all genes which have zero expression across the clusters (they will increase similarity)
query_average_nozeros=query_average[query_average['mean'] != 0]
#-drop mean column, since it has no further usefulness
query_average_nozeros=query_average_nozeros.drop('mean',axis=1)




#--------------------------------
# Format the Fantom5 Dataset
#--------------------------------

#- import the data file
background=pd.read_table(mapfile)


#! This is a half-way measure, need to think of something more practical
#- Perfom data preparation if mapfile is FANTOM5
if mapfile.endswith('nozeros_phase2_nopooled.expression'):
	#- import the metadata file
	metadata=pd.read_table(mapannfile)
	#- finding relative expression of each region
	#! create this functionality

	#- selecting gene regions based on peak 1
	gene_regions=[]
	name_list=[]
	gene_reg_association=[]
	for i in range(len(metadata)):
		try:
			peak,name=metadata.short_description[i].split('@')
			if peak=='p1':
				gene_regions.append(metadata['00Annotation'][i])
				name_list.append(name)
				gene_reg_association.append((metadata['00Annotation'][i],name))
			else:
				pass
		except:
			hold=metadata.short_description[i].split('@')
			all_items=[]
			for obj in hold:
				try:
					one,two=obj.split(',')
					all_items.append(one)
					all_items.append(two)
				except:
					all_items.append(obj)
			for pos in range(len(all_items)):
				if all_items[pos]=='p1':
					sm_name=all_items[pos+1]
					gene_regions.append(metadata['00Annotation'][i])
					name_list.append(sm_name)
					gene_reg_association.append((metadata['00Annotation'][i],sm_name))
				else:
					pass

	#- subset query dataset based on common genes names with genes which have peak1
	#query_subset=query_average[query_average['Gene_names'].isin(name_list)]          #- old version which keeps all genes, even with only zero values
	query_subset=query_average_nozeros[query_average_nozeros['Gene_names'].isin(name_list)]

	#- select common gene names and genomic regions combinations of fantom5 
	query_names=query_subset['Gene_names']
	common_gene=[]
	for i in range(len(gene_reg_association)):
		if gene_reg_association[i][1] in set(query_names):
			common_gene.append(gene_reg_association[i])

	#- make dataframe of gene-region dictionary to be used for quick name discovery (used later)
	common_genes=pd.DataFrame(common_gene)

	#- keep only gene regions to use for subselection of the fantom5 dataset
	common_gene_regions=[]
	for i in range(len(common_gene)):
		common_gene_regions.append(common_gene[i][0])

	#- keep gene regions which have gene equivalent in query 
	#- common_gene_association is the dictionary for the gene locations that are common between the two datasets
	#! this is redundant as it is the same as common_genes (choose which one to keep) !!!!
	common_gene_association=[]
	for i in range(len(gene_reg_association)):
		for jn in query_names:
			if gene_reg_association[i][1]==jn:
				common_gene_association.append(gene_reg_association[i])

	#- replace query gene names with fantom5 genomic regions
	for i in range(len(common_gene)):
		#for jl in query_subset['Unnamed: 0']:
		for jl in query_subset['Gene_names']:
			if jl==common_gene[i][1]:
				#query_subset['Unnamed: 0']=query_subset['Unnamed: 0'].replace(jl,common_gene[i][0])
				query_subset['Gene_names']=query_subset['Gene_names'].replace(jl,common_gene[i][0])

	#- create dataframes with proper indexes, sorted the same way for comparison
	x=pd.DataFrame()
	for i in range(len(common_gene_regions)):
		x=x.append(background[background['00Annotation']==common_gene[i][0]])         # ignore_index=True

	#! in final implementation will drop the addition of another variable keeping only x and subset 
	mapdf=x.sort_values('00Annotation', axis=0)
	mapdf = mapdf.reset_index(drop=True)
	#mapdf2=query_subset.sort_values('Unnamed: 0', axis=0)
	querydf3=query_subset.sort_values('Gene_names', axis=0)
	querydf3 = querydf3.reset_index(drop=True)
else:
	#- create list with gene names of map
	maplist=list(background[background.keys()[0]])
	#- create list with gene names of qiery
	querylist=list(query_average_nozeros[query_average_nozeros.keys()[0]])
	#- Subset both lists by each other's name set
	query_subset=query_average_nozeros[query_average_nozeros[query_average_nozeros.keys()[0]].isin(maplist)]
	map_subset=background[background[background.keys()[0]].isin(querylist)]
	#- create relevant df of map file
	mapdf = map_subset.sort_values(background.keys()[0], axis=0)
	mapdf = mapdf.reset_index(drop=True)

	#! this needs to be incorporatd in the beginning of the background file of this and every iteration
	mapdf = mapdf.drop(columns=['Illumina Probe ID'])

	#- create relevant df of query file
	querydf3 = query_subset.sort_values(query_average_nozeros.keys()[0], axis=0)
	querydf3 = querydf3.reset_index(drop=True)



#--------------------------------
# Calculate the correlation between query and map
#--------------------------------

#- send message to user
if slack_arg == True:
	#! add improvement to show names of files
	slack.chat.post_message(user,'Commencing correlation computation between Query and Map')

#- correlation matrix of the two datasets
#- create header row to be used for sanity check and in final df
header=[]
header.append('queryCluster')
for f_key in mapdf.keys():
	if f_key==mapdf.keys()[0]:
		pass
	else:
		header.append(f_key)

#- Run the correlation
corr_df=[]
for query_key in querydf3.keys():
	if query_key==querydf3.keys()[0]:
		pass
	else:
		#slack.chat.post_message(user,'Commencing correlation computation of query cluster: '+str(query_key))
		myline=[]
		myline.append(query_key)
		testline=[]
		testline.append('queryCluster')
		for fantom_key in mapdf.keys():
			if fantom_key==mapdf.keys()[0]:
				pass
			else:
				myline.append(sp.stats.spearmanr(mapdf[fantom_key],querydf3[query_key],nan_policy='omit')[0])
				testline.append(fantom_key)
		if testline==header:
			corr_df.append(myline)
		else:
			print('This line does not match Header, check manually ', query_key)

#- turn list of lists created above into a dataframe
df = pd.DataFrame(corr_df, columns=header)
#- set index
df = df.set_index('queryCluster')
#- save dataframe as csv file
#! improve naming convention
df.to_csv(outputdir+'query_map_spearcorr_matrix.csv')



#- set or compute mapfile correlation for emp p-value computation
if corr_args == True:
	corr_dat = mapdf.corr('spearman')
	#! possibly add the option to save it for future use
	corr_dat.to_csv(outputdir+'map_subset_correlation.csv')
else:
	corr_dat = pd.read_csv(config['File_paths']['corr_file'], index_col=0)


#- Create appropriate directory for data
resultsdir=os.path.join(outputdir+'results/')
check_dir(resultsdir)

#- Calculate emp pvalues
mydataframe=df.T
res_header=['QueryCluster','MapCluster','Spearman-corr','emp p-value','emp p-value Fraction','P(A^B)']
res_flist=[]
for i in mydataframe.keys():
	#if i=='ClusterMPs (4)_average' or i=='ClusterMPs (5)_average':
	res_flist=[]
	p_value_scores={}
	for j in mydataframe.index:
		w=mydataframe[i][j]
		thislist=sorted(corr_dat[j])
		del thislist[-1]
		p=empiricalp(w, thislist) # change to weight_list to comp dstn based on entire daatset, mydataframe[i] to find pvalue of the Henderson node, .iloc would be for position, not name
		#othlist=sorted(corrquery[i])
		#del othlist[-1]
		#p2=empiricalp(w,othlist)
		p2=empiricalp(w,sorted(mydataframe[i]))
		# Average Probability
		#pab=(p+p2)/2
		# Probability of Intersection
		pab=p*p2
		# Probability of Union
		#pab=p+p2-(p*p2)
		#pab=empiricalp_both(w,thislist,mydataframe[i])
		f=empiricalp_fract(w, thislist) # change to weight_list to comp dstn based on entire dataset, mydataframe[i] to find pvalue of the Henderson node, .iloc would be for position, not name
		p_value_scores[j]=p
		myline=[i,j,w,p,f,pab]
		#- identify top genes (for all associations will take VERY, VERY, VERY LOOONNGGG)
		#query_ranks = pd.Series(mapdf[j]).rank(ascending=False)
		#map_ranks = pd.Series(querydf3[i]).rank(ascending=False)
		#gene_list,gene_dict=find_genes(map_ranks,query_ranks,mapdf,querydf3,common_genes)
		#for gene_id in range(len(gene_list)):
			#- update header accordingly
			#res_header.append('gene_%s'%gene_id)
		#	myline.append(gene_list[gene_id])
		res_flist.append(myline)
	#- Create and sort emp p-value dataframe
	resdf = pd.DataFrame(res_flist, columns=res_header)
	ressort=resdf.sort_values(by=['QueryCluster','P(A^B)','Spearman-corr'],ascending=[False,True,False])
	#! improve naming convention
	ressort.to_csv(resultsdir+str(i)+'_new_emp_pval_full_sorted.txt',sep='\t',index=False)


#- select top 20 genes for top 5 emp p-vals across all cell types
k=5
newlines=[]
for filename in os.listdir(resultsdir):
	s=os.path.abspath(resultsdir)
	if filename.endswith("_new_emp_pval_full_sorted.txt"):
		mydf=pd.read_table(s+'/'+filename)
		nearest=mydf[:k]
		for index, row in nearest.iterrows():
			key1=row[nearest.keys()[0]]
			key2=row[nearest.keys()[1]]
			elem1 = row[nearest.keys()[2]]
			elem2 = row[nearest.keys()[3]]
			elem3 = row[nearest.keys()[4]]
			newline=[key1,key2,elem1,elem2,elem3]
			query_ranks = pd.Series(mapdf[key2]).rank(ascending=False)
			map_ranks = pd.Series(querydf3[key1]).rank(ascending=False)
			gene_list,gene_dict=find_genes(map_ranks,query_ranks,mapdf,querydf3,common_genes,maxrank=500)
			for gene_id in range(len(gene_list)):
				#- update header accordingly
				newline.append(gene_list[gene_id])
			newlines.append(newline)
#- create dataframe of result with gene names
results_df=pd.DataFrame(newlines)
#- save dataframe
#! improve naming convention
results_df.to_csv(resultsdir+'top_result_genes.txt',sep='\t',index=False)


#- Format result
#! improve functionality, this is a half way measure
if mapfile.endswith('nozeros_phase2_nopooled.expression'):
	for filename in os.listdir(resultsdir):
		s=os.path.abspath(resultsdir)
		if filename.endswith("_new_emp_pval_full_sorted.txt"):
			#- open unpolished file and save its lines for processing
			nf=open(s+'/'+filename)
			lines=[string.split(string.strip(x),'\t') for x in nf.readlines()]
			nf.close()
			#- clean names if they are part of the FANTOM5 CT
			savename=s+'/'+filename.split('emp')[0]+'emp_pval_full_sorted_clean.txt'
			fo=open(savename,'w')
			for line in lines:
				line[1] = convert_f5_samplename(convert_f5_samplename(line[1]))
				fo.write("{}\n".format('\t'.join(line)))
			fo.close()
		elif filename.endswith("top_result_genes.txt"):
			#- open unpolished file and save its lines for processing
			nf=open(s+'/'+filename)
			lines=[string.split(string.strip(x),'\t') for x in nf.readlines()]
			nf.close()
			#- clean names if they are part of the FANTOM5 CT
			savename=s+'/'+'top_result_genes_clean.txt'
			fo=open(savename,'w')
			for line in lines:
				line[1] = convert_f5_samplename(convert_f5_samplename(line[1]))
				fo.write("{}\n".format('\t'.join(line)))
			fo.close()

#- send message to user
if slack_arg == True:
	#! add improvement to show names of files
	slack.chat.post_message(user,'...Correlation computation between Query and Map completed...')
	slack.chat.post_message(user,'......Visualisation files creation initiated......')

#------------------------------
#- CREATE VISUALISATION FILES
#------------------------------

#! This part is new and might require improvement - further generalisation


#- Create appropriate directory for data
visualdir=os.path.join(outputdir+'visualisation/')
check_dir(visualdir)

#- if we are using FANTOM5 as a map #! aybe should generalise this statement throughout the script
if mapfile.endswith('nozeros_phase2_nopooled.expression'):
	#- find names of clusters from resulting files and make a list of them
	namelist=[]
	filelist=[]
	for filename in os.listdir(resultsdir):
		if filename.startswith('Cluster') and not filename.endswith('_clean.txt'):
			fname=filename.split('.')[0]
			namelist.append(fname)
			filelist.append(filename)
	#- from each file select the top 3 results and keep the relevant edgeweights
	#- also, select the name of the top fantom5 node to be used for x,y coord
	addlines=[]
	nodecoord=[]
	for file in filelist:
		f=open(resultsdir+file)
		lines=[string.split(string.strip(x),'\t') for x in f.readlines()]
		f.close()
		#- select top 3 rows (except header) == means that the format of the inputfiles is always the same
		for top in range(1,4):
			if top==1:
				nodecoord.append([lines[top][0],lines[top][1]])
				addline=[lines[top][0],lines[top][1],lines[top][5]]
				addlines.append(addline)
			else:
				addline=[lines[top][0],lines[top][1],lines[top][5]]
				addlines.append(addline)
	#- read in the coordinates from the relevant file
	ncoordfile=open(ncoord)
	ncoordlines=[string.split(string.strip(x),'\t') for x in ncoordfile.readlines()]
	ncoordfile.close()
	#- for each cluster, assign the proper coordinates
	newclines=[]
	for i in nodecoord:
		for cline in ncoordlines:
			if i[1]==cline[0]:
				#element='"{0}"'.format(i[0])
				newcoordline=[i[0], cline[1], cline[2], 1000]
				newclines.append(newcoordline)
	#- merge new and old coordinates
	for fcoord in newclines:
		ncoordlines.append(fcoord)
	#- save the coordinates files
	newfilename=visualdir+'query_coord.txt'
	o1=open(newfilename,'w')
	for thisline in ncoordlines:
		o1.write('\t'.join([str(x) for x in thisline])+'\n')
	o1.close()
	#- save the a2b file
	newfilename2=visualdir+'query_nodes.txt'
	o2=open(newfilename2,'w')
	for thislineagain in addlines:
		o2.write('\t'.join([str(x) for x in thislineagain])+'\n')
	o2.close()





#- Send user message that the analysis was completed, could improve message to contain filenames
if slack_arg == True:
	slack.chat.post_message(user,'ANALYSIS COMPLETE')

