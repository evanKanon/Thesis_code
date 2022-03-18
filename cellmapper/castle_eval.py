#- create network json file


import pandas as pd

results=pd.read_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/output/castle_comparison/GSE59114_DBA_stdt_13082019_CellMv0.1.0/results/top_result_genes.txt',sep='\t')
results=pd.read_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/output/castle_comparison/GSE59114_C57BL6_stdt_13082019_CellMv0.1.0/results/top_result_genes.txt',sep='\t')
results=pd.read_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/output/castle_comparison/goolam_stdt_13082019_CellMv0.1.0/results/top_result_genes.txt',sep='\t')
results=pd.read_csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/output/castle_comparison/deng_stdt_13082019_CellMv0.1.0/results/top_result_genes.txt',sep='\t')




query_list=results['0']
query_list=query_list.unique()

#- fraction uniform data into seperate sets based on similarity of disease
association_res=[]
assoc=0
disossiation_res=[]
disos=0
for cell in query_list:
	#- select subset with given celltype
	thisdf=results[results['0']==cell]
	#- sort based on how good the association between the two samples is
	thisdf=thisdf.sort_values('2',ascending=False)
	#- reset the index
	thisdf=thisdf.reset_index(drop=True)
	#- take top result
	mypair=thisdf.head(1)
	#- select the string with the first term
	mystring1=mypair['0'][0]
	#- replace double and single quotes if present
	mystring1=mystring1.replace('"', '')
	mystring1=mystring1.replace("'", '')
	#- select the string with the second term
	mystring2=mypair['1'][0]
	#- replace double and single quotes if present
	mystring2=mystring2.replace('"', '')
	mystring2=mystring2.replace("'", '')
	#- split the two terms into strings with their constitutent parts
	stringlist1=mystring1.split('_')
	stringlist2=mystring2.split('_')
	#- delete the string 'DBA' if present in the list
	for j in range(len(stringlist1)):
		if stringlist1[j]=='DBA':
			del stringlist1[j]
			break
	#- delete the term 'population' if present in the list 
	for j in range(len(stringlist1)):
		if stringlist1[j]=='population':
			del stringlist1[j]
			break
	#- delete the term if it is a digit
	for j in range(len(stringlist1)):
		if stringlist1[j].isdigit():
			del stringlist1[j]
			break
	for i in range(len(stringlist2)):
		if stringlist2[i]=='DBA':
			del stringlist2[i]
			break
	for i in range(len(stringlist2)):
		if stringlist2[i]=='population':
			del stringlist2[i]
			break
	for i in range(len(stringlist2)):
		if stringlist2[i].isdigit():
			del stringlist2[i]
			break
	mystring1=''.join(stringlist1).lower()
	mystring2=''.join(stringlist2).lower()
	stringlist1=mystring1.split(' ')
	stringlist2=mystring2.split(' ')
	mystring1=''.join(stringlist1).lower()
	mystring2=''.join(stringlist2).lower()
	stringlist1=mystring1.split('-')
	stringlist2=mystring2.split('-')
	mystring1=''.join(stringlist1).lower()
	mystring2=''.join(stringlist2).lower()
	if mystring1[:6]==mystring2[:6]:
		myasc=[mystring1,mystring2,'same']
		association_res.append(myasc)
		assoc+=1
	else:
		mydif=[mystring1,mystring2,'diff']
		disossiation_res.append(mydif)
		disos+=1


association_res=[]
assoc=0
disossiation_res=[]
disos=0
for cell in query_list:
	#- select subset with given celltype
	thisdf=results[results['0']==cell]
	#- sort based on how good the association between the two samples is
	thisdf=thisdf.sort_values('2',ascending=False)
	#- reset the index
	thisdf=thisdf.reset_index(drop=True)
	#- take top result
	mypair=thisdf.head(1)
	#- select the string with the first term
	mystring1=mypair['0'][0]
	#- replace double and single quotes if present
	mystring1=mystring1.replace('"', '')
	mystring1=mystring1.replace("'", '')
	#- select the string with the second term
	mystring2=mypair['1'][0]
	#- replace double and single quotes if present
	mystring2=mystring2.replace('"', '')
	mystring2=mystring2.replace("'", '')
	#- set string with the appropriate name
	if '2cell' in mystring1:
		mystring1='2cell'
	elif '4cell' in mystring1:
		mystring1='4cell'
	elif '8cell' in mystring1:
		mystring1='8cell'
	elif '16cell' in mystring1:
		mystring1='16cell'
	#- set string with the appropriate name
	if '2cell' in mystring2:
		mystring2='2cell'
	elif '4cell' in mystring2:
		mystring2='4cell'
	elif '8cell' in mystring2:
		mystring2='8cell'
	elif '16cell' in mystring2:
		mystring2='16cell'
	if mystring1[:6]==mystring2[:6]:
		myasc=[mystring1,mystring2,'same']
		association_res.append(myasc)
		assoc+=1
	else:
		mydif=[mystring1,mystring2,'diff']
		disossiation_res.append(mydif)
		disos+=1














