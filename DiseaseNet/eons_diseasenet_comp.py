

#compare EONS and DISEASEnet networks


import pandas as pd


def SEWD_similarity(mydataframe1,mydataframe2):
	''' Function which calculates Sum Edge Weight Distance similarity for dataframe matrix data '''
	#- subtract the two dataframes from each other
	diffdf=mydataframe1.sub(mydataframe2,fill_value=0)
	#- transfor each value to its absolute value
	diffdf=diffdf.transform(abs)
	#- sum them to calculate SEWD
	#- consider NAs as equal to zero
	sewd=diffdf.sum().sum()/2
	#sewd=diffdf.values.sum()/2
	#- calculate maximum possible difference
	dim=len(diffdf)
	maxdist=((dim*dim)-dim)/2    #- /2 is ommited because edge weights can take all values between -1 and 1, spanning 2 
	#- calculate similarity
	similarity=1-(sewd/maxdist)
	return similarity





kendall=pd.read_csv('/mnt/d/PhD/Results/Diseasenet_EONS/Diseasenetcomulative_GWAS_prepped.expression_kendall-9900.csv')

spearman=pd.read_csv('/mnt/d/PhD/Results/Diseasenet_EONS/Diseasenetcomulative_GWAS_prepped.expression_spearman-9900.csv')


mydnet2=pd.read_csv('/mnt/d/PhD/Results/DiseaseNet_v0.1.0_05102020_ldsc/network/result_df.csv')




#- clean up diseasenet name
def clean_d_name(string):
	new_string='Z'+string.split('/')[-1].split('.')[0]
	return new_string


mydnet2['p1']=mydnet2['p1'].apply(clean_d_name)
mydnet2['p2']=mydnet2['p2'].apply(clean_d_name)



#- turn disease net into matrix from ncol

name_list=mydnet2['p1'].unique()

df2=mydnet2.pivot(index='p1',columns='p2',values='rg')




































































































#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################




- add both name columns, find unique elements and create empty dataframe


#! add column selection feature if multiple columns available !!!!
def ncol_to_matrix(myfile,median_line=1):
	''' Function which converts ncol (a to b) file into matrix file '''
	#- read ncol file
	myncol=pd.read_csv(myfile)
	#- find element names
	name_list=myncol[myncol.keys()[0]]+myncol[myncol.keys()[1]]
	#- keep only unique element names and sort alphabetically
	name_list=sorted(name_list.unique())
	#- create empty pandas dataframe to store data
	df=pd.DataFrame(columns=name_list,index=name_list)
	#-- Add values to the corresponding cells of the matrix
	maxval=len(df)
	#- for each element in the empty dataframe find applicable value
	for i in range(len(df)):
		x=0
		while x!=maxval:
			#- identify name of items 
			a_it=df.keys()[i]
			b_it=df.index[x]
			value='nan'
			#- find corresponding value from ncol
			for i in range(len(myncol)):
				if myncol[myncol.keys()[0]][myncol.index[i]]==a_it and myncol[myncol.keys()[1]][myncol.index[i]]==b_it:
					value=myncol[myncol.keys()[2]][myncol.index[i]]
					break
				elif myncol[myncol.keys()[0]][myncol.index[i]]==b_it and myncol[myncol.keys()[1]][myncol.index[i]]==a_it:
					value=myncol[myncol.keys()[2]][myncol.index[i]]
					break
			if value=='nan':
				if a_it==b_it:
					value=median_line
			#if :
			#	value=
			#elif :
			#	value=
			#elif a_it==b_it:
			#	value=median_line
			#else:
			#	value='nan'
			#- add the value to the dataframe
			df[df.keys()[i]][df.index[x]]=value
			#- continue the loop
			x+=1
	#- return output
	return df