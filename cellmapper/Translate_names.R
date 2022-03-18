#Ortholog name selection script
#org.Mm.eg.db
#org.Hs.eg.db
#hom.Mm.inp.db
#org.Mm.egHGNC

options(stringsAsFactors = FALSE)

library(biomaRt)

#hend_data=read.csv('/Users/s1242130/Desktop/PhD/Projects/Singlecell/Differential_expression/query_cluster_sum_expression.csv')

murine_data=read.csv('/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/data/Murine_data/Henderson_mouse_liver_tissue_cleansed.csv')
#name_vector=murine_data['X']
name_vector2=murine_data[,1]

#- can't make it work yet
#res_genes=idConverter(name_vector,srcSpecies = "MUSMU",srcIDType = "GENENAME",destIDType = "GENENAME",destSpecies = "HOMSA")
#tes_genes=idConverter(name_vector,srcSpecies = "MUSMU",srcIDType = "GENENAME",destIDType = "GENENAME",destSpecies = "HOMSA")

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#- function which finds name of gene for particular input and appends human equivalent to specific dataframe column
f<-function(x){
  musname<-as.character(x[1])
  humname<-getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = musname , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  x$human_name<-as.character(humname[,2])
}

#- (dataset_to_apply_function, 1_rows/2_columns/c(1,2)_entire dataset, which function to apply)
murine_data$human_name<-apply(murine_data, 1, f)

write.csv(murine_data, file = '/mnt/ris-fas1a/linux_groups2/baillie_grp/singlecell/Data/Murine_data/human_comp_murine_data.csv', row.names = FALSE, quote = FALSE)
