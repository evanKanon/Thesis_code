#!/usr/bin/env Rscript

# RDS type data extraction

# import necessary libraries
library(scater)  # install from Bioconductor: source("https://bioconductor.org/biocLite.R"); biocLite("scater")
#library(xgboost) # install from CRAN: install.packages("xgboost")
library(igraph)  # install from CRAN: install.packages("igraph")
library(tools)

# read in arguments
args<-commandArgs(trailingOnly = TRUE)
# set up error message if no argument gets passed
if (length(args)==0) {stop("Input filepath is a required argument",call.=FALSE)}
# save file path in variable
#iname="/Users/s1242130/Downloads/goolam.rds"
iname=args[1]

# read in input
source = readRDS(iname)
# extract relevant data
ds1 = t(exprs(source))
# create dataframe variable of data to ensure writability
datfrm=as.data.frame(ds1)


# keep only file name
#fname=file_path_sans_ext(basename(iname))

# keep file path with name but without extension (in order to change extension of new file)
fpath=file_path_sans_ext(iname)

# create new file path and extension
nfpath=paste(fpath,'.csv',sep="")

# write data in CSV file format
write.csv(datfrm,nfpath,quote=FALSE)

#sourceCellTypes = as.factor(colData(source)[,"cell_type1"])










