### Merging multiple dataframe by rownames
## list of tables with the gene counts
files = list.files(pattern = '*.txt')
names <- gsub( "*.txt", "", files )
GeneSets = lapply(files, read.table,header=T,sep="\t")
names(GeneSets)=names
save(GeneSets,file="GeneSets_Disorders.RData")
