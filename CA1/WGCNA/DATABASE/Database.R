### Merging multiple dataframe by rownames
## list of tables with the gene counts
load("GeneSets_Disorders.RData")
load("GeneSets_SingleCell.RData")

file=list.files(pattern=".txt")
mod=read.table(file)

GeneSets$Mod=mod

myfiles=c(GeneSets,GeneSetsSingleCell)
res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, myfiles)

database=res[!(is.na(res$kWithin)),]

vec=c("Gene", "ModuleColor","kWithin","ASD", "ASD_SC","UP_CA1","DOWN_CA1","DGE_CA1","FMRP","ID","RBP","SYN","SZ_LOCI","SZ_Sklar","TF","ENDO","EPEND_ZS","EPITH_ZS","ASTRO","ASTRO_ZS","NEURO_CA1_ZS","NEURO_SA1_ZS","INTERN_ZS","NEURO","LAKE_EX","LAKE_IN","MICRO","MICROG_ZS","MYEL","OLIG","OLIG_ZS")
                  
database=database[,match(vec,names(database))]
 
## write out to XLSX
library(xlsx)
write.xlsx(database, file="CA1_DATABASE_WGCNA.xlsx",sheetName = "WGCNA_CA1",row.names=FALSE)
