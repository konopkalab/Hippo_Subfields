### Merging multiple dataframe by rownames
## list of tables with the gene counts
load("GeneSets_Disorders.RData")
load("GeneSets_SingleCell.RData")

file=list.files(pattern=".txt")
mod=read.table(file)
mod$Gene <- rownames(mod)
GeneSets$Mod=mod

myfiles=c(GeneSets,GeneSetsSingleCell)
res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, myfiles)

database=res[!(is.na(res$logFC)),]
database$Warning <- NULL
vec=c("Gene", "logFC","Pval","FDR","ASD", "ASD_SC", "ASD12", "ASD16","EPILEPSY","FMRP","ID","NEURO_DEV_RISK","RBP","SYN","SZ_LOCI","SZ_Sklar","SZ_Genes","TF","ENDO","EPEND_ZS","EPITH_ZS","ASTRO","ASTRO_ZS","NEURO_CA1_ZS","NEURO_SA1_ZS","INTERN_ZS","NEURO","LAKE_EX","LAKE_IN","MICRO","MICROG_ZS","MYEL","OLIG","OLIG_ZS")
                  
database=database[,match(vec,names(database))]
 
## write out to XLSX
openxlsx::write.xlsx(database, file = "CA1_DATABASE_DRUGS_DGE.xlsx", colNames = TRUE, borders = "columns")

