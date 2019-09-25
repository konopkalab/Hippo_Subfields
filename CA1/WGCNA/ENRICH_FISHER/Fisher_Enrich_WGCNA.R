library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
enableWGCNAThreads()

wgcna = list.files(pattern = 'ModuleOutput*')
tab=read.table(wgcna,sep="\t",header=T)
colnames(tab)=c("Gene","DEFINITION","Kwithin")
Genes=as.data.frame(table(tab$DEFINITION))

# Loop to make the overlap
# The loop goes along the length of the GeneSets lists
load("GeneSets_Disorders.RData")
#for(i in 1:length(GeneSets))
#{
#	GeneSets[[i]] <- GeneSets[[i]][GeneSets[[i]]$Gene %in% tab$Gene,]
#}
ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)

#
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
INT[[i]]$d <- 20268-Genes$Freq-nrow(GeneSets[[i]])
}

# sum(Genes$Freq)
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value), 
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat,file="Module_Enrich_Fisher.RData")

# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=12)
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0
pdf("Heatmap_GeneSets_DG_Fisher_coolmod_adj.pdf",width=7,height=4)
colnames(FisherAdj)=c("ASD","ASD_SC","CA1_DGE","CA1_DOWN","CA1_UP","FMRP","ID","RBP","SYN","SZ_Loci","FromerEtAl2016","TF")
colnames(FisherOR)=c("ASD","ASD_SC","CA1_DGE","CA1_DOWN","CA1_UP","FMRP","ID","RBP","SYN","SZ_Loci","FromerEtAl2016","TF")
coolmod=read.table("coolmod.txt")
FisherPf=FisherAdj[rownames(FisherAdj) %in% coolmod$V1,]
FisherORf=FisherOR[rownames(FisherOR) %in% coolmod$V1,]
df=-log10(FisherPf)
LabelMat = paste(signif(FisherORf, 2), "\n(",signif(FisherPf, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, xLabels = colnames(df), yLabels = rownames(df), colorLabels =FALSE,colors=colorRampPalette(c("white", "darkgrey"))(50),textMatrix=LabelMat, setStdMargins = FALSE, cex.text = 0.5, main = paste("CA1 Modules GeneSets Enrichment"))
dev.off()


# For Single cell
rm(list=ls())
library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
enableWGCNAThreads()

wgcna = list.files(pattern = 'ModuleOutput*')
tab=read.table(wgcna,sep="\t",header=T)
colnames(tab)=c("Gene","DEFINITION","Kwithin")
Genes=as.data.frame(table(tab$DEFINITION))

# Loop to make the overlap
# The loop goes along the length of the GeneSets lists
load("GeneSets_SingleCell.RData")
GeneSets <- GeneSetsSingleCell
#for(i in 1:length(GeneSets))
#{
#	GeneSets[[i]] <- GeneSets[[i]][GeneSets[[i]]$Gene %in% tab$Gene,]
#}
ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)

#
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
INT[[i]]$d <- 20268-Genes$Freq-nrow(GeneSets[[i]])
}


RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value), 
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat,file="Module_Enrich_Fisher.RData")

# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}

rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
FisherOR <- do.call(cbind,tmp)
}

rownames(FisherOR) <- rowNames
colnames(FisherOR) <- gsub("_NAMED","",colNames)


# Adjust FDR
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=16)
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- gsub("_NAMED","",colNames)

vec=c("ENDO","EPEND_ZS","EPITH_ZS","ASTRO","ASTRO_ZS","CA1_ZS","SA1_ZS","INTERN_ZS","NEURO","LAKE_EX","LAKE_IN","MICRO","MICRO_ZS","MYEL","OLIG","OLIG_ZS")

FisherAdj=FisherAdj[,match(vec,colnames(FisherAdj))]
FisherOR=FisherOR[,match(vec,colnames(FisherOR))]

# transform for visualization
FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0

pdf("Heatmap_CellType_DG_Fisher_coolmod.pdf",width=7,height=4)
coolmod=read.table("coolmod.txt")
FisherAdjF=FisherAdj[rownames(FisherAdj) %in% coolmod$V1,]
FisherORf=FisherOR[rownames(FisherOR) %in% coolmod$V1,]
df=-log10(FisherAdjF)
LabelMat = paste(signif(FisherORf, 2), "\n(",signif(FisherAdjF, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, xLabels = colnames(df), yLabels = rownames(df), colorLabels =FALSE,colors=colorRampPalette(c("white", "darkgrey"))(50),textMatrix=LabelMat, setStdMargins = FALSE, cex.text = 0.5, main = paste("CA1 Modules Cell Type Enrichment"))
dev.off()













