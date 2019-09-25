# WGCNA

library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(flashClust)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Loading Covariate
file <- list.files(pattern = "*_LMreg.txt")
tab <- read.table(file)
demo <- read.table("Demo_Full.txt",header=T,sep="\t",stringsAsFactors=TRUE)
suffix <- strsplit(file,"_")[[1]][1]
tmp <- demo[demo$Class == suffix,1:7]
rownames(tmp) <- tmp$Subject_ID
tmp$ID <- NULL
tmp$Subject_ID <- NULL
tmp$Diagnosis <- as.numeric(as.factor(tmp$Diagnosis))
tmp$Sex <- as.numeric(as.factor(tmp$Sex))
datTraits <- tmp

# Load tables 
colnames(tab) <- rownames(datTraits)
datExpr <- as.data.frame(t(tab));
names(datExpr) <- rownames(tab);
rownames(datExpr) <- names(tab);


## Powers analysis
powers = c(seq(2,30,2))
sft=pickSoftThreshold(datExpr,powerVector=powers,verbose = 5, blockSize= 18000, networkType = "signed") 
pdf("SoftThresholdingPower_signed.pdf")
par(mfrow = c(1,2), mar=c(5.1,5.1,4.1,2.1));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.5,col="red"); abline(h=0.8,col="blue");abline(h=0.9,col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


## Clustering
A=adjacency(t(datExpr),type="distance") 
k=as.numeric(apply(A,2,sum))-1
Z.k=scale(k)
thresholdZ.k=-5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average") 
traitColors=data.frame(numbers2colors(datTraits,signed=TRUE))
dimnames(traitColors)[[2]]=paste(colnames(datTraits),sep="") 
datColors=data.frame(outlierC=outlierColor,traitColors)
pdf("SampleClustering.pdf")
plotDendroAndColors(sampleTree,groupLabels=names(datTraits), colors=datColors,main="Sample dendrogram and trait heatmap")
dev.off()

######################################################################################################################
############################################ MODULE construction #####################################################
############################################    signed network   #####################################################
PWR=20
net = blockwiseModules(datExpr,corType="bicor", maxBlockSize = 18000, networkType="signed",power=PWR, minModuleSize=100,nThreads=15,
TOMType = "signed",TOMDenom = "mean",deepSplit=2,verbose=5,mergeCutHeight=0.15,detectCutHeight = 0.999,reassignThreshold = 1e-10,numericLabels=TRUE,saveTOMs=FALSE,pamStage=TRUE, pamRespectsDendro=TRUE)
moduleLabelsAutomatic=net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
MEsAutomatic=net$MEs
unique(moduleColorsAutomatic)
table(moduleColorsAutomatic)
write.table(moduleColorsAutomatic, "DG_colors.txt",sep="\t",quote=F)
save(net,file="DG_NetData.RData")

#KMEs
KMEs<-signedKME(datExpr, net$MEs,corFnc = "cor", corOptions = "use = 'p'")
kme=data.frame(rownames(tab), moduleColorsAutomatic, KMEs)
colnames(kme)[1]="Symbol"
rownames(kme)=NULL
write.table(kme,"KME_DG.txt",sep="\t",quote=F)

# Dat Trait Heatmap
moduleColorsIEGG=moduleColorsAutomatic
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColorsIEGG)$eigengenes
MEsIEGG = MEs0
MEsIEGG$MEgrey=NULL
modTraitCor= cor(MEsIEGG, datTraits,method="spearman")
write.table(modTraitCor,"modTraitCor_DG.txt",sep="\t",quote=F)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
write.table(modTraitP,"modTraitP_DG.txt",sep="\t",quote=F)
textMatrix = paste(signif(modTraitCor, 2), "\n(",signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor) 
par(mar = c(6, 8.5, 3, 3))
pdf("Heatmap_DatTraits.pdf",width=7,height=7)
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits), yLabels = names(MEsIEGG), ySymbols = names(MEsIEGG), colorLabels =FALSE,colors=blueWhiteRed(50),textMatrix=textMatrix, setStdMargins = FALSE, cex.text = 0.4, zlim = c(-1,1), main = paste("Module Association"))
dev.off()

#EigenGeneNet
MEList=moduleEigengenes(datExpr,colors=moduleColorsAutomatic,softPower = PWR,impute = TRUE)
MEs = MEList$eigengenes
MET=orderMEs(MEs)
pdf("EigengeneNetworks.pdf")
plotEigengeneNetworks(MET,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.5,xLabelsAngle=90)
dev.off()

#GGplotInput
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColorsAutomatic,softPower = PWR,impute = TRUE)$eigengenes
MEs0$Rows=colnames(tab)
MEs0$Class=paste("Class",datTraits$Diagnosis,sep="")
write.table(MEs0, "DG_Matrix_module_correlation.txt",sep="\t",quote=F)

# Adjacency matrix
Adj = adjacency(datExpr, power = PWR,type="signed",corFnc = "cor", corOptions = "use = 'p'")
moduleOutput <- data.frame(rownames(tab))
moduleOutput[,2]<- moduleColorsAutomatic
intraCon <- intramodularConnectivity(Adj, moduleColorsAutomatic)
moduleOutput[,3]<-intraCon$kWithin
colnames(moduleOutput) <- c("Gene", "ModuleColor", "kWithin")
write.table(moduleOutput, "ModuleOutput_DG.txt", sep="\t", quote=F)

# TO connectivity (you need the table as single gene list column in the directory)
TOM = TOMsimilarityFromExpr(datExpr, power= PWR,corType = "bicor",networkType="signed",TOMType="signed",TOMDenom = "mean",nThreads = 15,verbose = 5, indent = 0)
colnames(TOM)=rownames(TOM)=colnames(datExpr)
save(TOM,file="TOM_ourRNA_SIGNED.RData")
Connectivity=apply(TOM,1,sum)
save(Connectivity,file="Connectivity.RData")

# CytoScape output
dir.create("Cyto")
setwd("Cyto")
for(module in unique(moduleColorsAutomatic)){
inModule <- is.finite(match(moduleColorsAutomatic, module))
modTOM <- TOM[inModule, inModule]
cyt = exportNetworkToCytoscape(modTOM, edgeFile=paste("CytoEdge",paste(module,collapse="-"),".txt",sep=""), nodeFile=paste("CytoNode",paste(module,collapse="-"),".txt",sep=""), weighted = TRUE, threshold = 0, nodeAttr = moduleColorsAutomatic[inModule], nodeNames = names(datExpr)[inModule])
}

# WGCNA barplot output
library(ggplot2)
library(reshape2)
library(RColorBrewer)
df=melt(MEs0)
dir.create("Plot")
setwd("Plot")
df$Rows=factor(df$Rows, levels = colnames(tab))
df=df[!df$variable == "MEgrey", ]
# Function to make the barplot and saving as pdf according to the module name
cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
cols <- sample(cl, 7)
doPlot = function(sel_name) 
{
    df = subset(df, variable == sel_name)
    PLOT=ggplot(data=df, aes(x=Rows, y=value)) +
 geom_bar(aes(fill=Class),stat="identity",position = "dodge")+
 scale_y_continuous(limits=c(-1,+1))+
 theme_bw()+
 theme(strip.text.x = element_text(size=12, face="bold"),
 strip.background = element_rect(colour="black", fill="#CCCCFF"))+
 scale_fill_manual(values = cols)+
 theme(axis.title.x = element_blank(),
            axis.text.x  = element_text(face="bold", size=6,angle = 45, hjust = 1))+
 theme(axis.title.y = element_blank(),
            axis.text.y  = element_text(face="bold", size=6))
    print(PLOT)
    ggsave(sprintf("%s.pdf", sel_name))
 }

lapply(unique(df$variable), doPlot)




