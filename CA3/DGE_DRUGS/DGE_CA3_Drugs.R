
#' Load Libraries
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(pheatmap))
source("Utils.R")

file <- list.files(pattern = "log2CPM.txt")
Count <- read.table(file,header=T)
colnames(Count) <- gsub("X0","0-",colnames(Count))

demo <- read.table("Demo_DRUGS.txt",header=T,sep="\t")
suffix <- strsplit(file,"_")[[1]][1]
pd <- demo[demo$Class == suffix & demo$Diagnosis == "SZ",1:8]
vec <- c("0-869","0-870","0-871","0-872","0-873","0-882","0-883","0-884","0-855","0-886")
pd$Batch <- ifelse(pd$ID %in% vec,"Batch1","Batch2")
pd$Diagnosis <- NULL
Count <- Count[,match(pd$ID,colnames(Count))]
colnames(Count) <- pd$Subject_ID
filter=apply(Count, 1, function(x) (all(x[1:6] > 0) | all(x[7:13] > 0)))
dat <- Count[filter,]

# Format pheno for DESeq2
pd$ID <- NULL
rownames(pd) <- pd$Subject_ID
pd$Subject_ID <- NULL
pd <- pd[c(5,1,2,3,4,6)]


# Quantile Normalization
p <- normalize.quantiles(as.matrix(dat))
rownames(p) <- rownames(dat)
colnames(p) <- colnames(dat)
write.table(p, file = paste(suffix, "_QuantNorm", ".txt",sep = ""),sep="\t",quote=F)

# Variance explained
var <- VarExp(p,pd,5,FALSE)
pdf("CA3_Variance_Explained_DRUGS.pdf",width=4,height=4)
plotVarExp(var,"CA3 DRUGS VarExp")
dev.off()

#' Register cluster
cl <- makeCluster(3)
registerDoParallel(cl)
form <- ~ (1|Drugs) + (1|Sex) + Age + RIN + PMI
varPart <- fitExtractVarPartModel(p, form, pd)
vp <- varPart[order(varPart$Diagnosis, decreasing=TRUE),]

#' Plot Variance
plotPercentBars(vp[1:20,]) + 
theme_classic() +
ggsave("CA3_Top8genes_VarExp_DRUGS.pdf", width=6, height=8)

plotVarPart(vp) +
theme_classic() +
ggsave("CA3_AllGene_Explained_DRUGS.pdf", width=6, height=8)

# SVA calculation
mod = model.matrix(~Drugs+Sex+Age+PMI+RIN+Batch, data =pd)
mod0 <- model.matrix(~1, data=pd)
NSV <- num.sv(as.matrix(p),mod,method="be")
svaobj = sva(as.matrix(p),mod, n.sv=2,B=20)

#Clean with LM regression
svs=svaobj$sv
colnames(svs)=paste("SV",1:ncol(svs),sep="")
pd_sva=cbind(pd[c(-1)],svs)
avebeta.lm<-lapply(1:nrow(p), function(x){
  lm(unlist(p[x,])~.,data=pd_sva)
})
residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
adj.residuals<-residuals+matrix(apply(p, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(adj.residuals)<-rownames(p)
rownames(residuals) <- rownames(p)
write.table(adj.residuals, file = paste(suffix, "_QuantNorm_LMreg_DRUGS", ".txt",sep = ""),sep="\t",quote=F)

#' Add SVs to the Model Matrix
svs=svaobj$sv
colnames(svs)=paste("SV",1:ncol(svs),sep="")
TRAITSr=cbind(pd,svs)

Data=t(p)
output <- data.frame(matrix(nrow=ncol(Data), ncol=3, dimnames = list(colnames(Data), c("logFC", "Pval", "Warning"))))
output[,] <- NA

#" This is the formula I am using
print(formula(paste("Data[,i] ~ ", paste(colnames(TRAITSr),collapse = " + "))))

#' Looping on the Data
for (i in 1:ncol(Data)) {	
			Model=tryCatch(lm(as.formula(paste("Data[,i] ~ ", paste(colnames(TRAITSr),collapse = " + "))), data = TRAITSr),warning =  function(w) w)
				if (i %% 1000 == 0) {cat(paste("Done on gene ",i,"\n",sep=""))}
				if(typeof(Model) == "list"){
    				coefs = data.frame(coef(summary(Model)))
    				t_value = coefs["Drugs", "t.value"]
    				output[i,"Pval"] = 2 * (1 - pnorm(abs(t_value)))
    				output[i,"logFC"]= coefs["Drugs", "Estimate"]
    				} else {
    				output[i,"Warning"] = as.character(Model)
    				output[i, "logFC"] = 0
    				output[i,"Pval"] = 1
  			}
		}

output$FDR <- p.adjust(output$Pval,"BH")
sign <- output[output$FDR< 0.05,]
write.table(output, file = paste(suffix, "_LM_OUTPUT_DRUGS", ".txt",sep = ""),sep="\t",quote=F)
write.table(sign, file = paste(suffix, "_LM_OUTPUT_DRUGS_SIGN", ".txt",sep = ""),sep="\t",quote=F)

# Visualization DGE
mat=adj.residuals[rownames(adj.residuals) %in% rownames(sign),]
pdf("CA3_DRUGS_Heatmap.pdf",width=6,height=6)
pheatmap(mat,annotation=pd,cluster_cols=T,scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),clustering_distance_cols = "correlation",show_rownames = F,show_colnames = F,fontsize = 8,legend = TRUE,annotation_legend = TRUE)
dev.off()
