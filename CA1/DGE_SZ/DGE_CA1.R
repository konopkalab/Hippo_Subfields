
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

# Load the files
file <- list.files(pattern = "_log2CPM.txt")
Count <- read.table(file,header=T)
colnames(Count) <- gsub("X0","0-",colnames(Count))

demo <- read.table("Demo_Full.txt",header=T,sep="\t")
suffix <- strsplit(file,"_")[[1]][1]
pd <- demo[demo$Class == suffix,1:7]
Count <- Count[,match(pd$ID,colnames(Count))]
colnames(Count) <- pd$Subject_ID
filter=apply(Count, 1, function(x) (all(x[1:13] > 0) | all(x[14:26] > 0)))
dat <- Count[filter,]

# Format pheno for DESeq2
pd$ID <- NULL
rownames(pd) <- pd$Subject_ID
pd$Subject_ID <- NULL

# Quantile Normalization
p <- normalize.quantiles(as.matrix(dat))
rownames(p) <- rownames(dat)
colnames(p) <- colnames(dat)
write.table(p, file = paste(suffix, "_QuantNorm", ".txt",sep = ""),sep="\t",quote=F)


# Variance explained
var <- VarExp(p,pd,5,FALSE)
pdf("CA1_Variance_Explained.pdf",width=4,height=4)
plotVarExp(var,"CA1 VarExp")
dev.off()


#' Register cluster
cl <- makeCluster(3)
registerDoParallel(cl)
form <- ~ (1|Diagnosis) + (1|Sex) + Age + RIN + PMI
varPart <- fitExtractVarPartModel(p, form, pd)
vp <- varPart[order(varPart$Diagnosis, decreasing=TRUE),]

#' Plot Variance
plotPercentBars(vp[1:8,]) + 
theme_classic() +
ggsave("CA1_Top8genes_VarExp.pdf", width=6, height=5)

plotVarPart(vp) +
theme_classic() +
ggsave("CA1_AllGene_Explained.pdf", width=6, height=5)

# SVA calculation
mod = model.matrix(~Diagnosis+Sex+Age+PMI+RIN, data =pd)
mod0 <- model.matrix(~1, data=pd)
NSV <- num.sv(as.matrix(p),mod,method="be")
svaobj = sva(as.matrix(p),mod, n.sv=NSV,B=20)

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
write.table(adj.residuals, file = paste(suffix, "_QuantNorm_LMreg", ".txt",sep = ""),sep="\t",quote=F)

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
    				t_value = coefs["Diagnosis", "t.value"]
    				output[i,"Pval"] = 2 * (1 - pnorm(abs(t_value)))
    				output[i,"logFC"]= coefs["Diagnosis", "Estimate"]
    				} else {
    				output[i,"Warning"] = as.character(Model)
    				output[i, "logFC"] = 0
    				output[i,"Pval"] = 1
  			}
		}

output$FDR <- p.adjust(output$Pval,"BH")
sign <- output[output$FDR< 0.05,]
write.table(output, file = paste(suffix, "_LM_OUTPUT", ".txt",sep = ""),sep="\t",quote=F)
write.table(sign, file = paste(suffix, "_LM_OUTPUT_SIGN", ".txt",sep = ""),sep="\t",quote=F)

# Visualization DGE
mat=adj.residuals[rownames(adj.residuals) %in% rownames(sign),]
pdf("CA1_Heatmap.pdf",width=6,height=6)
pheatmap(mat,annotation=pd,cluster_cols=T,scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),clustering_distance_cols = "correlation",show_rownames = F,show_colnames = F,fontsize = 8,legend = TRUE,annotation_legend = TRUE)
dev.off()
