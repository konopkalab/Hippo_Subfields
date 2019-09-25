
library(GenomicRanges)
library(rtracklayer)
library(LOLA)
library(ggplot2)
library(ggrepel)

files=list.files(pattern="ModuleOutput_*")
myfiles=as.data.frame(lapply(files,read.table,sep="\t")[[1]])
coolmod=read.table("coolmod.txt")
mod=myfiles[myfiles$ModuleColor %in% coolmod$V1,]

totGTF <- import.bed(con="Homo_sapiens.GRCh37.87_Protein_Coding_GENEONLY.bed")

userSets <- list()
indbgSets <- list()
locResults <- list()
tab <- list()
for (i in unique(mod$ModuleColor))
{ 
userSets[[i]] <- totGTF[totGTF$name %in% mod[mod$ModuleColor == i,]$Gene]
indbgSets[[i]] <- totGTF[!(totGTF$name %in% mod[mod$ModuleColor == i,]$Gene)]
dbPath = system.file("extdata", "disorders3", package="LOLA")
regionDB = loadRegionDB(dbPath)
locResults[[i]] = runLOLA(userSets[[i]], indbgSets[[i]], regionDB, cores=3,minOverlap = 3)
tab[[i]]=as.data.frame(locResults[[i]])
tab[[i]]$qValue=NULL
write.table(tab[[i]], file = paste("GWAS_ENRICH_", i, ".txt",sep = ""),sep="\t",quote=F)
tab[[i]]$Class <- rep(i,nrow(tab[[1]]))
df <- do.call(rbind,tab)
df$adj=p.adjust(10^-df$pValueLog,"BH")
df$log=-log10(df$adj)
}


pdf("GWAS_ENRICH_COOLMOD.pdf",width=6,height=5,useDingbats=FALSE)
ggplot(df, aes(x = log, y = logOddsRatio, group=Class, size = logOddsRatio,label=gsub(".txt|.bed","",df$filename),color=gsub(".txt|.bed","",df$filename))) +
geom_point(shape=21) + 
geom_text_repel(size = 2, colour = "black") + 
theme_classic() + # clean up theme
theme(legend.position = "none") + # remove legend
xlab(expression(-log[10]("FDR"))) + # x-axis label
ylab(expression(log[2]("OR")))+
xlim(0,20)+
ylim(0,5)+
geom_vline(xintercept = 1.3, linetype="dotted", color = "red", size=1) + theme(legend.position="none")+
facet_wrap(~ Class,ncol=3,nrow=2)
dev.off()










