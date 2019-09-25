# XY plot cor diagnosis
library(ggplot2)
library(ggrepel)

files=list.files(pattern="modTrait*")
myfiles=lapply(files,read.table,sep="\t")
df=data.frame(Mod=gsub("ME","",rownames(myfiles[[1]])),Cor=myfiles[[1]]$Diagnosis,P=myfiles[[2]]$Diagnosis)
df$log=-log10(df$P)
df$abs=abs(df$Cor)

col=df$Mod


pdf("XY_PLOT_CorMod_Diagnosis_NeuN.pdf",width=5,height=5.5,useDingbats=FALSE)
ggplot(df, aes(x=Cor, y=log,size=abs)) +
geom_point(aes(fill=Mod),colour="black",pch=21)+
theme_classic()+
geom_hline(yintercept = 1.3,color="red",linetype ="dashed",alpha=0.5,size=1)+
scale_fill_manual(values=c(paste(dput(as.character(col)))))+
geom_text_repel(aes(Cor, log, label = df$Mod))+
xlim(-1,1)+
ylim(0,2.5)+
theme(legend.position="none")+
ggtitle("DG Schizophrenia Modueles")+
xlab("Pearson Correlation")+ 
ylab("-log10(p-value)")
dev.off()


