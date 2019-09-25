# Make Boxplot 
library(ggplot2)
library(reshape2)

tab=read.table("DG_Matrix_module_correlation.txt")
tab$Class=ifelse(tab$Class == "Class2","SZ","CT")
df=melt(tab)
df$variable=gsub("ME","",df$variable)

coolmod=read.table("coolmod.txt")
df=df[df$variable %in% coolmod$V1,]

pdf("Module_Eigengene_CoolMod_DG.pdf",width=4.5,height=2.5,useDingbats=FALSE)
ggplot(df, aes(x=variable, y=value, fill=Class,group=Class)) +
geom_pointrange(mapping = aes(x = variable, y = value,group=Class,colour=Class), position=position_dodge(width=0.5),
stat = "summary",
fun.ymin = function(z) {quantile(z,0.25)},
fun.ymax = function(z) {quantile(z,0.75)},
fun.y = median)+
theme_classic()+
geom_hline(yintercept = 0,linetype="dashed")+
scale_colour_manual(values=c("steelblue","gold"))+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
labs(title="DG Module EigenGene",x="", y = "Eigengene")
dev.off()

