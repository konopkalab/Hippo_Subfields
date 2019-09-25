# XY plot cor diagnosis
library(ggplot2)
library(ggrepel)
library(cowplot)

files=list.files(pattern="Module_*")
myfiles=lapply(files,read.table,sep="\t",header=T)

cool <- read.table("coolmod.txt")

df=myfiles[[1]][myfiles[[1]]$Rows %in% cool$V1,]
col=df$Rows
first <- ggplot(df, aes(x=SIZE, y=value,size=value)) +
geom_point(aes(fill=Rows),colour="black",pch=21)+
theme_classic()+
geom_hline(yintercept = 2,linetype="dashed",colour="red")+
geom_hline(yintercept = 10,linetype="dashed",colour="blue")+
scale_fill_manual(values=c(paste(dput(as.character(col)))))+
geom_text_repel(aes(SIZE, value, label = df$Rows))+
xlim(0,750)+
ylim(0,40)+
theme(legend.position="none")+
ggtitle("CA1 vs CA3")+
xlab("Module Size")+ 
ylab("Preservation Z-Score")


df=myfiles[[2]][myfiles[[2]]$Rows %in% cool$V1,]

col=df$Rows
second <- ggplot(df, aes(x=SIZE, y=value,size=value)) +
geom_point(aes(fill=Rows),colour="black",pch=21)+
theme_classic()+
geom_hline(yintercept = 2,linetype="dashed",colour="red")+
geom_hline(yintercept = 10,linetype="dashed",colour="blue")+
scale_fill_manual(values=c(paste(dput(as.character(col)))))+
geom_text_repel(aes(SIZE, value, label = df$Rows))+
xlim(0,750)+
ylim(0,40)+
theme(legend.position="none")+
ggtitle("CA1 vs DG")+
xlab("Module Size")+ 
ylab("Preservation Z-Score")

plot <- plot_grid(first, second,NULL,NULL,labels=c("A", "B","C","D"), ncol = 2,align = "v")
save_plot("ModulePreservation_DG.pdf", plot, ncol = 2,base_height=6,base_width=3.5)











