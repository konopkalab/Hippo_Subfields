library(made4)
library(pheatmap)

mat=read.table("DG_Matrix_module_correlation.txt")
mat$Rows=NULL
mat$Class=NULL

coolmod=read.table("coolmod.txt")
colnames(mat)=gsub("ME","",colnames(mat))

matF=mat[,colnames(mat)%in% coolmod$V1]

crosscor=cor(matF,method="pearson")
pdf("CrossCor_heatmap.pdf",width=5,height=5)
heatplot(crosscor,scaleKey=FALSE)
dev.off()