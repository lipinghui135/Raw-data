library(Seurat)
library(tidyverse)
library(sctransform)
library(harmony)
###############################################################################
mydata = SCTransform(mydata, vars.to.regress=c("nCount_RNA", "percent.mito", "S.Score", "G2M.Score"), verbose=FALSE)
mydata = RunPCA(mydata, verbose=FALSE)
mydata = RunHarmony(mydata, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(mydata)
mydata <- FindNeighbors(mydata, dims=1:20, reduction="harmony")
mydata <- RunUMAP(mydata, dims=1:20, reduction="harmony")
colors = c4a("brewer.set2", 5)
p = DimPlot(mydata, reduction="umap", group.by="Sample", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("UMAP_Sample.pdf", p, height=6, width=6)
######################################################### resolution=0.1
mydata <- FindClusters(mydata, resolution=0.1)
UMAPPlot(mydata, pt.size=1, label=T, cols=colors, label.size=5)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "./All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
############## 
VlnPlot(mydata, features=c("S100A9"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())

# MF1: EGFL6, NPNT, ELN, SNCG, CTXN1
# MF2: LPL, PIEZO2, CYGB, PROCR, EBF2, MARCKSL1, IFIT2, RFTN1, AGTR1, OAF
# MF3: CPB2, C5, ARG1, F2, ADH4, APOM, ALDH6A1, BAAT, ORM2, HMGCS2
# MF4: EFEMP1, CTHRC1, LUM, LOX, THBS2, FMOD, VCAN
# MF5: FCGR3A, GZMB, CD3D, TRAC

# 
cell_label = c("MF1", "MF2", "MF3", "MF4", "MF5")
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
p = UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)
ggsave("UMAP_subcluster.pdf", p, width=6, height=6)

genes = c("EGFL6", "NPNT", "ELN", "SNCG", "LPL", "PIEZO2", "CYGB", "AGTR1", "CPB2", "C5", "ARG1", "ADH4", "APOM", "EFEMP1", "CTHRC1", "LUM", "LOX", "FCGR3A", "GZMB", "CD3D", "TRAC")
DotPlot(mydata, features=genes)+coord_flip()+scale_color_distiller(palette="RdYlBu")+theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())


# marker基因的UMAP映射图
set = c("CD24", "CD19", "CCR7", "CD3D", "HBB", "CST3", "CCL5", "NKG7")
FeaturePlot(mydata, features=set, cols=c("snow", "red"), ncol=4)

genes = c("TNNT2", "CDH5", "DCN", "MYLK", "C1QA", "NRXN1")
VlnPlot(mydata, features=genes, pt.size=0, cols=colors1, ncol=3)+NoLegend()+theme(text=element_text(family="Times"), axis.title.x=element_blank())

#####################################################################  
bar = mydata@meta.data %>% group_by(Sample, cell_type) %>% count()
bar2 = bar %>% group_by(Diagnosis) %>% mutate(Percent=100*n/sum(n))
ggplot(data=bar2, aes(x=cell_type, y=percent, fill=Type))+
geom_bar(stat="identity", position=position_dodge())+
scale_fill_manual(values=c("skyblue", "orange"))+theme_classic()+geom_text(aes(label=percent), vjust=-0.2)+
theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())

#####################################################################  
Type_label = c("Control", "24h Treated", "72h Treated")
bar$Type = factor(bar$Type, levels=Type_label)
bar = bar %>% group_by(Type) %>% mutate(percent=100*n/sum(n))

ggplot(data=bar, aes(x=cell_type, y=percent, fill=Type))+
geom_bar(stat="identity", position=position_dodge())+
scale_fill_manual(values=colors)+theme_classic()+
ggtitle("Percent(%)")+
theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
ggsave("barplot_pair_number.pdf", p, width=6, height=6)


########################## 
library(ggplot2)
library(ggpubr)

df = read.table("./celltype_number_percent.txt", header=T, sep="\t")
ggplot(df, aes(x=reorder(cell_type, -percent, sum), y=percent, fill=Type))+
scale_fill_manual(values=c("orange1", "greenyellow"))+
geom_boxplot(outlier.size=0.1, width=0.3)+
theme_bw()+
stat_compare_means(aes(group=Type), label="p.signif", method="t.test")+
theme(axis.text.x=element_text(angle=15, hjust=1, face="bold", size=10), axis.text.y=element_text(face="bold", size=10), axis.title.x=element_blank())





