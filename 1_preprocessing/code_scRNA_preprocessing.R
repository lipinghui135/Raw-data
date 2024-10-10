library(Seurat)
library(tidyverse)
library(sctransform)
library(harmony)
###############################################################################
assays <- dir("./data/")
dir <- paste0("./data/", assays)
# 
samples_name = assays

# 
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir=dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])   
  if(T){    
    scRNAlist[[i]][["percent.mito"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-") 
  }
}
### 
names(scRNAlist) <- samples_name
# 
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<5000&percent.mito<25) # 
p = VlnPlot(scRNA, features=c("nFeature_RNA"), pt.size=0, cols=colors)
ggsave("nFeature_RNA_violin.pdf", p, width=6, height=6)
colors = c4a("carto.pastel", 8)
#########################################################################################################
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
scRNA <- CellCycleScoring(scRNA, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
scRNA = SCTransform(scRNA, method="glmGamPoi", vars.to.regress=c("nCount_RNA", "percent.mito", "S.Score", "G2M.Score"), verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
scRNA = RunHarmony(scRNA, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims=1:20, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:20, reduction="harmony")
p = DimPlot(scRNA, reduction="umap", group.by="Sample", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("UMAP_Sample.pdf", p, height=6, width=6)
######################################################### resolution=0.1
mydata <- FindClusters(scRNA, resolution=0.1)
UMAPPlot(mydata, pt.size=1, label=T, cols=colors, label.size=5)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "./All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
############## 
VlnPlot(mydata, features=c("S100A9"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())
# 0: Hepatocytes: HP, FGB, TTR, APOA2
# 1: Macrophages: C1QC, C1QA, MS4A7, CD163
# 2: Endothelial cells: FCN3, CRHBP
# 3: Proliferative NK/T cells: CD3D, NKG7, GNLY, TIGIT
# 4: Plasma B cells: CD79A, MZB1, JCHAIN
# 5: Myofibroblasts: ACTA2, TAGLN, MYLK, COL3A1, BGN
# 6: Dendritic cells: WDFY4, CLEC9A, XCR1
# 7: Proliferative hepatocytes: TOP2A, MKI67, CDK1, HP, FGB

# 
cell_label = c(
"Hepatocytes", "Macrophages", "Endothelial cells", "Proliferative NK/T cells", "Plasma B cells",
"Myofibroblasts", "Dendritic cells", "Proliferative hepatocytes"
)
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
p = UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)
ggsave("UMAP_subcluster.pdf", p, width=6, height=6)
genes = c("HP", "FGB", "TTR", "APOA2", "C1QA", "AIF1", "CD163", "FCN3", "CRHBP", "CD3D", "NKG7", "MZB1", "JCHAIN", "ACTA2", "TAGLN", "COL3A1", "BGN", "WDFY4", "CLEC9A", "TOP2A", "MKI67")
p = DotPlot(mydata, features=genes, cols=c("snow", "chartreuse4"))+coord_flip()+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())

# 
set = c("CD24", "CD19", "CCR7", "CD3D", "HBB", "CST3", "CCL5", "NKG7")
FeaturePlot(mydata, features=set, cols=c("snow", "red"), ncol=4)

genes = c("AKR1C3", "WNT7A", "FAM72B", "RERG", "IDO1", "HEY1")
VlnPlot(mydata, features=genes, pt.size=0, cols=colors, ncol=2)+NoLegend()+theme(axis.title.x=element_blank())


#####################################################################  
Type_label = c("Control", "24h Treated", "72h Treated")
bar$Type = factor(bar$Type, levels=Type_label)
bar = bar %>% group_by(Type) %>% mutate(percent=100*n/sum(n))

ggplot(data=bar2, aes(x=cell_type, y=percent, fill=Type))+
geom_bar(stat="identity", position=position_dodge())+
scale_fill_manual(values=c("skyblue", "orange"))+theme_classic()+geom_text(aes(label=percent), vjust=-0.2)+
theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())


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


	



