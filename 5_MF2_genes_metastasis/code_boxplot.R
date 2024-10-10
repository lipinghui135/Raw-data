library(tidyverse)
library(ggplot2)
library(gghalves)
library(cols4all)

df = read.table("./MF2_metastasis_gene_boxplot.txt", header=T, row.names=1, sep="\t")
mycol <- c4a('10', 6)     # 
df = subset(df, Gene=="CAMK2N1")
p = ggplot(data=df, aes(x=Type, y=Exp))+
geom_jitter(shape=21, fill='grey', color='grey', size=2, alpha=1, width=0.2)+
geom_violin(aes(fill=Type), color='grey', alpha=0.8, scale='width', linewidth=0.6, trim=TRUE)+
geom_boxplot(color='white',
             outlier.color='black',
             width=0.4, # 
             size=0.8, # 
             fill=NA) +
scale_fill_manual(values=mycol)+
scale_color_manual(values=mycol)+
theme_bw()+theme(axis.text=element_text(size=13, face="bold"), axis.title=element_blank())
ggsave("CAMK2N1_boxplot.pdf", p, width=6, height=6)













