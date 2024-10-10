library(ggplot2)

df = read.table("./pseudotime_hallmark_scatter.txt", header=T, row.names=1, sep="\t")

ggplot(df, aes(x=Pseudotime, y=HALLMARK_TGF_BETA_SIGNALING, color=Type))+
scale_color_manual(values=c("#87CEEB", "#FFA500"))+
geom_point(alpha=1, size=2)+
geom_smooth(method="lm", color="black")+
theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(face="bold", size=13), axis.text.y=element_text(face="bold", size=13), axis.title=element_text(face="bold", size=14), legend.position="none")







