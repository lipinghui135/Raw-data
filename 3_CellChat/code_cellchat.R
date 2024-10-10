library(CellChat)
#################### 
cellchat = createCellChat(object=mydata, group.by="cell_type", meta=mydata@meta.data)
groupSize = as.numeric(table(cellchat@idents))
CellChatDB = CellChatDB.human
#####  
> unique(CellChatDB$interaction$annotation)
[1] "Secreted Signaling"	"ECM-Receptor"	"Cell-Cell Contact"  #
CellChatDB.use <- subsetDB(CellChatDB, search="Cell-Cell Contact")
cellchat@DB <- CellChatDB.use
#####  
cellchat = subsetData(cellchat)
future::plan("multisession", workers=4)
## 
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)  #
## 
cellchat = projectData(cellchat, PPI.human)
## 
## 
#####################################################################################
###
### 
### 
### 
# 
cellchat = computeCommunProb(cellchat, raw.use=FALSE, population.size=TRUE)
# 
cellchat = filterCommunication(cellchat, min.cells=10)
#####  
cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)
#####################################################################################
#######  
## 

netVisual_bubble(cellchat, sources.use=c(2), targets.use=c(1))+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=0, hjust=0.5, face="bold", size=10), axis.text.y=element_text(face="bold", size=10))

# 

##############################################################################
TNBC = subset(scRNA, TNBC_type=="TNBC")
nonTNBC = subset(scRNA, TNBC_type=="non-TNBC")
chat_TNBC = createCellChat(object=TNBC, group.by="cell_type", meta=TNBC@meta.data)
chat_nonTNBC = createCellChat(object=nonTNBC, group.by="cell_type", meta=nonTNBC@meta.data)
CellChatDB.use <- subsetDB(CellChatDB, search="Secreted Signaling")
################################################################################################
chat_TNBC@DB <- CellChatDB.use
chat_TNBC = subsetData(chat_TNBC)
chat_TNBC = identifyOverExpressedGenes(chat_TNBC)
chat_TNBC = identifyOverExpressedInteractions(chat_TNBC)
chat_TNBC = projectData(chat_TNBC, PPI.human)
chat_TNBC = computeCommunProb(chat_TNBC, raw.use=FALSE, population.size=TRUE)
chat_TNBC = filterCommunication(chat_TNBC, min.cells=10)
chat_TNBC = computeCommunProbPathway(chat_TNBC)
chat_TNBC = aggregateNet(chat_TNBC)
chat_TNBC = netAnalysis_computeCentrality(chat_TNBC, slot.name="netP")
################################################################################################
chat_nonTNBC@DB <- CellChatDB.use
chat_nonTNBC = subsetData(chat_nonTNBC)
chat_nonTNBC = identifyOverExpressedGenes(chat_nonTNBC)
chat_nonTNBC = identifyOverExpressedInteractions(chat_nonTNBC)
chat_nonTNBC = projectData(chat_nonTNBC, PPI.human)
chat_nonTNBC = computeCommunProb(chat_nonTNBC, raw.use=FALSE, population.size=TRUE)
chat_nonTNBC = filterCommunication(chat_nonTNBC, min.cells=10)
chat_nonTNBC = computeCommunProbPathway(chat_nonTNBC)
chat_nonTNBC = aggregateNet(chat_nonTNBC)
chat_nonTNBC = netAnalysis_computeCentrality(chat_nonTNBC, slot.name="netP")
################################################################################################
chat_list = list(TNBC=chat_TNBC, `non-TNBC`=chat_nonTNBC)
cellchat = mergeCellChat(chat_list, add.names=names(chat_list), cell.prefix=TRUE)
rankNet(cellchat, mode="comparison", stacked=T, do.stat=TRUE)






