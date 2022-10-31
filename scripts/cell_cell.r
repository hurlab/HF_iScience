###load packages
library(CellChat)
library(Seurat)
Idents(rh.combined)<-"celltype"
HFD1M<-subset(rh.combined,Group=="HFD1M")
Ctrl1M<-subset(rh.combined,Group=="Ctrl1M")
HFD3M<-subset(rh.combined,Group=="HFD3M")
Ctrl3M<-subset(rh.combined,Group=="Ctrl3M")
#####
### Ctrl1M
data.input <- GetAssayData(Ctrl1M, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Ctrl1M)
meta <- data.frame(labels = labels, row.names = names(labels)) 
ct1M <- createCellChat(object = data.input,group.by = "labels",meta = meta)
levels(ct1M@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(ct1M@idents)) # number of cells in each cell group
ct1MDB <- CellChatDB.mouse # use ct1MDB.human if running on human data
ct1MDB.use <- CellChatDB.mouse
ct1M@DB<-ct1MDB.use
ct1M <- subsetData(ct1M) # subset the expression data of signaling genes for saving computation cost
ct1M <- identifyOverExpressedGenes(ct1M)
ct1M <- identifyOverExpressedInteractions(ct1M)
ct1M <- projectData(ct1M, PPI.mouse)
ct1M <- computeCommunProb(ct1M)
ct1M <- computeCommunProbPathway(ct1M)
ct1M <- aggregateNet(ct1M)
ct1M <- netAnalysis_computeCentrality(ct1M, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(ct1M, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
groupSize <- as.numeric(table(ct1M@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(ct1M@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(ct1M@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.print(pdf,file="Number of interactions.pdf")
mat <- ct1M@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.print(pdf,file="signaling_sent_cell_group.pdf")
pathways.show <- ct1M@netP$pathways
##
dir.create("ct1m")
setwd("ct1m")
vertex.receiver = seq(1,6)
for (i in 1:length(pathways.show)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual_aggregate(ct1M, signaling = pathways.show[i], vertex.receiver  = vertex.receiver, color.use = mycol )
  dev.print(pdf,file=paste("ctrl1M_",pathways.show[i],"aggregate.pdf",sep=""))
  par(mfrow=c(1,1))
  netVisual_aggregate(ct1M, signaling = pathways.show[i], layout = "circle",color.use = mycol)
  dev.print(pdf,file=paste("ctrl1M_",pathways.show[i],"circle.pdf",sep=""))
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  netAnalysis_signalingRole_network(ct1M, signaling = pathways.show[i], width = 8, height = 2.5, font.size = 10)
  dev.print(pdf,file=paste0(pathways.show[i], "_signalingRolen_ctrl1M.pdf"))
  gg<-netAnalysis_contribution(ct1M, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_L-R_contribution_ctrl1M.pdf"),plot = gg)
  # Visualize signaling roles of cell groups
  dev.off()
}
#####################
### HF1M
setwd("../hf1m")
data.input <- GetAssayData(HFD1M, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(HFD1M)
meta <- data.frame(labels = labels, row.names = names(labels)) 
HF1M <- createCellChat(object = data.input,group.by = "labels",meta = meta)
levels(HF1M@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(HF1M@idents)) # number of cells in each cell group
HF1MDB <- CellChatDB.mouse # use HF1MDB.human if running on human data
HF1MDB.use <- CellChatDB.mouse
HF1M@DB<-HF1MDB.use
HF1M <- subsetData(HF1M) # subset the expression data of signaling genes for saving computation cost
HF1M <- identifyOverExpressedGenes(HF1M)
HF1M <- identifyOverExpressedInteractions(HF1M)
HF1M <- projectData(HF1M, PPI.mouse)
HF1M <- computeCommunProb(HF1M)
HF1M <- computeCommunProbPathway(HF1M)
HF1M <- aggregateNet(HF1M)
HF1M <- netAnalysis_computeCentrality(HF1M, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(HF1M, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
groupSize <- as.numeric(table(HF1M@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(HF1M@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(HF1M@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.print(pdf,file="Number of interactions.pdf")
mat <- HF1M@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.print(pdf,file="signaling_sent_cell_group.pdf")
pathways.show <- HF1M@netP$pathways
##
vertex.receiver = seq(1,6)
for (i in 1:length(pathways.show)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual_aggregate(HF1M, signaling = pathways.show[i], vertex.receiver  = vertex.receiver, color.use = mycol )
  dev.print(pdf,file=paste("HF1M_",pathways.show[i],"aggregate.pdf",sep=""))
  par(mfrow=c(1,1))
  netVisual_aggregate(HF1M, signaling = pathways.show[i], layout = "circle",color.use = mycol)
  dev.print(pdf,file=paste("HF1M_",pathways.show[i],"circle.pdf",sep=""))
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  netAnalysis_signalingRole_network(HF1M, signaling = pathways.show[i], width = 8, height = 2.5, font.size = 10)
  dev.print(pdf,file=paste0(pathways.show[i], "_signalingRolen_HF1M.pdf"))
  gg<-netAnalysis_contribution(HF1M, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_L-R_contribution_HF1M.pdf"),plot = gg)
  # Visualize signaling roles of cell groups
  dev.off()
}
#####################
### ct3M
dir.create("../Ctrl3M")
setwd("../Ctrl3M/")
data.input <- GetAssayData(Ctrl3M, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Ctrl3M)
meta <- data.frame(labels = labels, row.names = names(labels)) 
ct3M <- createCellChat(object = data.input,group.by = "labels",meta = meta)
levels(ct3M@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(ct3M@idents)) # number of cells in each cell group
ct3MDB <- CellChatDB.mouse # use ct3MDB.human if running on human data
ct3MDB.use <- CellChatDB.mouse
ct3M@DB<-ct3MDB.use
ct3M <- subsetData(ct3M) # subset the expression data of signaling genes for saving computation cost
ct3M <- identifyOverExpressedGenes(ct3M)
ct3M <- identifyOverExpressedInteractions(ct3M)
ct3M <- projectData(ct3M, PPI.mouse)
ct3M <- computeCommunProb(ct3M)
ct3M <- computeCommunProbPathway(ct3M)
ct3M <- aggregateNet(ct3M)
ct3M <- netAnalysis_computeCentrality(ct3M, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(ct3M, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
#
groupSize <- as.numeric(table(ct3M@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(ct3M@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(ct3M@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.print(pdf,file="Number of interactions.pdf")
mat <- ct3M@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.print(pdf,file="signaling_sent_cell_group.pdf")
pathways.show <- ct3M@netP$pathways
##
vertex.receiver = seq(1,6)
for (i in 1:length(pathways.show)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual_aggregate(ct3M, signaling = pathways.show[i], vertex.receiver  = vertex.receiver, color.use = mycol )
  dev.print(pdf,file=paste("ct3M_",pathways.show[i],"aggregate.pdf",sep=""))
  par(mfrow=c(1,1))
  netVisual_aggregate(ct3M, signaling = pathways.show[i], layout = "circle",color.use = mycol)
  dev.print(pdf,file=paste("ct3M_",pathways.show[i],"circle.pdf",sep=""))
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  netAnalysis_signalingRole_network(ct3M, signaling = pathways.show[i], width = 8, height = 2.5, font.size = 10)
  dev.print(pdf,file=paste0(pathways.show[i], "_signalingRolen_ct3M.pdf"))
  gg<-netAnalysis_contribution(ct3M, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_L-R_contribution_ct3M.pdf"),plot = gg)
  # Visualize signaling roles of cell groups
  dev.off()
}
#####################
### HF3M
dir.create("../HF3M")
setwd("../HF3M/")
data.input <- GetAssayData(HFD3M, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(HFD3M)
meta <- data.frame(labels = labels, row.names = names(labels)) 
HF3M <- createCellChat(object = data.input,group.by = "labels",meta = meta)
levels(HF3M@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(HF3M@idents)) # number of cells in each cell group
HF3MDB <- CellChatDB.mouse # use HF3MDB.human if running on human data
HF3MDB.use <- CellChatDB.mouse
HF3M@DB<-HF3MDB.use
HF3M <- subsetData(HF3M) # subset the expression data of signaling genes for saving computation cost
HF3M <- identifyOverExpressedGenes(HF3M)
HF3M <- identifyOverExpressedInteractions(HF3M)
HF3M <- projectData(HF3M, PPI.mouse)
HF3M <- computeCommunProb(HF3M)
HF3M <- computeCommunProbPathway(HF3M)
HF3M <- aggregateNet(HF3M)
HF3M <- netAnalysis_computeCentrality(HF3M, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(HF3M, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
#
groupSize <- as.numeric(table(HF3M@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(HF3M@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(HF3M@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.print(pdf,file="Number of interactions.pdf")
mat <- HF3M@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.print(pdf,file="signaling_sent_cell_group.pdf")
pathways.show <- HF3M@netP$pathways
##
vertex.receiver = seq(1,6)
for (i in 1:length(pathways.show)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual_aggregate(HF3M, signaling = pathways.show[i], vertex.receiver  = vertex.receiver, color.use = mycol )
  dev.print(pdf,file=paste("HF3M_",pathways.show[i],"aggregate.pdf",sep=""))
  par(mfrow=c(1,1))
  netVisual_aggregate(HF3M, signaling = pathways.show[i], layout = "circle",color.use = mycol)
  dev.print(pdf,file=paste("HF3M_",pathways.show[i],"circle.pdf",sep=""))
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  netAnalysis_signalingRole_network(HF3M, signaling = pathways.show[i], width = 8, height = 2.5, font.size = 10)
  dev.print(pdf,file=paste0(pathways.show[i], "_signalingRolen_HF3M.pdf"))
  gg<-netAnalysis_contribution(HF3M, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_L-R_contribution_HF3M.pdf"),plot = gg)
  # Visualize signaling roles of cell groups
  dev.off()
}
#####################
####
cell1M<-mergeCellChat(list(ct1M, HF1M), add.names = c("SD1M","HFD1M"))
object.list1M <- list(SD1M=ct1M, HFD1M=HF1M)
cell3M<- mergeCellChat(list(ct3M,HF3M), add.names = c('SD3M','HFD3M'))
object.list3M <- list(SD3M=ct3M,HFD3M=HF3M)
################
cell1M <- computeNetSimilarityPairwise(cell1M, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cell1M <- netEmbedding(cell1M, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cell1M <- netClustering(cell1M, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cell1M, type = "functional", label.size = 3.5)
######
cell3M <- computeNetSimilarityPairwise(cell3M, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cell3M <- netEmbedding(cell3M, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cell3M <- netClustering(cell3M, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cell3M, type = "functional", label.size = 3.5)
rankSimilarity(cell1M, type = "functional")
dev.print(pdf,file="cell1M_distance.pdf")
rankSimilarity(cell3M, type = "functional")
dev.print(pdf,file="cell3M_distance.pdf")
gg1 <- rankNet(cell1M, mode = "comparison", stacked = T, do.stat = TRUE)
gg2<-rankNet(cell3M, mode = "comparison", stacked = T, do.stat = TRUE)
gg1+gg2
dev.print(pdf,file="information_flow.pdf")
##########