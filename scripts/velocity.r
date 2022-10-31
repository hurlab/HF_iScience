library(Seurat)
###load the packages
library(SeuratWrappers)
library(velocyto.R)
Idents(rh.combined)<-"celltype"
rhc <- subset(rh.combined,treat=="Ctrl")
rhh <- subset(rh.combined,treat=="HFD")
#####
#####run velocity 
rhc <- RunVelocity(object = rhc, deltaT = 2, kCells = 25, fit.quantile = 0.02,ncores = 50,spliced.average = 0.2,unspliced.average = 0.2)
ident.colors <- mycol
names(x = ident.colors) <- levels(x = rhc)
cell.colors <- ident.colors[Idents(object = rhc)]
names(x = cell.colors) <- colnames(x = rhc)
pvc<-show.velocity.on.embedding.cor(emb = Embeddings(object = rhc, reduction = "umap"), 
                                     vel = Tool(object = rhc, slot = "RunVelocity"), n = 100, scale = "sqrt",
                                     cell.colors = ac(x = cell.colors, alpha = 0.7), 
                                     cex = 0.4, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 75, arrow.lwd = 0.1, 
                                     do.par = FALSE, cell.border.alpha = 0.01)
## generate figures
show.velocity.on.embedding.cor(emb = Embeddings(object = rhc, reduction = "umap"), 
                               cell.colors = ac(x = cell.colors, alpha = 0.7),
                               vel = Tool(object = rhc, slot = "RunVelocity"), n = 100, scale = "sqrt",  
                               cex = 0.4, arrow.scale = 2, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 75, arrow.lwd = 0.35, 
                               do.par = FALSE, cell.border.alpha = 0.01,cc = pvc$cc,expression.scaling = T)
dev.print(pdf,file="rhc_velocity.pdf")
#####
rhh <- RunVelocity(object = rhh, deltaT = 2, kCells = 25, fit.quantile = 0.02,ncores = 50,spliced.average = 0.2,unspliced.average = 0.2)
ident.colors <- mycol
names(x = ident.colors) <- levels(x = rhh)
cell.colors <- ident.colors[Idents(object = rhh)]
names(x = cell.colors) <- colnames(x = rhh)
pvh<-show.velocity.on.embedding.cor(emb = Embeddings(object = rhh, reduction = "umap"), 
                                    vel = Tool(object = rhh, slot = "RunVelocity"), n = 100, scale = "sqrt",
                                    cell.colors = ac(x = cell.colors, alpha = 0.7), 
                                    cex = 0.4, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 75, arrow.lwd = 0.1, 
                                    do.par = FALSE, cell.border.alpha = 0.01)
##
show.velocity.on.embedding.cor(emb = Embeddings(object = rhh, reduction = "umap"), 
                               cell.colors = ac(x = cell.colors, alpha = 0.7),
                               vel = Tool(object = rhh, slot = "RunVelocity"), n = 100, scale = "sqrt",  
                               cex = 0.4, arrow.scale = 2, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 75, arrow.lwd = 0.35, 
                               do.par = FALSE, cell.border.alpha = 0.01,cc = pvh$cc,expression.scaling = T)

dev.print(pdf,file="rhh_velocity.pdf")