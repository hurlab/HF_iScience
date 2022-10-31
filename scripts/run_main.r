#loading packages
library(Matrix)
library(Seurat)
###read hash tag data
RH1<-readMM("RH1/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
RH1f = read.delim("RH1/outs/filtered_feature_bc_matrix/features.tsv.gz", 
                           header = FALSE,
                           stringsAsFactors = FALSE)
RH1b = read.delim("RH1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(RH1) = RH1b$V1
rownames(RH1) = RH1f$V2
###
RH11<-readMM("RH11/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
RH11f = read.delim("RH11/outs/filtered_feature_bc_matrix/features.tsv.gz", 
                  header = FALSE,
                  stringsAsFactors = FALSE)
RH11b = read.delim("RH11/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                  header = FALSE,
                  stringsAsFactors = FALSE)
colnames(RH11) = RH11b$V1
rownames(RH11) = RH11f$V2
###
RH2<-readMM("RH2/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
RH2f = read.delim("RH2/outs/filtered_feature_bc_matrix/features.tsv.gz", 
                  header = FALSE,
                  stringsAsFactors = FALSE)
RH2b = read.delim("RH2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                  header = FALSE,
                  stringsAsFactors = FALSE)
colnames(RH2) = RH2b$V1
rownames(RH2) = RH2f$V2
####
RH22<-readMM("RH22/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
RH22f = read.delim("RH22/outs/filtered_feature_bc_matrix/features.tsv.gz", 
                  header = FALSE,
                  stringsAsFactors = FALSE)
RH22b = read.delim("RH22/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                  header = FALSE,
                  stringsAsFactors = FALSE)
colnames(RH22) = RH22b$V1
rownames(RH22) = RH22f$V2
##
####
RH3<-readMM("RH31/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
RH3f = read.delim("RH31/outs/filtered_feature_bc_matrix/features.tsv.gz", 
                   header = FALSE,
                   stringsAsFactors = FALSE)
RH3b = read.delim("RH31/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                   header = FALSE,
                   stringsAsFactors = FALSE)
colnames(RH3) = RH3b$V1
rownames(RH3) = RH3f$V2
###
RH33<-readMM("RH32/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
RH33f = read.delim("RH32/outs/filtered_feature_bc_matrix/features.tsv.gz", 
                  header = FALSE,
                  stringsAsFactors = FALSE)
RH33b = read.delim("RH32/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                  header = FALSE,
                  stringsAsFactors = FALSE)
colnames(RH33) = RH33b$V1
rownames(RH33) = RH33f$V2
###
RH4<-readMM("RH41/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
RH4f = read.delim("RH41/outs/filtered_feature_bc_matrix/features.tsv.gz", 
                  header = FALSE,
                  stringsAsFactors = FALSE)
RH4b = read.delim("RH41/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                  header = FALSE,
                  stringsAsFactors = FALSE)
colnames(RH4) = RH4b$V1
rownames(RH4) = RH4f$V2
###
RH44<-readMM("RH42/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
RH44f = read.delim("RH42/outs/filtered_feature_bc_matrix/features.tsv.gz", 
                   header = FALSE,
                   stringsAsFactors = FALSE)
RH44b = read.delim("RH42/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                   header = FALSE,
                   stringsAsFactors = FALSE)
colnames(RH44) = RH44b$V1
rownames(RH44) = RH44f$V2
################################################
###
#tsb<-c("15831","15833","15835")
# 1a 1b 1c
#2a 2b 2c
## read and scale data
rh1 = CreateSeuratObject(RH1[grep('RH1',rownames(RH1),invert = T),])
rh1 = NormalizeData(rh1)
rh1 = FindVariableFeatures(rh1)
rh1 = ScaleData(rh1)
###
rh1[["TSB"]] = CreateAssayObject(counts = RH1[grep('RH1',rownames(RH1)),])
rh1 = NormalizeData(rh1, assay = "TSB", normalization.method = "CLR")
rh1 = ScaleData(rh1, assay = "TSB")
##
rh1 = HTODemux(rh1, assay = "TSB", positive.quantile =  0.98)
HTOHeatmap(rh1, assay = "TSB", ncells = 5000)
### check the hashtag status
table(rh1$TSB_classification.global)
####
rh11 = CreateSeuratObject(RH11[grep('RH1',rownames(RH11),invert = T),])
rh11 = NormalizeData(rh11)
rh11 = FindVariableFeatures(rh11)
rh11 = ScaleData(rh11)
rh11[["TSB"]] = CreateAssayObject(counts = RH11[grep('RH1',rownames(RH11)),])
rh11 = NormalizeData(rh11, assay = "TSB", normalization.method = "CLR")
rh11 = ScaleData(rh11, assay = "TSB")
##
rh11 = HTODemux(rh11, assay = "TSB", positive.quantile =  0.98)
HTOHeatmap(rh11, assay = "TSB", ncells = 5000,raster = F)
###
table(rh11$TSB_classification.global)
####
rh2 = CreateSeuratObject(RH2[grep('RH2',rownames(RH2),invert = T),])
rh2 = NormalizeData(rh2)
rh2 = FindVariableFeatures(rh2)
rh2 = ScaleData(rh2)
###
rh2[["TSB"]] = CreateAssayObject(counts = RH2[grep('RH2',rownames(RH2)),])
rh2 = NormalizeData(rh2, assay = "TSB", normalization.method = "CLR")
rh2 = ScaleData(rh2, assay = "TSB")
##
rh2 = HTODemux(rh2, assay = "TSB", positive.quantile =  0.98)
HTOHeatmap(rh2, assay = "TSB", ncells = 5000)
###
table(rh2$TSB_classification.global)
##
rh22 = CreateSeuratObject(RH22[grep('RH2',rownames(RH22),invert = T),])
rh22 = NormalizeData(rh22)
rh22 = FindVariableFeatures(rh22)
rh22 = ScaleData(rh22)
###
rh22[["TSB"]] = CreateAssayObject(counts = RH22[grep('RH2',rownames(RH22)),])
rh22 = NormalizeData(rh22, assay = "TSB", normalization.method = "CLR")
rh22 = ScaleData(rh22, assay = "TSB")
##
rh22 = HTODemux(rh22, assay = "TSB", positive.quantile =  0.98)
HTOHeatmap(rh22, assay = "TSB", ncells = 5000,raster =F )
###
table(rh22$TSB_classification.global)
########
###
rh3 = CreateSeuratObject(RH3[grep('RH1',rownames(RH3),invert = T),])
rh3 = NormalizeData(rh3)
rh3 = FindVariableFeatures(rh3)
rh3 = ScaleData(rh3)
###
rh3[["TSB"]] = CreateAssayObject(counts = RH3[grep('RH1',rownames(RH3)),])
rh3 = NormalizeData(rh3, assay = "TSB", normalization.method = "CLR")
rh3 = ScaleData(rh3, assay = "TSB")
##
rh3 = HTODemux(rh3, assay = "TSB", positive.quantile =  0.98)
HTOHeatmap(rh3, assay = "TSB", ncells = 5000)
###
table(rh3$TSB_classification.global)
########
rh33 = CreateSeuratObject(RH33[grep('RH2',rownames(RH33),invert = T),])
rh33 = NormalizeData(rh33)
rh33 = FindVariableFeatures(rh33)
rh33 = ScaleData(rh33)
###
rh33[["TSB"]] = CreateAssayObject(counts = RH33[grep('RH2',rownames(RH33)),])
rh33 = NormalizeData(rh33, assay = "TSB", normalization.method = "CLR")
rh33 = ScaleData(rh33, assay = "TSB")
##
rh33 = HTODemux(rh33, assay = "TSB", positive.quantile =  0.98)
HTOHeatmap(rh33, assay = "TSB", ncells = 5000)
###
table(rh33$TSB_classification.global)
################
rh4 = CreateSeuratObject(RH4[grep('RH1',rownames(RH4),invert = T),])
rh4 = NormalizeData(rh4)
rh4 = FindVariableFeatures(rh4)
rh4 = ScaleData(rh4)
###
rh4[["TSB"]] = CreateAssayObject(counts = RH4[grep('RH1',rownames(RH4)),])
rh4 = NormalizeData(rh4, assay = "TSB", normalization.method = "CLR")
rh4 = ScaleData(rh4, assay = "TSB")
##
rh4 = HTODemux(rh4, assay = "TSB", positive.quantile =  0.98)
HTOHeatmap(rh4, assay = "TSB", ncells = 5000)
###
table(rh4$TSB_classification.global)
########
rh44 = CreateSeuratObject(RH44[grep('RH2',rownames(RH44),invert = T),])
rh44 = NormalizeData(rh44)
rh44 = FindVariableFeatures(rh44)
rh44 = ScaleData(rh44)
###
rh44[["TSB"]] = CreateAssayObject(counts = RH44[grep('RH2',rownames(RH44)),])
rh44 = NormalizeData(rh44, assay = "TSB", normalization.method = "CLR")
rh44 = ScaleData(rh44, assay = "TSB")
##
rh44 = HTODemux(rh44, assay = "TSB", positive.quantile =  0.98)
HTOHeatmap(rh44, assay = "TSB", ncells = 5000)
###
table(rh44$TSB_classification.global)
###############################################################################
#########################################################
##generate figures
#############
Idents(rh1) <- "TSB_maxID"
RidgePlot(rh1, assay = "TSB", features = rownames(rh1[["TSB"]])[1:2], ncol = 2)
FeatureScatter(rh1, feature1 = "tsb_RH1a", feature2 = "tsb_RH1b")
dev.print(pdf,file="RH1a_b_scatter.pdf")
FeatureScatter(rh1, feature1 = "tsb_RH1b", feature2 = "tsb_RH1c")
dev.print(pdf,file="RH1b_c_scatter.pdf")
FeatureScatter(rh1, feature1 = "tsb_RH1a", feature2 = "tsb_RH1c")
dev.print(pdf,file="RH1a_c_scatter.pdf")
####
Idents(rh1) <- "TSB_classification.global"
VlnPlot(rh1, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.print(pdf,file="RH1a_nCount.pdf")
Idents(rh11) <- "TSB_classification.global"
####
Idents(rh2) <- "TSB_maxID"
RidgePlot(rh2, assay = "TSB", features = rownames(rh2[["TSB"]])[1:2], ncol = 2)
FeatureScatter(rh2, feature1 = "tsb_RH2a", feature2 = "tsb_RH2b")
dev.print(pdf,file="RH2a_b_scatter.pdf")
FeatureScatter(rh2, feature1 = "tsb_RH2b", feature2 = "tsb_RH2c")
dev.print(pdf,file="RH2b_c_scatter.pdf")
FeatureScatter(rh2, feature1 = "tsb_RH2a", feature2 = "tsb_RH2c")
dev.print(pdf,file="RH2a_c_scatter.pdf")
#### set idents
Idents(rh1) <- "TSB_classification.global"
Idents(rh11) <- "TSB_classification.global"
Idents(rh2) <- "TSB_classification.global"
Idents(rh22) <- "TSB_classification.global"
Idents(rh3) <- "TSB_classification.global"
Idents(rh33) <- "TSB_classification.global"
Idents(rh4) <- "TSB_classification.global"
Idents(rh44) <- "TSB_classification.global"
###############
##only select the singlet data and check the quality
rh1.singlet <- subset(rh1, idents = "Singlet")
rh1.singlet[["percent.mt"]]<-PercentageFeatureSet(rh1.singlet,pattern = "^mt")
VlnPlot(rh1.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.print(pdf,file="RH1_single_count_mt.pdf")
###
rh11.singlet <- subset(rh11, idents = "Singlet")
rh11.singlet[["percent.mt"]]<-PercentageFeatureSet(rh11.singlet,pattern = "^mt")
VlnPlot(rh11.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
####
rh2.singlet <- subset(rh2, idents = "Singlet")
rh2.singlet[["percent.mt"]]<-PercentageFeatureSet(rh2.singlet,pattern = "^mt")
VlnPlot(rh2.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.print(pdf,file="RH2_single_count_mt.pdf")
###
rh22.singlet <- subset(rh22, idents = "Singlet")
rh22.singlet[["percent.mt"]]<-PercentageFeatureSet(rh22.singlet,pattern = "^mt")
VlnPlot(rh22.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
###
rh3.singlet <- subset(rh3, idents = "Singlet")
rh3.singlet[["percent.mt"]]<-PercentageFeatureSet(rh3.singlet,pattern = "^mt")
VlnPlot(rh3.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.print(pdf,file="RH3_single_count_mt.pdf")
###
rh33.singlet <- subset(rh33, idents = "Singlet")
rh33.singlet[["percent.mt"]]<-PercentageFeatureSet(rh33.singlet,pattern = "^mt")
VlnPlot(rh33.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
###
rh4.singlet <- subset(rh4, idents = "Singlet")
rh4.singlet[["percent.mt"]]<-PercentageFeatureSet(rh4.singlet,pattern = "^mt")
VlnPlot(rh4.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.print(pdf,file="RH4_single_count_mt.pdf")
###
rh44.singlet <- subset(rh44, idents = "Singlet")
rh44.singlet[["percent.mt"]]<-PercentageFeatureSet(rh44.singlet,pattern = "^mt")
VlnPlot(rh44.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.print(pdf,file="RH44_single_count_mt.pdf")
#################################################
## filter data with nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt<25
rh1s <- subset(rh1.singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt<25)
rh11s <- subset(rh11.singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt<25 )
rh2s <- subset(rh2.singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt<25 )
rh22s <- subset(rh22.singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt<25 )
rh3s <- subset(rh3.singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt<25)
rh33s <- subset(rh33.singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt<25 )
rh4s <- subset(rh4.singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt<25 )
rh44s <- subset(rh44.singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt<25 )
##########################
rh1s <- FindVariableFeatures(rh1s, selection.method = "vst")
rh11s <- FindVariableFeatures(rh11s, selection.method = "vst")
rh2s<-FindVariableFeatures(rh2s,selection.method = "vst")
rh22s <- FindVariableFeatures(rh22s, selection.method = "vst")
rh3s<-FindVariableFeatures(rh3s,selection.method = "vst")
rh33s <- FindVariableFeatures(rh33s, selection.method = "vst")
rh4s<-FindVariableFeatures(rh4s,selection.method = "vst")
rh44s <- FindVariableFeatures(rh44s, selection.method = "vst")
####################
#### add group factor
rh1s$group<-paste0(rh1s$TSB_maxID,"1")
rh11s$group<-paste0(rh11s$TSB_maxID,"2")
rh2s$group<-paste0(rh2s$TSB_maxID,"1")
rh22s$group<-paste0(rh22s$TSB_maxID,"2")
rh3s$TSB_maxID<-sub('1','3',rh3s$TSB_maxID)
rh33s$TSB_maxID<-sub('2','3',rh33s$TSB_maxID)
rh4s$TSB_maxID<-sub('1','4',rh4s$TSB_maxID)
rh44s$TSB_maxID<-sub('2','4',rh44s$TSB_maxID)
rh3s$group<-paste0(rh3s$TSB_maxID,"1")
rh33s$group<-paste0(rh33s$TSB_maxID,"2")
rh4s$group<-paste0(rh4s$TSB_maxID,"1")
rh44s$group<-paste0(rh44s$TSB_maxID,"2")
#######
Idents(rh1s)<-"group"
Idents(rh11s)<-"group"
Idents(rh2s)<-"group"
Idents(rh22s)<-"group"
Idents(rh3s)<-"group"
Idents(rh33s)<-"group"
Idents(rh4s)<-"group"
Idents(rh44s)<-"group"
### integrate data
rh.anchors <- FindIntegrationAnchors(object.list = list(rh1s,rh11s,rh2s,rh22s,rh3s,rh33s,rh4s,rh44s), dims = 1:50)
rh.combined <- IntegrateData(anchorset = rh.anchors, dims = 1:50)
#################
DefaultAssay(rh.combined) <- "integrated"
all.genes <- rownames(rh.combined)
rh.combined <- ScaleData(rh.combined, verbose = TRUE,features = all.genes)
rh.combined <- RunPCA(rh.combined, npcs = 50, features = VariableFeatures(object = rh.combined))
#### find the best number of PCs
pct <- rh.combined[["pca"]]@stdev / sum(rh.combined[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
##################
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
######
pcs <- min(co1, co2)
pcs
######################
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") + xlab("Cumulative percent ")+ylab("Percent of variation")+
  theme_bw(base_size=15)+ theme(legend.position = "none")
dev.print(pdf,file="PCA_com_selection.pdf")
### 
###############################################
## UMAP
rh.combined <- RunUMAP(rh.combined, reduction = "pca", dims = 1:20)
rh.combined <- FindNeighbors(rh.combined, reduction = "pca", dims = 1:20)
rh.combined <- FindClusters(rh.combined, resolution = 0.6)
### generate figures
DimPlot(rh.combined,label=T,split.by = "group")
dev.print(pdf,file="rh_umap_all.pdf")
DimPlot(rh.combined,label=T)
###find markers
marker<-FindAllMarkers(rh.combined)
rownames(rhgroup)<-rhgroup$sequenceID
marker %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
##### assign cell type
clusters<-c('HMG1','HMG2','HMG3','Mac','UMG','IFLAMMG','PLFMG','Mac','Neuro','Mac','IFNMG','Neu','Mono','Mono')
names(clusters)<-levels(Idents(rh.combined))
rh.combined<-RenameIdents(rh.combined,clusters)
###
rh.combined$celltype<-Idents(rh.combined)
mycol<-c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C", "#7570B3","#E7298A","cyan3","violet","springgreen2" )
names(mycol)<-c("HMG1","HMG2","HMG3","Mac","UMG","IFLAMMG","PLFMG","Neuro","IFNMG","Neu","Mono")
DimPlot(rh.combined,label=T,cols=mycol)
dev.print(pdf,file="umap_new.pdf")
### add group information
### read the sample information
rhgroup<-read.csv("RH_group.csv")
group<-c("SD","SD","HFD","SD","HFD","HFD","SD","SD","HFD","SD","HFD","HFD")
names(group)<-c("RH1a1","RH1b1","RH1c1","RH2a1","RH2b1","RH2c1","RH1a2","RH1b2","RH1c2","RH2a2","RH2b2","RH2c2")
rhgroup$group<-sub('month','M',paste(rhgroup$Group,rhgroup$Time,sep=""))
rh.combined$treat<-rhgroup[rh.combined$group,"Group"]
rh.combined$Group<-rhgroup[rh.combined$group,"group"]
##################
### Find DEGs
DefaultAssay(rh.combined)<-"RNA"
rh.combined$test<-paste(rh.combined$Group,rh.combined$celltype,sep="_")
Idents(rh.combined)<-"test"
for(i in unique(rh.combined$celltype)){
  write.csv(FindMarkers(rh.combined,ident.1 = paste0("HFD1M_",i),ident.2 = paste0("Ctrl1M_",i),test.use = "DESeq2"),file=paste0("HFDvsSD1M_",i,"DEG.csv"))
  write.csv(FindMarkers(rh.combined,ident.1 = paste0("HFD3M_",i),ident.2 = paste0("Ctrl3M_",i),test.use = "DESeq2"),file=paste0("HFDvsSD3M_",i,"DEG.csv"))
}
filename<-list.files(pattern = "*DEG.csv")
deg<-lapply(filename, function(x)read.csv(x,row.names = 1))
names(deg)<-sub('DEG.csv','',filename)
###### functional enrichment analysis
sig<-lapply(deg, function(x)rownames(subset(x,p_val<0.01)))
### load the richR package 
## devtools::install_github('guokai8/richR')
library(richR)
mmko<-buildAnnot(species="mouse",anntype = "KEGG",builtin = F)
enr<-lapply(sig, function(x)richKEGG(x,mmko,builtin = F))
###1M figures
ggdot(enr$HFDvsSD1M_HMG1,usePadj = F)
dev.print(pdf,file="HFDvsSD1M_HMG1_KEGG.pdf")
ggdot(enr$HFDvsSD1M_HMG2,usePadj = F)
dev.print(pdf,file="HFDvsSD1M_HMG2_KEGG.pdf")
ggdot(enr$HFDvsSD1M_HMG3,usePadj = F)
dev.print(pdf,file="HFDvsSD1M_HMG3_KEGG.pdf")
ggdot(enr$HFDvsSD1M_IFLAMMG,usePadj = F)
dev.print(pdf,file="HFDvsSD1M_IFLAMMG_KEGG.pdf")
ggdot(enr$HFDvsSD1M_IFNMG,usePadj = F)
dev.print(pdf,file="HFDvsSD1M_IFN_KEGG.pdf")
ggdot(enr$HFDvsSD1M_Mac,usePadj = F)
dev.print(pdf,file="HFDvsSD1M_Mac_KEGG.pdf")
ggdot(enr$HFDvsSD1M_UMG,usePadj = F)
dev.print(pdf,file="HFDvsSD1M_UMG_KEGG.pdf")
ggdot(enr$HFDvsSD1M_Mono,usePadj = F)
dev.print(pdf,file="HFDvsSD1M_Mono_KEGG.pdf")
ggdot(enr$HFDvsSD1M_PLFMG,usePadj = F)
dev.print(pdf,file="HFDvsSD1M_PLFMG_KEGG.pdf")
ggdot(enr$HFDvsSD1M_Neu,usePadj = F)
ggdot(enr$HFDvsSD1M_Neuro,usePadj = F)
dev.print(pdf,file="HFDvsSD1M_Neuro_KEGG.pdf")
######### 3M figures
ggdot(enr$HFDvsSD3M_HMG1,usePadj = F)
dev.print(pdf,file="HFDvsSD3M_HMG1_KEGG.pdf")
ggdot(enr$HFDvsSD3M_HMG2,usePadj = F)
dev.print(pdf,file="HFDvsSD3M_HMG2_KEGG.pdf")
ggdot(enr$HFDvsSD3M_HMG3,usePadj = F)
dev.print(pdf,file="HFDvsSD3M_HMG3_KEGG.pdf")
ggdot(enr$HFDvsSD3M_IFLAMMG,usePadj = F)
dev.print(pdf,file="HFDvsSD3M_IFLAMMG_KEGG.pdf")
ggdot(enr$HFDvsSD3M_IFNMG,usePadj = F)
dev.print(pdf,file="HFDvsSD3M_IFN_KEGG.pdf")
ggdot(enr$HFDvsSD3M_Mac,usePadj = F)
dev.print(pdf,file="HFDvsSD3M_Mac_KEGG.pdf")
ggdot(enr$HFDvsSD3M_UMG,usePadj = F)
dev.print(pdf,file="HFDvsSD3M_UMG_KEGG.pdf")
ggdot(enr$HFDvsSD3M_Mono,usePadj = F)
dev.print(pdf,file="HFDvsSD3M_Mono_KEGG.pdf")
ggdot(enr$HFDvsSD3M_PLFMG,usePadj = F)
dev.print(pdf,file="HFDvsSD3M_PLFMG_KEGG.pdf")
ggdot(enr$HFDvsSD3M_Neu,usePadj = F)
dev.print(pdf,file="HFDvsSD3M_Neu_KEGG.pdf")
ggdot(enr$HFDvsSD3M_Neuro,usePadj = F)
### generate proportion figures
metadata<-rh.combined@meta.data
metadata%>%filter(treat=="Ctrl")%>%dplyr::select(treat,celltype)%>%group_by(treat,celltype)%>%
summarise(count=n())%>%mutate(cell=count/sum(count))%>%ggplot(aes(x=2,cell,fill=celltype))+geom_bar(width = 1, size = 0.1, color = "white", stat = "identity")+scale_fill_manual(values =mycol)+
  geom_text_repel(aes(label = paste0(round(100*cell,1), "%")),
                  position = position_stack(vjust = 0.5),max.overlaps = 100) + xlim(0.5,2.5)+
  coord_polar("y")+ guides(fill=guide_legend(ncol=2)) +  theme_classic() +labs(x = NULL, y = NULL, fill = NULL,
                                                                               title = "SD") +theme(axis.line = element_blank(),
                                                                                                    axis.text = element_blank(),
                                                                                                    axis.ticks = element_blank(),
                                                                                                    plot.title = element_text(hjust = 0.5, color = "#666666"))

metadata%>%filter(treat=="HFD")%>%dplyr::select(treat,celltype)%>%group_by(treat,celltype)%>%
summarise(count=n())%>%mutate(cell=count/sum(count))%>%ggplot(aes(x=2,cell,fill=celltype))+geom_bar(width = 1, size = 0.1, color = "white", stat = "identity")+scale_fill_manual(values =mycol)+
geom_text_repel(aes(label = paste0(round(100*cell,1), "%")),max.overlaps = 100,
position = position_stack(vjust = 0.5)) + xlim(0.5,2.5)+
coord_polar("y")+ guides(fill=guide_legend(ncol=2)) +  theme_classic() +labs(x = NULL, y = NULL, fill = NULL,
title = "HFD") +theme(axis.line = element_blank(),
axis.text = element_blank(),
axis.ticks = element_blank(),
plot.title = element_text(hjust = 0.5, color = "#666666"))
dev.print(pdf,file="HFD_pro.pdf")
##############################
ggdot(enr$GE,usePadj = F,top = Inf)
dev.print(pdf,file="GE_KEGG.pdf")
richR::ggdot(enr$HMG1,usePadj = F)
dev.print(pdf,file="HMG1_KEGG.pdf")
richR::ggdot(enr$HMG2,usePadj = F)
dev.print(pdf,file="HMG2_KEGG.pdf")
richR::ggdot(enr$HMG3,usePadj = F)
dev.print(pdf,file="HMG3_KEGG.pdf")
richR::ggdot(enr$HMG4,usePadj = F)
dev.print(pdf,file="HMG4_KEGG.pdf")
richR::ggdot(enr$HMG5,usePadj = F)
dev.print(pdf,file="HMG5_KEGG.pdf")
richR::ggdot(enr$IFLAMMG,usePadj = F)
dev.print(pdf,file="IFLAMMG_KEGG.pdf")
richR::ggdot(enr$IFNMG,usePadj = F)
dev.print(pdf,file="IFNMG_KEGG.pdf")
richR::ggdot(enr$Mac,usePadj = F)
dev.print(pdf,file="Mac_KEGG.pdf")
richR::ggdot(enr$Mono,usePadj = F)
dev.print(pdf,file="Mono_KEGG.pdf")
richR::ggdot(enr$Neu,usePadj = F)
dev.print(pdf,file="Neu_KEGG.pdf")
richR::ggdot(enr$Neuro,usePadj = F)
dev.print(pdf,file="Neuro_KEGG.pdf")
#########
### generate marker dot plot
genes<-c("Cx3cr1","Slc2a5","Rps11","Rps21","Ccl3","Ccl4","Spp1",
         "Ifit2","Ifit3","Mki67","Top2a","Ube2c","S100a8","S100a9","P2ry12", "Olfml3", "Tmem119", "P2ry13",
         "Nav2", "Nav3", "Macf1", "Tanc2",
         "Apoe","Mrc1","Ly6c2","Ccr2","Map2","Camk2a")

DotPlot(rh.combined,features = genes,cluster.idents = T)+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5))
################

DotPlot(rh.combined,features = genes,cluster.idents = T,cols = c("lightgrey","red"))+
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5))+coord_flip()
############


##################
### add loom information
### loom files were generated with velocyto
rh1loom<-ReadVelocity("../RH1/velocyto/RH1.loom")
rh11loom<-ReadVelocity("../RH11//velocyto/RH11.loom")
rh2loom<-ReadVelocity("../RH2/velocyto/RH2.loom")
rh22loom<-ReadVelocity("../RH22/velocyto/RH22.loom")
rh31loom<-ReadVelocity("../RH31/velocyto/RH31.loom")
rh32loom<-ReadVelocity("../RH32/velocyto/RH32.loom")
rh41loom<-ReadVelocity("../RH41/velocyto/RH41.loom")
rh42loom<-ReadVelocity("../RH42/velocyto/RH42.loom")

#####
saml<-list()
for(i in names(rh1loom)){
  colnames(rh1loom[[i]])<-sub('x','-1_1',sub('.*:','',colnames(rh1loom[[i]])))
  colnames(rh11loom[[i]])<-sub('x','-1_2',sub('.*:','',colnames(rh11loom[[i]])))
  colnames(rh2loom[[i]])<-sub('x','-1_3',sub('.*:','',colnames(rh2loom[[i]])))
  colnames(rh22loom[[i]])<-sub('x','-1_4',sub('.*:','',colnames(rh22loom[[i]])))
  colnames(rh31loom[[i]])<-sub('x','-1_5',sub('.*:','',colnames(rh31loom[[i]])))
  colnames(rh32loom[[i]])<-sub('x','-1_6',sub('.*:','',colnames(rh32loom[[i]])))
  colnames(rh41loom[[i]])<-sub('x','-1_7',sub('.*:','',colnames(rh41loom[[i]])))
  colnames(rh42loom[[i]])<-sub('x','-1_8',sub('.*:','',colnames(rh42loom[[i]])))
  saml[[i]]<-cbind(rh1loom[[i]],rh11loom[[i]],rh2loom[[i]],rh22loom[[i]],rh31loom[[i]],rh32loom[[i]],rh41loom[[i]],rh42loom[[i]])
}

##############################
for (i in names(x = saml)) {
  ### Store assay in a new variable
  assay <- saml[[i]]
  ### Subset to filtered cells in Seurat object
  assay <- assay[,colnames(rh.combined)]
  ### Add assay to Seurat object
  rh.combined[[i]] <- CreateAssayObject(counts = assay)
}
#############################################
### generate proportion figures with SD and HFD at 1M and 3M
metadata<-rh.combined@meta.data
metadata%>%filter(treat=="Ctrl",Time=="1month")%>%dplyr::select(treat,celltype)%>%group_by(treat,celltype)%>%
  summarise(count=n())%>%mutate(cell=count/sum(count))%>%ggplot(aes(x=2,cell,fill=celltype))+geom_bar(width = 1, size = 0.1, color = "white", stat = "identity")+scale_fill_manual(values =mycol)+
  geom_text_repel(aes(label = paste0(round(100*cell,1), "%")),
                  position = position_stack(vjust = 0.5),max.overlaps = 100) + xlim(0.5,2.5)+
  coord_polar("y")+ guides(fill=guide_legend(ncol=2)) +  theme_classic() +labs(x = NULL, y = NULL, fill = NULL,
                                                                               title = "SD1M") +theme(axis.line = element_blank(),
                                                                                                    axis.text = element_blank(),
                                                                                                    axis.ticks = element_blank(),
                                                                                                    plot.title = element_text(hjust = 0.5, color = "#666666"))


dev.print(pdf,file="SD1M_pro.pdf")
metadata%>%filter(treat=="Ctrl",Time=="3month")%>%dplyr::select(treat,celltype)%>%group_by(treat,celltype)%>%
  summarise(count=n())%>%mutate(cell=count/sum(count))%>%ggplot(aes(x=2,cell,fill=celltype))+geom_bar(width = 1, size = 0.1, color = "white", stat = "identity")+scale_fill_manual(values =mycol)+
  geom_text_repel(aes(label = paste0(round(100*cell,1), "%")),
                  position = position_stack(vjust = 0.5),max.overlaps = 100) + xlim(0.5,2.5)+
  coord_polar("y")+ guides(fill=guide_legend(ncol=2)) +  theme_classic() +labs(x = NULL, y = NULL, fill = NULL,
                                                                               title = "SD3M") +theme(axis.line = element_blank(),
                                                                                                    axis.text = element_blank(),
                                                                                                    axis.ticks = element_blank(),
                                                                                                    plot.title = element_text(hjust = 0.5, color = "#666666"))
dev.print(pdf,file="SD3M_pro.pdf")
#
metadata%>%filter(treat=="HFD",Time=="1month")%>%dplyr::select(treat,celltype)%>%group_by(treat,celltype)%>%
  summarise(count=n())%>%mutate(cell=count/sum(count))%>%ggplot(aes(x=2,cell,fill=celltype))+geom_bar(width = 1, size = 0.1, color = "white", stat = "identity")+scale_fill_manual(values =mycol)+
  geom_text_repel(aes(label = paste0(round(100*cell,1), "%")),max.overlaps = 100,
                  position = position_stack(vjust = 0.5)) + xlim(0.5,2.5)+
  coord_polar("y")+ guides(fill=guide_legend(ncol=2)) +  theme_classic() +labs(x = NULL, y = NULL, fill = NULL,
                                                                               title = "HFD1M") +theme(axis.line = element_blank(),
                                                                                                     axis.text = element_blank(),
                                                                                                     axis.ticks = element_blank(),
                                                                                                     plot.title = element_text(hjust = 0.5, color = "#666666"))

dev.print(pdf,file="HFD1M_pro.pdf")
#
metadata%>%filter(treat=="HFD",Time=="3month")%>%dplyr::select(treat,celltype)%>%group_by(treat,celltype)%>%
  summarise(count=n())%>%mutate(cell=count/sum(count))%>%ggplot(aes(x=2,cell,fill=celltype))+geom_bar(width = 1, size = 0.1, color = "white", stat = "identity")+scale_fill_manual(values =mycol)+
  geom_text_repel(aes(label = paste0(round(100*cell,1), "%")),max.overlaps = 100,
                  position = position_stack(vjust = 0.5)) + xlim(0.5,2.5)+
  coord_polar("y")+ guides(fill=guide_legend(ncol=2)) +  theme_classic() +labs(x = NULL, y = NULL, fill = NULL,
                                                                               title = "HFD3M") +theme(axis.line = element_blank(),
                                                                                                     axis.text = element_blank(),
                                                                                                     axis.ticks = element_blank(),
                                                                                                     plot.title = element_text(hjust = 0.5, color = "#666666"))
dev.print(pdf,file="HFD3M_pro.pdf")

