library(Seurat)
library(cowplot)
KO.data <- Read10X(data.dir = "D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Treg_KO/outs/filtered_feature_bc_matrix")
WT.data <- Read10X(data.dir = "D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/WT/outs/filtered_feature_bc_matrix")

KO.data <- as.data.frame(KO.data)
WT.data <- as.data.frame(WT.data)

for (i in 1:4105) {
  colnames(KO.data)[i] <- paste(colnames(KO.data)[i],"KO",i,sep = "-")  
}

for (i in 1:3905) {
  colnames(WT.data)[i] <- paste(colnames(WT.data)[i],"WT",i,sep = "-")  
}

KO.metadata<-data.frame(colnames(KO.data),rep("KO",4105))
WT.metadata<-data.frame(colnames(WT.data),rep("WT",3905))
colnames(KO.metadata)<-c("barcode","group")
colnames(WT.metadata)<-c("barcode","group")
rownames(KO.metadata)<-KO.metadata[,1]
rownames(WT.metadata)<-WT.metadata[,1]


KO <- CreateSeuratObject(counts = KO.data, project = "IMMUNE_KO", min.cells = 5,meta.data = KO.metadata)
KO$type <- "KO"
KO <- subset(KO, subset = nFeature_RNA > 500)
KO <- NormalizeData(KO, verbose = FALSE)
KO <- FindVariableFeatures(KO, selection.method = "vst", nfeatures = 2000)


WT <- CreateSeuratObject(counts = WT.data, project = "IMMUNE_WT", min.cells = 5,meta.data = WT.metadata)
WT$type <- "WT"
WT <- subset(WT, subset = nFeature_RNA > 500)
WT <- NormalizeData(WT, verbose = FALSE)
WT <- FindVariableFeatures(WT, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(KO, WT), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- ScaleData(immune.combined, verbose = FALSE,assay = "RNA")
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "type")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p1
p2
DimPlot(immune.combined, reduction = "umap", split.by = "type")

immune.combined <- RunTSNE(immune.combined, dims = 1:20)
DimPlot(immune.combined, reduction = "tsne",split.by = "type")
DimPlot(immune.combined, reduction = "tsne",group.by  = "type")
DimPlot(immune.combined, reduction = "tsne")
# immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
# write.table(immune.combined.markers_integrated,file = "c:/Users/xjmik/Desktop/immune.combined.markers_integrated.txt",sep = "\t")
# immune.combined.markers_RNA <- FindAllMarkers(immune.combined, only.pos = TRUE, assay = "RNA")
# immune.combined.markers_integrated <- FindAllMarkers(immune.combined, only.pos = TRUE, assay = "integrated")
features<-c("Lef1","Sell","Crem","Tox2","Tox","Zbed2","Lag3","Havcr2","Il2ra","S1pr1","Tnfrsf9","Ctla4","Layn","Stat1","Ifit1","Irf7","Ccr4")
DotPlot(immune.combined, features = features,assay = "RNA") + RotatedAxis()
new.cluster.ids <- c("TNFRSF9+", "TNFRSF9+", "TNFRSF9+", "S1pr1+", "S1pr1+", "TNFRSF9+","TNFRSF9+","ISG+","TNFRSF9+","ISG+","ISG+","TNFRSF9+","TNFRSF9+","TNFRSF9-","ISG+","TNFRSF9+")
                     
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
DimPlot(immune.combined, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
immune.combined@meta.data$celltype<-Idents(immune.combined)
immune.combined<-RunLDA(immune.combined,labels = immune.combined@meta.data$celltype,assay = "RNA",features = rownames(immune.combined))
immune.combined<-RunLDA(immune.combined,labels = immune.combined@meta.data$celltype,assay = "integrated",features = rownames(immune.combined),reduction.name = "LDA_integrated")
immune.combined<-RunTSNE(immune.combined,reduction = "lda",reduction.name = "lda_tsne",dims = 1:3)
immune.combined<-RunTSNE(immune.combined,reduction = "LDA_integrated",reduction.name = "lda_tsne_integrated",dims = 1:3)
