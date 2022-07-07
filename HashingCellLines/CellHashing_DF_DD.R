library(Seurat)
library(DoubletFinder)

# Read in UMI count matrix for RNA
hto12.umis <- readRDS("DoubletExperimentsData/CellHashingCellLines/hto12_umi_mtx.rds")

# Read in HTO count matrix
hto12.htos <- readRDS("DoubletExperimentsData/CellHashingCellLines/hto12_hto_mtx.rds")

# Select cell barcodes detected in both RNA and HTO
cells.use <- intersect(rownames(hto12.htos), colnames(hto12.umis))

# Create Seurat object and add HTO data
hto12 <- CreateSeuratObject(counts = hto12.umis[, cells.use], min.features = 300)
hto12[["HTO"]] <- CreateAssayObject(counts = t(x = hto12.htos[colnames(hto12), 1:12]))

# Normalize data
hto12 <- NormalizeData(hto12)
hto12 <- NormalizeData(hto12, assay = "HTO", normalization.method = "CLR")

hto12 <- HTODemux(hto12, assay = "HTO", positive.quantile = 0.99)

RidgePlot(hto12, assay = "HTO", features = c("HEK-A", "K562-B", "KG1-A", "THP1-C"), ncol = 2)

HTOHeatmap(hto12, assay = "HTO")

Idents(hto12) <- "HTO_classification.global"
VlnPlot(hto12, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Remove the negative cells
hto12 <- subset(hto12, idents = "Negative", invert = TRUE)


# Run PCA on most variable features
hto12 <- FindVariableFeatures(hto12, selection.method = "mean.var.plot")
hto12 <- ScaleData(hto12, features = VariableFeatures(hto12))

hto12 <- RunPCA(hto12)

hto12 <- RunTSNE(hto12, dims = 1:5, perplexity = 100)
DimPlot(hto12, reduction = "tsne") + NoLegend()

hto12 <- FindNeighbors(hto12,features = VariableFeatures(object = hto12), dims = 1:5)
hto12 <- FindClusters(hto12,resolution = 0.1)

hto12<- RunUMAP(hto12, dims = 1:5)
DimPlot(hto12, reduction = "umap",group.by = "HTO_classification.global")
FeaturePlot(hto12, features = c('nCount_RNA'),reduction = 'umap')

#write.csv(hto12[["RNA"]]@counts,'cellHashingCountsCellLines.csv')

#write.csv(hto12@meta.data,"CellHashingCellLines_meta.csv")

hto12.markers <- FindAllMarkers(hto12, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
hto12 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- hto12.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(hto12, features = top10$gene) + NoLegend()

sweep.res.hashtag <- paramSweep_v3(hto12, PCs = 1:10, sct = FALSE)
sweep.stats_hashtag <- summarizeSweep(sweep.res.hashtag, GT = FALSE)
bcmvn_hashtag <- find.pK(sweep.stats_hashtag)
pkValue<- as.double(levels(bcmvn_hashtag[which.max(bcmvn_hashtag$BCmetric),]$pK)[which.max(bcmvn_hashtag$BCmetric)])



homotypic.prop <- modelHomotypic(hto12@meta.data$seurat_clusters)         
nExp_poi <- round(0.20*hto12@assays$RNA@counts@Dim[2])
#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

hto12 <- doubletFinder_v3(hto12, PCs = 1:10, pN = 0.25, pK = pkValue, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

hto12[which(hto12@meta.data$HTO_classification.global == "Singlet" & hto12@meta.data$DF.classifications_0.25_0.06_1591 == "Doublet")] 

write.csv(hto12@meta.data,"CellHashingCellLines_metaDFResults.csv")

#write.csv(hto12@reductions$umap@cell.embeddings,"UMAP_embed.csv")


location = "DoubletExperiments/CellHashingCellLines/"
newFiles=Improved_Seurat_Pre_Process(hto12, num_genes=50, write_files=FALSE)

filename="CellLines_CellHashing2"
write.table(newFiles$newExpressionFile, paste0(location, filename, "_expression"), sep="\t")
write.table(newFiles$newFullExpressionFile, paste0(location, filename, "_fullExpression"), sep="\t")
write.table(newFiles$newGroupsFile, paste0(location, filename , "_groups"), sep="\t", col.names = F)

results=Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                           groupsFile=newFiles$newGroupsFile, 
                           filename=filename, 
                           location=location,
                           fullDataFile=NULL, 
                           removeCC=FALSE, 
                           species="hsa", 
                           rhop=1.1, 
                           write=TRUE, 
                           PMF=TRUE, 
                           useFull=FALSE, 
                           heatmap=FALSE,
                           centroids=TRUE,
                           num_doubs=100, 
                           only50=FALSE,
                           min_uniq=10)







