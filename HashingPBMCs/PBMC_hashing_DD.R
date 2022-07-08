library(Seurat)
library(DoubletDecon)
dir_data = "DoubletExperimentsData/CellHashing"
pbmc.umis <- readRDS(paste(dir_data,"/pbmc_umi_mtx.rds",sep=""))

pbmc.htos <- readRDS(paste(dir_data,"/pbmc_hto_mtx.rds",sep=""))

joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])
rownames(pbmc.htos)
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)

pbmc.hashtag <- NormalizeData(pbmc.hashtag)
# Find and scale variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))


# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

table(pbmc.hashtag$HTO_classification.global)

Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:2], ncol = 2)

FeatureScatter(pbmc.hashtag, feature1 = "hto_HTO-A", feature2 = "hto_HTO-B")

Idents(pbmc.hashtag) <- "HTO_classification.global"
VlnPlot(pbmc.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

pbmc.hashtag.subset <- subset(pbmc.hashtag, idents = "Negative", invert = TRUE)

plotCounts <- FeatureScatter(pbmc.hashtag.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotCounts

#pbmc.hashtag.subset <- subset(pbmc.hashtag.subset, subset = nFeature_RNA > 300 )

all.genes <- rownames(pbmc.hashtag.subset)
pbmc.hashtag.subset <- ScaleData(pbmc.hashtag.subset,features = all.genes, vars.to.regress = "nCount_RNA")

pbmc.hashtag.subset <- FindVariableFeatures(pbmc.hashtag.subset,
                                            selection.method = 'mean.var.plot',
                                            mean.cutoff = c(0.025,1 ),
                                            dispersion.cutoff = c(0.65,Inf))


pbmc.hashtag.subset <- RunPCA(pbmc.hashtag.subset, features = VariableFeatures(pbmc.hashtag.subset))

DimPlot(pbmc.hashtag.subset, reduction = "pca")

pbmc.hashtag.subset<- RunUMAP(pbmc.hashtag.subset, dims = 1:10)

DimPlot(pbmc.hashtag.subset, reduction = "umap")

pbmc.hashtag.subset <- FindNeighbors(pbmc.hashtag.subset,features = VariableFeatures(object = pbmc.hashtag.subset), dims = 1:10)
pbmc.hashtag.subset<- FindClusters(pbmc.hashtag.subset,resolution = 0.5)

DimPlot(pbmc.hashtag.subset, reduction = "umap")

pbmc.hashtag.subset.markers <- FindAllMarkers(pbmc.hashtag.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.hashtag.subset.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- pbmc.hashtag.subset.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc.hashtag.subset, features = top10$gene) + NoLegend()

FeaturePlot(pbmc.hashtag.subset, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))


location = "DoubletExperiments/CellHashing/DD_outputFilesNew"
newFiles=Improved_Seurat_Pre_Process(pbmc.hashtag.subset, num_genes=50, write_files=FALSE)

filename="PBMC_CellHashing"
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
                           rhop= 0.8, 
                           write=TRUE, 
                           PMF=TRUE, 
                           useFull=FALSE, 
                           heatmap=FALSE,
                           centroids=TRUE,
                           num_doubs=100, 
                           only50=FALSE,
                           min_uniq=4)
