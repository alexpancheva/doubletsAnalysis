library(Seurat)
library(dplyr)
library(ggplot2)
library(DoubletFinder)

dir_data = "/data/Alex/DoubletExperimentsData/CellHashing"
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
VlnPlot(pbmc.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)+ labs(y="nCount_RNA",title="Distribution of Counts")+ theme(plot.title = element_text(size = 12))          


pbmc.hashtag.subset <- subset(pbmc.hashtag, idents = "Negative", invert = TRUE)
plot1 <- FeatureScatter(pbmc.hashtag.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1


all.genes <- rownames(pbmc.hashtag.subset)
pbmc.hashtag.subset <- ScaleData(pbmc.hashtag.subset,features = all.genes, vars.to.regress = "nCount_RNA")

pbmc.hashtag.subset <- FindVariableFeatures(pbmc.hashtag.subset,
                                            selection.method = 'mean.var.plot',
                                            mean.cutoff = c(0.025,1 ),
                                            dispersion.cutoff = c(0.65,Inf))


pbmc.hashtag.subset <- RunPCA(pbmc.hashtag.subset, features = VariableFeatures(pbmc.hashtag.subset))

DimPlot(pbmc.hashtag.subset, reduction = "pca")

pbmc.hashtag.subset<- RunUMAP(pbmc.hashtag.subset, dims = 1:10)
DimPlot(pbmc.hashtag.subset, reduction = "umap",group.by = "HTO_classification.global")

pbmc.hashtag.subset <- FindNeighbors(pbmc.hashtag.subset,features = VariableFeatures(object = pbmc.hashtag.subset), dims = 1:10)
pbmc.hashtag.subset<- FindClusters(pbmc.hashtag.subset,resolution = 0.5)


pbmc.hashtag.subset.markers <- FindAllMarkers(pbmc.hashtag.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.hashtag.subset.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- pbmc.hashtag.subset.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)



sweep.res.hashtag <- paramSweep_v3(pbmc.hashtag.subset, PCs = 1:10, sct = FALSE)
sweep.stats_hashtag <- summarizeSweep(sweep.res.hashtag, GT = FALSE)
bcmvn_hashtag <- find.pK(sweep.stats_hashtag)
pkValue <- as.double(levels(bcmvn_hashtag[which.max(bcmvn_hashtag$BCmetric),]$pK)[which.max(bcmvn_hashtag$BCmetric)])


homotypic.prop <- modelHomotypic(pbmc.hashtag.subset@meta.data$seurat_clusters)           
nExp_poi <- round(0.20*pbmc.hashtag.subset@assays$RNA@counts@Dim[2])
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

pbmc.hashtag.subset <- doubletFinder_v3(pbmc.hashtag.subset, PCs = 1:10, pN = 0.25, pK = pkValue, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

DimPlot(pbmc.hashtag.subset, group.by="DF.classifications_0.25_0.005_2698", 
        reduction="umap")

DimPlot(pbmc.hashtag.subset, group.by="HTO_classification.global", 
        reduction="umap")

pbmc.hashtag.subset[which(pbmc.hashtag.subset@meta.data$HTO_classification.global == "Singlet" & pbmc.hashtag.subset@meta.data$DF.classifications_0.25_0.005_2698 == "Doublet")] 


write.csv(pbmc.hashtag.subset@meta.data,"cellHashingMetaData.csv")

