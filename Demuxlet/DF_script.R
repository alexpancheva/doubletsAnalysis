library(unixtools)
library(dplyr)
library(Seurat)
library(DoubletFinder)

raw_counts<-read.csv(file="DemuxletData/demuxlet_PBMCs/SeuratInput.csv",sep=",",header=TRUE,row.names =1)
head(raw_counts)

pbmc <- CreateSeuratObject(counts = raw_counts,min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc)

all.genes <- rownames(pbmc)

pbmc <- ScaleData(pbmc,features = all.genes,vars.to.regress = 'nCount_RNA')


plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

pbmc <- FindVariableFeatures(pbmc,selection.method = 'mean.var.plot',
                                      mean.cutoff = c(0.025,1 ),
                                      dispersion.cutoff = c(0.85,Inf))


pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
DimPlot(pbmc, reduction = "pca")

pbmc <- FindNeighbors(pbmc, features = VariableFeatures(pbmc),dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.4)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

FeaturePlot(pbmc, features = c('nCount_RNA'),reduction = 'umap')


#Running DoubletFinder
sweep.res.pmbc <- paramSweep_v3(pbmc, PCs = 1:10, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.pmbc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)


##pk is determined based on param sweep: highest peak of the sweep plot
pkValue <- as.double(levels(bcmvn_pbmc[which.max(bcmvn_pbmc$BCmetric),]$pK)[which.max(bcmvn_pbmc$BCmetric)])


#Testing how the doublet rate affects results
df_rates <- list(0.08,0.109,0.115,0.125)



for (k in df_rates) {
  print (k)
  nExp_poi <- round(k*pbmc@assays$RNA@counts@Dim[2]) 
  nExp_poi.adj <- round(nExp_poi*1)

  pbmc <- doubletFinder_v3(pbmc, PCs = 1:10, pN = 0.25, pK = pkValue, nExp=nExp_poi, reuse.pANN = FALSE, sct = FALSE)
}


DimPlot(pbmc, group.by="DF.classifications_0.25_0.005_1827", reduction="umap", pt.size=0.5)


write.csv(pbmc@meta.data,"DemuxletDataExperiments/metaDataDemuxletSeveralDFrates.csv")

