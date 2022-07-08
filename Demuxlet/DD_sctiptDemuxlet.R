library(dplyr)
library(Seurat)
library(DoubletDecon)

raw_counts<-read.csv(file="/DemuxletData/demuxlet_PBMCs/SeuratInput.csv",sep=",",header=TRUE,row.names =1)

head(raw_counts)

pbmc <- CreateSeuratObject(counts = raw_counts,min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc)

all.genes <- rownames(pbmc)

pbmc <- FindVariableFeatures(pbmc,selection.method = 'mean.var.plot')

pbmc <- ScaleData(pbmc,vars.to.regress = 'nCount_RNA')


pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
DimPlot(pbmc, reduction = "pca")

pbmc <- FindNeighbors(pbmc, features = VariableFeatures(pbmc),dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.2)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")


location = "DoubletExperiments/DemuxletDataExperiments/"
newFiles=Improved_Seurat_Pre_Process(pbmc, num_genes=50, write_files=FALSE)

filename="subSample_DEMUXLET"
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
                           rhop=0.8, 
                           write=TRUE, 
                           PMF=TRUE, 
                           useFull=FALSE, 
                           heatmap=FALSE,
                           centroids=TRUE,
                           num_doubs=100, 
                           only50=FALSE,
                           min_uniq=4)

FeaturePlot(pbmc, features = c('nCount_RNA'),reduction = 'umap')

