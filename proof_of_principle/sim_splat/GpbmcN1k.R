rm(list = ls())
library(Seurat)
library(ggplot2)
library(mclust)
library(plyr)
library(dplyr)    
library(splatter)
library(scater)
library(Rtsne)



# Load in the RNA UMI matrix

# Note that this dataset also contains ~5% of mouse cells, which we can use as negative controls
# for the protein measurements. For this reason, the gene expression matrix has HUMAN_ or MOUSE_
# appended to the beginning of each gene.
pbmc.rna <- as.sparse(read.csv(file = "input_data/pbmc_rna.csv", sep = ",", 
                               header = TRUE, row.names = 1))

# Load in the ADT UMI matrix
adt <- read.csv(file = "input_data/pbmc_adt.csv", sep = ",", 
                header = TRUE, row.names = 1)
pbmc.adt <- as.sparse(adt[1:49, ])

#-------------------------------------------------------------------------------------------------
#################################### Analyse RNA-seq ######################################
#-------------------------------------------------------------------------------------------------

pbmc <- CreateSeuratObject(counts = pbmc.rna)

# standard log-normalization
pbmc <- NormalizeData(pbmc)
# get data: m <- GetAssayData(pbmc, slot = 'data')
# write.table(t(as.matrix(m)), file = paste0("data/pbmc_process_rna_count.csv"), row.names = F, col.names = F, sep = ',' )


# choose ~1k variable features
pbmc <- FindVariableFeatures(pbmc)

# standard scaling (no regression)
pbmc <- ScaleData(pbmc)

# Run PCA, select 13 PCs for tSNE visualization and graph-based clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
ElbowPlot(pbmc, ndims = 50)


pbmc <- FindNeighbors(pbmc, dims = 1:25)
pbmc <- RunTSNE(pbmc, dims = 1:25, method = "FIt-SNE")
pbmc <- FindClusters(pbmc, resolution = 0.6)

DimPlot(pbmc, label = TRUE) + NoLegend()


# Find the markers that define each cluster, and use these to annotate the clusters, we use
# max.cells.per.ident to speed up the process
pbmc.rna.markers <- FindAllMarkers(pbmc, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)
top10 <- pbmc.rna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

new.cluster.ids <- c("CD8+ T","CD8+ T 2","Naive CD4+ T", "Memory CD4+ T", "No markers", 
                     "CD 14+ Monocypes","Monocypes","Nk","B", "FCGR3A+ Mono","CD8+ T 3","Memory CD4+ T 2","DC","Doublets?")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# plot
DefaultAssay(pbmc) <- "RNA"
DimPlot(pbmc, label = TRUE) + NoLegend()

# Before we recluster the data on ADT levels, we'll stash the RNA cluster IDs for later
pbmc[["rnaClusterID"]] <- Idents(pbmc)

#-------------------------------------------------------------------------------------------------
################ Add the protein expression levels to the Seurat object ####################
#-------------------------------------------------------------------------------------------------
# We will define an ADT assay, and store raw counts for it

# If you are interested in how these data are internally stored, you can check out the Assay
# class, which is defined in objects.R; note that all single-cell expression data, including RNA
# data, are still stored in Assay objects, and can also be accessed using GetAssayData
pbmc[["ADT"]] <- CreateAssayObject(counts = pbmc.adt)

# Now we can repeat the preprocessing (normalization and scaling) steps that we typically run
# with RNA, but modifying the 'assay' argument.  For CITE-seq data, we do not recommend typical
# LogNormalization. Instead, we use a centered log-ratio (CLR) normalization, computed
# independently for each feature.  This is a slightly improved procedure from the original
# publication, and we will release more advanced versions of CITE-seq normalizations soon.
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
pbmc <- ScaleData(pbmc, assay = "ADT")

# Visualize protein levels on RNA clusters
# in this plot, protein (ADT) levels are on top, and RNA levels are on the bottom
# FeaturePlot(pbmc, features = c("CD3", "CD11c", "CD8", "CD16", "CD3E", "ITGAX", "CD8A", 
#                                "FCGR3A"), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)

# Let's plot CD4 vs CD8 levels in T cells
# tcells <- subset(pbmc, idents = c("Naive CD4 T", "Memory CD4 T", "CD8 T"))
# FeatureScatter(tcells, feature1 = "adt_CD4", feature2 = "adt_CD8")


# Downsample the clusters to a maximum of 300 cells each (makes the heatmap easier to see for
# small clusters)
# pbmc.small <- subset(pbmc, downsample = 300)

# Find protein markers for all clusters, and draw a heatmap
adt.markers <- FindAllMarkers(pbmc, assay = "ADT", only.pos = TRUE)
DoHeatmap(pbmc, features = unique(adt.markers$gene), assay = "ADT", angle = 90) + NoLegend()

## cluster on protein
# Because we're going to be working with the ADT data extensively, we're going to switch the
# default assay to the 'CITE' assay.  This will cause all functions to use ADT data by default,
# rather than requiring us to specify it each time
DefaultAssay(pbmc) <- "ADT"
pbmc <- RunPCA(pbmc, features = rownames(pbmc), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
               verbose = FALSE)
# DimPlot(pbmc, reduction = "pca_adt")

# Since we only have 10 markers, instead of doing PCA, we'll just use a standard euclidean
# distance matrix here.  Also, this provides a good opportunity to demonstrate how to do
# visualization and clustering using a custom distance matrix in Seurat
adt.data <- GetAssayData(pbmc, slot = "data")
adt.dist <- dist(t(adt.data))

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein) levels.
pbmc[["tsne_adt"]] <- RunTSNE(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
pbmc[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
pbmc <- FindClusters(pbmc, resolution = 0.4, graph.name = "adt_snn")

# We can compare the RNA and protein clustering, and use this to annotate the protein clustering
# (we could also of course use FindMarkers)
clustering.table <- table(Idents(pbmc), pbmc$rnaClusterID)
clustering.table

new.cluster.ids <- c("CD8+ T & Nk", "CD4+ T", "CD8+ T", "Monocypes", "No markers", "CD4+ T 2", "Mono & CD8+", "B", "Unknown")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc[["adtClusterID"]] <- Idents(pbmc)


#### SPlatter
adt_labels <- pbmc[["adtClusterID"]]

subset_name <- "CD8+ T"
pbmc_cd8 <- pbmc[,adt_labels==subset_name]
# pbmc_cd8 <- pbmc


cd8_rna <- as.matrix(pbmc_cd8@assays[["RNA"]]@counts)
cd8_adt <- as.matrix(pbmc_cd8@assays[["ADT"]]@counts)

### Simulate data using Splatter.
## ADT
params <- splatEstimate(round(cd8_adt))
params <- setParam(params, "batchCells", 1000)
n_groups <- 5
deprob <- rep(0.08, n_groups)
de.facLoc <- 1.6
# DC: 2, NK: 20, B: 10, Mono: 8, T: 60
group.prob <- c(2, 8, 10, 20, 60)/100
sim <- splatSimulateGroups(params, group.prob = group.prob, de.prob = deprob, de.facLoc = de.facLoc,
                           verbose = FALSE)

## Analysis
sim <- logNormCounts(sim)
dat <- logcounts(sim)
dat[1:4,1:4] <- 300 # add few outliers
dataname <- "GpbmcN1k"
x <- as.factor(sim$Group)
levels(x) <- 1:length(levels(x))
labels <- as.numeric(x)

# write.table(t(dat),file = paste0("/data/hoan/multiomics/sim_splat/data/", dataname, "_adt.csv"),
#             sep = ",", row.names = F, col.names = F)
# write.table(labels,file = paste0("/data/hoan/multiomics/sim_splat/data/", dataname, "_labels.csv"),
#             sep = ",", row.names = F, col.names = F)

### RNA
params2 <- splatEstimate(cd8_rna)
params2 <- setParam(params2, "batchCells", 1000)
n_groups <- 5
deprob2 <- rep(0.02, n_groups)
de.facLoc2 <- c(0.8, 0.8, 0.8, 0.6, 0.6)
sim2 <- splatSimulateGroups(params2, group.prob = group.prob, de.prob = deprob2, de.facLoc = de.facLoc2,
                           verbose = FALSE)
## Analysis
sim2 <- logNormCounts(sim2)
dat2 <- logcounts(sim2)

# write.table(t(dat2),file = paste0("/data/hoan/multiomics/sim_splat/data/", dataname, "_rna.csv"),
#             sep = ",", row.names = F, col.names = F)




