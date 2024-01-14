#!/usr/bin/env Rscript

library(tidyverse)
library(ggrepel)
library(gridExtra)
library(Seurat)
library(edgeR)
library(GO.db)
options(bitmapType='cairo')

args = commandArgs(trailingOnly = TRUE)

#samplename <- 'CytAssist_11mm_FFPE_Human_Glioblastoma'
#datadirectory <- 'results_spaceranger/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma/outs'
#species <- 'Human'
#pca.variables <- 2000
#PC.number <- 30

samplename = args[1]
species = args[2]
datadirectory = args[3]
lower.limit = as.numeric(args[4])
upper.limit = as.numeric(args[5])
mito.limit = as.numeric(args[6])
pca.variables = as.numeric(args[7])
PC.number = as.numeric(args[8])
louvain.res = 0.8

target.dir <- paste0(getwd(), "/results_secondary_analysis/", samplename)
#samplemix <- strsplit(samplename, '_')[[1]][1]
#hashtag <- strsplit(samplename, '_')[[1]][2]


# Import data
master <- Load10X_Spatial(data.dir = datadirectory)


if ( species == "Human" ) {
	master[["percent.mt"]] <- PercentageFeatureSet(master, pattern = "^MT-")
	master[["percent.ribo"]] <- PercentageFeatureSet(master, pattern = "^RPL|^RPS")
	master[["percent.IG"]] <- PercentageFeatureSet(master, pattern = "^IGH[AEGVDJC]|^IGK[VDJC]|^IGL[VDJC]")
	master[["percent.TR"]] <- PercentageFeatureSet(master, pattern = "^TRA[VDJC]|^TRB[VDJC]|^TRD[VDJC]|^TRG[VDJC]")
	master[["percent.POLR"]] <- PercentageFeatureSet(master, pattern = "^POLR")
	master[["percent.HSP"]] <- PercentageFeatureSet(master, pattern = "^HSP")
	master[["percent.HIST"]] <- PercentageFeatureSet(master, pattern = "^HIST")
} else if ( species == "Mouse" ) {
	master[["percent.mt"]] <- PercentageFeatureSet(master, pattern = "^mt-")
	master[["percent.ribo"]] <- PercentageFeatureSet(master, pattern = "^Rpl|^Rps")
	master[["percent.IG"]] <- PercentageFeatureSet(master, pattern = "^Igh[aegvdjc]|^Igk[vdjc]|^Igl[vdjc]")
	master[["percent.TR"]] <- PercentageFeatureSet(master, pattern = "^Tra[vdjc]|^Trb[vdjc]|^Trd[vdjc]|^Trg[vdjc]")
	master[["percent.POLR"]] <- PercentageFeatureSet(master, pattern = "^Polr")
	master[["percent.HSP"]] <- PercentageFeatureSet(master, pattern = "^Hsp")
	master[["percent.HIST"]] <- PercentageFeatureSet(master, pattern = "^Hist")
}



# Pre-processing
remaining <- master
remaining <- SCTransform(remaining, assay = "Spatial")



# Running PCA
remaining <- RunPCA(remaining, assay = "SCT", reduction.name = "PCA_SCT", reduction.key = "PCASCT_")
write.csv(remaining@reductions$PCA_SCT@cell.embeddings, file = paste0(target.dir, "/pca_statistics/", samplename, "_PCA_SCT_cell_embeddings.csv"), row.names = TRUE)
write.csv(remaining@reductions$PCA_SCT@feature.loadings, file = paste0(target.dir, "/pca_statistics/", samplename, "_PCA_SCT_feature_loadings.csv"), row.names = TRUE)
write.csv(remaining@reductions$PCA_SCT@stdev, file = paste0(target.dir, "/pca_statistics/", samplename, "_PCA_SCT_standard_deviations.csv"), row.names = TRUE)



# Run t-SNE and UMAP
remaining <- RunTSNE(remaining, reduction = 'PCA_SCT', reduction.name = "TSNE_SCT", reduction.key = "TSNESCT_", dims = 1:PC.number)
write.csv(remaining@reductions$TSNE_SCT@cell.embeddings, file = paste0(target.dir, "/tables/", samplename, "_TSNE_SCT.csv"), row.names = TRUE)

remaining <- RunUMAP(remaining, reduction = 'PCA_SCT', reduction.name = "UMAP_SCT", reduction.key = "UMAPSCT_", dims = 1:PC.number)
write.csv(remaining@reductions$UMAP_SCT@cell.embeddings, file = paste0(target.dir, "/tables/", samplename, "_UMAP_SCT.csv"), row.names = TRUE)



## METHOD 1
# Clustering for RNA (SC-Trasformed)
sink(paste0(target.dir, "/clustering_SCT/", samplename, "_SCT_clustering_statistics.txt"))
remaining <- FindNeighbors(remaining, dims = 1:PC.number, reduction = 'PCA_SCT', graph.name = 'SNN_SCT')
remaining <- FindClusters(remaining, resolution = louvain.res, graph.name = 'SNN_SCT')
sink()

for ( clust in levels(remaining@meta.data[,paste0('SNN_SCT_res.', louvain.res)]) ) {
	write.csv(remaining@meta.data[remaining@meta.data[,paste0('SNN_SCT_res.', louvain.res)] == clust,], file = paste0(target.dir, "/clustering_SCT/", samplename, "_SCT_cluster_", clust, ".csv"))
}



# Identify marker GENES for RNA clusters
top.genes <- c()
for ( clust in levels(remaining@meta.data[,paste0('SNN_SCT_res.', louvain.res)]) ) {
	if ( summary(remaining@meta.data[,paste0('SNN_SCT_res.', louvain.res)])[[clust]] > 9 ) {
		markers <- FindMarkers(remaining, assay = 'SCT', ident.1 = as.numeric(clust), min.pct = 0.25, logfc.threshold = 0.1)
#	       markers <- markers[order(abs(markers$avg_log2FC)),]
		top.genes <- append(top.genes, dimnames(markers)[[1]][1])

		if ( dim(markers[markers["avg_log2FC"] > 0,])[1] > 0 ) {
			write.csv(markers[markers["avg_log2FC"] > 0,], file = paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "_positive.csv"))
		} else {
			sink(paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "_positive.csv"))
			cat("No positive avg_log2FC for this cluster.")
			sink()
		}

		if ( dim(markers[markers["avg_log2FC"] < 0,])[1] > 0 ) {
			write.csv(markers[markers["avg_log2FC"] < 0,], file = paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "_negative.csv"))
		} else {
			sink(paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "_negative.csv"))
			cat("No negative avg_log2FC for this cluster.")
			sink()

		}

	} else {
		top.genes <- append(top.genes, "unknown")
		sink(paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "_positive.csv"))
		cat("Cluster is too small to call.")
		sink()
		sink(paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "_negative.csv"))
		cat("Cluster is too small to call.")
		sink()

	}

}



# Pairwise comparisons of RNA clusters for RNA expression
clust.levels <- levels(remaining@meta.data[,paste0('SNN_SCT_res.', louvain.res)])
for ( clust in clust.levels ) {
	for ( clust.2 in clust.levels[! clust.levels %in% clust] ) {
		if ( summary(remaining@meta.data[,paste0('SNN_SCT_res.', louvain.res)])[[clust]] > 9 ) {
			markers <- FindMarkers(remaining, assay = 'SCT', ident.1 = as.numeric(clust), ident.2 = as.numeric(clust.2), min.pct = 0.25, logfc.threshold = 0.1)

			if ( dim(markers[markers["avg_log2FC"] > 0,])[1] > 0 ) {
				write.csv(markers[markers["avg_log2FC"] > 0,], file = paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "vs", clust.2, "_positive.csv"))
			} else {
				sink(paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "vs", clust.2, "_positive.csv"))
				cat("No positive avg_log2FC for this comparison.")
				sink()
			}

			if ( dim(markers[markers["avg_log2FC"] < 0,])[1] > 0 ) {
				write.csv(markers[markers["avg_log2FC"] < 0,], file = paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "vs", clust.2, "_negative.csv"))
			} else {
				sink(paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "vs", clust.2, "_negative.csv"))
				cat("No negative avg_log2FC for this comparison.")
				sink()
			}
		} else {
			sink(paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "vs", clust.2, "_positive.csv"))
			cat("Cluster is too small to call.")
			sink()
			sink(paste0(target.dir, "/diff_exp_RNA/", samplename, "_RNA_cluster_", clust, "vs", clust.2, "_negative.csv"))
			cat("Cluster is too small to call.")
			sink()
		}
	}
}



# Identify most and least variable genes in each cluster
for ( clust in levels(remaining@meta.data[,paste0('SNN_SCT_res.', louvain.res)]) ) {
	subpopulation <- subset(remaining, subset = seurat_clusters == clust)
	subpopulation <- FindVariableFeatures(subpopulation, nfeatures = pca.variables)
	top10 <- head(VariableFeatures(subpopulation), 10)
	top50 <- head(VariableFeatures(subpopulation), 50)
	bottom50 <- tail(VariableFeatures(subpopulation), 50)
	var.feat.plot <- VariableFeaturePlot(subpopulation)
	var.feat.plot <- LabelPoints(plot = var.feat.plot, points = top10, repel = TRUE)
	write.csv(as.data.frame(top50), file = paste0(target.dir, "/clustering_SCT/", samplename, "_RNA_cluster_", clust, "_top_50_variable_genes.csv"), row.names = FALSE)
	write.csv(as.data.frame(bottom50), file = paste0(target.dir, "/clustering_SCT/", samplename, "_RNA_cluster_", clust, "_bottom_50_variable_genes.csv"), row.names = FALSE)
	png(paste(target.dir, "/clustering_SCT/", samplename, "_RNA_cluster_", clust, "_variable_features.png", sep=""))
		print(var.feat.plot)
	dev.off()
}





## METHOD 2
#if ( species == "Human" ) {
#	features.to.use <- head(grep("^IGH[AEGVDJC]|^IGK[VDJC]|^IGL[VDJC]|^TRA[VDJC]|^TRB[VDJC]|^TRD[VDJC]|^TRG[VDJC]", VariableFeatures(remaining), value = TRUE, invert = TRUE), pca.variables)
#} else if ( species == "Mouse" ) {
#	features.to.use <- head(grep("^Igh[aegvdjc]|^Igk[vdjc]|^Igl[vdjc]|^Tra[vdjc]|^Trb[vdjc]|^Trd[vdjc]|^Trg[vdjc]", VariableFeatures(remaining), value = TRUE, invert = TRUE), pca.variables)
#}
#
#remaining <- FindSpatiallyVariableFeatures(remaining, assay = "SCT", features = features.to.use, selection.method = "moransi")






write.csv(remaining@meta.data, file = paste0(target.dir, "/tables/", samplename, "_remaining_final.csv"), row.names = TRUE)

saveRDS(remaining, file = paste0(target.dir, "/", samplename, "_Seurat_object.rds"))

unlink(paste0(getwd(), '/Rplots.pdf'))
