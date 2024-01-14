#!/usr/bin/env Rscript

library(png)
library(rjson)
library(tidyverse)
library(ggrepel)
options(bitmapType='cairo')

args = commandArgs(trailingOnly = TRUE)
samplename = args[1]
#samplename = 'CytAssist_11mm_FFPE_Human_Glioblastoma'
#lowres.image = '/SAN/colcc/Kordasti_RNA-Seq_Roberto/Z_spatial_trial/0_spatial/results_spaceranger/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma/outs/spatial/tissue_lowres_image.png'
louvain.res = 0.8
lowres.xdim = 600
lowres.ydim = 600

png.dim = 500
spng.xdim = 750
spng.ydim = 700
spng.dot = 1

target.dir <- paste0(getwd(), '/results_secondary_analysis/', samplename)


lowres.image <- paste0(getwd(), '/results_spaceranger/', samplename, '/', samplename, '/outs/spatial/tissue_lowres_image.png')
positions <- read.csv(paste0(getwd(), '/results_spaceranger/', samplename, '/', samplename, '/outs/spatial/tissue_positions.csv'), header = TRUE, sep = ',', row.names = 1)
positions <- positions[positions[,'in_tissue'] == 1,]
scale.factors <- fromJSON(file = paste0(getwd(), '/results_spaceranger/', samplename, '/', samplename, '/outs/spatial/scalefactors_json.json')) 
final.data <- read.csv(paste0(target.dir, '/tables/', samplename, '_remaining_final.csv'), header = TRUE, sep = ',', row.names = 1)



# RNA Clustering
RNA.column <- paste0('SNN_SCT_res.', louvain.res)
final.data[,RNA.column] <- factor(final.data[,RNA.column], levels = levels(factor(final.data[,RNA.column])))

png(paste0(target.dir, '/plots/', samplename, '_expression_count_RNA_clusters.png'), width = png.dim, height = png.dim)
final.data %>% ggplot(aes_string(x=RNA.column, y='nCount_Spatial', color=RNA.column)) +
	geom_jitter(size = 0.1, width = 0.2) +
	geom_violin(size = 0.5, alpha = 0) +
	theme_minimal() +
	xlab('RNA Cluster') +
	ylab('Gene Expression Count') +
	ggtitle(paste0(samplename, ': Gene Expression Count')) +
	theme(legend.position = 'none',
		plot.title = element_text(hjust = 0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = 'black'))
dev.off()


png(paste0(target.dir, '/plots/', samplename, '_feature_RNA_clusters.png'), width = png.dim, height = png.dim)
final.data %>% ggplot(aes_string(x=RNA.column, y='nFeature_Spatial', color=RNA.column)) +
	geom_jitter(size = 0.1, width = 0.2) +
	geom_violin(size = 0.5, alpha = 0) +
	theme_minimal() +
	xlab('RNA Cluster') +
	ylab('Feature Count') +
	ggtitle(paste0(samplename, ': Feature Count')) +
	theme(legend.position = 'none',
		plot.title = element_text(hjust = 0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = 'black'))
dev.off()


png(paste0(target.dir, '/plots/', samplename, '_mito_RNA_clusters.png'), width = png.dim, height = png.dim)
final.data %>% ggplot(aes_string(x=RNA.column, y='percent.mt', color=RNA.column)) +
	geom_jitter(size = 0.1, width = 0.2) +
	geom_violin(size = 0.5, alpha = 0) +
	theme_minimal() +
	xlab('RNA Cluster') +
	ylab('% Mitochondrial Content') +
	ggtitle(paste0(samplename, ': % Mitochondrial Content')) +
	theme(legend.position = 'none',
		plot.title = element_text(hjust = 0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = 'black'))
dev.off()


png(paste0(target.dir, '/plots/', samplename, '_ribo_RNA_clusters.png'), width = png.dim, height = png.dim)
final.data %>% ggplot(aes_string(x=RNA.column, y='percent.ribo', color=RNA.column)) +
	geom_jitter(size = 0.1, width = 0.2) +
	geom_violin(size = 0.5, alpha = 0) +
	theme_minimal() +
	xlab('RNA Cluster') +
	ylab('% Ribosomal Genes') +
	ggtitle(paste0(samplename, ': % Ribosomal Genes')) +
	theme(legend.position = 'none',
		plot.title = element_text(hjust = 0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = 'black'))
dev.off()


png(paste0(target.dir, '/plots/', samplename, '_POLR_RNA_clusters.png'), width = png.dim, height = png.dim)
final.data %>% ggplot(aes_string(x=RNA.column, y='percent.POLR', color=RNA.column)) +
	geom_jitter(size = 0.1, width = 0.2) +
	geom_violin(size = 0.5, alpha = 0) +
	theme_minimal() +
	xlab('RNA Cluster') +
	ylab('% RNA Polymerase Genes') +
	ggtitle(paste0(samplename, ': % RNA Polymerase Genes')) +
	theme(legend.position = 'none',
		plot.title = element_text(hjust = 0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = 'black'))
dev.off()


png(paste0(target.dir, '/plots/', samplename, '_HSP_RNA_clusters.png'), width = png.dim, height = png.dim)
final.data %>% ggplot(aes_string(x=RNA.column, y='percent.HSP', color=RNA.column)) +
	geom_jitter(size = 0.1, width = 0.2) +
	geom_violin(size = 0.5, alpha = 0) +
	theme_minimal() +
	xlab('RNA Cluster') +
	ylab('% HSP Genes') +
	ggtitle(paste0(samplename, ': % HSP Genes')) +
	theme(legend.position = 'none',
		plot.title = element_text(hjust = 0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = 'black'))
dev.off()



png(paste0(target.dir, '/plots/', samplename, '_HIST_RNA_clusters.png'), width = png.dim, height = png.dim)
final.data %>% ggplot(aes_string(x=RNA.column, y='percent.HIST', color=RNA.column)) +
	geom_jitter(size = 0.1, width = 0.2) +
	geom_violin(size = 0.5, alpha = 0) +
	theme_minimal() +
	xlab('RNA Cluster') +
	ylab('% HIST Genes') +
	ggtitle(paste0(samplename, ': % HIST Genes')) +
	theme(legend.position = 'none',
		plot.title = element_text(hjust = 0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = 'black'))
dev.off()



TSNE_RNA <- read.csv(paste0(target.dir, '/tables/', samplename, '_TSNE_SCT.csv'), header = TRUE, sep = ',', row.names = 1)
TSNE_RNA[,'Cluster'] <- final.data[rownames(TSNE_RNA),paste0('SNN_SCT_res.', louvain.res)]
png(paste0(target.dir, '/plots/', samplename, '_TSNE_RNA_RNA_clusters.png'), width = png.dim, height = png.dim)
TSNE_RNA %>% ggplot(aes(x=TSNESCT_1, y=TSNESCT_2, color=as.factor(Cluster))) +
	geom_point(size = 1) +
	theme_minimal() +
	xlab('t-SNE 1') +
	ylab('t-SNE 2') +
	ggtitle(paste0(samplename, ': t-SNE Plot')) +
	scale_color_discrete(name = 'Cluster') +
	theme(legend.position = 'right',
		plot.title = element_text(hjust=0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = 'black'))
dev.off()


UMAP_RNA <- read.csv(paste0(target.dir, '/tables/', samplename, '_UMAP_SCT.csv'), header = TRUE, sep = ',', row.names = 1)
UMAP_RNA[,'Cluster'] <- final.data[rownames(UMAP_RNA),paste0('SNN_SCT_res.', louvain.res)]
png(paste0(target.dir, '/plots/', samplename, '_UMAP_RNA_RNA_clusters.png'), width = png.dim, height = png.dim)
UMAP_RNA %>% ggplot(aes(x=UMAPSCT_1, y=UMAPSCT_2, color=as.factor(Cluster))) +
	geom_point(size = 1) +
	theme_minimal() +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	ggtitle(paste0(samplename, ': UMAP Plot')) +
	scale_color_discrete(name = 'Cluster') +
	theme(legend.position = 'right',
		plot.title = element_text(hjust=0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = 'black'))
dev.off()




# Spatial Results
grd <- expand.grid(1:lowres.xdim, 1:lowres.ydim)
img <- readPNG(lowres.image)
dim(img) <- c(600 * 600, 3)
img <- cbind(grd, img)
colnames(img) <- c("X", "Y", "R", "G", "B")
img[,'RGB'] <- rgb(img$R, img$G, img$B)

base.image <- ggplot() +
	geom_raster(data=img, aes(x = Y, y = -X, fill = RGB)) +
	scale_fill_identity() +
	theme_void()


# Spot Results
positions[,'x_axis'] <- positions[,'pxl_row_in_fullres'] * scale.factors$tissue_lowres_scalef
positions[,'y_axis'] <- positions[,'pxl_col_in_fullres'] * scale.factors$tissue_lowres_scalef
to.plot <- merge(positions, final.data, by = 0)
colnames(to.plot)[1] <- 'Cell'



png(paste0(target.dir, '/plots/', samplename, '_all_clusters.png'), width = spng.xdim, height = spng.ydim)
base.image + geom_point(data = to.plot, aes(x = y_axis, y = -x_axis, color = .data[[RNA.column]]), size = spng.dot) +
	ggtitle(paste0(samplename, ': Clustering Results')) +
	labs(color = 'Cluster') +
	theme(legend.position = 'right',
		plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
		legend.title=element_text(size = 15),
		legend.text=element_text(size = 15))
dev.off()


for ( i in levels(to.plot[,RNA.column]) ) {

	temp.plot <- to.plot[to.plot[,RNA.column] == i,]
	to.print <- base.image + geom_point(data = temp.plot, aes(x = y_axis, y = -x_axis), colour = 'red', size = spng.dot) +
		theme(legend.position = 'None',
			plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
			legend.title=element_text(size = 15),
			legend.text=element_text(size = 15))

	png(paste0(target.dir, '/plots/', samplename, '_cluster_', i, '.png'), width = lowres.xdim, height = lowres.xdim)
		print(to.print)
	dev.off()


}





unlink(paste0(getwd(), '/Rplots.pdf'))
