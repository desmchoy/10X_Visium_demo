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
#pca.variables <- 3000
#PC.number <- 30

samplename = args[1]
species = args[2]
datadirectory = args[3]
lower.limit = 200
upper.limit = 2500
mito.limit = 10
pca.variables = 3000

lowres.xdim = 600
lowres.ydim = 600

png.dim = 700
spng.xdim = 800
spng.ydim = 700
spng.dot = 1
plot.dim = 10
plot.dpi = 200

target.dir <- paste0(getwd(), "/results_preliminary_analysis/", samplename)



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



# Spatial Results
grd <- expand.grid(1:lowres.xdim, 1:lowres.ydim)
img <- master@images$slice1@image
dim(img) <- c(600 * 600, 3)
img <- cbind(grd, img)
colnames(img) <- c("X", "Y", "R", "G", "B")
img[,'RGB'] <- rgb(img$R, img$G, img$B)

base.image <- ggplot() +
	geom_raster(data=img, aes(x = Y, y = -X, fill = RGB)) +
	scale_fill_identity() +
	theme_void()



# Spot Results
to.plot <- data.frame(master@images$slice1@coordinates)
to.plot[,'x_axis'] <- to.plot[,'imagerow'] * master@images$slice1@scale.factors$lowres
to.plot[,'y_axis'] <- to.plot[,'imagecol'] * master@images$slice1@scale.factors$lowres
to.plot <- merge(to.plot, master@meta.data, by = 0)
colnames(to.plot)[1] <- 'Cell'



plot.overlay <- function ( column ) {

	plot <- base.image +
		geom_point(data = to.plot, aes(x=y_axis, y=-x_axis, color = .data[[column]])) +
		geom_point(size = 0.1) +

	return ( plot )

}





# Mitochondrial Content
mito.grouping <- cut(to.plot[,'percent.mt'], c(0, 5, 10, 15, 20, 100), include.lowest = TRUE)
to.plot[,'mito.group'] <- mito.grouping
mito.summary <- data.frame(table(mito.grouping))
mito.summary[,1] <- c('0 - 5%', '5 - 10%', '10 - 15%', '15 - 20%', '20 - 100%')
write.csv(mito.summary, file = paste0(target.dir, "/tables/", samplename, "_mito_summary.csv"), row.names = FALSE, col.names = FALSE)

mito.colours <- c('black', 'red', 'orange', 'green', 'blue')
png(paste0(target.dir, '/plots/', samplename, '_mito.png'), width = spng.xdim, height = spng.ydim)
base.image + geom_point(data = to.plot, aes(x = y_axis, y = -x_axis, color = mito.group), size = spng.dot) +
	ggtitle(paste0(samplename, ': Mitochondrial Content')) +
	scale_color_manual(name = "% MT content",
		values = mito.colours,
		labels = c("0 < % <= 5", "5 < % <= 10", "10 < % <= 15", "15 < % <= 20", "20 < % <= 100")) +
	theme(legend.position = 'right',
		plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
		legend.title=element_text(size = 15),
		legend.text=element_text(size = 15))
dev.off()

png(paste0(target.dir, '/plots/', samplename, '_violin_mito.png'), width = png.dim, height = png.dim)
to.plot %>% ggplot(aes(x = 'Sample', y = percent.mt)) +
	geom_jitter(aes(color = mito.group), width = 0.2, size = 0.1) +
	geom_violin(size = 0.5, color = 'red', alpha = 0) +
	theme_minimal() +
	xlab('') +
	ylab('% Mitochondrial Content') +
	ggtitle(paste0(samplename, ': ', '% Mitochondrial Content')) +
	scale_color_manual(name = "% MT content",
	values = mito.colours,
		labels = c("0 < % <= 5", "5 < % <= 10", "10 < % <= 15", "15 < % <= 20", "20 < % <= 100")) +
	theme(legend.position = 'right',
	plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = "black"),
		legend.title=element_text(size = 15),
		legend.text=element_text(size = 15))
dev.off()





# Feature Count
feature.grouping <- cut(to.plot[,'nFeature_Spatial'], c(0, 100, 200, 500, 2000, 2500, 3000, 4000, Inf), include.lowest = TRUE)
to.plot[,'feature.group'] <- feature.grouping
feature.summary <- data.frame(table(feature.grouping))
feature.summary[,1] <- c('0 - 100', '100 - 200', '200 - 500', '500 - 2000', '2000 - 2500', '2500 - 3000', '3000 - 4000', '> 4000')
write.csv(feature.summary, file = paste0(target.dir, "/tables/", samplename, "_feature_summary.csv"), row.names = FALSE, col.names = FALSE)


feature.colours <- c("blue", "blue3", "blueviolet", "black", "black", "rosybrown", "red3", "red")
png(paste0(target.dir, '/plots/', samplename, '_feature.png'), width = spng.xdim, height = spng.ydim)
base.image + geom_point(data = to.plot, aes(x = y_axis, y = -x_axis, color = feature.group), size = spng.dot) +
	ggtitle(paste0(samplename, ': Feature Count')) +
	scale_color_manual(name = "Number of Features",
		values = feature.colours,
		labels = c('0 - 100', '100 - 200', '200 - 500', '500 - 2000', '2000 - 2500', '2500 - 3000', '3000 - 4000', '> 4000')) +
	theme(legend.position = 'right',
		plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
		legend.title=element_text(size = 15),
		legend.text=element_text(size = 15))
dev.off()


png(paste0(target.dir, '/plots/', samplename, '_violin_feature.png'), width = png.dim, height = png.dim)
to.plot %>% ggplot(aes(x = 'Sample', y = nFeature_Spatial)) +
	geom_jitter(aes(color = feature.group), width = 0.2, size = 0.1) +
	geom_violin(size = 0.5, color = 'red', alpha = 0) +
	theme_minimal() +
	xlab('') +
	ylab('Feature Count') +
	ggtitle(paste0(samplename, ': ', 'Feature Count')) +
	scale_color_manual(name = "Number of Features",
	values = feature.colours,
		labels = c('0 - 100', '100 - 200', '200 - 500', '500 - 2000', '2000 - 2500', '2500 - 3000', '3000 - 4000', '> 4000')) +
	theme(legend.position = 'right',
	plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = "black"),
		legend.title=element_text(size = 15),
		legend.text=element_text(size = 15))
dev.off()



# Immune Contributions
IG.subset <- to.plot %>% filter(percent.IG > 0)
IG.grouping <- ifelse(to.plot[,'percent.IG'] > 0, 'Yes', 'No')
to.plot[,'IG.group'] <- IG.grouping
IG.summary <- data.frame(table(IG.grouping))
IG.summary <- IG.summary[order(IG.summary[,1]),]
write.csv(IG.summary, file = paste0(target.dir, "/tables/", samplename, "_IG_summary.csv"), row.names = FALSE, col.names = FALSE)


TR.subset <- to.plot %>% filter(percent.TR > 0)
TR.grouping <- ifelse(to.plot[,'percent.TR'] > 0, 'Yes', 'No')
to.plot[,'TR.group'] <- TR.grouping
TR.summary <- data.frame(table(TR.grouping))
TR.summary <- TR.summary[order(TR.summary[,1]),]
write.csv(TR.summary, file = paste0(target.dir, "/tables/", samplename, "_TR_summary.csv"), row.names = FALSE, col.names = FALSE)


plot_scan <- function ( to.highlight, colour.choice, plot.title ) {

	base.image + geom_point(data = to.plot, aes(x = y_axis, y = -x_axis), color = 'grey', size = spng.dot) +
		geom_point(data = to.highlight, aes(x = y_axis, y = -x_axis), color = colour.choice, size = spng.dot) +
		ggtitle(paste0(samplename, ': ', plot.title)) +
		theme(legend.position = 'none',
			plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))

}

plot_scatter <- function ( column, plot.title, y.label ) {
	to.plot %>% ggplot(aes(x = 'Sample', y = .data[[column]])) +
	geom_violin(size = 0.5, color = 'red', alpha = 0) +
	geom_jitter(width = 0.2, size = 0.1) +
	theme_minimal() +
	xlab('') +
	ylab(y.label) +
	ggtitle(paste0(samplename, ': ', plot.title)) +
	theme(legend.position = 'none',
		plot.title = element_text(hjust=0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = "black"))
}



png(paste0(target.dir, '/plots/', samplename, '_IG.png'), width = spng.xdim, height = spng.ydim)
plot_scan(IG.subset, 'red', 'IG expressing')
dev.off()
png(paste0(target.dir, '/plots/', samplename, '_TR.png'), width = spng.xdim, height = spng.ydim)
plot_scan(TR.subset, 'red', 'TR expressing')
dev.off()

png(paste0(target.dir, '/plots/', samplename, '_violin_IG.png'), width = png.dim, height = png.dim)
plot_scatter('percent.IG', 'IG expressing', '% IG genes')
dev.off()
png(paste0(target.dir, '/plots/', samplename, '_violin_TR.png'), width = png.dim, height = png.dim)
plot_scatter('percent.TR', 'TR expressing', '% TR genes')
dev.off()





# Genes QC
genes.to.plot <- c('CD2', 'CD3D', 'CD3E', 'CD3G', 'CD4', 'FOXP3', 'IKZF2', 'IL2RA', 'PTPRC', 'NELL2', 'ANK3')
rows.to.extract <- which(rownames(master@assays$Spatial@features@.Data) %in% genes.to.plot)

subset.frame <- as.data.frame(master@assays$Spatial@layers$counts[rows.to.extract,])
rownames(subset.frame) <- names(master@assays$Spatial@features@.Data[rows.to.extract,])
subset.frame <- as.data.frame(t(subset.frame))
subset.frame[,'Cell'] <- rownames(subset.frame)
write.csv(subset.frame, file = paste0(target.dir, "/tables/", samplename, "_selected_gene_counts.csv"), row.names = TRUE, col.names = NA)

png(paste0(target.dir, '/plots/', samplename, '_selected_gene_counts.png'), width = png.dim, height = png.dim)
subset.frame %>% pivot_longer(!Cell, names_to = 'Gene', values_to = 'Count') %>% ggplot(aes(x = Gene, y = Count)) +
	geom_jitter(aes(color = Gene), width = 0.2, size = 0.1) +
	geom_violin(size = 0.5, color = 'red', alpha = 0) +
	theme_minimal() +
	xlab('') +
	ylab('Normalised Gene Count') +
	ggtitle(paste0(samplename, ': ', 'Selected Gene Counts')) +
	theme(legend.position = 'none',
	plot.title = element_text(hjust=0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = "black"))
dev.off()





# Find Variable Features
master <- FindVariableFeatures(master, nfeatures = 5000)

var.genes <- master@assays$Spatial@meta.data
var.genes <- var.genes[rowSums(is.na(var.genes)) == 0,]
var.genes <- var.genes[order(var.genes[,'var.features.rank']),]
rownames(var.genes) <- NULL


top20 <- head(var.genes[,'var.features'], 20)
if ( species == "Human" ) {
	top20.no.immune <- head(grep("^IGH[AEGVDJC]|^IGK[VDJC]|^IGL[VDJC]|^TRA[VDJC]|^TRB[VDJC]|^TRD[VDJC]|^TRG[VDJC]", var.genes[,'var.features'], value = TRUE, invert = TRUE), 20)
} else if ( species == "Mouse" ) {
	top20.no.immune <- head(grep("^Igh[aegvdjc]|^Igk[vdjc]|^Igl[vdjc]|^Tra[vdjc]|^Trb[vdjc]|^Trd[vdjc]|^Trg[vdjc]", var.genes[,'var.features'], value = TRUE, invert = TRUE), 20)
}


top100 <- head(var.genes[,'var.features'], 100)
if ( species == "Human" ) {
	top100.no.immune <- head(grep("^IGH[AEGVDJC]|^IGK[VDJC]|^IGL[VDJC]|^TRA[VDJC]|^TRB[VDJC]|^TRD[VDJC]|^TRG[VDJC]", var.genes[,'var.features'], value = TRUE, invert = TRUE), 100)
} else if ( species == "Mouse" ) {
	top100.no.immune <- head(grep("^Igh[aegvdjc]|^Igk[vdjc]|^Igl[vdjc]|^Tra[vdjc]|^Trb[vdjc]|^Trd[vdjc]|^Trg[vdjc]", var.genes[,'var.features'], value = TRUE, invert = TRUE), 100)
}
write.csv(as.data.frame(top100), file = paste0(target.dir, "/tables/", samplename, "_top_100_variable_genes.csv"), row.names = FALSE)
write.csv(as.data.frame(top100.no.immune), file = paste0(target.dir, "/tables/", samplename, "_top_100_variable_genes_no_immune.csv"), row.names = FALSE)



temp.subset <- var.genes[var.genes[,'var.features'] %in% top20,]
png(paste0(target.dir, '/plots/', samplename, '_variable_genes.png'), width = png.dim, height = png.dim)
var.genes %>% ggplot(aes(x = vf_vst_counts_mean, y = vf_vst_counts_variance.standardized, label = var.features)) +
	geom_point(size = 0.1) +
	geom_point(data = temp.subset, color = 'red') +
	geom_text_repel(data = temp.subset, max.overlaps = 1000) +
	scale_x_continuous(trans = 'log10') +
	theme_minimal() +
	xlab('Average Expression') +
	ylab('Standardised Variance') +
	ggtitle(paste0(samplename, ': ', 'Most Variable Genes')) +
	theme(legend.position = 'none',
		plot.title = element_text(hjust = 0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = "black"))
dev.off()


temp.subset <- var.genes[var.genes[,'var.features'] %in% top20.no.immune,]
png(paste0(target.dir, '/plots/', samplename, '_variable_genes_no_immune.png'), width = png.dim, height = png.dim)
var.genes %>% ggplot(aes(x = vf_vst_counts_mean, y = vf_vst_counts_variance.standardized, label = var.features)) +
	geom_point(size = 0.1) +
	geom_point(data = temp.subset, color = 'red') +
	geom_text_repel(data = temp.subset, max.overlaps = 1000) +
	scale_x_continuous(trans = 'log10') +
	theme_minimal() +
	xlab('Average Expression') +
	ylab('Standardised Variance') +
	ggtitle(paste0(samplename, ': ', 'Most Variable Genes (excluding VDJ Genes)')) +
	theme(legend.position = 'none',
		plot.title = element_text(hjust = 0.5),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = "black"))
dev.off()


top100.GO <- goana(top100, species = 'Hs')
write.csv(top100.GO, file = paste0(target.dir, "/tables/", samplename, "_top_100_variable_genes_GO_terms.csv"), row.names = TRUE, col.names = NA)

top100.no.immune.GO <- goana(top100.no.immune, species = 'Hs')
write.csv(top100.no.immune.GO, file = paste0(target.dir, "/tables/", samplename, "_top_100_variable_genes_no_immune_GO_terms.csv"), row.names = TRUE, col.names = NA)



unlink(paste0(getwd(), '/Rplots.pdf'))
