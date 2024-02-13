library(DESeq2)
library(IHW)
library(ggplot2)

sampletable <- read.table("vul_WT_metadata.txt", header=T, sep="\t")
rownames(sampletable) <- sampletable$SampleName
count_matrix <- as.matrix(read.delim("vul_WT_gene_count.txt", header=T, sep="\t", row.names=1))
count_matrix<-round(count_matrix)
# Set the first column as row names
row.names(sampletable) <- sampletable$Name
count_matrix <- count_matrix[, grep("^(HBD2|Control)", colnames(count_matrix))]
sampletable <- sampletable[grep("^(HBD2|Control)", rownames(sampletable)), ]
se_star_matrix <- DESeqDataSetFromMatrix(countData = count_matrix,
                                         colData = sampletable,
                                         design = ~Status)
#tx2gene <- read.table("tx2gene.gencode.v29.csv", sep="\t", header=F)

keep <- rowSums(counts(se_star_matrix) >= 10) >= 1
dds <- se_star_matrix[keep,]
nrow(dds)

####variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
star_matrix_vsd_df <- data.frame(assay(vsd))
write.table(star_matrix_vsd_df, file="star_count_matrix_vsd.txt", sep = "\t")

####regularized-logarithm transformation or rlog
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
star_matrix_rld_df <- data.frame(assay(rld))
#####Plot Normalization
library("hexbin")
library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 5:6]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 5:6]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 5:6]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

####### Plot heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDists

library("pheatmap")
library("RColorBrewer")

jpeg("Heatmap_samples.jpg", height = 7, width = 7, units = 'in', res = 600)
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$Names, vsd$Status, sep = " - " )
colnames(sampleDistMatrix) <- paste( vsd$Names, vsd$Status, sep = " - " )
colors <- colorRampPalette( rev(brewer.pal(9, "BrBG")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
####### Plot heatmap using Poisson Distance 
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

jpeg("Heatmap_Poisson.jpg", height = 7, width = 7, units = 'in', res = 600)
samplePoisDistMatrix <- as.matrix( poisd$dd)
rownames(samplePoisDistMatrix) <- paste( vsd$Names, vsd$Status, sep=" - " )
colnames(samplePoisDistMatrix) <- paste( vsd$Names, vsd$Status, sep=" - " )
pheatmap(
  samplePoisDistMatrix,
  clustering_distance_rows = poisd$dd,
  clustering_distance_cols = poisd$dd,
  col = colorRampPalette(c("blue", "white", "red"))(100)  # Example color palette
)
dev.off()
###### PCA using vst-data
#plotPCA(vsd, intgroup = c("Status"))

jpeg("PCA_vst2.jpg", height = 7, width = 7, units = 'in', res = 600)

pcaData <- plotPCA(vsd, intgroup = c("Status"), returnData = TRUE)
#pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Status, shape = Status, label = name, group= Status)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + theme_bw() + 
  ggtitle("PCA with VST data") + geom_text()

dev.off()
###### PCA using rld-data
jpeg("PCA_rld.jpg", height = 7, width = 7, units = 'in', res = 600)

pcaData <- plotPCA(rld, intgroup = c("Status"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Status, shape = Status, label = name, group= Status)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + theme_bw() +
  ggtitle("PCA with RLD data") + geom_text()

dev.off()

#### Differential expression
dds <- DESeq(dds)

res <- results(dds, contrast=c("Status","HBD2","Control"), filterFun=ihw)
mcols(res, use.names = TRUE)
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
write.csv(res, file = "Control_HBD2_DE.csv")

sum(res$pvalue < 0.05, na.rm=TRUE)
sum(res$pvalue < 0.01, na.rm=TRUE)

sum(!is.na(res$pvalue))

sum(res$padj < 0.1, na.rm=TRUE)


resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])


topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Status"))
###################################
# Plot selected genes
###################################
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)


mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Status")])
pheatmap(mat)

######

library("IHW")
library("AnnotationDbi")
library("org.Mm.eg.db")

####volcano plot
library(EnhancedVolcano)
library(DESeq2)
library(IHW)
library(ggplot2)

sampletable <- read.table("vul_WT_metadata.txt", header=T, sep="\t")
rownames(sampletable) <- sampletable$SampleName
count_matrix <- as.matrix(read.delim("vul_WT_gene_count.txt", header=T, sep="\t", row.names=1))
count_matrix<-round(count_matrix)
# Set the first column as row names
row.names(sampletable) <- sampletable$Name
count_matrix <- count_matrix[, grep("^(HBD2|Control)", colnames(count_matrix))]
sampletable <- sampletable[grep("^(HBD2|Control)", rownames(sampletable)), ]
se_star_matrix <- DESeqDataSetFromMatrix(countData = count_matrix,
                                         colData = sampletable,
                                         design = ~Status)

keep <- rowSums(counts(se_star_matrix) >= 10) >= 1
dds <- se_star_matrix[keep,]
nrow(dds)

vsd <- vst(dds, blind = FALSE)
dds <- DESeq(dds)
res <- results(dds, contrast=c("Status","HBD2","Control"), filterFun=ihw)
jpeg("VolcanoPlot_HF_Chow_WSD.jpg", height = 20, width = 20, units = 'in', res = 600)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                pCutoff = 0.05,
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)
dev.off()

pdf(file="VolcanoPlot_HF_Chow_WSD.pdf", height = 20, width = 20)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                pCutoff = 0.05,
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)
dev.off()
