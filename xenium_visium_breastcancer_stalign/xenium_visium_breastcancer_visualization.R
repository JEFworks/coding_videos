## read in positions
xenium = read.csv('xenium_subset.csv.gz')
head(xenium)
xenium.pos <- data.frame(x = xenium$aligned_x, y = xenium$aligned_y)
rownames(xenium.pos) <- xenium$cell_id
plot(xenium.pos, pch='.')

visium = read.csv('visium_subset.csv.gz')
head(visium)
visium.pos <- data.frame(x = visium$aligned_x, y = visium$aligned_y)
rownames(visium.pos) <- visium$barcode
plot(visium.pos, pch=16)

## plot together
plot(xenium.pos, pch=".", col='blue')
points(visium.pos, pch=1, cex=2, col='orange')

## read in gene expression
library(Matrix)
xenium.gexp <- readMM('Xenium/outs/cell_feature_matrix/matrix.mtx.gz')
head(xenium.gexp)
xenium.barcodes <- read.csv('Xenium/outs/cell_feature_matrix/barcodes.tsv.gz', header=FALSE, sep="\t")
head(xenium.barcodes)
xenium.features <- read.csv('Xenium/outs/cell_feature_matrix/features.tsv.gz', header=FALSE, sep="\t")
head(xenium.features)
dim(xenium.gexp)
nrow(xenium.barcodes)
rownames(xenium.gexp) <- xenium.features$V2
colnames(xenium.gexp) <- xenium.barcodes$V1
head(xenium.gexp)

## subset
dim(xenium.pos)
xenium.gexp.subset <- xenium.gexp[, rownames(xenium.pos)]
dim(xenium.gexp.subset)

## plot
library(ggplot2)
df1 <- data.frame(xenium.pos, totgexp = log10(colSums(xenium.gexp.subset)+1))
p1 <- ggplot(df1, aes(x = x, y=y, col=totgexp)) + geom_point() + 
    scale_color_continuous(type="viridis") + theme_void()


visium.gexp <- readMM('Visium/filtered_feature_bc_matrix/matrix.mtx.gz')
head(visium.gexp)
visium.barcodes <- read.csv('Visium/filtered_feature_bc_matrix/barcodes.tsv.gz', header=FALSE, sep="\t")
head(visium.barcodes)
visium.features <- read.csv('Visium/filtered_feature_bc_matrix/features.tsv.gz', header=FALSE, sep="\t")
head(visium.features)
rownames(visium.gexp) <- visium.features$V2
colnames(visium.gexp) <- visium.barcodes$V1
head(visium.gexp)

dim(visium.pos)
visium.gexp.subset <- visium.gexp[, rownames(visium.pos)]
dim(visium.gexp.subset)

library(ggplot2)
df2 <- data.frame(visium.pos, totgexp = log10(colSums(visium.gexp.subset)+1))
p2 <- ggplot(df2, aes(x = x, y=y, col=totgexp)) + geom_point(size=5) + 
    scale_color_continuous(type="viridis") + theme_void()

library(gridExtra)
grid.arrange(p1, p2, nrow=2)

## look at shared gene
genes.shared <- intersect(rownames(visium.gexp.subset), rownames(xenium.gexp.subset))
length(genes.shared)
genes.shared

df1 <- data.frame(xenium.pos, totgexp = log10(xenium.gexp.subset['EGFR',]+1))
p1 <- ggplot(df1, aes(x = x, y=y, col=totgexp)) + geom_point() + 
    scale_color_continuous(type="viridis") + theme_void()
df2 <- data.frame(visium.pos, totgexp = log10(visium.gexp.subset['EGFR',]+1))
p2 <- ggplot(df2, aes(x = x, y=y, col=totgexp)) + geom_point(size=5) + 
    scale_color_continuous(type="viridis") + theme_void()
grid.arrange(p1, p2, nrow=2)

dim(visium.gexp.subset)
dim(xenium.gexp.subset)