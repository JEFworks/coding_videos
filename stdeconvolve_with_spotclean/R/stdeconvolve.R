library(STdeconvolve)
library(SpotClean)

#mbrain_raw <- read10xRaw("data/raw_feature_bc_matrix/")
#dim(mbrain_raw)
#head(mbrain_raw)

library(Matrix)
countmat <- Matrix::readMM("data/raw_feature_bc_matrix/matrix.mtx.gz")
dim(countmat)
barcode.info <- read.csv("data/raw_feature_bc_matrix/barcodes.tsv.gz", header=FALSE)
head(barcode.info)
barcode <- barcode.info[,1]
head(barcode)
gene.info <- read.csv("data/raw_feature_bc_matrix/features.tsv.gz", header=FALSE, sep="\t")
head(gene.info)
gene <- gene.info[,2]
dim(countmat)
length(barcode)
length(gene)
colnames(countmat) <- barcode
rownames(countmat) <- make.unique(gene)

pos.info <- read.csv('data/spatial/tissue_positions_list.csv', header=FALSE)
dim(pos.info)
head(pos.info)
pos <- pos.info[,c(5,6)]
rownames(pos) <- pos.info[,1]
head(pos)
plot(pos)

## fix order
pos <- pos[colnames(countmat),]

## fix position
pos <- pos[, c(2,1)]
pos[,2] <- -pos[,2] 
pos[,2] <- pos[,2] - min(pos[,2])
head(pos)

library(Matrix)
libsize <- colSums(countmat)
par(mfrow=c(1,1))
MERINGUE::plotEmbedding(pos, col=libsize)

g <- 'Slc17a7'
g <- 'Gabbr2'
g <- 'Gad2'
gexp <- countmat[g,]
MERINGUE::plotEmbedding(pos, col=gexp[rownames(pos)], main=g)

###### run STdeconvolve
## remove pixels with too few genes
counts <- cleanCounts(countmat, 
                      min.lib.size = 10^3.75,
                      min.reads = 2000)
dim(counts)
par(mfrow=c(1,1))
MERINGUE::plotEmbedding(pos[colnames(counts),], col=colSums(counts))
g %in% rownames(counts)
## feature select for genes
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
dim(corpus)
## choose optimal number of cell-types
ldas <- fitLDA(t(as.matrix(corpus)), Ks = c(25))
## get best model results
optLDA <- optimalModel(models = ldas, opt = "25")
## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta
head(deconProp)
dim(deconProp)
dim(pos)
dim(pos[rownames(deconProp),])
## visualize deconvolved cell-type proportions
ppos <- pos[rownames(deconProp),]
colnames(ppos) <- c('x', 'y')
vizAllTopics(deconProp, ppos, r=100, lwd = 0)	 

## save
save(deconProp, deconGexp, ppos, file="results1.RData")
