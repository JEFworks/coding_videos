data <- read.csv('genomic_data_visualization_2023_course_prep/charmander.csv.gz', row.names=1)
head(data)
pos <- data[, 1:2]
area <- data[,3]
gexp <- data[,4:ncol(data)]
gexp[1:5,1:5]
dim(gexp)

## plot cell density
plot(pos, pch=".")
## distribution of gene expression
hist(log10(colSums(gexp)+1))
## distribution of cell expression
hist(log10(rowSums(gexp)+1), xlim=c(0,3))
## color cells by gene expression
head(sort(colSums(gexp), decreasing=TRUE))
g <- "ERBB2"
g <- "LUM"
hist(log10(gexp[, g]+1))
gcol <- colorRampPalette(c('lightgrey', 'red'))(100)[round((gexp[, g]/max(gexp[, g]))*99+1)]
plot(pos, pch=".", col=gcol, main=g)

## clean up 
good.cells <- rownames(gexp)[rowSums(gexp) >= 3]
length(good.cells)
gexp <- gexp[good.cells,]
pos <- pos[good.cells,]

## pca
pcs <- prcomp(log10(gexp+1))
plot(pcs$sdev[1:30], type="l")
head(pcs$x[1:5,])
g <- "ERBB2"
gcol <- colorRampPalette(c('lightgrey', 'red'))(100)[round((gexp[, g]/max(gexp[, g]))*99+1)]
plot(pcs$x[,1:2], pch=".", col=gcol, main=g)

libsize <- rowSums(gexp)
gcol <- colorRampPalette(c('lightgrey', 'red'))(100)[round((libsize/max(libsize))*99+1)]
plot(pcs$x[,1:2], pch=".", col=gcol, main='libsize')
## note to self: explore normalization

## tsne
library(Rtsne)
emb <- Rtsne::Rtsne(pcs$x[,1:10])
#pcs$x[,1:10][duplicated(pcs$x[,1:10]),]
#table(rowSums(gexp[names(which(duplicated(pcs$x[,1:10]))),]))
emb <- emb$Y
rownames(emb) <- rownames(gexp)
head(emb)
plot(emb, pch=".")
g <- "ERBB2"
g <- 'LUM'
gcol <- colorRampPalette(c('lightgrey', 'red'))(100)[round((gexp[, g]/max(gexp[, g]))*99+1)]
plot(emb, pch=".", col=gcol, main=g)

## clustering
k <- 5
com <- kmeans(pcs$x[,1:10], centers = k)
head(com$cluster)
ccol <- rainbow(k)[com$cluster]
par(mfrow=c(1,2))
plot(emb, pch=".", col=ccol, main='kmeans: 10')
plot(pos, pch=".", col=ccol, main='kmeans: 10')

## differential expression
cells1 <- names(which(com$cluster == 5))
plot(pos[cells1,], pch=".")
cells2 <- names(which(com$cluster != 5))
plot(pos[cells2,], pch=".")

#g <- 'ERBB2'
g <- 'LUM'
wilcox.test(gexp[cells1, g], gexp[cells2,g], alternative = 'greater')
boxplot(log10(gexp[cells1, g]+1))
boxplot(log10(gexp[cells2, g]+1))



