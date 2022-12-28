# read in data
?read.csv
data <- read.csv('genomic_data_visualization_2023_course_prep/charmander.csv.gz', row.names=1)
head(data)

pos <- data[,1:2]
area <- data[,3]
gexp <- data[,4:ncol(data)]
head(pos)
head(area)
names(area) <- rownames(data)
head(gexp)

# prelim look
?plot
plot(pos, pch=".")

g <- "MS4A1"
head(gexp[,g])
boxplot(gexp[,g])
hist(log10(gexp[,g]+1))

cellwithcd20 <- rownames(gexp)[which(gexp[,g] > 0)]
cellwithoutcd20 <- rownames(gexp)[which(gexp[,g] == 0)]
range(gexp[cellwithcd20, g])
range(gexp[cellwithoutcd20, g])

plot(pos[cellwithcd20,], pch=".", col='red')
points(pos[cellwithoutcd20,], pch=".", col='grey')

g2 <- "CXCR4"
hist(log10(gexp[,g2]+1))

cellwithg <- rownames(gexp)[which(gexp[,g2] > 0)]
cellwithoutg <- rownames(gexp)[which(gexp[,g2] == 0)]
range(gexp[cellwithg, g2])
range(gexp[cellwithoutg, g2])
plot(pos[cellwithg,], pch=".", col='red')
points(pos[cellwithoutg,], pch=".", col='grey')

plot(gexp[,g], gexp[,g2], pch=16)

dim(gexp)

# QC
head(gexp)
hist(log10(colSums(gexp>0)+1))

hist(log10(rowSums(gexp)+1))
good.cells <- names(which(rowSums(gexp) > 10))
plot(pos, pch=".")
points(pos[good.cells,], pch=".", col='green')

pos <- pos[good.cells,]
gexp <- gexp[good.cells,]
dim(pos)
dim(gexp)

# principal components analysis
?prcomp
pcs <- prcomp(log10(gexp+1))
names(pcs)
head(pcs$x)
plot(pcs$x[,1:2], pch=".")
plot(pcs$x[,3:4], pch=".")
plot(pcs$x[,9:10], pch=".")

plot(pcs$x[,1:2][cellwithcd20,], pch=".", col='red')
points(pcs$x[,1:2][cellwithoutcd20,], pch=".", col='grey')

head(pcs$rotation[g,])

plot(pcs$sdev[1:30], type='l')

range(pcs$x[,1])
mycolors <- colorRampPalette(c('blue', 'red'))(100)
x <- pcs$x[,1] - min(pcs$x[,1])
y <- x/max(x)*99+1
pccol <- mycolors[floor(y)]

plot(pos, col=pccol, pch=".")
plot(pcs$x[,1:2], col=pccol, pch=".")

# tSNE
library(Rtsne)
?Rtsne::Rtsne
emb <- Rtsne::Rtsne(pcs$x[,1:10])

dups <- which(duplicated(pcs$x[,1:10]))
table(rowSums(gexp[names(dups), ]))

names(emb)
embb <- emb$Y
rownames(embb) <- rownames(gexp)

plot(embb, pch=".")

head(sort(colSums(gexp), decreasing=TRUE))
g <- "LUM"
cellwithcd20 <- rownames(gexp)[which(gexp[,g] > 0)]
cellwithoutcd20 <- rownames(gexp)[which(gexp[,g] == 0)]
plot(embb[cellwithcd20,], pch=".", col='red')
points(embb[cellwithoutcd20,], pch=".", col='grey')

# kmeans
?kmeans
com <- kmeans(pcs$x[,1:10], centers =10)
names(com)
clusters <- com$cluster

ccol <- rainbow(10)[clusters]
plot(embb, pch=".", col=ccol)
plot(pos, pch=".", col=ccol)

# differential expression analysis
cellsofinterest <- names(which(clusters == 5))
othercells <- names(which(clusters != 5))

plot(embb, pch=".")
points(embb[cellsofinterest,], pch=".", col='green')

g <- 'LUM'
boxplot(gexp[cellsofinterest, g])
boxplot(gexp[othercells,g])
?wilcox.test
wilcox.test(gexp[cellsofinterest, g], gexp[othercells,g], 'greater')
wilcox.test(gexp[cellsofinterest, g], gexp[othercells,g], 'less')





