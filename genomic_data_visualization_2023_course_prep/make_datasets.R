## read in data

file <- '~/Downloads/Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.csv.gz'
cell.info <- read.csv(file)
head(cell.info)
dim(cell.info)

file <- '~/Downloads/Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5'
library(rhdf5)
h5ls(file)

data <- h5read(file,'matrix')
names(data)

length(data$barcodes)
length(data$data)
length(data$features$name)
head(data$features$name)
length(data$indices)
data$shape

library(Matrix)
?Matrix::sparseMatrix
cd <- sparseMatrix(i = data$indices + 1,
                   p = data$indptr,
                   x = data$data,
                   dims = data$shape)
rownames(cd) <- data$features$name
colnames(cd) <- data$barcodes
cd[1:5,1:5]

## sanity checks
pos <- cell.info[, 2:3]
rownames(pos) <- cell.info[,1]
head(pos)
plot(pos, pch=".")

table(colnames(cd) %in% rownames(pos))

rownames(cd)
head(sort(rowSums(cd), decreasing=TRUE))
g <- "ERBB2"
ggexp <- cd[g, ]
range(ggexp)
gcol <- colorRampPalette(c('lightgrey', 'red'))(100)[floor(ggexp/max(ggexp)*99) + 1]
names(gcol) <- names(ggexp)
head(gcol)
plot(pos, col=gcol, main=g, pch=".")

######## make smaller datasets for students

## downsample
vi <- sample(rownames(pos), nrow(pos)/4)
plot(pos[vi,], col=gcol[vi], main=g, pch=16, cex=0.5)
pos <- pos[vi,]

plot(pos, pch=".")

## make 4 smaller datasets
vi1 <- pos[,1] > 0 & pos[,1] < 3000 & pos[,2] > 0 & pos[,2] < 2000
table(vi1)
points(pos[vi1,], col='red', pch=".")

vi2 <- pos[,1] > 0 & pos[,1] < 3000 & pos[,2] > 2000 & pos[,2] < 4000
table(vi2)
points(pos[vi2,], col='green', pch=".")

vi3 <- pos[,1] > 4000 & pos[,1] < 8000 & pos[,2] > 0 & pos[,2] < 3000
table(vi3)
points(pos[vi3,], col='blue', pch=".")

vi4 <- pos[,1] > 4000 & pos[,1] < 8000 & pos[,2] > 3000 & pos[,2] < 6000
table(vi4)
points(pos[vi4,], col='magenta', pch=".")


########## write out some csv files
rownames(cd)
good.genes <- rownames(cd)[!grepl('NegControl|antisense|BLANK', rownames(cd))]
good.genes
  
area <- cell.info$nucleus_area
names(area) <- cell.info$cell_id

dim(pos[vi1,])
dim(cd[, rownames(pos[vi1,])])
dataset1 <- cbind(pos[vi1,], area = area[rownames(pos[vi1,])], t(as.matrix(cd[, rownames(pos[vi1,])]))[, good.genes])
write.csv(dataset1, file='pikachu.csv')

dim(pos[vi2,])
dim(cd[, rownames(pos[vi2,])])
dataset2 <- cbind(pos[vi2,], area = area[rownames(pos[vi2,])], t(as.matrix(cd[, rownames(pos[vi2,])]))[, good.genes])
write.csv(dataset2, file='bulbasaur.csv')

dim(pos[vi3,])
dim(cd[, rownames(pos[vi3,])])
dataset3 <- cbind(pos[vi3,], area = area[rownames(pos[vi3,])], t(as.matrix(cd[, rownames(pos[vi3,])]))[, good.genes])
write.csv(dataset3, file='charmander.csv')

dim(pos[vi4,])
dim(cd[, rownames(pos[vi4,])])
dataset4 <- cbind(pos[vi4,], area = area[rownames(pos[vi4,])], t(as.matrix(cd[, rownames(pos[vi4,])]))[, good.genes])
write.csv(dataset4, file='squirtle.csv')







