load('results1.RData')
load('results2.RData')

head(deconGexp)
head(deconGexp2)

dim(deconGexp)
dim(deconGexp2)

genes.shared <- intersect(colnames(deconGexp), colnames(deconGexp2))
length(genes.shared)

cmat <- cor(t(deconGexp[, genes.shared]), t(deconGexp2[, genes.shared]))
dim(cmat)
heatmap(cmat)

spot.shared <- intersect(rownames(ppos), rownames(ppos2))

vizTopic(theta = deconProp[spot.shared,], pos = ppos[spot.shared,], 
         topic = "21", plotTitle = "STdeconvolve vanilla: topic 21",
         size = 1.5, stroke = 0.1, alpha = 0.5,
         low = "white",
         high = "red")

vizTopic(theta = deconProp2[spot.shared,], pos = ppos2[spot.shared,], 
         topic = "1", plotTitle = "STdeconvolve Spotclean: topic 1",
         size = 1.5, stroke = 0.1, alpha = 0.5,
         low = "white",
         high = "red")

x <- deconProp[spot.shared,21]
y <- deconProp2[spot.shared,1]
length(x)
length(y)
plot(x,y, pch=".")
cor.test(as.numeric(x),as.numeric(y))


