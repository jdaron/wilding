library(ggplot2)
library(reshape2)
library(viridis)
library(grid)
library(gtable)
library(RColorBrewer)
library(strucchange)
library(plyr)
library(scatterpie)
library(ggrepel)
library(geosphere)
library(measurements)
library(grid)
library("gridExtra")
library(ggpubr)

## 1. PCA
# info
infowi <- read.table("~/bioInf/wilding/ressources/wilding.samples.meta.txt", h=T)
infowi = infowi[,c("id", "population")]
infoag <- read.table("~/bioInf/wilding/ressources/ag1000g.samples.meta.txt", h=T, sep="\t")
colnames(infoag)[1] = "id"
infoag = infoag[,c("id", "population")]

info = rbind(infoag, infowi)

# colors
colag <- read.table("~/bioInf/wilding/ressources/ag1000g.samples.colors.txt", h=T, sep="\t")
colag$shape = "ag"
colwi <- read.table("~/bioInf/wilding/ressources/wilding.samples.colors.txt", h=T, sep="\t")
colwi$shape = "wi"

cols = rbind(colag, colwi)
cols$order = 1:dim(cols)[1]

group.colors = as.character(cols$colors[cols$order])
names(group.colors) = factor(cols$population, levels=cols$population[cols$order])
group.shape = c(rep(21, dim(colag)[1]), rep(24, dim(colwi)[1]))
names(group.shape) = factor(cols$population, levels=cols$population[cols$order])

prefix = "~/bioInf/wilding/popstructure/pca/wilding.biallelic.3.pca"
#prefix = "~/bioInf/wilding/popstructure/pca/wilding.biallelic.3.pca.eigenval.txt"
inF = paste(prefix, ".eigenval.txt", sep="")
eigenval <- data.frame(read.table(inF, header=FALSE, skip=0, sep=" "))
eigenval = eigenval$V1*100
eigenval = round(eigenval, digits = 2)
df = data.frame(eval = eigenval, coord = 1:length(eigenval))

p1 = ggplot(df, aes(x=coord, y=eval)) +
  geom_bar(stat = "identity", color = "cornsilk4", fill = "cornsilk4", width = 0.6) +
  theme_classic() +
  ylab("% of explained variance") +
  xlab("Principal Components")
p1

inF = paste(prefix, ".eigenvec.txt", sep="")
eigenvec <- data.frame(read.table(inF, header=FALSE, skip=0, sep="\t"))
eigenvec = eigenvec[,1:7]
colnames(eigenvec) = c("id", "EV1", "EV2", "EV3", "EV4", "EV5", "EV6")
eigenvec = cbind(samples, eigenvec)
head(eigenvec)
dim(eigenvec)

merge = c()
merge = merge(eigenvec, info, by="id")
head(merge)

mat = data.frame()
mat <- data.frame(
  id = merge$id,
  pop = merge$population,
#  location = merge$location,
#  site = merge$site,
  EV1 = merge$EV1,    # the first eigenvector
  EV2 = merge$EV2,    # the second eigenvector
  EV3 = merge$EV3,    # the second eigenvector
  EV4 = merge$EV4,    # the second eigenvector
  stringsAsFactors = FALSE)

mat$pop = factor(mat$pop, levels=cols$population[cols$order])

xlab = paste("EV1", round(eigenval[1], digits = 2), sep = " ")
ylab = paste("EV2", round(eigenval[2], digits = 2), sep = " ")

p2 = ggplot(mat, aes(x=EV1, y=EV2, color=pop, fill=pop, shape=pop)) +
  geom_point(size=2) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = group.shape) +
  xlab(xlab) +
  ylab(ylab) +
  theme_classic() +
  theme(legend.position="none")
p2

xlab = paste("EV1", round(eigenval[1], digits = 2), sep = " ")
ylab = paste("EV3", round(eigenval[3], digits = 2), sep = " ")

p3 = ggplot(mat, aes(x=EV1, y=EV3, color=pop, fill=pop, shape=pop)) +
  geom_point(size=2) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = group.shape) +
  xlab(xlab) +
  ylab(ylab) +
  theme_classic() +
  theme(legend.position="none")

xlab = paste("EV2", round(eigenval[2], digits = 2), sep = " ")
ylab = paste("EV3", round(eigenval[3], digits = 2), sep = " ")

p4 = ggplot(mat, aes(x=EV2, y=EV3, color=pop, fill=pop, shape=pop)) +
  geom_point(size=2) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = group.shape) +
  xlab(xlab) +
  ylab(ylab) +
  theme_classic() + 
  guides(color=guide_legend(ncol=1))
# theme(legend.position="none")

# to extract legend and print it
leg <- get_legend(p4)
pl = as_ggplot(leg)


pdf("~/bioInf/wilding/popstructure/pca/wilding.biallelic.3.pca.pdf", width = conv_unit(150, "mm", "inch"), height = conv_unit(120, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2, top = "PCA 857,565 SNPs pruned")
dev.off()

pdf("~/bioInf/wilding/popstructure/pca/ag1000g_phase2.wilding.merged.biallelic.3.pca.pdf", width = conv_unit(150, "mm", "inch"), height = conv_unit(120, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p1, p2, p3, pl, ncol = 2, nrow = 2, top = "PCA 1M SNPs pruned")
dev.off()
