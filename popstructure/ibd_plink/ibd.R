library(ggplot2)
library(grid)
library("gridExtra")
library(ggplot2)
library(reshape2)
library(viridis)
library(grid)
library(gtable)
library(RColorBrewer)
library(strucchange)
library(plyr)
library(ggrepel)
library(geosphere)
library("gridExtra")
library(measurements)
library(ggpubr)
library(tidyr)

#### IBD for AG1000G
info <- read.table("~/bioInf/wilding/ressources/ag1000g.samples.meta.txt", h=T, sep="\t")
head(info)
colnames(info)[1] = "FID1"

colors <- read.table("~/bioInf/wilding/ressources/ag1000g.samples.colors.txt", h=T, sep = "\t")
head(colors)

### AG1000G
mat = read.table("~/bioInf/wilding/popstructure/ibd/ag1000g.phase2.ar1.pass.biallelic.allPop.ibd.genome.txt", header=TRUE, skip=0, sep="\t")
head(mat)

mat = merge(mat, info, by="FID1")
mat$population = factor(mat$population, levels = colors$population)

pop.colors = as.character(colors$colors)
names(pop.colors) = colors$population

p1 = ggplot(mat, aes(x=population, y=PI_HAT, color=population, fill=population)) +
  geom_jitter(width = 0.2) +
#  geom_hline(yintercept = 0.7, linetype = "dotted") +
  scale_color_manual(values=pop.colors) +
  scale_fill_manual(values=pop.colors) +
  ylim(c(-0.01,1)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1))
p1

pdf("~/bioInf/wilding/popstructure/ibd/ag1000g.ibd.pdf", width = conv_unit(150, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
#grid.arrange(p1, ncol = 1, nrow = 1, top = "IBD intra pop AG1000G")
grid.arrange(p1, ncol = 1, nrow = 1, top = textGrob("IBD intra pop AG1000G", gp=gpar(fontsize=8,font=8)))
dev.off()

#### IBD for WILDING 
mat = read.table("~/bioInf/wilding/popstructure/ibd/wilding.chr3.ibd.genome.txt", header=TRUE, skip=0, sep="\t")
head(mat)

info <- read.table("~/bioInf/wilding/ressources/wilding.samples.meta.txt", h=T, sep="\t")
head(info)
info = info[,c(1,7)]
colnames(info)[1] = "FID1"
colnames(info)[2] = "pop_FID1"
mat = merge(mat, info, by="FID1")
head(mat)
colnames(info)[1] = "FID2"
colnames(info)[2] = "pop_FID2"
mat = merge(mat, info, by="FID2")
head(mat)

mat$compLabel = NA
mat$compLabel[which(mat$pop_FID1==mat$pop_FID2)] = as.character(mat$pop_FID1[which(mat$pop_FID1==mat$pop_FID2)])
mat$compLabel[which( (mat$pop_FID1=="LBVdom" & mat$pop_FID2=="LPdom") | (mat$pop_FID1=="LPdom" & mat$pop_FID2=="LBVdom") )] = "LBVdom-LPdom"
mat$compLabel[which( (mat$pop_FID1=="LBVdom" & mat$pop_FID2=="LPfor") | (mat$pop_FID1=="LPfor" & mat$pop_FID2=="LBVfor") )] = "LBVdom-LPfor"
mat$compLabel[which( (mat$pop_FID1=="LPdom" & mat$pop_FID2=="LPfor") | (mat$pop_FID1=="LPfor" & mat$pop_FID2=="LPdom") )] = "LPdom-LPfor"
mat$compLabel = factor(mat$compLabel, levels = unique(mat$compLabel)[c(1,3,6,2,4,5)])

p2 = ggplot(mat, aes(x=compLabel, y=PI_HAT, color=compLabel, fill=compLabel)) +
  geom_jitter(width = 0.2) +
#  geom_hline(yintercept = 0.7, linetype = "dotted") +
#  scale_color_manual(values=pop.colors) +
#  scale_fill_manual(values=pop.colors) +
  theme_classic() +
  ylim(c(-0.01,1)) +
  xlab("") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1))
p2

pdf("~/bioInf/wilding/popstructure/ibd/wilding.ibd.pdf", width = conv_unit(150, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p2, ncol = 1, nrow = 1, top = "IBD intra and inter pop Wilding")
dev.off()

mat[which(mat$PI_HAT>0.4),]

#### IBD for WILDING and URBANO
mat = read.table("~/bioInf/wilding/popstructure/ibd/wilding_urbano.chr3.ibd.genome.txt", header=TRUE, skip=0, sep="\t")
head(mat)

info <- read.table("~/bioInf/wilding/ressources/wilding_urbano.samples.meta.txt", h=T, sep="\t")
info = info[,c(1,2)]
head(info)
colnames(info)[1] = "FID1"
colnames(info)[2] = "pop_FID1"
mat = merge(mat, info, by="FID1")
head(mat)
colnames(info)[1] = "FID2"
colnames(info)[2] = "pop_FID2"
mat = merge(mat, info, by="FID2")
head(mat)

mat = mat[which(mat$pop_FID1==mat$pop_FID2),]
mat$txt = NA
mat$txt[which(mat$PI_HAT>0.2)] = paste(mat$FID1[which(mat$PI_HAT>0.2)], mat$FID2[which(mat$PI_HAT>0.2)], sep = "_") 

p3 = ggplot(mat, aes(x=pop_FID1, y=PI_HAT, color=pop_FID1, fill=pop_FID1, label=txt)) +
  geom_jitter(width = 0.2) +
  geom_text_repel(color="black", arrow = arrow(length = unit(0.02, "npc"))) +
  #  geom_hline(yintercept = 0.7, linetype = "dotted") +
  #  scale_color_manual(values=pop.colors) +
  #  scale_fill_manual(values=pop.colors) +
  theme_classic() +
  ylim(c(-0.01,1)) +
  xlab("") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1))

pdf("~/bioInf/wilding/popstructure/ibd/wilding_urbano.ibd.pdf", width = conv_unit(150, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p3, ncol = 1, nrow = 1, top = "IBD Wilding & Urbano")
dev.off()

