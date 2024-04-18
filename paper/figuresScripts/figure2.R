library(reshape)
library(ggplot2)
library(ggthemes)
library(tidyverse)
library("factoextra")
library("FactoMineR")
library(ggrepel)
library(plyr)
library(reshape2)
library(scatterpie)
library(sp)
library(RColorBrewer)
library(sf)
library(grid)
library(gtable)
library("gridExtra")
library(measurements)
library(scales)

setwd("~/bioInf/wilding/github/wilding/paper/input/figure2/")

##########################
### Figure 2:
##########################

### Info
infowi <- read.table("wilding.samples.meta.txt", h=T)
infowi = infowi[,c("id", "population")]
infowi$grp = "wil"
infoag <- read.table("ag1000g.samples.meta.txt", h=T, sep="\t")
colnames(infoag)[1] = "id"
infoag = infoag[,c("id", "population")]
infoag$grp = "ag"

info = rbind(infoag, infowi)

### colors
# colors
colag <- read.table("ag1000g.samples.colors.txt", h=T, sep="\t")
colag$shape = "ag"
colwi <- read.table("wilding.samples.colors.txt", h=T, sep="\t")
colwi$shape = "wi"

cols = rbind(colag, colwi)
cols$order = 1:dim(cols)[1]

group.colors = as.character(cols$colors[cols$order])
names(group.colors) = factor(cols$population, levels=cols$population[cols$order])

group.shape = c(rep(21, dim(colag)[1]), rep(24, dim(colwi)[1]))
names(group.shape) = factor(cols$population, levels=cols$population[cols$order])

group.alpha = c(rep(0.75, dim(colag)[1]), rep(1, dim(colwi)[1]))
names(group.alpha) = factor(cols$population, levels=cols$population[cols$order])

group.size = c(rep(0.8, dim(colag)[1]), rep(1.2, dim(colwi)[1]))
names(group.size) = factor(cols$population, levels=cols$population[cols$order])

order=c("AOcol", "LBVwil", "LPdom", "LPfor", "CIcol", "GHcol", "GNcol", "BFcol")

### Figure 2A: nucleotide diversity (pi)
###===============
mat <- read.table("wilding_ag1000g.statWindows.tab")
colnames(mat) = c("population", "chrom", "start", "stop", "nbase", "counts", "pi", "tajimaD", "inb_coef")
head(mat)

mat = mat[which(mat$population %in% order),]
mat$population = factor(mat$population, levels=order)

p1 = ggplot(mat, aes(y=population, x=pi, color=population, alpha=population)) +
  geom_boxplot(outlier.shape = NA, width=0.5) +
  coord_cartesian(xlim=c(0, 0.03)) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_alpha_manual(values = group.alpha) +
  #  facet_wrap(~grp)+ 
  ylab("") +
  xlab("Nucleotide Diversity (pi)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=8, angle = 0, vjust = 1, hjust=1),
        legend.position="none")
p1

### Figure 2B: tajimaD
###===============
mat <- read.table("wilding_ag1000g.statWindows.tab")
colnames(mat) = c("population", "chrom", "start", "stop", "nbase", "counts", "pi", "tajimaD", "inb_coef")
head(mat)

mat = mat[which(mat$population %in% order),]
mat$population = factor(mat$population, levels=order)

p2 = ggplot(mat, aes(y=population, x=tajimaD ,color=population, alpha=population)) +
  geom_boxplot(outlier.shape = NA, width=0.5) +
  geom_vline(xintercept=0, linetype="dashed") +
  #  coord_cartesian(ylim=c(0, 0.015)) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_alpha_manual(values = group.alpha) +
  ylab("") +
  xlab("Tajima's D") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=8, angle = 0, vjust = 1, hjust=1),
    legend.position="none"
  )
p2

### Figure 2C: LDdecay
###===============
mat <- read.table("wilding_ag1000g.lddecay.txt", header = T)
sub = mat[which(mat$population=="LPdom" | mat$population=="LPfor" | mat$population=="LBVwil"),]
mat = mat[which(mat$population %in% order),]
mat = mat[-which(mat$population=="GNcol"),]

mat$population = factor(mat$population, levels=order)
sub$population = factor(sub$population, levels=order)

p3 = ggplot() +
  geom_line(data = mat, aes(x=distance, y=ld, color=population ,fill=population, alpha=population, size=population)) +
#  geom_line(data = sub, aes(x=distance, y=ld, color=population ,fill=population, alpha=population, size=population)) +
  #  facet_wrap(~phase) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=group.colors[order]) +
  scale_fill_manual(values=group.colors[order]) +
  scale_alpha_manual(values = group.alpha[order]) +
  scale_size_manual(values = group.size[order]) +
  xlab("Physical Distance (log10)") +
  ylab("Linkage Disequilibrium") +
  theme_classic() +
  theme(
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) +
  guides(color=guide_legend(ncol=1))

p3

### Figure 2D: Folded site frequency spectrum
###===============
mat <- read.table("wilding_ag1000g.folded_sfs_scaled.tab", header = F)
colnames(mat) = c("population", "x", "scaled_y")
mat.d = ddply(mat, .(population), plyr::summarize, sum=sum(scaled_y))

mat = mat[-which(mat$x==0),]
mat = merge(mat, mat.d, by="population")

mat$sum[which(mat$population%in%c("LBVwil", "LPdom", "LPfor"))] = mat.d$sum[which(mat.d$population=="AOcol")]

mat = mat[which(mat$population %in% order),]
mat$population = factor(mat$population, levels=order)
mat = mat[-which(mat$population=="GNcol"),]

p4 = ggplot(mat, aes(x=x, y=scaled_y/sum, color=population ,fill=population, alpha=population, size=population)) +
  geom_line() +
  #  facet_wrap(~phase) +
  #  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=group.colors[order]) +
  scale_fill_manual(values=group.colors[order]) +
  scale_alpha_manual(values = group.alpha[order]) +
  scale_size_manual(values = group.size[order]) +
  xlab("Minor allele frequency") +
  ylab("SNP density") +
  theme_classic() +
  theme(
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) +
  guides(color=guide_legend(ncol=1))

p4

### Figure 2E: ROH Run of homozygoties
###===============
froh <- read.table("~/Dropbox/ag1000g_wilding.froh.txt")
colnames(froh) = c("population", "id", "chrom", "froh")
froh.s = ddply(froh, .(id), plyr::summarize, froh=mean(froh))
head(froh)
head(froh.s)

roh <- read.table("~/Dropbox/ag1000g_wilding.roh.txt")
colnames(roh) = c("population", "start", "stop", "length", "is_marginal", "id", "chrom")
roh = roh[which(roh$length>100000),]
roh = roh[which(roh$is_marginal=="False"),]
roh.s = ddply(roh, .(id), plyr::summarize, count=length(id), sum=sum(stop-start+1))
head(roh)
head(roh.s)
dim(roh.s)

mat = merge(froh.s, roh.s, by=c("id"))
mat = merge(mat, info,by="id")
head(mat)

group.size = c(rep(1, dim(colag)[1]), rep(2, dim(colwi)[1]))
names(group.size) = factor(cols$population, levels=cols$population[cols$order])

p5 = ggplot(mat, aes(x=froh, y=count, color=population, fill=population, shape=population, size=population)) +
  geom_jitter(stroke = 0) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = group.shape) +
  scale_size_manual(values = group.size) +
  xlab("froh") +
  ylab("count roh") +
  theme_classic() +
  theme_classic() +
  theme(
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p5

### Figure 2F: 
###===============

# count vs sum IBD
mat <- read.table("wilding_ag1000g.ibdseq.txt")
colnames(mat) = c("id", "haploidx", "id2", "haploidx2", "chrom", "start", "stop", "lod")
head(mat)
mat$id_id2 = paste(mat$id, mat$id2, sep="_")
mat$len = mat$stop-mat$start+1
head(mat)

mat = merge(mat, info, by="id")
mat = mat[which(mat$population %in% order),]

mat.s = ddply(mat, .(population, grp, id_id2), plyr::summarize, count=length(id_id2), sum=sum(stop-start+1))
head(mat.s)

group.size = c(rep(0.3, dim(colag)[1]), rep(0.7, dim(colwi)[1]))
names(group.size) = factor(cols$population, levels=cols$population[cols$order])

p6 = ggplot(mat.s, aes(x=sum, y=count, color=population, fill=population, shape=population, size=population)) +
  geom_point(stroke = 0) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = group.shape) +
  scale_alpha_manual(values = group.alpha) +
  scale_size_manual(values = group.size) +
  xlab("sum IBD") +
  ylab("count IBD") +
  theme_classic() +
  theme(
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

p6

##### plot all figures
pdf("wilding.ag1000g.statDesc.pdf", width = conv_unit(200, "mm", "inch"), height = conv_unit(120, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p1, p2, p5, p3, p4, p6, ncol = 3, nrow = 2)
dev.off()
