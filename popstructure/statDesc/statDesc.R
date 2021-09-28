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
library(scales)

### Info
infowi <- read.table("~/bioInf/wilding/ressources/wilding.samples.meta.txt", h=T)
infowi = infowi[,c("id", "population")]
infowi$grp = "wil"
infour <- read.table("~/bioInf/wilding/ressources/wilding_urbano.samples.meta.txt", h=T, sep = "\t")
infour = infour[,c("id", "population")]
infour$grp = "urb"
infoag <- read.table("~/bioInf/wilding/ressources/ag1000g.samples.meta.txt", h=T, sep="\t")
colnames(infoag)[1] = "id"
infoag = infoag[,c("id", "population")]
infoag$grp = "ag"

info = rbind(infoag, infour, infowi)

### colors
# colors
colag <- read.table("~/bioInf/wilding/ressources/ag1000g.samples.colors.txt", h=T, sep="\t")
colag$shape = "ag"
colwi <- read.table("~/bioInf/wilding/ressources/wilding.samples.colors.txt", h=T, sep="\t")
colwi$shape = "wi"
colur <- read.table("~/bioInf/wilding/ressources/urbano.samples.colors.txt", h=T, sep="\t")
colur$shape = "ur"

cols = rbind(colag, colur, colwi)
cols$order = 1:dim(cols)[1]

group.colors = as.character(cols$colors[cols$order])
names(group.colors) = factor(cols$population, levels=cols$population[cols$order])
group.shape = c(rep(21, dim(colag)[1]), rep(22, dim(colur)[1]), rep(24, dim(colwi)[1]))
names(group.shape) = factor(cols$population, levels=cols$population[cols$order])

group.alpha = c(rep(0.5, dim(colag)[1]), rep(1, dim(colur)[1]), rep(1, dim(colwi)[1]))
names(group.alpha) = factor(cols$population, levels=cols$population[cols$order])

group.size = c(rep(0.8, dim(colag)[1]), rep(1.2, dim(colur)[1]), rep(1.2, dim(colwi)[1]))
names(group.size) = factor(cols$population, levels=cols$population[cols$order])

### stat by windows
mat <- read.table("~/bioInf/wilding/popstructure/statDesc/pi_tajima/wilding_urbano_ag1000g.statWindows.tab")
colnames(mat) = c("population", "chrom", "start", "stop", "nbase", "counts", "pi", "tajimaD", "inb_coef")
head(mat)
length(levels(mat$population))
length(levels(info$population))
mat$population = factor(mat$population, levels=names(group.colors))
head(mat)
head(info)

length(levels(mat$population))
length(group.alpha)

# Pi
p1 = ggplot(mat, aes(x=population, y=pi, color="b", fill=population, alpha=population)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim=c(0, 0.03)) +
  scale_color_manual(values="black") +
  scale_fill_manual(values=group.colors) +
  scale_alpha_manual(values = group.alpha) +
#  facet_wrap(~grp)+ 
  xlab("pop") +
  ylab("pi") +
  theme_classic() +
  theme(axis.text.x = element_text(size=8, angle = 45, vjust = 1, hjust=1),
        legend.position="none")
p1

# tajimaD
p2 = ggplot(mat, aes(x=population, y=tajimaD, color="b" ,fill=population, alpha=population)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept=0, linetype="dashed") +
#  coord_cartesian(ylim=c(0, 0.015)) +
  scale_color_manual(values="black") +
  scale_fill_manual(values=group.colors) +
  scale_alpha_manual(values = group.alpha) +
  xlab("pop") +
  ylab("tajima D") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=8, angle = 45, vjust = 1, hjust=1),
    legend.position="none"
    )
p2

pdf("~/bioInf/wilding/popstructure/statDesc/pi_tajima/wilding_urbano_ag1000g.pi.tajima.pdf", width = conv_unit(170, "mm", "inch"), height = conv_unit(90, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p1, p2, ncol = 2, nrow = 1, top = "Pi / tajima D")
dev.off()

### ldDecay
mat <- read.table("~/bioInf/wilding/popstructure/statDesc/ldDecay/wilding_urbano_ag1000g.lddecay.txt", header = T)
sub = mat[which(mat$population=="LPdom" | mat$population=="LPfor" | mat$population=="LBVwil" | mat$population=="LBVurb" | mat$population=="DLAurb" | mat$population=="BZVurb"),]
head(sub)
head(mat)

mat$population = factor(mat$population, levels=names(group.colors))

p3 = ggplot() +
  geom_line(data = mat, aes(x=distance, y=ld, color=population ,fill=population, alpha=population, size=population)) +
  geom_line(data = sub, aes(x=distance, y=ld, color=population ,fill=population, alpha=population, size=population)) +
#  facet_wrap(~phase) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_alpha_manual(values = group.alpha) +
  scale_size_manual(values = group.size) +
  xlab("Physical Distance (log10)") +
  ylab("linkage disequilibrium") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) +
  guides(color=guide_legend(ncol=2))
p3

pdf("~/bioInf/wilding/popstructure/statDesc/ldDecay/wilding_urbano_ag1000g.lddecay.pdf", width = conv_unit(170, "mm", "inch"), height = conv_unit(90, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p3, ncol = 1, nrow = 1, top = "LD decay")
dev.off()

### folded Site Frequency Spetrum
mat <- read.table("~/bioInf/wilding/popstructure/statDesc/sfs/wilding_urbano_ag1000g.sfs.txt", header = F)
colnames(mat) = c("population", "x", "scaled_y")
#mat$phase = "phase2"
#mat$phase[which(mat$population=="LPdom" | mat$population=="LPfor" | mat$population=="LBVdom")] = "wilding"
mat = mat[-which(mat$x==0),]

mat$population = factor(mat$population, levels=names(group.colors))

p4 = ggplot(mat, aes(x=x, y=scaled_y, color=population ,fill=population, alpha=population, size=population)) +
  geom_line() +
  #  facet_wrap(~phase) +
#  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_alpha_manual(values = group.alpha) +
  scale_size_manual(values = group.size) +
  xlab("Minor allele frequency") +
  ylab("SNP density") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) +
  guides(color=guide_legend(ncol=2))
p4

pdf("~/bioInf/wilding/popstructure/statDesc/sfs//wilding_urbano_ag1000g.folded_sfs_scaled.pdf", width = conv_unit(170, "mm", "inch"), height = conv_unit(90, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p4, ncol = 1, nrow = 1, top = "Folded Site Frequency Spetrum")
dev.off()

### IBD count/sum
mat <- read.table("~/bioInf/wilding/popstructure/statDesc/ibdseq/wilding_urbano.ag1000g.ibdseq.txt")
colnames(mat) = c("id", "haploidx", "id2", "haploidx2", "chrom", "start", "stop", "lod")
head(mat)
mat$id_id2 = paste(mat$id, mat$id2, sep="_")
mat$len = mat$stop-mat$start+1
head(mat)

mat = merge(mat, info, by="id")

# violin plot
p5 = ggplot(mat, aes(x=population, y=len, color=population, fill=population, alpha=population)) +
  geom_violin(color="black") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_alpha_manual(values=group.alpha) +
  xlab("Population") +
  ylab("IBD track length (bp)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
    legend.position="none") 
p5

pdf("~/bioInf/wilding/popstructure/statDesc/ibdseq//wilding_urbano_ag1000g.IBD_violin.pdf", width = conv_unit(150, "mm", "inch"), height = conv_unit(90, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p5, ncol = 1, nrow = 1, top = "")
dev.off()

# count vs sum IBD
mat.s = ddply(mat, .(population, grp, id_id2), plyr::summarize, count=length(id_id2), sum=sum(stop-start+1))
head(mat.s)

group.size = c(rep(0.3, dim(colag)[1]), rep(0.7, dim(colur)[1]), rep(0.7, dim(colwi)[1]))
names(group.size) = factor(cols$population, levels=cols$population[cols$order])

p6 = ggplot(mat.s, aes(x=sum, y=count, color=population, fill=population, shape=population, alpha=population, size=population)) +
  geom_point(stroke = 0) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = group.shape) +
  scale_alpha_manual(values = group.alpha) +
  scale_size_manual(values = group.size) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
#  facet_wrap(~population)+ 
  xlab("sum IBD") +
  ylab("count IBD") +
  theme_classic() +
  guides(color=guide_legend(ncol=2))

p6

pdf("~/bioInf/wilding/popstructure/statDesc/ibdseq//wilding_urbano_ag1000g.IBD_count_sum.pdf", width = conv_unit(150, "mm", "inch"), height = conv_unit(90, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p6, ncol = 1, nrow = 1, top = "")
dev.off()

### Comparison method ibdseq vs hapibd ; ibdne
mat = read.table("~/bioInf/wilding/popstructure/statDesc/ibd/comparison_ibdtrack_hapibd_ibdseq.txt", h=F, fill = TRUE, dec = ".")
colnames(mat) = c("method", "pop", "id", "haploidx", "id2", "haploidx2", "chrom", "start", "stop", "lod")
head(mat)
mat$id_id2 = paste(mat$id, mat$id2, sep="_")
mat$len = mat$stop-mat$start+1
head(mat)

mat.s = ddply(mat, .(method, id_id2), plyr::summarize, count=length(id_id2), sum=sum(stop-start+1))
head(mat.s)

p1 = ggplot(mat.s, aes(x=sum, y=count, color=method, fill=method)) +
  geom_jitter(size=2) +
  facet_wrap(~method) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  #  facet_wrap(~population)+ 
  xlab("sum IBD") +
  ylab("count IBD") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

pdf("~/bioInf/wilding/popstructure/statDesc/ibd/comparison_ibdtrack_hapibd_ibdseq.pdf", width = conv_unit(210, "mm", "inch"), height = conv_unit(120, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p1, ncol = 1, nrow = 1, top = "Comparison IBD detection methods (KE)")
dev.off()

### ROH: Run of Homozygoties
roh <- read.table("~/bioInf/wilding/popstructure/statDesc/roh/wilding_ag1000g.roh.tab")
colnames(roh) = c("population", "start", "stop", "length", "is_marginal", "id", "chrom")
head(roh)
roh = roh[-which(roh$length<100000),]
roh.s = ddply(roh, .(id), plyr::summarize, count=length(id))
head(roh.s)

froh <- read.table("~/bioInf/wilding/popstructure/statDesc/roh/wilding_ag1000g.froh.tab")
colnames(froh) = c("population", "id", "chrom", "froh")
head(froh)
froh.s = ddply(froh, .(id, population), plyr::summarize, froh=sum(froh)/2)

mat = merge(roh.s, froh.s, by="id")
head(mat)

mat$population = factor(mat$population, levels=names(group.colors))

group.size = c(rep(1, dim(colag)[1]), rep(2, dim(colur)[1]), rep(2, dim(colwi)[1]))
names(group.size) = factor(cols$population, levels=cols$population[cols$order])

p7 = ggplot(mat, aes(x=froh, y=count, color=population, fill=population, shape=population, alpha=population, size=population)) +
  geom_jitter(stroke = 0) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = group.shape) +
  scale_alpha_manual(values = group.alpha) +
  scale_size_manual(values = group.size) +
  #  facet_wrap(~population)+ 
  xlab("froh") +
  ylab("count roh") +
  theme_classic() +
  guides(color=guide_legend(ncol=2))
p7

pdf("~/bioInf/wilding/popstructure/statDesc/roh/wilding_ag1000g.countroh_froh.pdf", width = conv_unit(170, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p7, ncol = 1, nrow = 1, top = "")
dev.off()

