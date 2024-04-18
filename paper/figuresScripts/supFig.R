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
library(grid)
library(gtable)
library("gridExtra")
library(measurements)
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
library(grid)
library(measurements)
library(ggpubr)
library(dplyr)
library(eulerr)

setwd("~/bioInf/wilding/github/wilding/paper/input/supFig/")

##########################
### Supplementary Figures:
##########################

####################################
### SupFig1 Mapping and SNP density:
####################################

# load heterochromatine
heterochrom = read.table("Anopheles_gambiae.AgamP4.dna.heterochromatin.bed", h=F, fill = TRUE, dec = ".")
colnames(heterochrom) = c("chrom", "start", "stop")
heterochrom$chrom = factor(heterochrom$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
head(heterochrom)

# chrom size
chromSize <- read.table("Anopheles_gambiae.AgamP4.dna.chr.genome", h=F)
colnames(chromSize) = c("chr", "len")
k = c("2R", "2L", "3R", "3L", "X")
chromSize = chromSize[which(chromSize$chr %in% k),]
chromSize$proportion = chromSize$len/mean(chromSize$len)

### 1. Genome Mean Coverage
colwi <- read.table("wilding.samples.colors.txt", h=T, sep="\t")

group.colors = as.character(colwi$colors[colwi$order])
names(group.colors) = factor(colwi$population, levels=colwi$population[colwi$order])

info <- read.table("~/bioInf/wilding/github/wilding/ressources/wilding.samples.meta.txt", h=T)
head(info)
info$id = factor(info$id, levels = info$id[order(info$seq_depth)])

p1 = ggplot(info, aes(x=id, y=seq_depth, color=population, fill=population)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 14) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_classic() +
  ylab("Mean Coverage") +
  xlab("Individuals") +
  theme(
    legend.position="none",
    axis.text.x=element_blank(),
#    axis.ticks.x=element_blank(),
#    axis.text.x = element_text(size = 6, angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank())

p1

mean(info$seq_depth[which(info$seq_depth>14)])
sd(info$seq_depth[which(info$seq_depth>14)])
min(info$seq_depth[which(info$seq_depth>14)])
max(info$seq_depth[which(info$seq_depth>14)])

length(info$seq_depth[which(info$seq_depth<14)])

### 2. Sexe typing:
mat = read.table("wilding.coveragePerChr.tab", header=F, skip=0, sep="\t")
colnames(mat) = c("id", "chr", "chr_size", "bp_map", "chr_mean_cov", "md_cov")

mat = merge(mat, info, by="id")
mat = mat[-which(mat$chr=="Y_unplaced" | mat$chr=="Mt" | mat$chr=="UNKN"),]
head(mat)

p2 = ggplot(mat, aes(x=chr, y=log10(chr_mean_cov/seq_depth), color=population, fill=population)) +
  geom_jitter(width = 0.2) +
  theme_classic() +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  ylab("Ratio chr depth / total depth (log10)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position="none")

p2

# subselect male
sub_chrX = mat[which(mat$chr=="X"),]
t = sub_chrX[which(log10(sub_chrX$chr_mean_cov/sub_chrX$seq_depth) < -0.1),]
t$id

# write
tmp = mat[which(mat$chr=="X"),]
head(mat)
write.table(tmp, "~/bioInf/wilding/tmp", quote = F, sep = "\t")

### 3. SNP density
chromSize <- read.table("Anopheles_gambiae.AgamP4.dna.chr.genome", h=F)
colnames(chromSize) = c("chr", "len")
k = c("2R", "2L", "3R", "3L", "X")
chromSize = chromSize[which(chromSize$chr %in% k),]
chromSize$proportion = chromSize$len/mean(chromSize$len)

mat = read.table("wilding.unifiedGenotyper.cov14x.passQC.w200kb.s20kb.snpDensity.txt", h=T, fill = TRUE, dec = ".")
mat$pos = (mat$stop+mat$start)/2
head(mat)

mat.m = melt(mat, id.vars = c("chrom", "pos", "countsAccessPos"), measure.vars = c("singleton", "LBVwil", "LPdom", "LPfor", "shared"))
head(mat.m)
mat.m$chrom = factor(mat.m$chrom, levels = c("2R", "2L", "3R", "3L", "X"))

addcolors = c("bisque1", "gainsboro")
names(addcolors) = c("singleton", "shared")

group.colors = c(group.colors, addcolors)
group.colors = group.colors[which(names(group.colors) %in% levels(mat.m$variable))]

plot = ggplot() +
  geom_area(data=mat.m, aes(x=pos/1000000, y=value/countsAccessPos, color=variable, fill=variable)) +
  geom_rect(data=heterochrom, aes(xmin=start/1000000, xmax=stop/1000000, ymin=-Inf, ymax=Inf), alpha = 1, fill = 'grey') +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  facet_wrap(~chrom, ncol = 5, scales = "free_x") +
  ylab("SNP density") +
  xlab("position Mb") +
  theme_classic() +
  theme(
    #    strip.text.x = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

plot

gt3 = ggplot_gtable(ggplot_build(plot))

# take a look at the grob object's layout
#gtable::gtable_show_layout(gt)

gt3$widths[5] = gt3$widths[5]*chromSize$proportion[which(chromSize$chr=="2R")]
gt3$widths[9] = gt3$widths[9]*chromSize$proportion[which(chromSize$chr=="2L")]
gt3$widths[13] = gt3$widths[13]*chromSize$proportion[which(chromSize$chr=="3R")]
gt3$widths[17] = gt3$widths[17]*chromSize$proportion[which(chromSize$chr=="3L")]
gt3$widths[21] = gt3$widths[21]*chromSize$proportion[which(chromSize$chr=="X")]

gt1 = ggplot_gtable(ggplot_build(p1))
gt2 = ggplot_gtable(ggplot_build(p2))

lay <- rbind(c(1,2),
             c(3,3))

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/wilding.genomeCoverage_sexeDiscrimination_SNPdensity.pdf", width = conv_unit(190, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(gt1, gt2, gt3, layout_matrix = lay)
dev.off()

####################################
### SupFig2: IBD kinship analysis
####################################
mat = read.table("wilding.chr3.ibd.genome.txt", header=TRUE, skip=0, sep="\t")
head(mat)

info <- read.table("wilding.samples.meta.txt", h=T, sep="\t")
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
mat$compLabel[which( (mat$pop_FID1=="LBVwil" & mat$pop_FID2=="LPdom") | (mat$pop_FID1=="LPdom" & mat$pop_FID2=="LBVwil") )] = "LBVwil-LPdom"
mat$compLabel[which( (mat$pop_FID1=="LBVwil" & mat$pop_FID2=="LPfor") | (mat$pop_FID1=="LPfor" & mat$pop_FID2=="LBVwil") )] = "LBVwil-LPfor"
mat$compLabel[which( (mat$pop_FID1=="LPdom" & mat$pop_FID2=="LPfor") | (mat$pop_FID1=="LPfor" & mat$pop_FID2=="LPdom") )] = "LPdom-LPfor"
mat$compLabel = factor(mat$compLabel, levels = unique(mat$compLabel)[c(1,3,6,2,4,5)])

cols = c("#1a8c8d","darkgoldenrod", "#71b53f", "#009e60", "#fcd116", "#3a75c4")
names(cols) = c("LPdom-LPfor", "LBVwil-LPdom", "LBVwil-LPfor", "LPfor", "LBVwil", "LPdom")

p = ggplot(mat, aes(x=compLabel, y=PI_HAT, color=compLabel, fill=compLabel)) +
  geom_jitter(width = 0.2) +
  geom_hline(yintercept = 0.7, linetype = "dotted") +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme_classic() +
  ylim(c(-0.01,1)) +
  xlab("") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1))
p

pdf("~/Dropbox/professional/05_wilding/paper/figures/supFig/wilding.ibd.pdf", width = conv_unit(150, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p, ncol = 1, nrow = 1, top = textGrob("IBD intra and inter pop Wilding", gp=gpar(fontsize=8,font=8)))
dev.off()

mat[which(mat$PI_HAT>0.4),]


####################################
### SupFig2: PCA
####################################

### A: Eigen barplot values PCA all samples
prefix = "ag1000g_phase2_col.wilding.merged.biallelic.3.pca"
inF = paste(prefix, ".eigenval.txt", sep="")
eigenval <- data.frame(read.table(inF, header=FALSE, skip=0, sep=" "))
sum(eigenval)
eigenval = eigenval$V1*100
eigenval = round(eigenval, digits = 2)
df = data.frame(eval = eigenval, coord = 1:length(eigenval))

p1 = ggplot(df, aes(x=coord, y=eval)) +
  geom_bar(stat = "identity", color = "cornsilk4", fill = "cornsilk4", width = 0.6) +
  theme_classic() +
  ylab("% of explained variance") +
  xlab("Principal Components")
p1

ggsave("~/Dropbox/professional/05_wilding/paper/figures/supFig/pca.ag1000col_wilding.eigenval.pdf", units = "mm", width = 95, height = 75, useDingbats=FALSE) # export PDF 7x9

### B: PCA and Eigen barplot values PCA all wilding samples only

prefix = "wilding.biallelic.3.pca"
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

#eigenvec = cbind(samples, eigenvec)
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
  EV5 = merge$EV5,    # the second eigenvector
  EV6 = merge$EV6,    # the second eigenvector
  stringsAsFactors = FALSE
)

xlab = paste("EV1 (", round(eigenval[1], digits = 2), "%)", sep = "")
ylab = paste("EV2 (", round(eigenval[2], digits = 2), "%)", sep = "")

size = 2
alpha = 1
stroke = 0

p2 = ggplot() +
  geom_point(data=mat, aes(x=EV1, y=EV2, color=pop, fill=pop, shape=pop), size=size, alpha=alpha, stroke = stroke) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_shape_manual(values = c(24, 24, 24)) +
  xlab(xlab) +
  ylab(ylab) +
  theme_classic() +
  theme(legend.position="none")
# +  guides(color=guide_legend(ncol=2))

p2

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/wilding.pca.pdf", width = conv_unit(200, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p2, p1, ncol = 2, nrow = 1, top = "")
dev.off()


###############################
### Admixture cross validation:
###############################
cross <- read.table("ag1000gCol_wilding_crossValidationError.data.txt")
colnames(cross) = c("k", "error")
cross$k = factor(cross$k)
head(cross)

cross.avg = ddply(cross, .(k), plyr::summarize, avg=mean(error))
cross.avg$k = factor(cross.avg$k)

ggplot() +
  geom_jitter(data = cross, aes(x=k, y=error), position=position_jitter(width=.1, height=0), size=2) +
  geom_boxplot(data = cross, aes(x=k, y=error)) +
  geom_point(data = cross.avg, aes(x=k, y=avg), color="red", size=2) +
  theme_classic() +
  xlab("K") +
  ylab("Cross-validation Error")

ggsave("~/Dropbox/professional/wilding/paper/figures/supFig/wilding.ag1000g.col.crossValidationError.pdf", units = "mm", width = 100, height = 90, useDingbats=FALSE) # export PDF 7x9

###############################
### Admixture barplot:
###############################

admix = data.frame()
for(run in 2:10){
  run = 2
  print(run)
  inF = paste("ag1000gCol_wilding.admixture.K", run, ".txt", sep = "")
  ad = read.table(inF)
  colnames(ad)[1] = c("popId")
  colnames(ad)[2:dim(ad)[2]] = paste(rep("k", dim(ad)[2]-1), 1:(dim(ad)[2]-1), sep="")

  # 2.2 Plot bargraph admixture
  # 2.2.1 Reorder individu within each population to make barplot pretty
  ad$ind = 1:dim(ad)[1]
  reorder = c() # reorder barplot
  for(i in unique(ad$popId)){
    tmp = ad[which(ad$popId==i),]
    if(dim(tmp)[1]>1){
      reorder = c(reorder, tmp$ind[do.call(order, rev(as.list(tmp[,2:(dim(tmp)[2]-1)])))] )
    } else{
      reorder = c(reorder, tmp$ind)
    }
  }
  length(reorder)
  length(ad$popId)
  ad = ad[match(reorder, ad$ind),]
  ad$newOrder = 1:dim(ad)[1]
  head(ad)
  
  # 2.2.2 Add spacer betwen each pop
  spacerSize = 5
  tmp = ad
  for(i in unique(tmp$popId)){
    tmp$newOrder[-which(tmp$popId<i)] = tmp$newOrder[-which(tmp$popId<i)] + spacerSize
  }
  ad = tmp
  
  # 2.3 Make dataframe for population label
  labs = matrix(data=NA, nrow = length(unique(ad$popId)), ncol = 3)
  for(i in unique(ad$popId)){
    md = median(ad$newOrder[which(ad$popId==i)])
    labs[i, ] = c(md, i, -0.1)
  }
  colnames(labs) = c("xcoord", "popId", "ycoord")
  labs = data.frame(labs)
  labs = merge(labs, coord, by="popId")
  
  # 2.4 Plot bargraph
  ad.m = melt(ad, id.vars =c("newOrder", "popId"), measure.vars = colnames(tmp[,2:(dim(tmp)[2]-2)]))
  colnames(ad.m) = c("newOrder", "popId", "k", "p")
  
  cols = colsFrame[1:length(grep("k[0-9]", colnames(ad)))] # select colors
  names(cols) = colnames(ad)[grep("k[0-9]", colnames(ad))]
  head(ad.m)
  ad.m$run = run 
  admix = rbind(admix, ad.m)
}

colsFrame = getPalette(10)
cols = colsFrame[1:length(grep("k[0-9]", levels(admix$k)))] # select colors
names(cols) = levels(admix$k)
cols

ggplot() +
  geom_bar(admix, mapping = aes(x=newOrder, y=p, color=k, fill=k), stat = "identity", position = "stack") +
#  geom_text(data = labs, aes(x=xcoord, y=-0.05, label=popName), angle = 45, size=3, hjust=1 ) +
  ylim(c(-0.3,1.1)) +
  facet_wrap(~run, ncol = 1, strip.position="left") +
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) +
  theme(
    legend.position="none",
    axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("~/Dropbox/professional/wilding/paper/figures/supFig/wilding.ag1000g.col.admixturek2-10.pdf", units = "mm", width = 100, height = 90, useDingbats=FALSE) # export PDF 7x9

###############################
### Test if pops are panmictic:
###############################
null = read.table("wilding.nullDistri.txt", h=T, dec = ".")
head(null)

null.m = melt(null)
colnames(null.m) = c("popcomp", "chi2")
head(null.m)

chi2Obs = data.frame(
  chi = c(22779.07788979063, 18457850.560944583, 30468962.87874197),
  pvalue = c(0.45, 0, 0),
  popcomp = c("LPdom_LPfor", "LBVwil_LPdom", "LBVwil_LPfor"),
  x = c(150000, 1.5e7, 3e7),
  y = c(240, 140, 140)
)

cols = c("#1a8c8d","darkgoldenrod", "#71b53f")
names(cols) = c("LPdom_LPfor", "LBVwil_LPdom", "LBVwil_LPfor")

ggplot() +
  geom_histogram(data=null.m, aes(x=chi2, color=popcomp, fill=popcomp), binwidth=10000) +
  geom_vline(data=chi2Obs, aes(xintercept=chi), color="black") +
  geom_text(data = chi2Obs, aes(x = x, y=y, label=round(chi,4))) +
  facet_wrap(~popcomp, scales = "free", nrow = 3, ncol = 1) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  xlab("chi2 values") +
  ylab("count") +
  theme_classic() +
  theme(
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) 

ggsave("~/Dropbox/professional/wilding/paper/figures/supFig/testPopsPanmictic.pdf", units = "mm", width = 80, height = 160, useDingbats=FALSE) # export PDF 7x9


###############################
### FST matrix
###############################

### Load matrix
mat <- read.table("ag1000g_wilding.fst.matrix.txt", h=F, sep="\t")
#colnames(mat) = c('pop1', 'pop2', 'fst', 'zscore')
colnames(mat) = c('pop1', 'pop2', 'fst', 'se', 'zscore')
mat$pvalue = pnorm(q=mat$zscore, lower.tail=FALSE)
head(mat)

mat$pop1 = factor(mat$pop1, levels = c("BFcol", "GNcol", "GHcol", "CIcol", "LBVwil", "LPdom", "LPfor", "AOcol"))
mat$pop2 = factor(mat$pop2, levels = c("BFcol", "GNcol", "GHcol", "CIcol", "LBVwil", "LPdom", "LPfor", "AOcol"))

# Heatmap fst 
p1 = ggplot(data = mat, aes(pop2, pop1, fill= fst)) + 
  geom_tile(color="black") +
  geom_text(aes(label=round(fst, 2)), size=3) +
  scale_fill_gradient(low="blue", high="red") +
  coord_flip() +
  theme_classic() 
#+ theme(legend.position="none")
p1

mat$pop1 = factor(mat$pop1, levels = c("BFcol", "GNcol", "GHcol", "CIcol", "LBVwil", "LPdom", "LPfor", "AOcol"))

p2 = ggplot(data = mat, aes(y=pop1, x=pop2)) + 
  geom_tile(fill="white", color="black") +
  geom_text(aes(label=round(pvalue, 3)), size=3) +
  theme_classic() 
p2

pdf("~/Dropbox/professional/05_wilding/paper/figures/supFig/fst_matrix.pdf", width = conv_unit(200, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p1, p2, ncol = 2, nrow = 1, top = "")
dev.off()


###############################
### Stairplot2
###############################
order=c("AOcol", "LBVwil", "LPdom", "LPfor", "CIcol", "GHcol", "BFcol")

mat = read.table("wilding_ag1000g.plot.final.summary.txt", h=F, fill = TRUE, dec = ".")
colnames(mat) = c("pop", "mutation_per_site", "n_estimation", "theta_per_site_median", "theta_per_site_2.5", "theta_per_site_97.5", "year", "Ne_median", "Ne_2.5", "Ne_97.5", "Ne_12.5", "Ne_87.5")
head(mat)

mat = mat[which(mat$pop%in%order),]
mat$pop = factor(mat$pop, levels=order)

stwp = ggplot() +
  geom_line(data=mat, aes(x=year, y=Ne_median, color=pop), size=1) +
  geom_ribbon(data=mat, aes(x=year, ymin=Ne_2.5, ymax=Ne_97.5, color=pop, fill=pop), size=0, alpha=0.4) +
  facet_wrap(~pop) +
  scale_x_continuous(limits=c(1e3,1e5)) +
  scale_y_log10(limits=c(1e4,1e8), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = group.colors) +
  scale_fill_manual(values = group.colors) +
  ylab("Ne") +
  xlab("Years before present") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

stwp

ggsave("~/Dropbox/professional/wilding/paper/figures/supFig/wilding.ag1000g.stairwayPlot.pdf", units = "mm", width = 180, height = 120, useDingbats=FALSE) # export PDF 7x9

########################
### Find dadi best model
########################
mat = read.table("wilding_3_LPdom_LBVdom.dadi2D8Model.txt", h=F, dec = ".")
colnames(mat) = c("rep", "model", "optm", "stats", "value")
head(mat)

mat = mat[which(mat$optm=="BFGS"),]
mat.c = dcast(mat, rep+model~stats)
head(mat.c)

cols = brewer.pal(8, "Paired")
mat.c = mat.c[order(mat.c$optm_ll),]
mat.c$order = 1:dim(mat.c)[1]

pdadi = ggplot() +
  geom_bar(data=mat.c, aes(x=order, y=AIC, color=model, fill=model, width = 0.1), stat = "identity") +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme_classic()

pdadi

ggsave("~/Dropbox/professional/wilding/paper/figures/supFig/wilding.dadi.findBestModel.pdf", units = "mm", width = 100, height = 80, useDingbats=FALSE) # export PDF 7x9

########################
### Selection fst dxy
########################

group.colors = c("#00A585","darkorange2", "darkorchid3")
names(group.colors) = c("LPdom_LPfor", "LPdom_LBVwil", "LPfor_LBVwil")

mat_pop2 = read.table("wilding.1000snp.fst_dxy.scan.tab", h=T, fill = TRUE, dec = ".")
colnames(mat_pop2) = c("comp_pop", "chrom", "start", "stop", "fst", "dxy",  "xpehh_mean", "xpehh_md", "xpehh_std", "xpehh_max", "xpehh_min")
mat_pop2$pos = (mat_pop2$stop+mat_pop2$start)/2
mat_pop2$fst[which(mat_pop2$fst<0)] = 0
head(mat_pop2)

mat_pop2.m = melt(mat_pop2, id.vars = c("chrom", "pos", "comp_pop"), measure.vars = c("fst", "dxy"))
head(mat_pop2.m)
colnames(mat_pop2.m) = c("chrom", "pos", "comp_pop", "stat", "value")
mat_pop2.m = mat_pop2.m[!is.na(mat_pop2.m$value),] # remove NaN
mat_pop2.m$chrom = factor(mat_pop2.m$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
mat_pop2.m$comp_pop = factor(mat_pop2.m$comp_pop, levels = c("LPdom_LBVwil", "LPfor_LBVwil", "LPdom_LPfor"))

# keep only pos in euchro
head(mat_pop2.m)
dim(mat_pop2.m)
mat_pop2.m.f = data.frame()
for(line in 1:dim(euchrom)[1]){
  sub = mat_pop2.m[which(mat_pop2.m$chrom==euchrom$chr[line] & mat_pop2.m$pos>euchrom$start[line] & mat_pop2.m$pos<euchrom$stop[line]),]
  mat_pop2.m.f = rbind(mat_pop2.m.f, sub)
}
dim(mat_pop2.m.f)
head(mat_pop2.m.f)

### set up the y axis bounderies
fstr = range(mat_pop2.m.f[which(mat_pop2.m.f$stat=="fst"),]$value)
dxyr = range(mat_pop2.m.f[which(mat_pop2.m.f$stat=="dxy"),]$value)

dummy <- data.frame(
  chrom= rep(levels(mat_pop2.m.f$chrom), each=2),
  stat= c(rep("fst",5),rep("dxy",5)),
  min = c(rep(0,5), rep(0,5)),
  max = c(rep(roundUpNice(fstr[2]),5), rep(roundUpNice(dxyr[2]),5))
)

dummy
dummy.m = melt(data=dummy, id.vars = c("chrom", "stat"), measure.vars = c("min", "max"))
dummy.m$pos = 1
dummy.m$value[which(dummy.m$stat=="dxy" & dummy.m$variable=="max")] = 0.03

p1 = ggplot() +
  geom_rect(data=heterochrom, aes(xmin=start/1000000, xmax=stop/1000000, ymin=0, ymax=0.01), alpha = 1, fill = 'gray') +
  geom_line(data=mat_pop2.m.f, aes(x=pos/1000000, y=value, color=comp_pop), size=0.5) +
  geom_vline(data=rgenes, aes(xintercept=start/1000000), linetype="dotted") +
  geom_point(data=rgenes, aes(x=start/1000000, y=Inf), size=1.5, color="black", fill="black", shape=25) +
  geom_text_repel(data=rgenes, aes(x=(start+100000)/1000000, y=Inf, label=name), size= 3) +
  geom_blank(data=dummy.m, aes(x=pos, y=value)) +
  facet_grid(stat+comp_pop~chrom, scales="free") +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_classic() +
  theme(
#    strip.text.x = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p1

gt = ggplot_gtable(ggplot_build(p1))
gt$widths[5] = gt$widths[5]*chromSize$proportion[which(chromSize$chr=="2R")]
gt$widths[7] = gt$widths[7]*chromSize$proportion[which(chromSize$chr=="2L")]
gt$widths[9] = gt$widths[9]*chromSize$proportion[which(chromSize$chr=="3R")]
gt$widths[11] = gt$widths[11]*chromSize$proportion[which(chromSize$chr=="3L")]
gt$widths[13] = gt$widths[13]*chromSize$proportion[which(chromSize$chr=="X")]

y = arrangeGrob(gt, ncol = 1, nrow = 1)
grid.draw(y)

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/gscan_fst.pdf", width = conv_unit(180, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
y = arrangeGrob(gt, ncol = 1)
grid.draw(y)
dev.off()

########################
### Selection rehh rsb
########################
areaAroundRGenes = data.frame(
  chrom = c("3R", "2R", "2L", "2L"),
  start = c( 27599344, 27492278, 24399104, 1358158),
  stop = c(29599344, 29492278, 26399104, 4431617),
  name = c("Gste", "Cyp6p", "Gaba", "Vgsc")
)

chromSize <- read.table("Anopheles_gambiae.AgamP4.dna.chr.genome", h=F)
colnames(chromSize) = c("chr", "len")
k = c("2R", "2L", "3R", "3L", "X")
chromSize = chromSize[which(chromSize$chr %in% k),]
chromSize$proportion = chromSize$len/mean(chromSize$len)

euchrom = read.table("Anopheles_gambiae.AgamP4.dna.euchromatin.bed", h=F, fill = TRUE, dec = ".")
colnames(euchrom) = c("chrom", "start", "stop")
euchrom$chrom = factor(euchrom$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
head(euchrom)

### load rehh by SNP
rehh_df = data.frame()
for(p in c("LPdom.LBVwil", "LPfor.LBVwil", "LPdom.LPfor")){
  rehh_df_pop = data.frame()
  for(c in c("2L", "2R", "3L", "3R", "X")){
    path = paste(p, ".", c, ".gscann_rehh_bySNP.txt.gz", sep = "")
    zz=gzfile(path,'rt')
    mat=read.csv(zz,header=F,sep = "\t")
    head(mat)
    colnames(mat) = c("chrom", "pos", "rsb", "rsb_pvalue", "xpehh", "xpehh_pvalue", "comp_pop")
    mat = mat[complete.cases(mat), ]
    
    rehh_df_pop = rbind(rehh_df_pop, mat)
    
  }
  # keep only pos in euchro
  rehh_df_pop.f = data.frame()
  for(line in 1:dim(euchrom)[1]){
    sub = rehh_df_pop[which(rehh_df_pop$chrom==euchrom$chr[line] & rehh_df_pop$pos>euchrom$start[line] & rehh_df_pop$pos<euchrom$stop[line]),]
    rehh_df_pop.f = rbind(rehh_df_pop.f, sub)
  }
  
  # reformat df 
  rsb = rehh_df_pop.f[,c("chrom", "pos", "rsb", "rsb_pvalue", "comp_pop")]
  colnames(rsb) = c("chrom", "pos", "value", "pvalue", "comp_pop")
  rsb$stat= "rsb"
  xpehh = rehh_df_pop.f[,c("chrom", "pos", "xpehh", "xpehh_pvalue", "comp_pop")]
  colnames(xpehh) = c("chrom", "pos", "value", "pvalue", "comp_pop")
  xpehh$stat= "xpehh"
  rehh_df_pop.f.m = rbind(rsb, xpehh)
  
  ### pvalue transformation FDR and bonferonni
  rehh_df_pop.f.m$pvalue_back_trans = 10**(-rehh_df_pop.f.m$pvalue)  # transform back to p-values
  
  ### Perform p-value adjustments and add to data frame
  rehh_df_pop.f.m$Bonferroni = NA
  rehh_df_pop.f.m$fdr = NA
  for(s in unique(rehh_df_pop.f.m$stat)){
    sub = rehh_df_pop.f.m[which(rehh_df_pop.f.m$stat==s),]
    sub$Bonferroni = p.adjust(sub$pvalue_back_trans, method = "bonferroni")
    sub$fdr = p.adjust(sub$pvalue_back_trans, method = "BH")
    rehh_df_pop.f.m$Bonferroni[which(rehh_df_pop.f.m$stat==s)] = sub$Bonferroni
    rehh_df_pop.f.m$fdr[which(rehh_df_pop.f.m$stat==s)] = sub$fdr
  }
  
  # add color categories
  rehh_df_pop.f.m$colcat = NA
  rehh_df_pop.f.m$colcat[which(rehh_df_pop.f.m$fdr>=1e-4 & rehh_df_pop.f.m$value<0)] = ntile(rehh_df_pop.f.m$pvalue[which(rehh_df_pop.f.m$fdr>=1e-4 & rehh_df_pop.f.m$value<0)], 5)*-1
  rehh_df_pop.f.m$colcat[which(rehh_df_pop.f.m$fdr>=1e-4 & rehh_df_pop.f.m$value>=0)] = ntile(rehh_df_pop.f.m$pvalue[which(rehh_df_pop.f.m$fdr>=1e-4 & rehh_df_pop.f.m$value>=0)], 5)
  rehh_df_pop.f.m$colcat[which(rehh_df_pop.f.m$fdr<1e-4  & rehh_df_pop.f.m$value<0)] = (ntile(rehh_df_pop.f.m$pvalue[which(rehh_df_pop.f.m$fdr<1e-4 & rehh_df_pop.f.m$value<0)], 5)+5)*-1
  rehh_df_pop.f.m$colcat[which(rehh_df_pop.f.m$fdr<1e-4  & rehh_df_pop.f.m$value>=0)] = ntile(rehh_df_pop.f.m$pvalue[which(rehh_df_pop.f.m$fdr<1e-4 & rehh_df_pop.f.m$value>=0)], 5)+5
  rehh_df_pop.f.m$colcat = factor(rehh_df_pop.f.m$colcat)
  
  # discard a fraction of nosig SNP
  #  nb = length(which(rehh_df_pop.f.m$value > -3 & rehh_df_pop.f.m$value < 3 & rehh_df_pop.f.m$fdr>1e-4))
  #  nosig = rehh_df_pop.f.m[which(rehh_df_pop.f.m$value > -3 & rehh_df_pop.f.m$value < 3 & rehh_df_pop.f.m$fdr>1e-4)[seq(1, nb, 20)],]
  #  sig = rehh_df_pop.f.m[-which(rehh_df_pop.f.m$value > -3 & rehh_df_pop.f.m$value < 3 & rehh_df_pop.f.m$fdr>1e-4),]
  #  rehh.f = rbind(nosig, sig)
  rehh.f = rehh_df_pop.f.m
  
  rehh_df = rbind(rehh_df, rehh.f)
}
head(rehh_df)

### select only RSB value
rehh_df_rsb = rehh_df[which(rehh_df$stat=="rsb"),]

# set up the y axis bounderies
rehh_r = max(range(rehh_df_rsb$value))

dummy <- data.frame(
  chrom= rep(levels(rehh_df_rsb$chrom), each=2),
  stat= c(rep("rsb",5)),
  min = c(rep(-1*rehh_r[1],5)),
  max = c(rep(rehh_r[1],5))
)

dummy
dummy.m = melt(data=dummy, id.vars = c("chrom", "stat"), measure.vars = c("min", "max"))
dummy.m$pos = 1

### set colors
colfunc = colorRampPalette(c("gray", "royalblue4"))
colnosigneg = colfunc(10)[1:5]
names(colnosigneg) = seq(-1,-5)

colfunc = colorRampPalette(c("gray", "red4"))
colnosigpos = colfunc(10)[1:5]
names(colnosigpos) = seq(1,5)
plot(rep(1,10),col=colfunc(10),  pch=19,cex=2)

colfunc = colorRampPalette(c("mediumblue", "mediumblue"))
colsigneg = colfunc(5)
names(colsigneg) = seq(-6,-10)

colfunc = colorRampPalette(c("red", "red"))
colsigpos = colfunc(5)
names(colsigpos) = seq(6,10)

group.colors = c(colnosigneg, colnosigpos, colsigneg, colsigpos)

### plot
sig = rehh_df_rsb[which(rehh_df_rsb$fdr<1e-4),]
nosig = rehh_df_rsb[which(rehh_df_rsb$fdr>=1e-4),]
sig$comp_pop = factor(sig$comp_pop, levels = c("LPdom_LBVwil", "LPfor_LBVwil", "LPdom_LPfor"))
nosig$chrom = factor(nosig$chrom, levels = c("2R", "2L", "3R", "3L", "X"))

plot = ggplot() +
  #  geom_rect(data=heterochrom, aes(xmin=start/1000000, xmax=stop/1000000, ymin=-0.1, ymax=0.1), alpha = 1, fill = 'black') +
  #  geom_rect(data=areaAroundRGenes, aes(xmin=start/1000000, xmax=stop/1000000, ymin=-Inf, ymax=Inf), alpha = 1, fill = 'gray') +
  geom_point(data=nosig, aes(x=pos/1000000, y=value, color=colcat, fill=colcat), size=0.2) +
  geom_point(data=sig, aes(x=pos/1000000, y=value, color=colcat, fill=colcat), size=0.5) +
  #  geom_vline(data=rgenes, aes(xintercept=start/1000000), linetype="dotted") +
  #  geom_point(data=dat_text, aes(x=start/1000000, y=Inf), size=1.5, color="black", fill="black", shape=25) +
  #  geom_text_repel(data=dat_text, aes(x=(start+100000)/1000000, y=Inf, label=name), size= 3) +
  geom_blank(data=dummy.m, aes(x=pos, y=value)) +
  facet_grid(comp_pop~chrom, scales="free") +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_blank())

# theme_classic() +
#     theme(
#     strip.text.x = element_blank(),
#     axis.title.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     axis.title.y=element_blank(),
#     legend.position = "none",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     strip.background = element_blank()
#   )

plot

gt = ggplot_gtable(ggplot_build(plot))
gt$widths[5] = gt$widths[5]*chromSize$proportion[which(chromSize$chr=="2R")]
gt$widths[7] = gt$widths[7]*chromSize$proportion[which(chromSize$chr=="2L")]
gt$widths[9] = gt$widths[9]*chromSize$proportion[which(chromSize$chr=="3R")]
gt$widths[11] = gt$widths[11]*chromSize$proportion[which(chromSize$chr=="3L")]
gt$widths[13] = gt$widths[13]*chromSize$proportion[which(chromSize$chr=="X")]

# pdf("~/Dropbox/professional/wilding/paper/figures/figure4/gscan_xpehh.LPdom_LBVwil.pdf", width = conv_unit(190, "mm", "inch"), height = conv_unit(30, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
# y = arrangeGrob(gt, ncol = 1)
# grid.draw(y)
# dev.off()

tiff("~/Dropbox/professional/wilding/paper/figures/supFig/wilding.rehhMaxGap20kb.rsb.tiff", width = 3600, height = 2000, res = 300, compression = 'lzw')
y = arrangeGrob(gt)
grid.draw(y)
dev.off()

#######################
### Sup Figure XX: proportion on gscan rehh
#######################
rehh_df$rg ="out"
for(i in 1:dim(areaAroundRGenes)[1]){
  rehh_df$rg[which(rehh_df$chrom%in%areaAroundRGenes$chrom[i] & rehh_df$pos>areaAroundRGenes$start[i] & rehh_df$pos<areaAroundRGenes$stop[i])] = "in"
}
head(rehh_df)

rehh_df.t = ddply(rehh_df, .(comp_pop, stat), plyr::summarize, total=length(comp_pop))
sig = rehh_df[which(rehh_df$fdr<1e-4),]
sig$sign = "+"
sig$sign[which(sig$value<0)] = "-"
head(sig)
sig.d = ddply(sig, .(comp_pop, stat, rg, sign), plyr::summarize, count=length(comp_pop))

rehh.prop = merge(sig.d, rehh_df.t, by=c("comp_pop", "stat") )
rehh.prop$comp_pop = factor(rehh.prop$comp_pop, levels=c("LPdom_LBVwil", "LPfor_LBVwil", "LPdom_LPfor"))
rehh.prop$rg = factor(rehh.prop$rg, levels=c("in", "out"))
rehh.prop$stat = factor(rehh.prop$stat, levels=c("xpehh", "rsb"))
rehh.prop$variables = factor(paste(rehh.prop$rg, rehh.prop$sign, rehh.prop$comp_pop, sep = "_"))

# set colors
group.colors = c("#fcd116", NA, "#fcd116", NA, "#fcd116", "#009e60", "#fcd116", "#3a75c4", "#3a75c4", "#009e60")
names(group.colors) = levels(rehh.prop$variables)

# order factor
rehh.prop$variables = factor(rehh.prop$variables, levels = c("out_+_LPdom_LPfor", "in_-_LPdom_LPfor", "in_+_LPdom_LPfor", "out_-_LPdom_LPfor", "out_+_LPdom_LBVwil", "out_+_LPfor_LBVwil", "out_-_LPdom_LBVwil", "out_-_LPfor_LBVwil", "in_-_LPdom_LBVwil", "in_-_LPfor_LBVwil"))

# plot proportion and localization
plot = ggplot() +
  geom_bar(data = rehh.prop, aes(x=comp_pop, y=count/total, color=variables, fill=variables), stat = "identity", position ="stack", width = 0.8) +
  theme_classic() +
  facet_wrap(~stat) +
  ylab("Proportion of significant SNPs") +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme(
    #    axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=0.99),
    legend.position="none",
    axis.text.y = element_text(size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
plot

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/rehh_proportion.pdf", width = conv_unit(140, "mm", "inch"), height = conv_unit(80, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
plot
dev.off()

# chisquare test

rehh.prop$nosig=rehh.prop$total-rehh.prop$count

test.df = ddply(rehh.prop, .(comp_pop, stat, total), plyr::summarize, sig=sum(count))
test.df$nosig=test.df$total-test.df$sig
test.df.t = t(test.df[which(test.df$stat=="xpehh"), c("sig", "nosig")])
colnames(test.df.t) = c("LPdom_LBVwil", "LPfor_LBVwil", "LPdom_LPfor")

chisq.test(test.df.t)

#######################
### Figure 4X: sig SNP annotation
#######################
zz=gzfile('wilding.annot.snpEff.short.txt.gz','rt')
annotsnpeff=read.csv(zz,header=F,sep = "\t")
head(annotsnpeff)
colnames(annotsnpeff) = c("chrom", "pos", "Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type")

# select significant SNPs
sig = rehh_df[which(rehh_df$fdr<1e-4),]
sig$sign = 1
sig$sign[which(sig$value<0)] = -1
dim(sig)
head(sig)

# merge with SNP annotation
sigAnnot = merge(sig, annotsnpeff, by=c("chrom", "pos"))
sigAnnot$Annotation = sub('()_variant', '\\1', sigAnnot$Annotation)
head(sigAnnot)

# by pop
sigAnnot$pop1 = sub("(.*)_.*", "\\1", sigAnnot$comp_pop)
sigAnnot$pop2 = sub(".*_(.*)", "\\1", sigAnnot$comp_pop)
sigAnnot$pop = NA
sigAnnot$pop[which(sigAnnot$value>0)] = sigAnnot$pop1[which(sigAnnot$value>0)]
sigAnnot$pop[which(sigAnnot$value<0)] = sigAnnot$pop2[which(sigAnnot$value<0)]

# get proportion of sig SNP per annotation categories
sigAnnot.d = ddply(sigAnnot, .(comp_pop, Annotation), plyr::summarize, count=length(chrom))
sigAnnot.t = ddply(sigAnnot, .(comp_pop), plyr::summarize, total=length(chrom))
sigAnnot.m = merge(sigAnnot.d, sigAnnot.t, by=c("comp_pop"))

# sub set to main interesting annotation
sigAnnot.m$Annotation = sub('()_variant', '\\1', sigAnnot.m$Annotation)

df = data.frame(
  cat = c("coding", "coding", "coding", "coding", "coding", "coding", "noncoding", "noncoding"),
  Annotation = c("5_prime_UTR", "3_prime_UTR", "downstream_gene", "upstream_gene", "synonymous", "missense", "intron", "intergenic_region")
)

sigAnnot.m = merge(sigAnnot.m, df, by="Annotation")
sigAnnot.m$Annotation = factor(sigAnnot.m$Annotation, levels = df$Annotation)

group.colors = c("#fcd116", "#3a75c4", "#009e60")
names(group.colors) = c("LBVwil", "LPdom", "LPfor")

p = ggplot() +
  geom_bar(data = sigAnnot.m[-grep('UTR', sigAnnot.m$Annotation),], aes(x=Annotation, y=count/total, color=comp_pop, fill=comp_pop), stat = "identity", position ="dodge", width = 0.8) +
  theme_classic() +
  #  scale_color_manual(values=group.colors) +
  #  scale_fill_manual(values=group.colors) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=0.99)
  ) +
  theme(
    strip.text.x = element_blank(),
#    axis.text.x = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/wilding.areaDiff.SNPeff.xpehh.pdf", width = conv_unit(100, "mm", "inch"), height = conv_unit(60, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p)
dev.off()

### test chi2
sigAnnot.m
sigAnnot.d = dcast(sigAnnot.m, comp_pop~Annotation, value.var = "count")
row.names(sigAnnot.d) = sigAnnot.d$comp_pop

khi_test = chisq.test(sigAnnot.d[,-1])
khi_test

khi_test = chisq.test(sigAnnot.d[-1,-1])
khi_test

########################
### Diploshic confusion matrix
########################
conf = data.frame()
for(p in c("LBVwil", "LPdom", 'LPfor')){
  for(i in 1:10){
    for(j in 1:10){
      prefix = paste(p, i, j, sep = "_")
      inF = paste0(prefix, ".confusionFile.txt")
      print(inF)
      mat <- read.table(inF, h=F)
      colnames(mat) = c("Hard", "Hard-linked", "Soft", "Soft-linked", "Neutral")
      mat$row = c("Hard", "Hard-linked", "Soft", "Soft-linked", "Neutral")
      mat$prefix = prefix
      mat = melt(mat)
      mat$pop = p
      colnames(mat) = c("true_label", "prefix", "predict_label", "value", "pop")
      conf = rbind(conf, mat)
    }
  }
}
head(conf)

conf.avg = ddply(conf, .(pop, true_label, predict_label), plyr::summarize, N=length(value), mean=mean(value), sd=sd(value), ci = 1.960*(sd/sqrt(N)) )

conf.avg$true_label = factor(conf.avg$true_label, levels = c("Hard", "Hard-linked", "Soft", "Soft-linked", "Neutral"))
conf.avg$predict_label = factor(conf.avg$predict_label, levels = c("Hard", "Hard-linked", "Soft", "Soft-linked", "Neutral"))

group.colors = c("#fcd116", "#3a75c4", "#009e60")
names(group.colors) = c("LBVwil", "LPdom", "LPfor")

p1 = ggplot(data = conf.avg, aes(x=predict_label, y=mean, color=pop, fill=pop)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), color="black", width=.2, position=position_dodge(.9)) +
  theme_classic() +
  facet_wrap(~true_label, ncol=5) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=0.99),
    axis.text.y = element_text(size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) +
  xlab("Predicted label") +
  ylab("Proportion of windows asign to each predicted label") +
  ggtitle("True label")

p1
pdf("~/Dropbox/professional/wilding/paper/figures/supFig/wilding.barplot.confusionMatrix.pdf", width = conv_unit(180, "mm", "inch"), height = conv_unit(90, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p1, ncol = 1, nrow = 1, top = "")
dev.off()

########################
### Diploshic goodness of fit
########################

### load diploshic predictions
zz=gzfile('wilding.preds.100reps.txt.gz','rt')
mat=read.csv(zz,header=F,sep = "\t")
head(mat)
colnames(mat) = c('chrom', 'classifiedWinStart', 'classifiedWinEnd', 'bigWinRange', 'predClass', 'prob_neutral', 'prob_linkedSoft', 'prob_linkedHard', 'prob_soft', 'prob_hard', 'pop', 'rep')
vec = c("neutral", "linkedSoft", "soft", "linkedHard", "hard")
mat$predClass = factor(mat$predClass, levels=vec)
head(mat)

zz=gzfile('wilding.allpops.diploid.noscaled.w110000.win5.fvec.txt.gz','rt')
emp <- read.csv(zz,header=T,sep = "\t")
dim(emp)
head(emp)

### For each pop
plot_list = list()
for(p in c("LBVwil", "LPdom", "LPfor")){
  inF = paste0("~/bioInf/wilding/selection/diploshic/goodfit/", p, ".train_1.allpredClass.fvec.gz")
  sim <- read.table(inF, h=T)
  dim(sim)
  sim$id = 1:dim(sim)[1]
  sim.m <- melt(sim, id.vars = c("id", "predClass"))
  sim.m$statdesc = sub("(.*)_.*","\\1",sim.m$variable)
  sim.m$win = sub(".*_(.*)","\\1",sim.m$variable)
  sim.m$dt = "simulation"
  sim.m = sim.m[which(sim.m$win=="win5"),]
  sim.m$predClass = factor(sim.m$predClass, levels=c(levels(sim.m$predClass), "neutral"))
  sim.m$predClass = as.character(sim.m$predClass)
  sim.m$predClass[which(sim.m$predClass=="neut")] = "neutral"
  sim.m$predClass = factor(sim.m$predClass, levels = vec)
  sim.m$prob_neutral = 0
  sim.m = sim.m[,c(1,2,8,3,4,5,6,7)]
  head(sim.m)
  
  sub_emp = emp[which(emp$pop==p),]
  head(sub_emp)
  
  matsub = mat[which(mat$pop==p & mat$rep=="1_1"),]
  emp.merge = merge(sub_emp, matsub, by=c("chrom", "classifiedWinStart", "classifiedWinEnd", "bigWinRange", "pop"))
  head(emp.merge)
  emp.merge = emp.merge[,6:19]
  head(emp.merge)
  
  emp.merge$id = 1:dim(emp.merge)[1]
  head(emp.merge)
  emp.m <- melt(emp.merge, id.vars = c("id", "predClass", "prob_neutral"))
  emp.m$statdesc = sub("(.*)_.*","\\1",emp.m$variable)
  emp.m$win = sub(".*_(.*)","\\1",emp.m$variable)
  emp.m$dt = "empirical"
  head(emp.m)
  
  # merge empircal and simu
  df = rbind(sim.m, emp.m)
  head(df)
  
  # distribution of stat desc simu and emp 
  df_prop = data.frame()
  for(d in unique(df$dt)){
    for(p in unique(df$predClass)){
      for(s in unique(df$statdesc)){
        sub = df[which(df$predClass==p & df$statdesc==s & df$dt==d),]
        h = hist(sub$value, breaks = 10, plot=F)
        tmp = data.frame(breaks=h$breaks[2:length(h$breaks)], prop=h$counts/sum(h$counts))
        tmp$predClass = p 
        tmp$dt = d
        tmp$statdesc = s
        df_prop = rbind(df_prop, tmp)
      }
    }
  }
  vec = c("neutral", "linkedSoft", "soft", "linkedHard", "hard")
  df_prop$predClass = factor(df_prop$predClass, levels=vec)
  
  p1 = ggplot(data = df_prop, aes(x=breaks, y=prop, color=dt, fill=dt)) +  
    geom_line(size=0.8) + 
#    geom_bar(stat = "identity", position = "dodge") + 
    facet_grid(predClass~statdesc, scales = "free") +
    scale_color_manual(values=c("mediumvioletred", "mediumseagreen")) +
    theme_classic() +
    theme(
      legend.position="none",
      axis.text.y = element_text(size=8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank()
    ) 
  p1
  
  gt = ggplot_gtable(ggplot_build(p1))
  plot_list[[length(plot_list)+1]] = gt
  
  # PCA
  df.d = dcast(df, id+dt+predClass+prob_neutral+win~statdesc, value.var = "value")
  df.d$id = 1:dim(df.d)
  head(df.d)
  info = df.d[,c(1,2,3,4)]
  info$prob_neutral_class = "neut proba<0.01"
  info$prob_neutral_class[which(info$prob_neutral>0.01 & info$predClass!="neutral")] = "neut proba>0.01"
  head(info)
  
  df.d.pca = df.d[,c(6:dim(df.d)[2])]
  head(df.d.pca)
  
  res.pca <- PCA(df.d.pca,  graph = FALSE, scale.unit = T)
  
  p2 = fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
  gt = ggplot_gtable(ggplot_build(p2))
  plot_list[[length(plot_list)+1]] = gt
  
  p3 = fviz_pca_var(res.pca, axes = c(1, 2), col.var="contrib")+
    scale_color_gradient2(low="white", mid="blue",
                          high="red", midpoint=1.25, space ="Lab") +
    theme_minimal()
  gt = ggplot_gtable(ggplot_build(p3))
  plot_list[[length(plot_list)+1]] = gt
  
  eigenvec = data.frame(res.pca$ind$coord)
  eigenvec$id = row.names(eigenvec)
  eigenvec = merge(eigenvec, info, by="id")
  head(eigenvec)
  
  cols = c("#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C", "gray40")
  names(cols) = c("linkedSoft", "soft", "linkedHard", "hard", "neutral")
  eigenvec$predClass = factor(eigenvec$predClass, levels=c("neutral", "neutral_proba", "linkedSoft", "soft", "linkedHard", "hard"))
  eigenvec$grp = paste(eigenvec$dt, eigenvec$prob_neutral_class, sep = "_")
  
  p4 = ggplot() +
    geom_point(data=eigenvec, aes(x=Dim.1, y=Dim.2, color=predClass, fill=predClass, shape=prob_neutral_class)) +
    scale_color_manual(values=cols) +
    scale_fill_manual(values=cols) +
    scale_shape_manual(values=c(20,3)) +
    #  xlab(xlab) +
    #  ylab(ylab) +
    theme_classic() +
    facet_grid(grp~predClass) +
    theme(
      legend.position="none",
      axis.text.y = element_text(size=8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank()
    )
  gt = ggplot_gtable(ggplot_build(p4))
  plot_list[[length(plot_list)+1]] = gt
}

length(plot_list)

lay <- rbind(c(1,2,3))


pdf("~/Dropbox/professional/wilding/paper/figures/supFig/diploshicGoodfitStatDesc.pdf", width = conv_unit(400, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE)
y = arrangeGrob(plot_list[[1]], ncol = 1)
grid.draw(y)
dev.off()

lay <- rbind(c(1,2),
             c(3,3),
             c(3,3))

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/LBVwil.diploshicGoodfitPCA.pdf", width = conv_unit(180, "mm", "inch"), height = conv_unit(140, "mm", "inch"), useDingbats=FALSE)
y = arrangeGrob(plot_list[[2]], plot_list[[3]], plot_list[[4]], layout_matrix = lay)
grid.draw(y)
dev.off()

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/diploshicGoodfitEigenVal.pdf", width = conv_unit(180, "mm", "inch"), height = conv_unit(60, "mm", "inch"), useDingbats=FALSE)
y = arrangeGrob(plot_list[[2]], plot_list[[6]], plot_list[[10]], layout_matrix = lay)
grid.draw(y)
dev.off()

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/diploshicGoodfitVariables.pdf", width = conv_unit(240, "mm", "inch"), height = conv_unit(80, "mm", "inch"), useDingbats=FALSE)
y = arrangeGrob(plot_list[[3]], plot_list[[7]], plot_list[[11]], layout_matrix = lay)
grid.draw(y)
dev.off()

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/diploshicGoodfitEigenVec.pdf", width = conv_unit(180, "mm", "inch"), height = conv_unit(60, "mm", "inch"), useDingbats=FALSE)
y = arrangeGrob(plot_list[[4]], plot_list[[8]], plot_list[[12]], layout_matrix = lay)
grid.draw(y)
dev.off()


########################
### Diploshic rgenes examples
########################
# set up colors
cols = c("#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C", "gray90")
names(cols) = c("linkedSoft", "soft", "linkedHard", "hard", "neutral")

# load diploshic predictions
zz=gzfile('~/bioInf/wilding/selection/diploshic/wilding.preds.100reps.txt.gz','rt')
mat=read.csv(zz,header=F,sep = "\t")
head(mat)
colnames(mat) = c('chrom', 'classifiedWinStart', 'classifiedWinEnd', 'bigWinRange', 'predClass', 'prob_neutral', 'prob_linkedSoft', 'prob_linkedHard', 'prob_soft', 'prob_hard', 'pop', 'rep')
vec = c("neutral", "linkedSoft", "soft", "linkedHard", "hard")
mat$predClass = factor(mat$predClass, levels=vec)
head(mat)

### Keep only windows in present in euchromatine
euchrom = read.table("~/bioInf/wilding/ressources/Anopheles_gambiae.AgamP4.dna.euchromatin.bed", h=F, fill = TRUE, dec = ".")
colnames(euchrom) = c("chrom", "start", "stop")
euchrom$chrom = factor(euchrom$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
head(euchrom)

mat_f = data.frame()
for(i in 1:dim(euchrom)[1]){
  tmp = mat[which(mat$chrom==euchrom$chrom[i] & mat$classifiedWinStart>=euchrom$start[i] & mat$classifiedWinEnd<=euchrom$stop[i]),]
  mat_f = rbind(mat_f, tmp)
}
head(mat_f)

### look at example of knwon hard and soft sweep
head(mat_f)

plot_ex <- function(chr, st, sp, p) {
  mat_f_neu = mat_f
  mat_f_neu$predClass[which(mat_f_neu$prob_neutral>p)] = "neutral"
  range = 500000
  
  submat = mat_f_neu[which(mat_f_neu$chrom==chr & mat_f_neu$classifiedWinStart>(st-range) & mat_f_neu$classifiedWinEnd<(sp+range)),]
  submat.d = ddply(submat, .(chrom, classifiedWinStart, classifiedWinEnd, pop, predClass), plyr::summarize, count=length(predClass))
  submat.d$variable="classification"
  
  plot <- ggplot() +
    geom_bar(data=submat.d, aes(x=( (classifiedWinStart+classifiedWinEnd)/2)/100000, y=count, color=predClass, fill=predClass), stat = "identity", position = "stack") +
    facet_wrap(~pop, ncol = 3, scales = "free_y") +
    scale_color_manual(values=cols) +
    scale_fill_manual(values=cols) +
    geom_vline(data=targetRegions[i,], aes(xintercept=st/100000), size=0.8, linetype="dotted") +
    geom_vline(data=targetRegions[i,], aes(xintercept=sp/100000), size=0.8, linetype="dotted") +
    xlab("Pos") +
    ylab("# reps") +
    theme_bw() +
    theme(
      legend.position="none",
      axis.text.y = element_text(size=8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank()
    ) +
    ggtitle(paste0("Neutral proba<", p))
  return(plot)
}

chr = rgenes$chrom[which(rgenes$name=="Gste")]
st = rgenes$start[which(rgenes$name=="Gste")]
sp = rgenes$stop[which(rgenes$name=="Gste")]

p1 = plot_ex(chr, st, sp, 1)
p2 = plot_ex(chr, st, sp, 0.1)
p3 = plot_ex(chr, st, sp, 0.01)
p4 = plot_ex(chr, st, sp, 0.001)
y1 = grid.arrange(p1, p2, p3, p4, ncol = 1, nrow = 4, top = "Gste")

chr = rgenes$chrom[which(rgenes$name=="Cyp6p")]
st = rgenes$start[which(rgenes$name=="Cyp6p")]
sp = rgenes$stop[which(rgenes$name=="Cyp6p")]

p1 = plot_ex(chr, st, sp, 1)
p2 = plot_ex(chr, st, sp, 0.1)
p3 = plot_ex(chr, st, sp, 0.01)
p4 = plot_ex(chr, st, sp, 0.001)
y2 = grid.arrange(p1, p2, p3, p4, ncol = 1, nrow = 4, top = "Cyp6p")

chr = rgenes$chrom[which(rgenes$name=="Gaba")]
st = rgenes$start[which(rgenes$name=="Gaba")]
sp = rgenes$stop[which(rgenes$name=="Gaba")]

p1 = plot_ex(chr, st, sp, 1)
p2 = plot_ex(chr, st, sp, 0.1)
p3 = plot_ex(chr, st, sp, 0.01)
p4 = plot_ex(chr, st, sp, 0.001)
y3 = grid.arrange(p1, p2, p3, p4, ncol = 1, nrow = 4, top = "Gaba")

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/rgenesExample.pdf", width = conv_unit(200, "mm", "inch"), height = conv_unit(200, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(y1, ncol = 1, nrow = 1)
dev.off()

########################
### gscan XPEHH/RSB LPdom vs LPdom and LPfor vs LPfor
########################

### chrom size
chromSize <- read.table("~/bioInf/wilding/selection/Anopheles_gambiae.AgamP4.dna.chr.genome", h=F)
colnames(chromSize) = c("chr", "len")
k = c("2R", "2L", "3R", "3L", "X")
chromSize = chromSize[which(chromSize$chr %in% k),]
chromSize$proportion = chromSize$len/mean(chromSize$len)

### load annotation: resistance genes, heterochromatine
heterochrom = read.table("~/bioInf/wilding/ressources/Anopheles_gambiae.AgamP4.dna.heterochromatin.bed", h=F, fill = TRUE, dec = ".")
colnames(heterochrom) = c("chrom", "start", "stop")
heterochrom$chrom = factor(heterochrom$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
head(heterochrom)

euchrom = read.table("~/bioInf/wilding/ressources/Anopheles_gambiae.AgamP4.dna.euchromatin.bed", h=F, fill = TRUE, dec = ".")
colnames(euchrom) = c("chrom", "start", "stop")
euchrom$chrom = factor(euchrom$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
head(euchrom)

rgenes = read.table("~/bioInf/wilding/ressources/Anopheles_gambiae.AgamP4.dna.insecticide_genes.4plot.bed", h=F, fill = TRUE, dec = ".")
colnames(rgenes) = c("chrom", "start", "stop", "name", "geneid")
rgenes$chrom = factor(rgenes$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
head(rgenes)

areaAroundRGenes = data.frame(
  chrom = c("3R", "2R", "2L"),
  start = c(28298871, 28391415, 24363652),
  stop = c(28899816, 28593141, 26434556),
  name = c("Gste", "Cyp6p", "Gaba")
)

rehh_df = data.frame()
for(p in c("LPdom_1.LPdom_2", "LPfor_1.LPfor_2")){
  #for(cp in c("Gado_Zembe")){
  rehh_df_pop = data.frame()
  for(c in c("2L", "2R", "3L", "3R", "X")){
    path = paste("~/bioInf/wilding/selection/2pop/", p, ".", c, ".gscann_rehh_bySNP.txt.gz", sep = "")
    zz=gzfile(path,'rt')
    mat=read.csv(zz,header=F,sep = "\t")
    head(mat)
    colnames(mat) = c("chrom", "pos", "rsb", "rsb_pvalue", "xpehh", "xpehh_pvalue", "comp_pop")
    mat = mat[complete.cases(mat), ]
    
    rehh_df_pop = rbind(rehh_df_pop, mat)
    
  }
  # keep only pos in euchro
  rehh_df_pop.f = data.frame()
  for(line in 1:dim(euchrom)[1]){
    sub = rehh_df_pop[which(rehh_df_pop$chrom==euchrom$chr[line] & rehh_df_pop$pos>euchrom$start[line] & rehh_df_pop$pos<euchrom$stop[line]),]
    rehh_df_pop.f = rbind(rehh_df_pop.f, sub)
  }
  
  # reformat df 
  rsb = rehh_df_pop.f[,c("chrom", "pos", "rsb", "rsb_pvalue", "comp_pop")]
  colnames(rsb) = c("chrom", "pos", "value", "pvalue", "comp_pop")
  rsb$stat= "rsb"
  xpehh = rehh_df_pop.f[,c("chrom", "pos", "xpehh", "xpehh_pvalue", "comp_pop")]
  colnames(xpehh) = c("chrom", "pos", "value", "pvalue", "comp_pop")
  xpehh$stat= "xpehh"
  rehh_df_pop.f.m = rbind(rsb, xpehh)
  
  # discard a fraction of nosig SNP
  nb = length(which(rehh_df_pop.f.m$pvalue<=4))
  nosig = rehh_df_pop.f.m[which(rehh_df_pop.f.m$pvalue<=4)[seq(1, nb, 20)],]
  sig = rehh_df_pop.f.m[which(rehh_df_pop.f.m$pvalue>4),]
  dim(nosig)
  dim(sig)
  rehh.f = rbind(nosig, sig)
  head(rehh.f)
  
  # add color categories
  rehh.f$colcat = NA
  rehh.f$colcat[which(rehh.f$pvalue<=4 & rehh.f$value<0)] = ntile(rehh.f$pvalue[which(rehh.f$pvalue<=4 & rehh.f$value<0)], 5)*-1
  rehh.f$colcat[which(rehh.f$pvalue<=4 & rehh.f$value>=0)] = ntile(rehh.f$pvalue[which(rehh.f$pvalue<=4 & rehh.f$value>=0)], 5)
  rehh.f$colcat[which(rehh.f$pvalue>4  & rehh.f$value<0)] = (ntile(rehh.f$pvalue[which(rehh.f$pvalue>4 & rehh.f$value<0)], 5)+5)*-1
  rehh.f$colcat[which(rehh.f$pvalue>4  & rehh.f$value>=0)] = ntile(rehh.f$pvalue[which(rehh.f$pvalue>4 & rehh.f$value>=0)], 5)+5
  rehh.f$colcat = factor(rehh.f$colcat)
  head(rehh.f)
  
  rehh_df = rbind(rehh_df, rehh.f)
}

head(rehh_df)

unique(rehh_df$colcat)

### set colors
group.colors = c()

colfunc = colorRampPalette(c("gray", "royalblue4"))
colnosigneg = colfunc(10)[1:5]
names(colnosigneg) = seq(-1,-5)

colfunc = colorRampPalette(c("gray", "red4"))
colnosigpos = colfunc(10)[1:5]
names(colnosigpos) = seq(1,5)

colfunc = colorRampPalette(c("mediumblue", "mediumblue"))
colsigneg = colfunc(5)
names(colsigneg) = seq(-6,-10)

colfunc = colorRampPalette(c("red", "red"))
colsigpos = colfunc(5)
names(colsigpos) = seq(6,10)

group.colors = c(colnosigneg, colnosigpos, colsigneg, colsigpos)

### plot
pop2_plot_list = list()
for(cp in c("LPdom_1_LPdom_2", "LPfor_1_LPfor_2")){
  sig = rehh_df[which(rehh_df$comp_pop==cp & rehh_df$pvalue>4),]
  nosig = rehh_df[which(rehh_df$comp_pop==cp & rehh_df$pvalue<=4),]
  
  dat_text = rgenes
  dat_text$stat = sig$stat[1]
  
  plot = ggplot() +
    geom_rect(data=heterochrom, aes(xmin=start/1000000, xmax=stop/1000000, ymin=-Inf, ymax=Inf), alpha = 1, fill = 'gray85') +
    geom_rect(data=areaAroundRGenes, aes(xmin=start/1000000, xmax=stop/1000000, ymin=-Inf, ymax=Inf), alpha = 1, fill = 'gray85') +
    geom_point(data=sig, aes(x=pos/1000000, y=value, color=colcat, fill=colcat), size=0.2) +
    geom_point(data=nosig, aes(x=pos/1000000, y=value, color=colcat, fill=colcat), size=0.2) +
    #    geom_vline(data=rgenes, aes(xintercept=start/1000000), linetype="dotted") +
    geom_point(data=dat_text, aes(x=start/1000000, y=Inf), size=1.5, color="black", fill="black", shape=25) +
    geom_text_repel(data=dat_text, aes(x=(start+100000)/1000000, y=Inf, label=name), size= 3) +
    facet_grid(stat~chrom, scales="free") +
    scale_color_manual(values=group.colors) +
    scale_fill_manual(values=group.colors) +
    theme_classic() +
    theme(
      #    strip.text.x = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank()
    ) +
    ggtitle(cp)
  plot
  
  gt = ggplot_gtable(ggplot_build(plot))
  gt$widths[5] = gt$widths[5]*chromSize$proportion[which(chromSize$chr=="2R")]
  gt$widths[7] = gt$widths[7]*chromSize$proportion[which(chromSize$chr=="2L")]
  gt$widths[9] = gt$widths[9]*chromSize$proportion[which(chromSize$chr=="3R")]
  gt$widths[11] = gt$widths[11]*chromSize$proportion[which(chromSize$chr=="3L")]
  gt$widths[13] = gt$widths[13]*chromSize$proportion[which(chromSize$chr=="X")]
  pop2_plot_list[[length(pop2_plot_list)+1]] = gt
}


pdf("~/Dropbox/professional/wilding/paper/figures/supFig/mock_xpehh_rsb.pdf", width = conv_unit(190, "mm", "inch"), height = conv_unit(150, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
y = arrangeGrob(pop2_plot_list[[1]], pop2_plot_list[[2]], ncol = 1, nrow = 2)
grid.draw(y)
dev.off()

########################
### proportion sig SNP XPEHH/RSB
########################

mat = read.table("~/bioInf/wilding/selection/2pop/wilding.rehh_ctrl_analysis.txt", h=F, fill = TRUE, dec = ".")
colnames(mat) = c("chrom", "stat", "sigsnp", "snp", "comp_pop", "rep")
mat$comp_pop = factor(mat$comp_pop, levels=c("LPdom_1_LPdom_2", "LPfor_1_LPfor_2", "LPdom_LPfor"))
head(mat)
group.colors = c("chartreuse3", "#3a75c4", "#009e60")
names(group.colors) = c("LPdom_LPfor", "LPdom_1_LPdom_2", "LPfor_1_LPfor_2")

mat.d = ddply(mat, .(stat, comp_pop, rep), plyr::summarize, sigsnp=sum(sigsnp), snp=sum(snp))

ggplot() +
  geom_jitter(data=mat, aes(x=comp_pop, y=sigsnp/snp, color=comp_pop, fill=comp_pop), width = 0.2) +
  facet_grid(stat~chrom) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p

########################
### proportion sig SNP XPEHH/RSB
########################

euchrom = read.table("~/bioInf/wilding/ressources/Anopheles_gambiae.AgamP4.dna.euchromatin.bed", h=F, fill = TRUE, dec = ".")
colnames(euchrom) = c("chrom", "start", "stop")
euchrom$chrom = factor(euchrom$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
head(euchrom)

rehh_df = data.frame()
for(p in c("LPdom_1.LPdom_2", "LPfor_1.LPfor_2", "LPdom.LPfor")){
  #for(cp in c("Gado_Zembe")){
  rehh_df_pop = data.frame()
  for(c in c("2L", "2R", "3L", "3R", "X")){
    path = paste("~/bioInf/wilding/selection/2pop/", p, ".", c, ".gscann_rehh_bySNP.txt.gz", sep = "")
    zz=gzfile(path,'rt')
    mat=read.csv(zz,header=F,sep = "\t")
    head(mat)
    colnames(mat) = c("chrom", "pos", "rsb", "rsb_pvalue", "xpehh", "xpehh_pvalue", "comp_pop")
    mat = mat[complete.cases(mat), ]
    
    rehh_df_pop = rbind(rehh_df_pop, mat)
    
  }
  # keep only pos in euchro
  rehh_df_pop.f = data.frame()
  for(line in 1:dim(euchrom)[1]){
    sub = rehh_df_pop[which(rehh_df_pop$chrom==euchrom$chr[line] & rehh_df_pop$pos>euchrom$start[line] & rehh_df_pop$pos<euchrom$stop[line]),]
    rehh_df_pop.f = rbind(rehh_df_pop.f, sub)
  }
  
  # reformat df 
  rsb = rehh_df_pop.f[,c("chrom", "pos", "rsb", "rsb_pvalue", "comp_pop")]
  colnames(rsb) = c("chrom", "pos", "value", "pvalue", "comp_pop")
  rsb$stat= "rsb"
  xpehh = rehh_df_pop.f[,c("chrom", "pos", "xpehh", "xpehh_pvalue", "comp_pop")]
  colnames(xpehh) = c("chrom", "pos", "value", "pvalue", "comp_pop")
  xpehh$stat= "xpehh"
  rehh_df_pop.f.m = rbind(rsb, xpehh)
  
  rehh_df = rbind(rehh_df, rehh_df_pop.f.m)
}

head(rehh_df)
unique(rehh_df$comp_pop)

total = ddply(rehh_df, .(stat, comp_pop), plyr::summarize, count=length(chrom))
total$var = "total"
sig = ddply(rehh_df[which(rehh_df$pvalue>4),], .(stat, comp_pop), plyr::summarize, count=length(chrom))
sig$var = "sig"

df = dcast(rbind(total, sig),stat+comp_pop~var, value.var="count")

group.colors = c("chartreuse3", "#3a75c4", "#009e60")
names(group.colors) = c("LPdom_LPfor", "LPdom_1_LPdom_2", "LPfor_1_LPfor_2")

p = ggplot() +
  geom_bar(data=df, aes(x=comp_pop, y=sig/total, color=comp_pop, fill=comp_pop), stat = "identity", position = "dodge", width = 0.75) +
  facet_wrap(~stat) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/proportionSigSNP.pdf", width = conv_unit(120, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.draw(p)
dev.off()

#####################
##### supfig X: good fit
#####################

### load diploshic predictions
zz=gzfile('~/bioInf/wilding/selection/diploshic/wilding.preds.100reps.txt.gz','rt')
mat=read.csv(zz,header=F,sep = "\t")
head(mat)
colnames(mat) = c('chrom', 'classifiedWinStart', 'classifiedWinEnd', 'bigWinRange', 'predClass', 'prob_neutral', 'prob_linkedSoft', 'prob_linkedHard', 'prob_soft', 'prob_hard', 'pop', 'rep')
vec = c("neutral", "linkedSoft", "soft", "linkedHard", "hard")
mat$predClass = factor(mat$predClass, levels=vec)
head(mat)

### For LBVwil
sim <- read.table("~/bioInf/wilding/selection/diploshic/goodfit//LBVwil.train_1.allpredClass.fvec.gz", h=T)
dim(sim)
sim$id = 1:dim(sim)[1]
sim.m <- melt(sim, id.vars = c("id", "predClass"))
sim.m$statdesc = sub("(.*)_.*","\\1",sim.m$variable)
sim.m$win = sub(".*_(.*)","\\1",sim.m$variable)
sim.m$dt = "simulation"
sim.m = sim.m[which(sim.m$win=="win5"),]
sim.m$predClass = factor(sim.m$predClass, levels=c(levels(sim.m$predClass), "neutral"))
sim.m$predClass = as.character(sim.m$predClass)
sim.m$predClass[which(sim.m$predClass=="neut")] = "neutral"
sim.m$predClass = factor(sim.m$predClass, levels = vec)
sim.m$prob_neutral = 0
sim.m = sim.m[,c(1,2,8,3,4,5,6,7)]
head(sim.m)

zz=gzfile('~/bioInf/wilding/selection/diploshic/wilding.allpops.diploid.noscaled.w110000.win5.fvec.txt.gz','rt')
emp <- read.csv(zz,header=T,sep = "\t")
dim(emp)
head(emp)
emp = emp[which(emp$pop=="LBVwil"),]
head(emp)

matsub = mat[which(mat$pop=="LBVwil" & mat$rep=="1_1"),]
emp.merge = merge(emp, matsub, by=c("chrom", "classifiedWinStart", "classifiedWinEnd", "bigWinRange", "pop"))
head(emp.merge)
emp = emp.merge[,6:19]
head(emp)

emp$id = 1:dim(emp)[1]
head(emp)
emp.m <- melt(emp, id.vars = c("id", "predClass", "prob_neutral"))
emp.m$statdesc = sub("(.*)_.*","\\1",emp.m$variable)
emp.m$win = sub(".*_(.*)","\\1",emp.m$variable)
emp.m$dt = "empirical"
head(emp.m)

# distribution
df = rbind(sim.m, emp.m)
tail(df)
# df = df[which(df$win=="win5"),]

df_prop = data.frame()
for(d in unique(df$dt)){
  for(p in unique(df$predClass)){
    for(s in unique(df$statdesc)){
      sub = df[which(df$predClass==p & df$statdesc==s & df$dt==d),]
      h = hist(sub$value, breaks = 10, plot=F)
      tmp = data.frame(breaks=h$breaks[2:length(h$breaks)], prop=h$counts/sum(h$counts))
      tmp$predClass = p 
      tmp$dt = d
      tmp$statdesc = s
      df_prop = rbind(df_prop, tmp)
    }
  }
}
vec = c("neutral", "linkedSoft", "soft", "linkedHard", "hard")
df_prop$predClass = factor(df_prop$predClass, levels=vec)

p1 = ggplot(data = df_prop, aes(x=breaks, y=prop, color=dt, fill=dt)) +  
  geom_line(size=0.8) + 
  facet_grid(predClass~statdesc, scales = "free") +
  ylab("Proportion of windows") +
  xlab("Descriptive statistics value") +
  theme_classic() +
  theme(
#    axis.text.x = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p1

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/LBVwil.simulation_vs_empirical.pdf", width = conv_unit(190, "mm", "inch"), height = conv_unit(75, "mm", "inch"), useDingbats=FALSE)
grid.arrange(p1)
dev.off()

### For LPdom
sim <- read.table("~/bioInf/wilding/selection/diploshic/goodfit/LPdom.train_1.allpredClass.fvec.gz", h=T)
dim(sim)
sim$id = 1:dim(sim)[1]
sim.m <- melt(sim, id.vars = c("id", "predClass"))
sim.m$statdesc = sub("(.*)_.*","\\1",sim.m$variable)
sim.m$win = sub(".*_(.*)","\\1",sim.m$variable)
sim.m$dt = "simulation"
sim.m = sim.m[which(sim.m$win=="win5"),]
sim.m$predClass = factor(sim.m$predClass, levels=c(levels(sim.m$predClass), "neutral"))
sim.m$predClass = as.character(sim.m$predClass)
sim.m$predClass[which(sim.m$predClass=="neut")] = "neutral"
sim.m$predClass = factor(sim.m$predClass, levels = vec)
sim.m$prob_neutral = 0
sim.m = sim.m[,c(1,2,8,3,4,5,6,7)]
head(sim.m)

zz=gzfile('~/bioInf/wilding/selection/diploshic/wilding.allpops.diploid.noscaled.w110000.win5.fvec.txt.gz','rt')
emp <- read.csv(zz,header=T,sep = "\t")
dim(emp)
head(emp)
emp = emp[which(emp$pop=="LPdom"),]
head(emp)

matsub = mat[which(mat$pop=="LPdom" & mat$rep=="1_1"),]
emp.merge = merge(emp, matsub, by=c("chrom", "classifiedWinStart", "classifiedWinEnd", "bigWinRange", "pop"))
head(emp.merge)
emp = emp.merge[,6:19]
head(emp)

emp$id = 1:dim(emp)[1]
head(emp)
emp.m <- melt(emp, id.vars = c("id", "predClass", "prob_neutral"))
emp.m$statdesc = sub("(.*)_.*","\\1",emp.m$variable)
emp.m$win = sub(".*_(.*)","\\1",emp.m$variable)
emp.m$dt = "empirical"
head(emp.m)

# distribution
df = rbind(sim.m, emp.m)
tail(df)
# df = df[which(df$win=="win5"),]

df_prop = data.frame()
for(d in unique(df$dt)){
  for(p in unique(df$predClass)){
    for(s in unique(df$statdesc)){
      sub = df[which(df$predClass==p & df$statdesc==s & df$dt==d),]
      h = hist(sub$value, breaks = 10, plot=F)
      tmp = data.frame(breaks=h$breaks[2:length(h$breaks)], prop=h$counts/sum(h$counts))
      tmp$predClass = p 
      tmp$dt = d
      tmp$statdesc = s
      df_prop = rbind(df_prop, tmp)
    }
  }
}
vec = c("neutral", "linkedSoft", "soft", "linkedHard", "hard")
df_prop$predClass = factor(df_prop$predClass, levels=vec)

p1 = ggplot(data = df_prop, aes(x=breaks, y=prop, color=dt, fill=dt)) +  
  geom_line(size=0.8) + 
  facet_grid(predClass~statdesc, scales = "free") +
  ylab("Proportion of windows") +
  xlab("Descriptive statistics value") +
  theme_classic() +
  theme(
#    axis.text.x = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p1

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/LPdom.simulation_vs_empirical.pdf", width = conv_unit(190, "mm", "inch"), height = conv_unit(75, "mm", "inch"), useDingbats=FALSE)
grid.arrange(p1)
dev.off()

### For LPfor
sim <- read.table("~/bioInf/wilding/selection/diploshic/goodfit/LPfor.train_1.allpredClass.fvec.gz", h=T)
dim(sim)
sim$id = 1:dim(sim)[1]
sim.m <- melt(sim, id.vars = c("id", "predClass"))
sim.m$statdesc = sub("(.*)_.*","\\1",sim.m$variable)
sim.m$win = sub(".*_(.*)","\\1",sim.m$variable)
sim.m$dt = "simulation"
sim.m = sim.m[which(sim.m$win=="win5"),]
sim.m$predClass = factor(sim.m$predClass, levels=c(levels(sim.m$predClass), "neutral"))
sim.m$predClass = as.character(sim.m$predClass)
sim.m$predClass[which(sim.m$predClass=="neut")] = "neutral"
sim.m$predClass = factor(sim.m$predClass, levels = vec)
sim.m$prob_neutral = 0
sim.m = sim.m[,c(1,2,8,3,4,5,6,7)]
head(sim.m)

zz=gzfile('~/bioInf/wilding/selection/diploshic/wilding.allpops.diploid.noscaled.w110000.win5.fvec.txt.gz','rt')
emp <- read.csv(zz,header=T,sep = "\t")
dim(emp)
head(emp)
emp = emp[which(emp$pop=="LPfor"),]
head(emp)

matsub = mat[which(mat$pop=="LPfor" & mat$rep=="1_1"),]
emp.merge = merge(emp, matsub, by=c("chrom", "classifiedWinStart", "classifiedWinEnd", "bigWinRange", "pop"))
head(emp.merge)
emp = emp.merge[,6:19]
head(emp)

emp$id = 1:dim(emp)[1]
head(emp)
emp.m <- melt(emp, id.vars = c("id", "predClass", "prob_neutral"))
emp.m$statdesc = sub("(.*)_.*","\\1",emp.m$variable)
emp.m$win = sub(".*_(.*)","\\1",emp.m$variable)
emp.m$dt = "empirical"
head(emp.m)

# distribution
df = rbind(sim.m, emp.m)
tail(df)
# df = df[which(df$win=="win5"),]

df_prop = data.frame()
for(d in unique(df$dt)){
  for(p in unique(df$predClass)){
    for(s in unique(df$statdesc)){
      sub = df[which(df$predClass==p & df$statdesc==s & df$dt==d),]
      h = hist(sub$value, breaks = 10, plot=F)
      tmp = data.frame(breaks=h$breaks[2:length(h$breaks)], prop=h$counts/sum(h$counts))
      tmp$predClass = p 
      tmp$dt = d
      tmp$statdesc = s
      df_prop = rbind(df_prop, tmp)
    }
  }
}
vec = c("neutral", "linkedSoft", "soft", "linkedHard", "hard")
df_prop$predClass = factor(df_prop$predClass, levels=vec)

p1 = ggplot(data = df_prop, aes(x=breaks, y=prop, color=dt, fill=dt)) +  
  geom_line(size=0.8) + 
  facet_grid(predClass~statdesc, scales = "free") +
  ylab("Proportion of windows") +
  xlab("Descriptive statistics value") +
  theme_classic() +
  theme(
#    axis.text.x = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p1

pdf("~/Dropbox/professional/wilding/paper/figures/supFig/LPfor.simulation_vs_empirical.pdf", width = conv_unit(190, "mm", "inch"), height = conv_unit(75, "mm", "inch"), useDingbats=FALSE)
grid.arrange(p1)
dev.off()
