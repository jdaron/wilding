library(reshape)
library(ggplot2)
library(ggthemes)
#library(tidyverse)
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
library(scales)

setwd("~/bioInf/wilding/github/wilding/paper/input/figure3/")

##########################
### Figure 3: Demographic history
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

order=c("LBVwil", "LPdom", "LPfor")

### Figure 3A-B: stairplot2
###===============
mat = read.table("wilding_ag1000g.plot.final.summary.txt", h=F, fill = TRUE, dec = ".")
colnames(mat) = c("pop", "mutation_per_site", "n_estimation", "theta_per_site_median", "theta_per_site_2.5", "theta_per_site_97.5", "year", "Ne_median", "Ne_2.5", "Ne_97.5", "Ne_12.5", "Ne_87.5")
head(mat)

sub = mat[which(mat$pop%in%order),]
sub$pop = factor(sub$pop, levels=order)

stwp1 = ggplot() +
  geom_line(data=sub, aes(x=year, y=Ne_median, color=pop, alpha=pop), size=1) +
  geom_ribbon(data=sub, aes(x=year, ymin=Ne_2.5, ymax=Ne_97.5, color=pop, fill=pop), size=0, alpha=0.4) +
  #  facet_wrap(~pop) +
  #  scale_x_log10(limits=c(100,100000) ,breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
#  scale_x_continuous(limits=c(1e3,1e5)) +
  scale_y_log10(limits=c(1e4,1e8), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = group.colors) +
  scale_fill_manual(values = group.colors) +
  scale_alpha_manual(values = group.alpha) +
  ylab("Ne") +
  xlab("Years before present") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) +
  annotation_logticks(size = 0.2, sides = "l")

stwp1

sub = mat[-which(mat$pop%in%order),]
sub = sub[grep('col', sub$pop),]
sub$pop = factor(sub$pop, levels=levels(sub$pop)[grep('col',levels(sub$pop))])

sub$year[which(sub$year>100000)] = 100000

stwp2 = ggplot() +
  geom_line(data=sub, aes(x=year, y=Ne_median, color=pop, alpha=pop), size=1) +
  geom_ribbon(data=sub, aes(x=year, ymin=Ne_2.5, ymax=Ne_97.5, color=pop, fill=pop), size=0, alpha=0.4) +
  #  facet_wrap(~pop) +
  #  scale_x_log10(limits=c(100,100000) ,breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(limits=c(1e3,1e5)) +
  scale_y_log10(limits=c(1e4,1e8), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = group.colors) +
  scale_fill_manual(values = group.colors) +
  scale_alpha_manual(values = group.alpha) +
  ylab("Ne") +
  xlab("Years before present") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) +
  annotation_logticks(size = 0.2, sides = "l")

stwp2

pdf("~/Dropbox/professional/wilding/paper/figures/figure3/wilding_ag1000col.demographic_history.pdf", width = conv_unit(200, "mm", "inch"), height = conv_unit(70, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(stwp1, stwp2, ncol = 2, nrow = 1)
dev.off()
