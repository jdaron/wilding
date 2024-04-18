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
library(grid)
library(measurements)
library(ggpubr)
library(dplyr)
library(eulerr)

setwd("~/bioInf/wilding/github/wilding/paper/input/figure4/")

### function
roundUpNice <- function(x, nice=c(1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.2,2.3,2.4,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

### chrom size
chromSize <- read.table("Anopheles_gambiae.AgamP4.dna.chr.genome", h=F)
colnames(chromSize) = c("chr", "len")
k = c("2R", "2L", "3R", "3L", "X")
chromSize = chromSize[which(chromSize$chr %in% k),]
chromSize$proportion = chromSize$len/mean(chromSize$len)

### load annotation: resistance genes, heterochromatine
heterochrom = read.table("Anopheles_gambiae.AgamP4.dna.heterochromatin.bed", h=F, fill = TRUE, dec = ".")
colnames(heterochrom) = c("chrom", "start", "stop")
heterochrom$chrom = factor(heterochrom$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
head(heterochrom)

euchrom = read.table("Anopheles_gambiae.AgamP4.dna.euchromatin.bed", h=F, fill = TRUE, dec = ".")
colnames(euchrom) = c("chrom", "start", "stop")
euchrom$chrom = factor(euchrom$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
head(euchrom)

rgenes = read.table("Anopheles_gambiae.AgamP4.dna.insecticide_genes.4plot.bed", h=F, fill = TRUE, dec = ".")
colnames(rgenes) = c("chrom", "start", "stop", "name", "geneid")
rgenes$chrom = factor(rgenes$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
head(rgenes)

# annotate rgene regions
areaAroundRGenes = rgenes[which(rgenes$name%in%c("Gste", "Cyp6p", "Gaba")),]
areaAroundRGenes$start=areaAroundRGenes$start-1e6
areaAroundRGenes$start=areaAroundRGenes$stop+1e6

areaAroundRGenes = data.frame(
  chrom = c("3R", "2R", "2L"),
  start = c(28298871, 28391415, 24363652),
  stop = c(28899816, 28593141, 26434556),
  name = c("Gste", "Cyp6p", "Gaba")
)

areaAroundRGenes = data.frame(
  chrom = c("3R", "2R", "2L", "2L"),
  start = c( 27599344, 27492278, 24399104, 1358158),
  stop = c(29599344, 29492278, 26399104, 4431617),
  name = c("Gste", "Cyp6p", "Gaba", "Vgsc")
)

((areaAroundRGenes$start+areaAroundRGenes$stop)/2)+1e6

#######################
### Figure 4A: gscan XPEHH
#######################

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

rm(rehh_mat)
rm(rehh_xpehh)
rm(rehh_rsb)
rm(rehh.f)
rm(zz)
gc()

### select only XPEHH value
rehh_df_xpehh = rehh_df[which(rehh_df$stat=="xpehh"),]

# set up the y axis bounderies
rehh_r = max(range(rehh_df_xpehh$value))

dummy <- data.frame(
  chrom= rep(levels(rehh_df_xpehh$chrom), each=2),
  stat= c(rep("xpehh",5)),
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
dat_text = rgenes
dat_text$stat = "xpehh"

sig = rehh_df_xpehh[which(rehh_df_xpehh$fdr<1e-4),]
nosig = rehh_df_xpehh[which(rehh_df_xpehh$fdr>=1e-4),]
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

plot

gt = ggplot_gtable(ggplot_build(plot))
gt$widths[5] = gt$widths[5]*chromSize$proportion[which(chromSize$chr=="2R")]
gt$widths[7] = gt$widths[7]*chromSize$proportion[which(chromSize$chr=="2L")]
gt$widths[9] = gt$widths[9]*chromSize$proportion[which(chromSize$chr=="3R")]
gt$widths[11] = gt$widths[11]*chromSize$proportion[which(chromSize$chr=="3L")]
gt$widths[13] = gt$widths[13]*chromSize$proportion[which(chromSize$chr=="X")]

png("wilding.rehhMaxGap20kb.xpehh.png", width = 900, height = 500) # export PDF 7x9
y = arrangeGrob(gt)
grid.draw(y)
dev.off()

tiff("wilding.rehhMaxGap20kb.xpehh.tiff", width = 3600, height = 2000, res = 300, compression = 'lzw')
y = arrangeGrob(gt)
grid.draw(y)
dev.off()

##### plot pvalue
# transform colors
dim(nosig)
nosig$cols = "nosig"
sig$cols = NA
sig$cols[which(sig$comp_pop=="LPdom_LBVwil" & sig$value>0)] = "LPdom"
sig$cols[which(sig$comp_pop=="LPdom_LBVwil" & sig$value<0)] = "LBVwil"

sig$cols[which(sig$comp_pop=="LPfor_LBVwil" & sig$value>0)] = "LPfor"
sig$cols[which(sig$comp_pop=="LPfor_LBVwil" & sig$value<0)] = "LBVwil"

sig$cols[which(sig$comp_pop=="LPdom_LPfor" & sig$value>0)] = "LPdom"
sig$cols[which(sig$comp_pop=="LPdom_LPfor" & sig$value<0)] = "LPfor"

# colors set up 
group.colors = c("gray", "#009e60", "#fcd116", "#3a75c4")
names(group.colors) = c("nosig", "LPfor", "LBVwil", "LPdom")

# dummy
dummy <- data.frame(
  chrom= rep(levels(rehh_df_xpehh$chrom)),
  pos = 1,
  stat= c(rep("fdr",5)),
  min = c(rep(0,5)),
  max = c(rep(max(-1*log10(sig$fdr)),5))
)
dummy.m = melt(data=dummy, id.vars = c("chrom", "stat"), measure.vars = c("min", "max"))
dummy.m$pos = 1

th = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank()
)

plot = ggplot() +
  geom_point(data=nosig, aes(x=pos/1000000, y=-1*log10(fdr), color=cols, fill=cols), size=0.2) +
  geom_point(data=sig, aes(x=pos/1000000, y=-1*log10(fdr), color=cols, fill=cols), size=0.5) +
  geom_blank(data=dummy.m, aes(x=pos, y=value)) +
  facet_grid(comp_pop~chrom, scales="free") +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_classic() + th +
#  theme_void() +
  theme(
    legend.position = "none", 
    strip.text = element_blank())

gt = ggplot_gtable(ggplot_build(plot))
gt$widths[5] = gt$widths[5]*chromSize$proportion[which(chromSize$chr=="2R")]
gt$widths[7] = gt$widths[7]*chromSize$proportion[which(chromSize$chr=="2L")]
gt$widths[9] = gt$widths[9]*chromSize$proportion[which(chromSize$chr=="3R")]
gt$widths[11] = gt$widths[11]*chromSize$proportion[which(chromSize$chr=="3L")]
gt$widths[13] = gt$widths[13]*chromSize$proportion[which(chromSize$chr=="X")]

png("wilding.rehhMaxGap20kb.xpehhPvalue.png", width = 900, height = 500) # export PDF 7x9
y = arrangeGrob(gt)
grid.draw(y)
dev.off()

tiff("wilding.rehhMaxGap20kb.xpehhPvalue.tiff", width = 3600, height = 2000, res = 300, compression = 'lzw')
y = arrangeGrob(gt)
grid.draw(y)
dev.off()

pdf("gscan_xpehh.wilding.pdf", width = conv_unit(190, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
y = arrangeGrob(gt, ncol = 1)
grid.draw(y)
dev.off()

#######################
### Supplemntary Figure 11: proportion on gscan rehh
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

pdf("rehh_proportion.pdf", width = conv_unit(140, "mm", "inch"), height = conv_unit(50, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
plot
dev.off()

#######################
### Supplementary Figure 11 Significant SNPs annotation
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
dim(sigAnnot)

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
    axis.text.x = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p

pdf("wilding.areaDiff.SNPeff.pdf", width = conv_unit(100, "mm", "inch"), height = conv_unit(60, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
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

#######################
### Generate supplementary table
#######################
sigAnnot <- read.table("wilding.sigSNP.rehh.annot.txt", h=T)
head(sigAnnot)
dim(sigAnnot)
#sigAnnot = sigAnnot[which(sigAnnot$stat=="xpehh"),]
sigAnnot[which(sigAnnot$chrom=="3L" & sigAnnot$pos>41275000 & sigAnnot$pos<41500000),]

sigAnnot[which(sigAnnot$pos==3439855),]

candAreas <- read.table("~/bioInf/wilding/function/candidate_regions.rehh.txt", h=F)
colnames(candAreas) = c("chrom", "start", "end", "nbSNP", "comp_pop", "pop")
head(candAreas)

out = data.frame()
for(i in 1:dim(candAreas)[1]){
  tmp = sigAnnot[which(candAreas$chrom[i]==sigAnnot$chrom & candAreas$comp_pop[i]==sigAnnot$comp_pop & candAreas$pop[i]==sigAnnot$pop & candAreas$start[i] <= sigAnnot$pos & candAreas$end[i] >= sigAnnot$pos),]
  tmp = cbind(candAreas[i,], tmp)
  out = rbind(out, tmp)
}

out.d = dcast(out,chrom+start+end+nbSNP+comp_pop+pop+Gene_ID~Annotation)
head(out.d)

write.table(out.d, "~/bioInf/wilding/function/supTableFunction.txt", quote = F, sep = "\t")

m = min(out$fdr[which(out$chrom=="3R" & out$start==18575000 & out$end==19025000)])

out[which(out$chrom=="3R" & out$start==18575000 & out$end==19025000 & out$fdr==m),]

#######################
### Figure 4B: diploshic
#######################
# set up colors
cols = c("#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C", "gray90")
names(cols) = c("linkedSoft", "soft", "linkedHard", "hard", "neutral")

group.colors = c("#fcd116", "#3a75c4", "#009e60")
names(group.colors) = c("LBVwil", "LPdom", "LPfor")

# load diploshic predictions
zz=gzfile('wilding.w110000.sameModelLP.100reps.txt.gz','rt')
mat=read.csv(zz,header=F,sep = "\t")
colnames(mat) = c('chrom', 'classifiedWinStart', 'classifiedWinEnd', 'bigWinRange', 'predClass', 'prob_neutral', 'prob_linkedSoft', 'prob_linkedHard', 'prob_soft', 'prob_hard', 'pop', 'rep')
vec = c("neutral", "linkedSoft", "soft", "linkedHard", "hard")
mat$predClass = factor(mat$predClass, levels=vec)
head(mat)

# Keep only windows in present in euchromatine
euchrom = read.table("~/bioInf/wilding/ressources/Anopheles_gambiae.AgamP4.dna.euchromatin.bed", h=F, fill = TRUE, dec = ".")
colnames(euchrom) = c("chrom", "start", "stop")
euchrom$chrom = factor(euchrom$chrom, levels = c("2R", "2L", "3R", "3L", "X"))
head(euchrom)
mat$chrom

mat_f = data.frame()
for(i in 1:dim(euchrom)[1]){
  tmp = mat[which(mat$chrom==euchrom$chrom[i] & mat$classifiedWinStart>=euchrom$start[i] & mat$classifiedWinEnd<=euchrom$stop[i]),]
  mat_f = rbind(mat_f, tmp)
}
head(mat_f)

# stat desc
head(mat_f)
mat_f$newPredClass = mat_f$predClass
stat = ddply(mat_f, .(pop, newPredClass, rep), plyr::summarize, count=length(pop))
stat = ddply(stat, .(pop, newPredClass), plyr::summarize, avg=mean(count), sd=sd(count))
stat

total = ddply(mat_f, .(pop, newPredClass, rep), plyr::summarize, count=length(pop))
total = total[-which(total$newPredClass=="neutral"),]
total = ddply(total, .(pop, rep), plyr::summarize, avg=mean(count), sd=sd(count))
ddply(total, .(pop), plyr::summarize, avg=mean(avg), sd=sd(avg))
head(total)

# filtration based on probability
mat_f$newPredClass = mat_f$predClass
mat_f$newPredClass[which(mat_f$predClass!="neutral" & mat_f$prob_neutral>0.01)] = "neutral"
head(mat_f)

# filtration based on majority predClass
mat_f.d = ddply(mat_f, .(chrom, classifiedWinStart, classifiedWinEnd, pop, newPredClass), plyr::summarize, nbRep=length(pop))
head(mat_f.d)

#write.table(mat_f.d, "wilding.w110000.postprocessing.txt", quote = F, sep = "\t")

mat.count = ddply(mat_f.d[which(mat_f.d$nbRep>=50),], .(pop, newPredClass), plyr::summarize, count=length(pop))
colnames(mat.count) = c("pop", "predClass", "count")
mat.total = ddply(mat[which(mat$rep=="1_1"),], .(pop), plyr::summarize, total=length(pop))

mat.merge = merge(mat.count, mat.total, by=c("pop"))

p = ggplot() +
  geom_bar(data = mat.merge[which(mat.merge$predClass!="neutral"),], aes(x=predClass, y=count/total, color=pop, fill=pop), stat = "identity", position ="dodge") +
  facet_wrap(~predClass, scales = "free", ncol = 2) +
  theme_classic() +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme(
    strip.text = element_blank(),
    axis.title.y=element_blank(),
    strip.background = element_blank(),
    legend.position = "none")

pdf("wilding.diploshic.barplot.pdf", width = conv_unit(90, "mm", "inch"), height = conv_unit(70, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p)
dev.off()


