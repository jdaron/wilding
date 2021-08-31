library(ggplot2)
library(measurements)
library("gridExtra")

infow <- read.table("~/bioInf/wilding/ressources/wilding.samples.meta.txt", h=T)
infow$grp = "wilding"
head(infow)
infoj <- read.table("~/bioInf/wilding/ressources/jesus.samples.meta.txt", h=T, sep="\t")
head(infoj)
infoj$grp = "jesus"

info = rbind(infow[,c("id", "population", "mean_coverage", "grp")], infoj[,c("id", "population", "mean_coverage", "grp")])

### plot mean coverage per individual
mat = info
head(mat)
mat$id = factor(mat$id, levels = mat$id[order(mat$mean_coverage)])

p1 = ggplot(mat, aes(x=id, y=mean_coverage, color=population, fill=population)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 14) +
  theme_classic() +
  facet_wrap(~grp, scales = "free_x") +
  ylab("Mean Coverage") +
  xlab("Individuals") +
#  scale_fill_brewer(palette = "Set1") +
#  scale_color_brewer(palette = "Set1") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank())

p1
  
pdf("~/bioInf/wilding/fastqBamInfo/bamInfo/wilding_jesus_mean_coverage.pdf", width = conv_unit(180, "mm", "inch"), height = conv_unit(80, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p1, ncol = 1, nrow = 1, top = "Genome mean coverage")
dev.off()

### plot reads depth coverage per sample and per chromosome 
matw = read.table("~/bioInf/wilding/fastqBamInfo/bamInfo/wilding/wilding.coveragePerChr.tab", header=F, skip=0, sep="\t")
colnames(matw) = c("id", "chr", "chr_size", "bp_map", "mean_cov", "md_cov")

matj = read.table("~/bioInf/wilding/fastqBamInfo/bamInfo/jesus/jesus.coveragePerChr.tab", header=F, skip=0, sep="\t")
colnames(matj) = c("id", "chr", "chr_size", "bp_map", "mean_cov", "md_cov")
head(matj)

mat = rbind(matw, matj)

mat = merge(mat, info, by="id")
mat = mat[-which(mat$chr=="Y_unplaced" | mat$chr=="Mt"),]
head(mat)

p2 = ggplot(mat, aes(x=chr, y=log10(mean_cov/mean_coverage), color=population, fill=population)) +
  geom_jitter(width = 0.2) +
  theme_classic() +
  facet_wrap(~grp) +
  ylab("Ratio chr depth / total depth (log10)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
  
p2

pdf("~/bioInf/wilding/fastqBamInfo/bamInfo/wilding_jesus.chr_coverageDepth.pdf", width = conv_unit(180, "mm", "inch"), height = conv_unit(80, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p2, ncol = 1, nrow = 1, top = "Coverage depth per Chr and per Samples")
dev.off()

sub_chrX = mat[which(mat$chr=="X"),]
t = sub_chrX[which(log10(sub_chrX$mean_cov/sub_chrX$mean_coverage) < -0.1),]

write.csv(t, "~/bioInf/wilding/fastqBamInfo/bamInfo/tmp", quote = F)
