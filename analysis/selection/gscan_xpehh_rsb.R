#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#### R script to calculate XPEHH and RSB using rehh library
#===============================================================================
#   Author: Josquin Daron
#
#   File: gscan_xpehh_rsb.R
#   Date: 07-03-2022
#   Version: 0.1
#
#   Usage:
#	 R --vanilla --args 'vcf pop1' 'vcf pop2' 'pop1 name' 'pop2 name' 'output name and path' < gscan_xpehh_rsb.R 
#
#===============================================================================
options(digits=22)

###################
# SCAN  HAPLOTYPE #
###################
#Package
library(vcfR)
library(rehh)
library(readr)
library(plyr)
library(reshape2)

path_to_vcf_pop1 = args[1]
path_to_vcf_pop2 = args[2]
pop1 = args[3]
pop2 = args[4]
out = args[5]

#path_to_vcf_pop1 = "~/bioInf/wilding/selection/LBVwil.2L.sub.vcf.gz"
#path_to_vcf_pop2 = "~/bioInf/wilding/selection/LPdom.2L.sub.vcf.gz"
#pop1 = "LBV"
#pop2 = "LPdom"
#out = "~/bioInf/wilding/selection/test"

### slide function
slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- median(data[spots[i]:(spots[i]+window)], na.rm = T)
  }
  return(result)
}

file <- paste0(path_to_vcf_pop1)
hh <- data2haplohh(hap_file = file, polarize_vcf = FALSE, vcf_reader = "vcfR", allele_coding = "01")
#scan_pop1 <- scan_hh(hh, phased = T, polarized = F, maxgap=5000)
scan_pop1 <- scan_hh(hh, phased = T, polarized = F)
id <- paste0("scan_", pop1)
assign(id, scan_pop1)

file <- paste0(path_to_vcf_pop2)
hh <- data2haplohh(hap_file = file, polarize_vcf = FALSE, vcf_reader = "vcfR", allele_coding = "01")
#scan_pop2 <- scan_hh(hh, phased = T, polarized = F, maxgap=5000)
scan_pop2 <- scan_hh(hh, phased = T, polarized = F)
id <- paste0("scan_", pop2)
assign(id, scan_pop2)

#rsb
rsb <- ines2rsb(scan_pop1 = scan_pop1,
                        scan_pop2 = scan_pop2,
                        popname1 = pop1,
                        popname2 = pop2)

#xpehh
xpehh <- ies2xpehh(scan_pop1 = scan_pop1,
                scan_pop2 = scan_pop2,
                popname1 = pop1,
                popname2 = pop2)

# print stat by SNP
mat = cbind(rsb, xpehh)
mat = mat[,c(1,2,3,4,7,8)]
colnames(mat) = c("chrom", "pos", "rsb", "pval_rsb", "xpehh", "pval_xpehh")
mat = mat[!(is.na(mat[,3]) & is.na(mat[,4])),]
out_snp = paste0(out, "_bySNP.txt")
write.table(mat, out_snp, quote = F, sep = "\t", row.names = F)

# print stat by window
pos_w = slideFunct(mat[,2], 1000, 1000)
rsb_w = slideFunct(mat[,3], 1000, 1000)
xpehh_w = slideFunct(mat[,5], 1000, 1000)

mat_w = data.frame(cbind(rep(as.character(mat$CHR[1]), length(pos_w)), pos_w, rsb_w, xpehh_w))
out_w = paste0(out, "_byWindows.txt")
write.table(mat_w, out_w, quote = F, sep = "\t", row.names = F)

### reformat
rsb = mat[,c("chr", "pos", "rsb", "pval_rsb")]
colnames(rsb) = c("chrom", "pos", "value", "pvalue")
rsb$stat = "rsb"
xpehh = mat[,c("chr", "pos", "rsb", "pval_rsb")]
colnames(xpehh) = c("chrom", "pos", "value", "pvalue")
xpehh$stat = "xpehh"
rehh_df = rbind(rsb, xpehh)

### Perform p-value adjustments and add to data frame
head(rehh_df)
rehh_df$Bonferroni = NA
rehh_df$fdr = NA
rehh_df$pvalue_back_trans = 10**(-rehh_df$pvalue)  # transform back to p-values
for(s in unique(rehh_df$stat)){
  sub = rehh_df[which(rehh_df$stat==s),]
  sub$Bonferroni = p.adjust(sub$pvalue_back_trans, method = "bonferroni")
  sub$fdr = p.adjust(sub$pvalue_back_trans, method = "BH")
  rehh_df$Bonferroni[which(rehh_df$stat==s)] = sub$Bonferroni
  rehh_df$fdr[which(rehh_df$stat==s)] = sub$fdr
}

# print summary stat
total = ddply(rehh_df, .(chrom, stat), plyr::summarize, count=length(chrom))
total$var = "total"
sigpval = ddply(rehh_df[which(rehh_df$pvalue>4),], .(chrom, stat), plyr::summarize, count=length(chrom))
sigpval$var = "pvalue"
sigfdr = ddply(rehh_df[which(rehh_df$fdr<0.5),], .(chrom, stat), plyr::summarize, count=length(chrom))
sigfdr$var = "fdr"
sigbonf = ddply(rehh_df[which(rehh_df$Bonferroni<0.5),], .(chrom, stat), plyr::summarize, count=length(chrom))
sigbonf$var = "Bonferroni"

df = dcast(rbind(total, sigpval, sigfdr, sigbonf),chrom+stat~var, value.var="count")
out_w = paste0(out, "_summaryStat.txt")
write.table(mat_w, out_w, quote = F, sep = "\t", row.names = F)

