#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#### R script to reorder admix output
#===============================================================================
#   Author: Josquin Daron
#
#   File: orderAdmix.R
#   Date: 19-12-2018
#   Version: 0.1
#
#   Usage:
#	 R --vanilla --args $PWD/'admix' $PWD/'ind ordered as in the admix' $PWD/'ind new order' < orderAdmix.R 
#===============================================================================
options(digits=22)

ex = sub('.*\\.(.*)$', '\\1', args[1])
prefix = sub('(.*)\\..*$', '\\1', args[1])

filename = paste(prefix, ex, sep=".")
tbl=read.table(filename, header = F)

ind=read.table(args[2])
colnames(ind) = c("id")
rownames(tbl) = ind$id

ind.o=read.table(args[3])
colnames(ind.o) = c("id")

mat = tbl[match(ind.o$id, rownames(tbl)), ]
mat <- format(mat, scientific = FALSE)

outputFile = paste(prefix, ".ordered.", ex, sep="")
print(outputFile)
write.table(mat, file = outputFile, row.names = FALSE, col.names = FALSE, quote = FALSE)


