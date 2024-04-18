#!/usr/bin/env python3

import os
from os import path
import re
import argparse
import sys
import random
random.seed(42)
import time
import numpy as np
np.random.seed(42)
import h5py
import matplotlib.pyplot as plt
# %matplotlib inline
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import pandas as pd
import zarr
import allel; print('scikit-allel', allel.__version__, file=sys.stderr)
# import bcolz
# check which version is installed

def getROH(gt_pop_f, callset, chrom, pos):
	# 4. FROH/ROHcount (<100kb)
	genome_size = pos[len(pos)-1]
	samples_path = '/'+chrom+'/samples'

	total_df_roh = pd.DataFrame()
	for i in range(0,gt_pop_f.shape[1]):
		gv = gt_pop_f[:,i]
		df_roh, froh = allel.roh_mhmm(gv, pos, contig_size=genome_size)
		df_roh = df_roh.assign(id=callset[samples_path][i])
		total_df_roh = pd.concat([total_df_roh, df_roh], axis=0)
	return(total_df_roh)

def readBedFile(bedFilePath):
	bed = {}
	k = 0
	if bedFilePath:
		with open(bedFilePath, 'r') as bedFile:
			for row in bedFile:
				line = row.strip().split("\t")
				if line[0] not in bed:
					bed[line[0]] = {}
					k = 0
				# bed[line[0]][k] = []
				bed[line[0]][k] = [int(line[1]), int(line[2])]
				k += 1
	return(bed)

def selectInd(callset, chrom, keep):
	gt_path = '/'+chrom+'/calldata/GT'
	samples_path = '/'+chrom+'/samples'

	gt = allel.GenotypeChunkedArray(callset[gt_path])

	# coorespondance of samples order between snp and meta file
	panel = pd.read_csv(keep, sep='\t', header=None)
	samples_list = list(callset[samples_path])
	samples_callset_index = [samples_list.index(s) for s in panel[0]]

	# create a genotype array of the selected samples
	gt_pop = gt.take(samples_callset_index, axis=1)

	# d = {'keep':list(panel.iloc[:, 0]),'index':samples_callset_index}
	# df = pd.DataFrame(d)
	return(gt_pop, samples_callset_index)

if __name__=="__main__": # main
	pgrmName = os.path.basename(__file__)

	parser=argparse.ArgumentParser(
		prog=pgrmName,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='''Progam will get the ROH (runs of homozygosity).
		using scikit allel.
		''',
		epilog="""...""")

	parser.add_argument("-a", "--access", dest="access", metavar="FILE",
		help="Genome accessibility in h5 format")

	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	parser.add_argument("-m", "--method", dest="method", type=str, default="mhmm",
		help="Method use to detect ROH [mhmm|poisson].")

	requiredNamed = parser.add_argument_group('required arguments')

	requiredNamed.add_argument("--snp", dest="snp", metavar="FILE", required=True,
		help="Input file, VCF (.vcf or .vcf.gz) or zarr (.zarr)")

	requiredNamed.add_argument("-o", "--out", dest="out", metavar="STR", required=True,
		help="Output file prefix.")

	requiredNamed.add_argument("--keep", dest="keep", metavar="FILE", required=True,
		help="List of individual to keep")

	requiredNamed.add_argument("-b", "--bed", dest="bed", metavar="FILE", required=True,
		help="Sorted bed file of position to look for")

	## if no args then return help message
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)

	args = parser.parse_args()

	### Check option
	if path.exists(args.snp):
		print('--> reading zarr repository:', args.snp, file=sys.stderr)
	else:
		print('ERROR: --snp zarr file doesn\'t exist.\n', file=sys.stderr)
		parser.parse_args(['-h'])
		exit()
	if path.exists(args.keep):
			print('--> reading keep file:', args.keep, file=sys.stderr)
	else:
		print('ERROR: --keep keep file doesn\'t exist.\n', file=sys.stderr)
		parser.parse_args(['-h'])
		exit()
	if args.access:
		if path.exists(args.access):
			print('--> reading accessibility file:', args.access, file=sys.stderr)
		else:
			print('ERROR: --access accessibility file doesn\'t exist.\n', file=sys.stderr)
			parser.parse_args(['-h'])
			exit()
	if args.bed:
		if path.exists(args.bed):
			print('--> reading bed file:', args.bed, file=sys.stderr)
		else:
			print('ERROR: --bed file doesn\'t exist.\n', file=sys.stderr)
			parser.parse_args(['-h'])
			exit()
	#### end check option

	panel = pd.read_csv(args.keep, sep='\t', header=None)
	panel = list(panel[0])

	bed = readBedFile(args.bed)
	total_df_roh = pd.DataFrame()
	total_df_froh = pd.DataFrame()

	for chrom in bed:
		callset = zarr.open_group(args.snp)
		gt_pop, samples_callset_index = selectInd(callset, chrom, args.keep)
		pos_path = '/'+chrom+'/variants/POS'
		samples_path = '/'+chrom+'/samples'
		pos = allel.SortedIndex(callset[pos_path])

                # accessibility
		if args.access:
			accessibility = h5py.File(args.access, mode='r')
			is_accessible = accessibility[chrom]['is_accessible'][:]

		for feat in bed[chrom]:

			# roh
			loc_region = pos.locate_range(bed[chrom][feat][0], bed[chrom][feat][1])
			gt_sub = allel.GenotypeArray(gt_pop[loc_region]) # change gt_pop to gt_pop_f
			pos_sub = pos[loc_region]
			genome_size = pos_sub[len(pos_sub)-1]

			for i in range(0,gt_sub.shape[1]):
				gv = gt_sub[:,i]
				# df_roh, froh = allel.roh_mhmm(gv, pos_sub, contig_size=genome_size)
				if args.method=="mhmm":
					df_roh, froh = allel.roh_mhmm(gv, pos_sub, contig_size=genome_size, is_accessible=is_accessible)

				elif args.method=="poisson":
					df_roh, froh = allel.roh_poissonhmm(gv, pos_sub, contig_size=genome_size, is_accessible=is_accessible)

				df_roh = df_roh.assign(id=callset[samples_path][samples_callset_index[i]], chrom=chrom)
				total_df_roh = pd.concat([total_df_roh, df_roh], axis=0)

				df_froh = pd.DataFrame({'id':callset[samples_path][samples_callset_index[i]], 'chrom':chrom, 'froh':froh}, index=[0])
				total_df_froh = pd.concat([total_df_froh, df_froh], axis=0)

	## Print dataframe
	txt_fn = args.out+".roh.txt"
	total_df_roh.to_csv(txt_fn, sep='\t', index=False,  float_format='%.8f')
	txt_fn = args.out+".froh.txt"
	total_df_froh.to_csv(txt_fn, sep='\t', index=False,  float_format='%.8f')
