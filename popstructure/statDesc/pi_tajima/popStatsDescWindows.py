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
import bcolz
# check which version is installed

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
	return(gt_pop)

def getMafDistribution(gt_pop_f):
	ac = gt_pop_f.count_alleles()
	return(ac)

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

def getInbreedCoeffWindows(gt_pop, pos, windows):
	## inbreed coeff
	inb_coef = []
	for i in range(0,windows.shape[0]):
		if counts[i]>0:
			loc_region = pos.locate_range(windows[i][0], windows[i][1])
			gt_sub = allel.GenotypeArray(gt_pop[loc_region]) # change gt_pop to gt_pop_f
			if gt_sub.shape[0]>0:
				ac = gt_sub.count_alleles()
				flt = ((ac.max_allele()[:] == 1) & ac.is_segregating()[:])
				gt_sub_f = gt_sub.compress(flt, axis=0)
				if gt_sub_f.shape[0]>0:
					# inbreeding coef (F)
					inb_coef.append(np.nanmean(allel.inbreeding_coefficient(gt_sub_f)))
				else:
					inb_coef.append(0)
			else:
				inb_coef.append(0)
		else:
			inb_coef.append(0)
	return(inb_coef)


if __name__=="__main__": # main
	pgrmName = os.path.basename(__file__)

	parser=argparse.ArgumentParser(
		prog=pgrmName,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='''Progam is made to calculate summary statistics of a given population
		using scikit allel.
		''',
		epilog="""...""")

	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	parser.add_argument("--maf", dest="maf", type=float, default=0.00,
		help="Minor allele frequency (default 0.00)")

	parser.add_argument("-w", "--windows", dest="windows", type=int, default=20000,
		help="Window size in bp (default 20kb)")

	parser.add_argument("-a", "--access", dest="access", metavar="FILE", required=False,
		help="Accessibility in h5 format")

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

	statWindowsGenome = ""
	acGenome = ""
	bed = readBedFile(args.bed)
	data = []

	for chrom in bed:
		callset = zarr.open_group(args.snp)
		gt_pop = selectInd(callset, chrom, args.keep)
		pos_path = '/'+chrom+'/variants/POS'
		pos = allel.SortedIndex(callset[pos_path])
		wsize = 20000
		ac = gt_pop.count_alleles()

		# accessibility
		accessibility = h5py.File(args.access, mode='r')
		is_accessible = accessibility[chrom]['is_accessible'][:]

		for feat in bed[chrom]:

			# pi nucleotide diversity
			pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size=20000, start=bed[chrom][feat][0], stop=bed[chrom][feat][1], is_accessible=is_accessible)

			# tajimaD
			D, windows, counts = allel.windowed_tajima_d(pos, ac, size=20000, start=bed[chrom][feat][0], stop=bed[chrom][feat][1])

			# inbcoef
			inb_coef = getInbreedCoeffWindows(gt_pop, pos, windows) # by windows
			allel.inbreeding_coefficient(gt_sub_f)

			# save to df
			start = windows[:, 0]
			stop = windows[:, 1]
			chrarr = np.array([chrom] * len(start))

			featData = np.stack([chrarr, start, stop, n_bases, counts, pi, D, np.array(inb_coef)]).T

			if len(data)==0:
				data = featData
			else : 
				data = np.vstack((data, featData))

			## here could check if same length for all array
			# df = pd.DataFrame({'chrom':chrarr, 'start':start, 'stop':stop, 'n_bases':n_bases, 'counts':counts, 'pi':pi, 'tajimaD':D, 'inb_coef':np.array(inb_coef)})
			# dfGenome = pd.concat([dfGenome, df], ignore_index=True)
	## Print dataframe

	df = pd.DataFrame(data, columns=['chrom', 'start', 'stop', 'n_bases', 'counts', 'pi', 'tajimaD', 'inb_coef'])
	txt_fn = args.out+".statWindows.txt"
	df.to_csv(txt_fn, sep='\t', index=False,  float_format='%.8f')
