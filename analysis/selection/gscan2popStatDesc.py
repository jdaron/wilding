#!/usr/bin/env python3

import os
from os import path
import re
import argparse
import sys
import random
random.seed(42)
import time
import h5py
import numpy as np
np.random.seed(42)
import pandas as pd
import zarr
import allel; print('scikit-allel', allel.__version__, file=sys.stderr)
import bcolz

def jackknife(x, func):
    """Jackknife estimate of the estimator func"""
    n = len(x)
    idx = np.arange(n)
    return np.sum(func(x[idx!=i]) for i in range(n))/float(n)

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
				bed[line[0]][k] = [int(line[1])+1, int(line[2])]
				k += 1
	return(bed)

if __name__=="__main__": # main
	pgrmName = os.path.basename(__file__)

	parser=argparse.ArgumentParser(
		prog=pgrmName,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='''Progam has been created to get inter pop statistics (FST & DXY) between 2 populations.
		''',
		epilog="""...""")

	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	parser.add_argument("--maf", dest="maf", type=float, default=0.00, 
		help="Minor allele frequency of the union of both pop (default 0.00)")

	parser.add_argument("-w", "--windows", dest="windows", type=int, metavar="int", default=1000,
		help="Window size in nb. of SNPbps")

	parser.add_argument("-s", "--step", dest="step",  type=int, metavar="INT", default=None,
		help="Window step size in nb. of snp (default None)")

	parser.add_argument("-a", "--access", dest="access", metavar="FILE",
		help="Genome accessibility in h5 format")

	# parser.add_argument("-m", "--map", dest="recom_map", metavar="FILE",
		# help="Recombination maps")


	requiredNamed = parser.add_argument_group('required arguments')

	requiredNamed.add_argument("-z", "--zarr", dest="zarr", metavar="FILE", required=True,
		help="Input file in zarr format.")

	requiredNamed.add_argument("-o", "--out", dest="out", metavar="STR", required=True,
		help="Output file prefix.")

	requiredNamed.add_argument("--pop1", dest="pop1", metavar="FILE", required=True,
		help="List of individuals from pop1")

	requiredNamed.add_argument("--pop2", dest="pop2", metavar="FILE", required=True,
		help="List of individuals from pop2")

	requiredNamed.add_argument("-b", "--bed", dest="bed", metavar="FILE", required=True,
		help="Sorted bed file of position to look for")

	# requiredNamed.add_argument("-b", "--bed", dest="bed", metavar="FILE", required=True,
	# 	help="Sorted bed file of position to look for")

	## if no args then return help message
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)

	args = parser.parse_args()

	### Check option
	if path.exists(args.zarr):
		print('--> reading zarr repository:', args.zarr, file=sys.stderr)
	else:
		print('ERROR: --snp zarr file doesn\'t exist.\n', file=sys.stderr)
		parser.parse_args(['-h'])
		exit()
	if path.exists(args.pop1):
			print('--> reading --pop1 file:', args.pop1, file=sys.stderr)
	else:
		print('ERROR: --pop1 keep file doesn\'t exist.\n', file=sys.stderr)
		parser.parse_args(['-h'])
		exit()
	if path.exists(args.pop2):
			print('--> reading --pop2 file:', args.pop2, file=sys.stderr)
	else:
		print('ERROR: --pop2 keep file doesn\'t exist.\n', file=sys.stderr)
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
	### End of check option

	# read bed file
	bed = readBedFile(args.bed)

	# accessibility
	accessibility = []
	if args.access:
		accessibility = h5py.File(args.access, mode='r')

	### create callset
	callset = zarr.open_group(args.zarr, mode='r')

	### get chromosome id in zarr dir
	# directory_contents = os.listdir(args.zarr)

	output_fst_dxy = []
	output_xpehh = []

	for chrom in bed:

		# get Genotype array
		gt_path = '/'+chrom+'/calldata/GT'
		gt = allel.GenotypeChunkedArray(callset[gt_path])
		pos_path = '/'+chrom+'/variants/POS'
		pos_all = allel.SortedIndex(callset[pos_path])

		# read pop id
		samples_path = '/'+chrom+'/samples'
		samples_list = list(callset[samples_path])

		pop1 = pd.read_csv(args.pop1, sep='\t', header=None)
		pop2 = pd.read_csv(args.pop2, sep='\t', header=None)

		samples_callset_index1 = [samples_list.index(s) for s in pop1[0]]
		samples_callset_index2 = [samples_list.index(s) for s in pop2[0]]
		subpops = {"pop1": list(samples_callset_index1), "pop2": list(samples_callset_index2)}
		pop1_idx = subpops["pop1"]
		pop2_idx = subpops["pop2"]

		# accessibility
		if args.access:
			accessibility = h5py.File(args.access, mode='r')
			is_accessible = accessibility[chrom]['is_accessible'][:]

		for feat in bed[chrom]:

			# get Genotype arrary for the feature
			# loc_region = pos_all.locate_range(bed[chrom][feat][0], bed[chrom][feat][1])
			# gt_feat = allel.GenotypeArray(gt[loc_region])
			flt = (pos_all>bed[chrom][feat][0]) & (pos_all<bed[chrom][feat][1])
			gt_feat = gt.compress(flt, axis=0)
			pos_feat = pos_all.compress(flt)

			# allele counts
			acs = gt_feat.count_alleles_subpops(subpops)

			# filter out variants that aren’t segregating in the union of our two populations.
			# Let’s also filter out multiallelic variants for simplicity.

			acu = allel.AlleleCountsArray(acs["pop1"][:] + acs["pop2"][:])
			flt = acu.is_segregating() & (acu.max_allele() == 1 & (1-acu.to_frequencies().max(axis=1) > args.maf) )
			print('INFO: keep segregating sites and biallelic sites then retaining:', np.count_nonzero(flt), 'SNPs' , file=sys.stderr)

			pos = pos_feat.compress(flt)
			ac1 = allel.AlleleCountsArray(acs["pop1"].compress(flt, axis=0)[:, :2])
			ac2 = allel.AlleleCountsArray(acs["pop2"].compress(flt, axis=0)[:, :2])
			gt_feat_flt = gt_feat.compress(flt, axis=0)

			# get haplotype array for each pop
			gt_h1 = gt_feat_flt.take(subpops["pop1"], axis=1)
			h1 = gt_h1.to_haplotypes()
			gt_h2 = gt_feat_flt.take(subpops["pop2"], axis=1)
			h2 = gt_h2.to_haplotypes()

			# get windows pos
			st = allel.moving_statistic(pos, statistic=lambda v: v[0], size=args.windows, start=0, stop=None, step=args.step)
			sp = allel.moving_statistic(pos, statistic=lambda v: v[-1], size=args.windows, start=0, stop=None, step=args.step)
			w_pos = np.vstack((st, sp)).T

			# get fst 
			fst_windows = allel.moving_weir_cockerham_fst(gt_feat_flt, subpops=[pop1_idx, pop2_idx], size=args.windows, step=args.step, max_allele=1)

			# get Dxy
			dxy = []
			if args.access is not None:
				is_accessible = accessibility[chrom]['is_accessible'][:]
				dxy, windows, n_bases, counts = allel.windowed_divergence(pos, ac1, ac2, windows=w_pos, is_accessible=is_accessible)
			else:
				dxy, windows, n_bases, counts = allel.windowed_divergence(pos, ac1, ac2, windows=w_pos)

			# get xp-ehh
			st_xpehh = []
			if args.access:
				is_accessible = accessibility[chrom]['is_accessible'][:]
				un_xpehh = allel.xpehh(h1, h2, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None)
				st_xpehh = allel.standardize(un_xpehh)
				st_xpehh[~np.isnan(st_xpehh)]
			else:
				un_xpehh = allel.xpehh(h1, h2, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None)
				st_xpehh = allel.standardize(un_xpehh)
				st_xpehh[~np.isnan(st_xpehh)]

			st_xpehh_w_mean = allel.moving_statistic(st_xpehh, np.mean, size=args.windows, step=args.step, start=0)
			st_xpehh_w_median = allel.moving_statistic(st_xpehh, statistic=lambda values: np.percentile(values, 50), size=args.windows, step=args.step, start=0)
			st_xpehh_w_std = allel.moving_statistic(st_xpehh, np.std, size=args.windows, step=args.step, start=0)
			st_xpehh_w_max = allel.moving_statistic(st_xpehh, np.max, size=args.windows, step=args.step, start=0)
			st_xpehh_w_min = allel.moving_statistic(st_xpehh, np.min, size=args.windows, step=args.step, start=0)

			# get all stat together
			windowsStat = np.stack([np.array([chrom] * len(st)), st, sp, counts, np.round_(fst_windows, 8), np.round_(dxy, 8), np.round_(st_xpehh_w_mean, 8), np.round_(st_xpehh_w_median, 8), np.round_(st_xpehh_w_std, 8), np.round_(st_xpehh_w_min, 8), np.round_(st_xpehh_w_max, 8)]).T
			xpehh = np.stack([np.array([chrom] * len(pos)), pos, np.round_(st_xpehh, 3)]).T

			if len(output_fst_dxy)==0:
				output_fst_dxy = windowsStat
				output_xpehh = xpehh
			else : 
				output_fst_dxy = np.vstack((output_fst_dxy, windowsStat))
				output_xpehh = np.vstack((output_xpehh, xpehh))

			# get fst avg by jackniff
			# https://people.duke.edu/~ccc14/sta-663/ResamplingAndMonteCarloSimulations.html

			# Jackknife estimate of standard deviation
			print(jackknife(fst_windows, np.std))

			# Jackknife estimate of mean
			print(jackknife(fst_windows, np.mean))
			# fst_avg, se, vb, vj = allel.average_weir_cockerham_fst(gt.flt, subpops=[pop1_idx, pop2_idx], blen=1000, max_allele=1)
			# chrom_cutoff = [chrom, round(fst_avg, 8), round(se, 8)]

			# if len(output_cutoff)==0:
			# 	output_cutoff = chrom_cutoff
			# else : 
			# 	output_cutoff = np.vstack((output_cutoff, chrom_cutoff))

	df = pd.DataFrame(output_fst_dxy, columns=['chrom', 'start', 'stop', 'count', 'fst', 'dxy', 'xpehh_mean', 'st_xpehh_w_median', 'st_xpehh_w_std', 'xpehh_max', 'xpehh_min'])
	txt_fn = args.out+".fst_dxy.scan.txt"
	df.to_csv(txt_fn, sep='\t', index=False,  float_format='%.8f')

	df = pd.DataFrame(output_xpehh, columns=['chrom', 'pos', 'xpehh'])
	txt_fn = args.out+".xpehh.txt"
	df.to_csv(txt_fn, sep='\t', index=False,  float_format='%.8f')

	exit()
