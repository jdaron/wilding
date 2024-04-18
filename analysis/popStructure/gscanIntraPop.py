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

# def getInbreedCoeffWindows(gt_pop, pos, windows, counts):
# 	## inbreed coeff
# 	inb_coef = []
# 	for i in range(0,windows.shape[0]):
# 		if counts[i]>0:
# 			loc_region = pos.locate_range(windows[i][0], windows[i][1])
# 			gt_sub = allel.GenotypeArray(gt_pop[loc_region]) # change gt_pop to gt_pop_f
# 			if gt_sub.shape[0]>0:
# 				ac = gt_sub.count_alleles()
# 				flt = ((ac.max_allele()[:] == 1) & ac.is_segregating()[:])
# 				gt_sub_f = gt_sub.compress(flt, axis=0)
# 				if gt_sub_f.shape[0]>0:
# 					# inbreeding coef (F)
# 					inb_coef.append(np.nanmean(allel.inbreeding_coefficient(gt_sub_f)))
# 				else:
# 					inb_coef.append(0)
# 			else:
# 				inb_coef.append(0)
# 		else:
# 			inb_coef.append(0)
# 	return(inb_coef)

# def getGarud_hWindows(gt_pop, pos, windows, counts):
# 	## inbreed coeff
# 	garud_h = []
# 	for i in range(0,windows.shape[0]):
# 		if counts[i]>0:
# 			loc_region = pos.locate_range(windows[i][0], windows[i][1])
# 			gt_sub = allel.GenotypeArray(gt_pop[loc_region]) # change gt_pop to gt_pop_f
# 			if gt_sub.shape[0]>0: # no variants
# 				ac = gt_sub.count_alleles()
# 				flt = ((ac.max_allele()[:] == 1) & ac.is_segregating()[:])
# 				gt_sub_f = gt_sub.compress(flt, axis=0)
# 				if gt_sub_f.shape[0]>0: # no segregating variants
# 					h = gt_sub_f.to_haplotypes()
# 					h1, h12, h123, h1_2 = allel.garud_h(h)

# 					garud_h.append([h1, h12, h123, h1_2])
# 				else:
# 					garud_h.append([0,0,0,0])
# 			else:
# 				garud_h.append([0,0,0,0])
# 		else:
# 			garud_h.append([0,0,0,0])
# 	return(garud_h)

if __name__=="__main__": # main
	pgrmName = os.path.basename(__file__)

	parser=argparse.ArgumentParser(
		prog=pgrmName,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='''Progam is made to calculate summary statistics of a given population
		using scikit allel. Some stat need phased data to be correctly calculated (H12, ihs).
		''',
		epilog="""...""")

	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	parser.add_argument("--maf", dest="maf", type=float, metavar="flot", default=0.00,
		help="Minor allele frequency (default 0.00)")

	parser.add_argument("-w", "--windows", dest="windows", type=int, metavar="INT", default=20000,
		help="Window size in bp (default 20kb)")

	parser.add_argument("-s", "--step", dest="step",  type=int, metavar="INT", default=None,
		help="Window step size in bp or snp (default None)")

	parser.add_argument("-wt", "--windType", dest="windType", type=str, metavar="STR", default="bp",
		help="Type of window to make, either bp or snp [bp|snp].")

	parser.add_argument("-a", "--access", dest="access", metavar="FILE", help="Accessibility in h5 format")

	requiredNamed = parser.add_argument_group('required arguments')

	requiredNamed.add_argument("-z", "--zarr", dest="snp", metavar="FILE", required=True,
		help="Input zarr folder.")

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
	if args.bed:
		if path.exists(args.bed):
			print('--> reading bed file:', args.bed, file=sys.stderr)
		else:
			print('ERROR: --bed file doesn\'t exist.\n', file=sys.stderr)
			parser.parse_args(['-h'])
			exit()
	#### end check option

	bed = readBedFile(args.bed)
	statWindowsGenome = ""
	acGenome = ""
	wsize = args.windows
	data = []
	dataInb = []

	for chrom in bed:
		callset = zarr.open_group(args.snp)
		gt_pop = selectInd(callset, chrom, args.keep)
		pos_path = '/'+chrom+'/variants/POS'
		pos_all = allel.SortedIndex(callset[pos_path])

		# filter invariant sites
		ac = gt_pop.count_alleles()
		flt = ac.is_segregating() & ac.max_allele() == 1
		pos = pos_all.compress(flt)
		gt_pop = gt_pop.compress(flt, axis=0)
		ac = gt_pop.count_alleles()

		# accessibility
		accessibility = ""
		is_accessible = ""
		if args.access:
			accessibility = h5py.File(args.access, mode='r')
			is_accessible = accessibility[chrom]['is_accessible'][:]

		for feat in bed[chrom]:
			windows = []
			if args.windType == "bp":
				mask = np.isnan(pos)
				x, windows, counts = allel.windowed_statistic(pos, mask, statistic=np.count_nonzero, size=wsize, start=bed[chrom][feat][0], stop=bed[chrom][feat][1], step=args.step)
				windows = windows[counts>0] # remove windows with 0 SNP
			elif args.windType == "snp":
				st = allel.moving_statistic(pos, statistic=lambda v: v[0], size=wsize, start=0, stop=None, step=args.step)
				sp = allel.moving_statistic(pos, statistic=lambda v: v[-1], size=wsize, start=0, stop=None, step=args.step)
				idx = (st>bed[chrom][feat][0]) & (st<bed[chrom][feat][1])
				st = st[idx]
				sp = sp[idx]
				windows =  np.vstack((st, sp)).T

			if args.access:
				# pi nucleotide diversity
				pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, windows=windows, is_accessible=is_accessible)
				# theta watterson
				theta_hat_w, windows, n_bases, counts = allel.windowed_watterson_theta(pos, ac, windows=windows, is_accessible=is_accessible)
			else:
				pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size=wsize, start=bed[chrom][feat][0], stop=bed[chrom][feat][1])
				theta_hat_w, windows, n_bases, counts = allel.windowed_watterson_theta(pos, ac, windows=windows)

			# tajimaD
			# D, windows, counts = allel.windowed_tajima_d(pos, ac, size=wsize, start=bed[chrom][feat][0], stop=bed[chrom][feat][1])
			D, windows, counts = allel.windowed_tajima_d(pos, ac, windows=windows)

			# inbcoef
			inb_coef, windows, counts = allel.windowed_statistic(pos, gt_pop, statistic=allel.inbreeding_coefficient, windows=windows)
			inb_coef_avg = []
			for i in range(0,len(inb_coef)):
				inb_coef_avg.append(np.nanmean(inb_coef[i]))
			inb_coef_avg = np.array(inb_coef_avg)

			# inbcoef by SNP
			inb_coef = allel.inbreeding_coefficient(gt_pop)
			chrarr = np.array([chrom] * len(pos))
			inb_coef_gs = np.stack([chrarr, pos, np.round_(inb_coef,8)]).T
			if len(dataInb)==0:
				dataInb = inb_coef_gs
			else : 
				dataInb = np.vstack((dataInb, inb_coef_gs))

			# garud_h
			h = gt_pop.to_haplotypes()
			garud_h, windows, counts = allel.windowed_statistic(pos, h, statistic=allel.garud_h, windows=windows)

			# iHS
			# ihs = allel.ihs(h, pos, map_pos=None, min_ehh=0.05, min_maf=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=is_accessible, use_threads=True)
			# mask = np.isnan(ihs)
			# ihs_avg, windows, counts = allel.windowed_statistic(pos[~mask], ihs[~mask], statistic=np.mean, windows=windows)

			# save to df
			start = windows[:, 0]
			stop = windows[:, 1]
			chrarr = np.array([chrom] * len(start))

			featData = np.stack([chrarr, start, stop, n_bases, counts, np.round_(pi, 8), np.round_(theta_hat_w, 8), np.round_(D, 8), np.round_(inb_coef_avg,8)]).T
			featData = np.hstack((featData, np.round_(np.array(garud_h),8)))

			if len(data)==0:
				data = featData
			else : 
				data = np.vstack((data, featData))

			## here could check if same length for all array
			# df = pd.DataFrame({'chrom':chrarr, 'start':start, 'stop':stop, 'n_bases':n_bases, 'counts':counts, 'pi':pi, 'tajimaD':D, 'inb_coef':np.array(inb_coef)})
			# dfGenome = pd.concat([dfGenome, df], ignore_index=True)
	## Print dataframe

	df = pd.DataFrame(data, columns=['chrom', 'start', 'stop', 'n_bases', 'counts', 'pi', 'theta_wat', 'tajimaD', 'inb_coef', 'h1', 'h12', 'h123', 'h2.1'])
	txt_fn = args.out+"."+str(args.windows)+"w."+str(args.step)+"s."+"statWindows.txt"
	df.to_csv(txt_fn, sep='\t', index=False,  float_format='%.8f')

	df = pd.DataFrame(dataInb, columns=['chrom', 'pos', 'inbreeding_coefficient'])
	txt_fn = args.out+".statSNPs.txt"
	df.to_csv(txt_fn, sep='\t', index=False,  float_format='%.8f')
