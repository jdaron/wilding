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
# import bcolz

if __name__=="__main__": # main
	pgrmName = os.path.basename(__file__)

	parser=argparse.ArgumentParser(
		prog=pgrmName,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='''Progam is made to get the SNP density for 3 different populations and output private and shared SNPs.''',
		epilog="""...""")

	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	parser.add_argument("--maf", dest="maf", type=float, default=0.00, 
		help="Minor allele frequency of the union of both pop (default 0.00)")

	parser.add_argument("-w", "--windows", dest="windows", type=int, metavar="int", default=1000,
		help="Window size in bp")

	parser.add_argument("-s", "--step", dest="step", type=int, metavar="int", default=None,
		help="Window step size in bp")

	parser.add_argument("-a", "--access", dest="access", metavar="FILE", required=False,
		help="Accessibility in h5 format")

	requiredNamed = parser.add_argument_group('required arguments')

	requiredNamed.add_argument("-z", "--zarr", dest="zarr", metavar="FILE", required=True,
		help="Input file in zarr format.")

	requiredNamed.add_argument("--pop1", dest="pop1", metavar="FILE", required=True,
		help="List of individuals from pop1")

	requiredNamed.add_argument("--pop2", dest="pop2", metavar="FILE", required=True,
		help="List of individuals from pop2")

	requiredNamed.add_argument("--pop3", dest="pop3", metavar="FILE", required=True,
		help="List of individuals from pop3")

	requiredNamed.add_argument("-o", "--out", dest="out", metavar="STR", required=True,
		help="Output file prefix.")

	# requiredNamed.add_argument("-b", "--bed", dest="bed", metavar="FILE", required=True,
	# 	help="Sorted bed file of position to look for")

	## if no args then return help message
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)

	args = parser.parse_args()

	### Check option
	### End of check option

	### accessibility
	accessibility = h5py.File(args.access, mode='r')

	### read pop file
	pop1 = pd.read_csv(args.pop1, sep='\t', header=None)
	pop2 = pd.read_csv(args.pop2, sep='\t', header=None)
	pop3 = pd.read_csv(args.pop3, sep='\t', header=None)

	callset = zarr.open_group(args.zarr, mode='r')

	genomeDensity = []
	directory_contents = os.listdir(args.zarr)
	for chrom in directory_contents:
		if os.path.isdir(os.path.join(args.zarr,chrom)):	

			# read pop id
			samples_path = '/'+chrom+'/samples'
			samples_list = list(callset[samples_path])

			samples_callset_index1 = [samples_list.index(s) for s in pop1[0]]
			samples_callset_index2 = [samples_list.index(s) for s in pop2[0]]
			samples_callset_index3 = [samples_list.index(s) for s in pop3[0]]

			subpops = {"pop1": list(samples_callset_index1), "pop2": list(samples_callset_index2), "pop3": list(samples_callset_index3)}

			# get GT and POS
			gt_path = '/'+chrom+'/calldata/GT'
			pos_path = '/'+chrom+'/variants/POS'
			gt = allel.GenotypeChunkedArray(callset[gt_path])
			pos = allel.SortedIndex(callset[pos_path])

			# filter out variants that aren’t segregating in the union of our two populations.
			# Let’s also filter out multiallelic variants for simplicity.
			acs = gt.count_alleles_subpops(subpops)
			acu = allel.AlleleCountsArray(acs["pop1"][:] + acs["pop2"][:] + acs["pop3"][:])
			flt = acu.is_segregating() & (acu.max_allele() == 1) & (1-acu.to_frequencies().max(axis=1) > args.maf)
			pos_flt = pos.compress(flt)
			gt_flt = gt.compress(flt, axis=0)
			print('INFO: keep segregating sites and biallelic sites then retaining:', np.count_nonzero(flt), 'SNPs' , file=sys.stderr)

			### create windows
			is_accessible = accessibility[chrom]['is_accessible'][:]
			pos_accessible = np.nonzero(is_accessible)[0] + 1
			values = np.isnan(pos_accessible)
			x, windows, countsAccessPos = allel.windowed_statistic(pos_accessible, values, statistic=np.count_nonzero, size=args.windows, start=1, stop=None, step=args.step)
			windows = windows[countsAccessPos>0] # remove windows with 0 SNP
			countsAccessPos = countsAccessPos[countsAccessPos>0] # remove windows with 0 SNP

			# values = np.isnan(pos_flt)
			# x, windows, countsSNPs = allel.windowed_statistic(pos_flt, values, statistic=np.count_nonzero, size=args.windows, start=1, stop=None, step=args.step)
			# windows = windows[countsSNPs>0] # remove windows with 0 SNP
			# countsSNPs = countsSNPs[countsSNPs>0] # remove windows with 0 SNP

			# find singleton
			acs = gt_flt.count_alleles_subpops(subpops)
			acu = allel.AlleleCountsArray(acs["pop1"][:] + acs["pop2"][:] + acs["pop3"][:])
			singleton_mask = (acu[:, :2].min(axis=1) == 1)
			singletonw, windows, counts = allel.windowed_statistic(pos_flt, singleton_mask, statistic=np.count_nonzero, windows=windows)
			flt = (acu[:, :2].min(axis=1) > 1)
			gt_flt = gt_flt.compress(flt, axis=0)
			pos_flt = pos_flt.compress(flt)

			# allele counts
			ac1 = gt_flt.count_alleles(subpop=subpops["pop1"])
			ac2 = gt_flt.count_alleles(subpop=subpops["pop2"])
			ac3 = gt_flt.count_alleles(subpop=subpops["pop3"])

			loc_private_alleles = allel.locate_private_alleles(ac1, ac2, ac3)

			### get shared variants
			loc_private_variants = np.any(loc_private_alleles, axis=1)
			shared = np.array([not elem for elem in loc_private_variants])
			pos_priv_variants = pos_flt.compress(loc_private_variants)

			### get private variant pop1
			ac1 = np.array(ac1)
			priv_mask_pop1 = ac1[loc_private_alleles]>0

			### get private variant pop2
			ac2 = np.array(ac2)
			priv_mask_pop2 = ac2[loc_private_alleles]>0

			### get private variant pop1
			ac3 = np.array(ac3)
			priv_mask_pop3 = ac3[loc_private_alleles]>0

			priv_pop1w, windows, counts = allel.windowed_statistic(pos_priv_variants, priv_mask_pop1, statistic=np.count_nonzero, windows=windows)
			priv_pop1w[np.isnan(priv_pop1w)] = 0
			priv_pop2w, windows, counts = allel.windowed_statistic(pos_priv_variants, priv_mask_pop2, statistic=np.count_nonzero, windows=windows)
			priv_pop2w[np.isnan(priv_pop2w)] = 0
			priv_pop3w, windows, counts = allel.windowed_statistic(pos_priv_variants, priv_mask_pop3, statistic=np.count_nonzero, windows=windows)
			priv_pop3w[np.isnan(priv_pop3w)] = 0

			sharedw, windows, counts = allel.windowed_statistic(pos_flt, shared, statistic=np.count_nonzero, windows=windows)
			sharedw[np.isnan(sharedw)] = 0

			chrDensity = np.stack([np.array([chrom] * len(counts)), windows[:,0], windows[:,1], countsAccessPos, singletonw, priv_pop1w, priv_pop2w, priv_pop3w, sharedw]).T

			if len(genomeDensity)==0:
				genomeDensity = chrDensity
			else : 
				genomeDensity = np.vstack((genomeDensity, chrDensity))

	df = pd.DataFrame(genomeDensity, columns=['chrom', 'start', 'stop', 'countsAccessPos', 'singleton', 'priv_1', 'priv_2', 'priv_3', 'shared'])
	txt_fn = args.out+".snpDensity.txt"
	df.to_csv(txt_fn, sep='\t', index=False,  float_format='%.8f')

