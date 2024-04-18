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
import bcolz
import pandas as pd
import zarr
import allel; print('scikit-allel', allel.__version__, file=sys.stderr)
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
				bed[line[0]][k] = [int(line[1]), int(line[2])]
				k += 1
	return(bed)

def get_ac(snp, bed, keep, keep2):
	bed = readBedFile(bed)
	data = []

	for chrom in bed:
		callset = zarr.open_group(snp)
		pos_path = '/'+chrom+'/variants/POS'
		pos = allel.SortedIndex(callset[pos_path])

		### create dic of samples index
		samples_path = '/'+chrom+'/samples'
		panel = pd.read_csv(keep, sep='\t', header=None)
		samples_list = list(callset[samples_path])
		samples_callset_index = [samples_list.index(s) for s in panel[0]]
	
		if keep2 is None:
			subpops = {"k1": list(samples_callset_index)}
		elif keep2 is not None:
			panel = pd.read_csv(keep2, sep='\t', header=None)
			samples_callset_index2 = [samples_list.index(s) for s in panel[0]]
			subpops = {"k1": list(samples_callset_index), "k2": list(samples_callset_index2)}

		### get GT
		gt_path = '/'+chrom+'/calldata/GT'
		gt_pop = allel.GenotypeChunkedArray(callset[gt_path])

		for feat in bed[chrom]:
			### get sub GT
			loc_region = pos.locate_range(bed[chrom][feat][0], bed[chrom][feat][1])
			gt_sub = allel.GenotypeArray(gt_pop[loc_region]) # change gt_pop to gt_pop_f

			print(gt_sub.shape)
			miss = gt_sub.is_missing()[:]
			gt_sub_miss = gt_sub[miss]
			print(gt_sub_miss)
			print(gt_sub_miss.shape)
			exit()

			### get allel count
			ac_subpops = gt_sub.count_alleles_subpops(subpops, max_allele=1)

			if len(data)==0:
				data = ac_subpops
			else :
				ds = [data, ac_subpops]
				dtmp = {}
				for key in data.keys():
				  dtmp[key] = np.concatenate(list(dtmp[key] for dtmp in ds))
				  data = dtmp
	return data

	## old
	# for chrom in bed:
	# 	callset = zarr.open_group(snp)
	# 	pos_path = '/'+chrom+'/variants/POS'
	# 	pos = allel.SortedIndex(callset[pos_path])

	# 	gt_pop = selectInd(callset, chrom, keep)

	# 	for feat in bed[chrom]:
	# 		loc_region = pos.locate_range(bed[chrom][feat][0], bed[chrom][feat][1])
	# 		gt_sub = allel.GenotypeArray(gt_pop[loc_region]) # change gt_pop to gt_pop_f

	# 		# filter variant if needed
	# 		ac = gt_sub.count_alleles()
	# 		flt = ( (1-ac.to_frequencies()[:].max(axis=1) > maf) & (ac.max_allele()[:] == 1) & ac.is_segregating()[:] )
	# 		gt_sub_f = gt_sub.compress(flt, axis=0)
	# 		ac = gt_sub_f.count_alleles()

	# 		if len(data)==0:
	# 			data = ac
	# 		else : 
	# 			data = np.vstack((data, ac))
	# return data


if __name__=="__main__": # main
	pgrmName = os.path.basename(__file__)

	parser=argparse.ArgumentParser(
		prog=pgrmName,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='''Progam will get the SFS site frequency spectrum for a given population.
		''',
		epilog="""...""")

	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	parser.add_argument("--scaled", dest="scaled", action='store_true',
		help="Activate scaled method.")

	parser.add_argument("--unfolded", dest="unfolded", action='store_true',
		help="Activate unfolded method.")

	parser.add_argument("-k2", "--keep2", dest="keep2", metavar="FILE",
		help="Second list of individual to keep to create a join sfs")

	requiredNamed = parser.add_argument_group('required arguments')

	requiredNamed.add_argument("-z", "--zarr", dest="snp", metavar="FILE", required=True,
		help="Input file zarr")

	requiredNamed.add_argument("-o", "--out", dest="out", metavar="STR", required=True,
		help="Output file prefix.")

	requiredNamed.add_argument("-k", "--keep", dest="keep", metavar="FILE", required=True,
		help="List of individual to keep")

	requiredNamed.add_argument("-b", "--bed", dest="bed", metavar="FILE", required=True,
		help="Sorted bed file")

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

	if args.keep2:
		if path.exists(args.keep2):
			print('--> reading keep2 repository:', args.keep2, file=sys.stderr)
		else:
			print('ERROR: --keep2 file doesn\'t exist.\n', file=sys.stderr)
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

	dic_ac = get_ac(args.snp, args.bed, args.keep, args.keep2)

	if "k2" not in dic_ac:
		if args.unfolded==True:
			if args.scaled==True:
				sfs = np.round(allel.sfs_scaled(dic_ac["k1"][:,1]), decimals=4)
			elif args.scaled==False:
				sfs = allel.sfs(dic_ac["k1"][:,1])
		elif args.unfolded==False:
			if args.scaled==True:
				sfs = np.round(allel.sfs_folded_scaled(dic_ac["k1"]), decimals=4)
			elif args.scaled==False:
				sfs = np.round(allel.sfs_folded(dic_ac["k1"]), decimals=4)

		x = np.arange(sfs.shape[0])
		x_norm = np.round(x / dic_ac["k1"].sum(axis=1).max(), decimals=4)
		sfs = np.vstack((x, x_norm, sfs)).T

		df = pd.DataFrame(sfs, columns=['x', 'x_norm', 'y'])
		txt_fn = args.out+".sfs.txt"
		df.to_csv(txt_fn, sep='\t', index=False,  float_format='%.8f')
	else:
		### print join sfs
		if args.scaled==True:
			jsfs = allel.joint_sfs_scaled(dic_ac['k1'][:, 1], dic_ac['k2'][:, 1])
			print("here")
		elif args.scaled==False:
			jsfs = allel.joint_sfs(dic_ac['k1'][:, 1], dic_ac['k2'][:, 1])

		df = pd.DataFrame(jsfs)
		txt_fn = args.out+".jsfs.txt"
		df.to_csv(txt_fn, sep='\t', index=False,  float_format='%.8f', header=False)
		x = np.arange(jsfs.shape[0])
		x_norm = np.round(x / dic_ac["k1"].sum(axis=1).max(), decimals=4)
		y = np.arange(jsfs.shape[1])
		y_norm = np.round(y / dic_ac["k2"].sum(axis=1).max(), decimals=4)


	# path_sfs = args.out+".sfs_folded_scaled.txt"
	# np.savetxt(path_sfs, sfs, newline="\n", delimiter="\t", fmt='%s', header="minor_allele_frequency\tscaled_site_frequency")

