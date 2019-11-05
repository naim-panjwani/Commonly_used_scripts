#!/usr/local/bin/python

import numpy as np
import pandas as pd
import gzip
import subprocess
import scipy.stats as stats
import argparse
import os

parser = argparse.ArgumentParser(description='Subset gct file to the desired columns given')
parser.add_argument('expression_gct', help='GCT file with expression')
parser.add_argument('expression_reads', help='GCT file with read counts')
parser.add_argument('columns_to_extract', help='File with the column names to extract; make sure to also include "Name" and "Description"')
parser.add_argument('outfilename', help="Desired output filename")
args = parser.parse_args()

#expression_gct="../ExpressionFiles/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1//GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
#columns_to_extract="output//02-RNAseq_Ileum_subset/gct_columns_to_extract_181126"
#outfilename="output//02-RNAseq_Ileum_subset/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_Ileum.gz"

print('expression_gct '+args.expression_gct)
print('expression_reads '+args.expression_reads)
print('columns_to_extract '+args.columns_to_extract)
print('outfilename '+args.outfilename)

print("Loading expression file")
gct = pd.read_csv(args.expression_gct, sep="\t", compression="gzip", skiprows=2)
gct_reads = pd.read_csv(args.expression_reads, sep="\t", compression="gzip", skiprows=2)
columns_to_subset = pd.read_csv(args.columns_to_extract, sep="\t", header=None)
columns_to_subset = columns_to_subset.iloc[:,0]

#print("Columns to subset:")
#print(columns_to_subset.head())

#print("gct columns")
#print(gct.columns)

subsetted_gct = gct[[i for i in gct.columns if i in list(columns_to_subset)]]
subsetted_reads = gct_reads[[i for i in gct_reads.columns if i in list(columns_to_subset)]]
subsetted_gct.to_csv(args.outfilename, sep="\t", index=False)
subsetted_reads.to_csv(args.outfilename+'_reads', sep="\t", index=False)
subprocess.check_call('bgzip -f '+args.outfilename, shell=True, executable='/bin/bash')
#subprocess.check_call('tabix -f '+args.outfilename+'.gz', shell=True, executable='/bin/bash')
subprocess.check_call('bgzip -f '+args.outfilename+'_reads', shell=True, executable='/bin/bash')
#subprocess.check_call('tabix -f '+args.outfilename+'_reads.gz', shell=True, executable='/bin/bash')
num_samples=subsetted_gct.shape[1]
print(f'Subsetting done. {num_samples} were extracted')

