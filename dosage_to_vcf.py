# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 10:08:11 2020

@author: Naim
Description: This script takes a dosage.gz file as formatted by the GWAS2 Consortium imputation 
and converts it to a simple vcf file with rounded genotypes
"""

import pandas as pd
import numpy as np
import argparse
import subprocess, os, sys
import gzip
from tqdm import tqdm
from datetime import datetime


########################################
# HELPER FUNCTIONS
########################################
def getchrom(dosage_filename_):
    with gzip.open(dosage_filename_, 'r') as f_in:
        chrom = f_in.readline().decode().split('\t')[0]
        return chrom

def getID(phenotype):
    fids = [x.replace('_','-') for x in list(phenotype['FID'])]
    iids = [x.replace('_','-') for x in list(phenotype['IID'])]
    fid_iid = [fids[i]+'_'+iids[i] for i in np.arange(len(iids))]
    return fid_iid

def writeHeader(chrom_, fidiids, outfilename):
    with open(outfilename, 'w') as f_out:
        header = '##fileformat=VCF4.2\n'
        header += '##filedate=' + datetime.now().strftime('%Y%m%d') + '\n'
        header += '##source=dosage\n'
        header += '##contig=<ID=' + chrom_ + ',length=' + str(2**31 - 3) + '\n'
        header += '##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">\n'
        header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
        header += '\t'.join(fidiids) + '\n'
        f_out.write(header)

def dosageToGT(dosage):
    GT = './.'
    try:
        dosage = int(float(dosage))
        if dosage == 0:
            GT = '0/0'
        elif dosage == 1:
            GT = '0/1'
        elif dosage == 2:
            GT = '1/1'
        else:
            GT = './.'
    except ValueError:
        return GT
    return GT
    


########################################
# MAIN
########################################
if __name__=='__main__':
    
    # Read in arguments:
    parser = argparse.ArgumentParser(description="Tool to format dosage format to VCF")
    parser.add_argument('dosage_filename', help='The dosage filename that needs reformatting')
    parser.add_argument('phenotype_filename', help='The phenotype filename in csv format')
    parser.add_argument('output_filename', help='VCF output filename')
    args = parser.parse_args()
    
    # For testing:
#    dosage_filename = 'test.dosage.gz'
#    outvcf_filename = 'test.vcf.gz'
#    dosage_filename = dosage_filename.replace('.gz','').replace('.dosage','')+'.dosage.gz'
#    phenotype_filename = 'phenotype.gwas1.csv'
#    outvcf_filename = outvcf_filename.replace('.gz','').replace('.vcf','')+'.vcf'
#    logfilename = outvcf_filename.replace('.vcf','')+'.log'

    # Parse read arguments:
    dosage_filename = args.dosage_filename.replace('.gz','').replace('.dosage','')+'.dosage.gz'
    phenotype_filename = args.phenotype_filename
    outvcf_filename = args.output_filename.replace('.gz','').replace('.vcf','')+'.vcf'
    logfilename = outvcf_filename.replace('.vcf','')+'.log'

    # Output logfile in parallel
    old_stdout = sys.stdout
    log_file = open(logfilename, "w")
    sys.stdout = log_file
    
    # Log specific command details:
    print(datetime.now().strftime('%c'))
    print('args.dosage_filename: ' + dosage_filename)
    print('args.phenotypefilename: ' + phenotype_filename)
    print('args.output_filename: ' + outvcf_filename)
    
    
    print('Creating header')
    numHeaderLines=0
    phenotype = pd.read_csv(phenotype_filename)
    chrom = getchrom(dosage_filename)
    fidiid = getID(phenotype)
    writeHeader(chrom, fidiid, outvcf_filename)
    
    print('Reading in dosage file ' + dosage_filename)
    dosage = pd.read_csv(dosage_filename, skiprows=numHeaderLines, compression='gzip', encoding='utf-8', sep="\t", header=None)
    chromcol = dosage.iloc[:,0]
    idcol = dosage.iloc[:,1]
    poscol = dosage.iloc[:,2]
    refcol = dosage.iloc[:,3]
    altcol = dosage.iloc[:,4]
    dosagemat = np.around(dosage.iloc[:,5:])
    GTmat = []
    for row in np.arange(dosagemat.shape[0]):
        currow = [dosageToGT(x) for x in dosagemat.iloc[row,:]]
        GTmat.append(currow)
    GTdf = pd.DataFrame(GTmat)
    qualcol = pd.Series(np.repeat('.', len(chromcol)))
    filtercol = pd.Series(np.repeat('.', len(chromcol)))
    infocol = pd.Series(np.repeat('PR', len(chromcol)))
    formatcol = pd.Series(np.repeat('GT', len(chromcol)))
    
    print('Building VCF file')
    newvcf = pd.concat([chromcol, poscol, idcol, refcol, altcol, qualcol, filtercol, infocol, formatcol, GTdf ], axis = 1)
    print('Writing to disk')
    newvcf.to_csv(outvcf_filename, index=False, sep="\t", encoding='utf-8', mode='a', header=False)
    print('VCF written')
    print('Compressing')
    subprocess.run(args=['bgzip', outvcf_filename])
    print('Tabix indexing')
    subprocess.run(args=['tabix', '-p', 'vcf', outvcf_filename+'.gz'])
    print('Done')
    print(datetime.now().strftime('%c'))
    sys.stdout = old_stdout
    log_file.close()
    
    
    
    