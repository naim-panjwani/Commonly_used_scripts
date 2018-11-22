# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 12:01:17 2018

@author: naim
"""

import pandas as pd
import argparse
import subprocess

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Fix the BED file to collapse several gene row entries to one')
    parser.add_argument('bedfile', help='The BED file that requires fixing')
    args = parser.parse_args()

    print('args.befile ' + args.bedfile)    
    
    outname = ""
    if args.bedfile.endswith('.bed.gz'):
        outname = args.bedfile.replace('.bed.gz','_fixed.bed')
    elif args.bedfile.endswith('.gz'):
        outname = args.bedfile.replace('.gz','_fixed')
    elif args.bedfile.endswith('.bed'):
        outname = args.bedfile.replace('.bed','_fixed.bed')
    else:
        outname= str(args.bedfile + "_fixed")
    
    print('Reading file')
    bedfile = pd.DataFrame()
    if args.bedfile.endswith('gz'):
        bedfile = pd.read_csv(args.bedfile, compression='gzip', sep="\t")
    else:
        bedfile = pd.read_csv(args.bedfile, sep="\t")
    
    print('Grouping same gene_id\'s')
    bed_grp = bedfile.groupby('gene_id')
    bedstart = bed_grp['start'].min().reset_index()
    bedend = bed_grp['end'].max().reset_index()
    bedchr = bed_grp['#chr'].first().reset_index()
    
    newbed = pd.merge(bedchr,bedstart,on='gene_id')
    newbed = pd.merge(newbed,bedend,on='gene_id')
    
    print('Re-building the new bed file')
    for col in bedfile.columns[4:]:
        newbedcol = bed_grp[col].first().reset_index()
        newbed = pd.merge(newbed, newbedcol, on='gene_id')
    
    print('Re-indexing')
    cols = ['#chr','start','end','gene_id']
    cols.extend([i for i in newbed.columns[4:]])
    newbed = newbed.reindex(columns=cols)
    
    print('Removing \'chr\' from the chromosome field')
    newbed['#chr'] = [newbed['#chr'][i].replace('chr','') for i in range(len(newbed['#chr']))]
    
    print('Sorting')
    newbed = newbed.sort_values(by=['#chr','start'])
    
    print('Writing to file ' + outname)
    newbed.to_csv(outname, index=False, sep="\t")
    subprocess.check_call('bgzip -f '+outname, shell=True, executable='/bin/bash')
    subprocess.check_call('tabix -f '+outname+'.gz', shell=True, executable='/bin/bash')

    print('Done')