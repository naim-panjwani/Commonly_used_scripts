import pandas as pd
import subprocess
import argparse
import gzip

parser = argparse.ArgumentParser(description="Fixes missing value format for chrX file generated from PLINK")
parser.add_argument('vcf_filename', help='The VCF filename that needs reformatting')
parser.add_argument('outfilename', help="Output filename")
args = parser.parse_args()
vcf_filename = args.vcf_filename
outfilename = args.outfilename

def getNumHeaderLines(vcf_filename, num_lines_to_check = 1000):
    num_header_lines = 0
    num_lines_checked = 0
    with gzip.open(vcf_filename, 'rb') as f:
        nextline = f.readline().decode('utf-8')
        num_lines_checked += 1
        while nextline[0:2] == "##" and num_lines_checked <= num_lines_to_check:
            num_header_lines += 1
            nextline = f.readline().decode('utf-8')
            num_lines_checked += 1
    return num_header_lines

def writeHeader(infile, outfilename):
    with gzip.open(infile, 'rb') as f_in:
        with open(outfilename, 'wb') as f_out:
            nextline = f_in.readline().decode('utf-8')
            while nextline[0:1] == "#":
                f_out.write(bytes(nextline, 'UTF-8'))
                nextline = f_in.readline().decode('utf-8')


afile=vcf_filename
outfile=outfilename.replace('.gz','')
print('Reading VCF')
vcf = pd.read_csv(afile, sep="\t", skiprows=getNumHeaderLines(afile))
print('Fixing . for ./.')
newvcf = pd.concat([ vcf.iloc[:,0:9], vcf.iloc[:,9:].replace(".","./.") ], axis=1)
print('Saving to file')
writeHeader(afile, outfile)
newvcf.to_csv(outfile, sep="\t", index=False, encoding='utf-8', mode='a', header=False)
print('Compressing')
compressvcf = subprocess.run(args=['bgzip', outfile], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
if compressvcf.returncode != 0:
  print('Could not compress vcf')
print('Done')

