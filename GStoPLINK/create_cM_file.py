#!/usr/bin/python
# Small script to create Step3_cM file
# Syntax: python create_cM_file.py -i <inputfile> -o <outputfile>
# Syntax2: python create_cM_file.py -i <number of snps> -o <outputfile>
# Example: python create_cM_file.py -i 750530 -o Step3_cM

import sys, getopt

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print 'test.py -i <inputfile> -o <outputfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'test.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg

   numsnps = int(inputfile)+1
   with open(outputfile,'w') as f:
      for i in range(1,numsnps):
         f.write('0\n')
   f.close()


if __name__ == "__main__":
   main(sys.argv[1:])
