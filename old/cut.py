#!/usr/local/bin/python3.3

import sys, getopt, gzip

def main():
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError as err:
      print(str(err))
      print('cut.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('cut.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      else:
        assert False, "unhandled option"
   print('Input file is ', inputfile)
   print('Output file is ', outputfile)
   with gzip.open(inputfile,'rb') as f_in:
     with open(outputfile, 'wb') as f_out:
       for line in f_in:
         line = line.decode('utf-8')
         if(line[0:2] != "##"):
           columns = line.split('\t')
           columns = columns[0:8]
           newline = '\t'.join(columns)
           newline += '\n'
           f_out.write(bytes(newline, 'UTF-8'))

if __name__ == '__main__':
  main()

