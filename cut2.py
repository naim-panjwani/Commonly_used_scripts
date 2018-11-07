#!/usr/local/bin/python3.3

import sys, getopt, gzip

def usage():
   print('cut2.py -i <inputfile> -o <outputfile> --columnsfile=<columnsfile> --colindices=<True/False> --headerLine=<int>')
   print('colindices is optional and defaults to False')

def getHeader(inputfile, headerLine):
   header=[]
   currLine=0
   with gzip.open(inputfile,'rb') as f_in:
      for line in f_in:
         currLine += 1
         if currLine == int(headerLine):
            line = line.decode('utf-8')
            header = line.split('\t')
            f_in.close()
            break
   return(header)

def main():
   inputfile = ''
   outputfile = ''
   columnsfile=''
   colindices=False
   acceptableBoolStrings = ['True','False','T','F']
   headerLine=1
   try:
      opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile=", "columnsfile=", "colindices=", "headerLine="])
   except getopt.GetoptError as err:
      print(str(err))
      usage()
      sys.exit(2)
   print('opts: ', opts)
   print('args: ', args)

   for opt, arg in opts:
      if opt == '-h':
         usage()
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("--columnsfile"):
         columnsfile = arg
      elif opt in ("--colindices"):
         if arg.title() not in acceptableBoolStrings:
            print("colindices must be either 'True' or 'False' or 'T' or 'F'")
            sys.exit()
         else:
            if arg.title() in ['True', 'T']:
               colindices=True
      elif opt in ("--headerLine"):
         headerLine = int(arg)
      else:
        assert False, "unhandled option"

   print('Input file is ', inputfile)
   print('Output file is ', outputfile)
  
   print('Reading in column names') 
   columnnames = []
   with open(columnsfile, 'rb') as f_in:
     for line in f_in:
        columnnames.append(line.decode('utf-8').split('\t')[0]) # will only take first column values and ignore any other columns in the file 

   print('Determining column indices to extract')
   columnIndices = []
   columnsUnavailable = []
   if colindices:
      columnIndices = [(int(x)-1) for x in columnnames]
   else:
      header = getHeader(inputfile, headerLine)
      columnIndices = [header.index(x) for x in columnnames if x in header]

   print('Writing the header to file')
   currLine=0
   header = getHeader(inputfile, headerLine)
   outHeader = [header[i] for i in columnIndices]
   newline = '\t'.join(outHeader)
   newline += '\n'

   print('Extracting the desired columns line by line')
   with gzip.open(inputfile,'rb') as f_in:
     with gzip.open(outputfile, 'wb') as f_out:
       f_out.write(bytes(newline, 'UTF-8')) # writing the header
       for line in f_in:
         currLine += 1
         line = line.decode('utf-8')
         if(currLine > headerLine):
            columns = line.split('\t')
            columns_to_keep = [columns[i] for i in columnIndices]
            newline = '\t'.join(columns_to_keep)
            newline += '\n'
            f_out.write(bytes(newline, 'UTF-8')) 

if __name__ == '__main__':
  main()

