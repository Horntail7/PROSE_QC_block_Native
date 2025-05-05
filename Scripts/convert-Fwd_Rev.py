import math
import numpy as np
import argparse
import os.path
from Bio.Seq import Seq

parser = argparse.ArgumentParser()

parser.add_argument('--inp', '--input_file', type=str, default='Errors_Single_Mol.dat',
                    help='Full path to the input file (Errors_Single_Mol.dat)')
parser.add_argument('--ncol', '--ncolums', type=int, default='4',
                        help='No. of columns in input (4)')
parser.add_argument('--icol', '--icolum', type=int, default='3',
                        help='Insert column number (3)')
parser.add_argument('--pcol', '--pcolum', type=int, default='4',
                        help='Primer column number (4)')
parser.add_argument('--acol', '--acolum', type=int, default='2',
                        help='Insert column number (2)')
parser.add_argument('--bcol', '--bcolum', type=int, default='1',
                        help='Insert column number (1)')


args = parser.parse_args()

infile=open(args.inp, "r") 
lines=infile.readlines()


outfile=open("inserts_Rev.out","w")

for line in lines:
    words = line.split()
    if len(words) == args.ncol:
        aseq=Seq(words[args.acol-1]).reverse_complement()
        bseq=Seq(words[args.bcol-1]).reverse_complement()
        iseq=Seq(words[args.icol-1]).reverse_complement()
        pseq=Seq(words[args.pcol-1]).reverse_complement()
            
        outfile.write("%s\t%s\t%s\t%s\n"%(pseq,iseq,aseq,bseq))

outfile.close()



    
