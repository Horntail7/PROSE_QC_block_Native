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
parser.add_argument('--ncut', '--ncutoff', type=int, default='10',
                        help='No of copies (10)')
parser.add_argument('--spol', '--strand_polarity', type=str, default='Fwd',
                    help='Polarity of the strand (Fwd or Rev)')


args = parser.parse_args()

infile=open(args.inp, "r") 
lines=infile.readlines()


inserts={}
umis={}
copies={}
count=0

for line in lines:
    words = line.split()
    if len(words) == args.ncol:
        if len(words[args.pcol-1]) == 7 and len(words[args.icol-1]) < 10 and len(words[args.acol-1]) == 13 \
           and len(words[args.bcol-1]) == 16:
            temp = words[args.acol-1][:4]+words[args.acol-1][-4:]+words[args.bcol-1][:4]+words[args.bcol-1][-4:]
            if temp in umis.keys():
                if words[args.pcol-1] not in umis[temp]:
                    inserts[temp].append(words[args.icol-1])
                    umis[temp].append(words[args.pcol-1])
                    copies[temp].append(1)
                else:
                    copies[temp][umis[temp].index(words[args.pcol-1])]+=1
            else:
                inserts[temp]=[words[args.icol-1]]
                umis[temp]=[words[args.pcol-1]]
                copies[temp]=[1]

                
print(len(inserts))
inserts_sorted=dict(sorted(inserts.items(), key=lambda item: len(item[1]), reverse=True))
umis_sorted=dict(sorted(umis.items(), key=lambda item: len(item[1]), reverse=True))
copies_sorted=dict(sorted(copies.items(), key=lambda item: len(item[1]), reverse=True))

count=0
outfile=open("inserts_per_template.out","w")
if args.spol == "Fwd":
    for key,value in inserts_sorted.items():
        if len(value) > args.ncut-1:
            count+=1
            outfile.write("%s\t"%Seq(key).reverse_complement())
            for val in value:
                outfile.write("%s\t"%Seq(val).reverse_complement())
            outfile.write("\n")
elif args.spol == "Rev":
    for key,value in inserts_sorted.items():
        if len(value) > args.ncut-1:
            count+=1
            outfile.write("%s\t"%key)
            for val in value:
                outfile.write("%s\t"%val)
            outfile.write("\n")
        
#print("Templates with atleast 10 copies = ", count)
outfile.close()

#for key in inserts_sorted.keys():
#    print(key,len(umis_sorted[key]), umis_sorted[key])
#    print(key,len(inserts_sorted[key]), inserts_sorted[key])
#    print(key,len(copies_sorted[key]), copies_sorted[key])


    
