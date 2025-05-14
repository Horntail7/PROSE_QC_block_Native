import math
import numpy as np
import argparse
import os.path
import random


def main():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--inp', '--input_file', type=str, default='inserts.out',
                        help='Full path to the input file (inserts.out)')
    parser.add_argument('--ncol', '--ncolums', type=int, default='3',
                        help='No. of columns in input (3)')
    parser.add_argument('--icol', '--icolums', type=int, default='1',
                        help='Insert column number (1)')
    parser.add_argument('--Nsample', '--Num_of_samples', type=int, default='10000',
                        help='No. of samples (10000)')

    
    args = parser.parse_args()
    
    infile=open(args.inp, "r") 
    lines=infile.readlines()

    random_sample = random.sample(lines, args.Nsample)
    for line in random_sample:
        words=line.split()
        if len(words) == args.ncol:
            print(len(words[args.icol-1]))

            
if __name__ == "__main__":
    main()

    
