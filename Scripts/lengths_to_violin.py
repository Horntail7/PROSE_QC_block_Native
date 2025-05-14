import math
import numpy as np
import argparse
import os.path
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def main():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--inp', '--input_file', type=str, default='lengths.dat',
                        help='Full path to the input file (lengths.dat)')
    parser.add_argument('--out', '--output_file', type=str, default='Lengths_violin.png',
                        help='Full path to the output file (Lengths_violin.png)')
    
    args = parser.parse_args()
    
    infile=open(args.inp, "r") 
    lines=infile.readlines()

    barcodes=lines[0].split()

    group=[]
    value=[]
    
    for line in lines[1:]:
        words=line.split()
        for i in range(len(words)):
            if int(words[i]) < 15:
                group.append(barcodes[i])
                value.append(int(words[i]))


    df = pd.DataFrame({"Barcodes": group,"Lengths": value})

    sns.violinplot(x="Barcodes", y="Lengths", data=df)

    # Group data and add mean ± std annotations
    grouped = df.groupby("Barcodes")["Lengths"]

    for i, (day, values) in enumerate(grouped):
        mean = values.mean()
        std = values.std()
        plt.text(
            i+0.1,               # x position
            mean + std + 1,  # y position
            f"μ={mean:.1f}\nσ={std:.1f}",
            ha='center', va='bottom', fontsize=10, color='black'
        )

    
    plt.savefig(args.out)
    plt.show()
    
if __name__ == "__main__":
    main()

    
