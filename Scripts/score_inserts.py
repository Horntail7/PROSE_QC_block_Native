import math
import numpy as np
import argparse
import os.path
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import csv

MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_OPEN_SCORE = -1
GAP_EXTEND_SCORE = -1
DEGENERATE_MATCH_SCORE = 0

def create_substitution_matrix():
    matrix = substitution_matrices.load("NUC.4.4")
    for base in "ATCG":
        matrix[base, base] = MATCH_SCORE
        for other_base in "ATCG":
            if base != other_base:
                matrix[base, other_base] = MISMATCH_SCORE
        matrix[base, "N"] = DEGENERATE_MATCH_SCORE
        matrix["N", base] = DEGENERATE_MATCH_SCORE
    matrix["N", "N"] = DEGENERATE_MATCH_SCORE
    return matrix

def score(seq, template, matrix):
    max_possible_score = len(template) * MATCH_SCORE

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = GAP_OPEN_SCORE
    aligner.extend_gap_score = GAP_EXTEND_SCORE
    
    alignments = aligner.align(seq, template)
    best_alignment = alignments[0]

    return best_alignment.score/max_possible_score


def main():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--inp', '--input_file', type=str, default='inserts_per_template.out',
                        help='Full path to the input file (inserts_per_template.out)')

    args = parser.parse_args()
    
    infile=open(args.inp, "r") 
    lines=infile.readlines()

    insert_templates=['TGTACAA','TGTGCAA','TGTTCAA','TGTCCAA','TGTCAA']
    sub_matrix=create_substitution_matrix()

    outfile=open("scores.csv","w")
    outfile.write("UMI,sA,sG,sT,sC,nS,nR\n")
    inserts={}
    for line in lines:
        words=line.split()
        outfile.write("%s,"%words[0])
        inserts[words[0]]=words[1:]
        ncopies=len(inserts[words[0]])
        scores={}
        for temp in insert_templates:
            scores[temp]=0
        scores['rest']=0
        for ins in inserts[words[0]]:
            flag=0
            for temp in insert_templates:
                temp_score=score(ins,temp,sub_matrix)
                if temp_score > 0.99:
                    scores[temp]+=1
                    flag=1
                    break
            if flag == 0:
                scores['rest']+=1

        for key in scores.keys():
            scores[key]/=ncopies

        for temp in insert_templates:
            outfile.write("%6.2f,"%scores[temp])
        outfile.write("%6.2f\n"%scores['rest'])

if __name__ == "__main__":
    main()

    
