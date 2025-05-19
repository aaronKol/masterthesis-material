# pass as many MSA of fasta format as you want. MSA_{sample}.fa
# there has to be the corresponding (same_named) MSA_{sample}.fa in the same directory (will be linked in summary)

# the debug prints are kept, since they can be greatly used for understanding how a MSA is scored. uncomment every print statement in calcDistance() for that.

import argparse
import csv
import os
from Bio import AlignIO

GAP_OPEN_COST = 1
GAP_EXT_COST = 0
MISMATCH = 1
IGNORE_AMBIG_DIFFS = True

CHARS_GAP = {'-', '.', 'n', 'N'}
CODES = [
    {'A', 'R', 'W', 'M', 'D', 'H', 'V'},  # codes including A
    {'C', 'S', 'Y', 'M', 'B', 'H', 'V'},  # codes including C
    {'G', 'R', 'K', 'S', 'B', 'D', 'V'},  # codes including G
    {'T', 'K', 'W', 'Y', 'B', 'D', 'H'}   # codes including T
]


def calcDistance(f):
    if os.path.getsize(f) == 0:  # empty MSA
        return -2, False
    print(f)
    distance = 0
    continousGapCost = 0
    inMatchingRegion = False
    inGapExtenstionSeqOne = False
    inGapExtenstionSeqTwo = False
    containsAmbig = False
    align = list(AlignIO.parse(f, "fasta"))
    if len(align[0]) != 2: # only one sequenz
        return -1, False
    seq1 = align[0][0]
    seq2 = align[0][1]
    for col in range(align[0].get_alignment_length()):

        if (seq1[col] in CHARS_GAP or seq2[col] in CHARS_GAP) and not inMatchingRegion:  # remove front overhang
            #print(seq1[col], seq2[col], "Front overhang")
            continue
        
        inMatchingRegion = True

        if inGapExtenstionSeqOne and not (seq1[col] in CHARS_GAP):
            distance += continousGapCost
            #print(f"gap ends s1. ADD: {continousGapCost} [{distance}]")
            continousGapCost = 0
            inGapExtenstionSeqOne = False
        if inGapExtenstionSeqTwo and not (seq2[col] in CHARS_GAP):
            #distance += continousGapCost
            #print(f"gap ends s2. ADD: {continousGapCost} [{distance}]")
            continousGapCost = 0
            inGapExtenstionSeqTwo = False

        if seq1[col] != seq2[col]:
            # count gaps from back and substract them

            if seq1[col] in CHARS_GAP and seq2[col] in CHARS_GAP:
                #print(seq1[col], seq2[col], f"BOTH GAPS [{distance}]")   
                continue

            if (not inGapExtenstionSeqOne) and (seq1[col] in CHARS_GAP):
                continousGapCost += GAP_OPEN_COST
                inGapExtenstionSeqOne = True
                #print(seq1[col], seq2[col], f"OPEN GAP s1 ({GAP_OPEN_COST}) [{distance}]")
                continue
            elif inGapExtenstionSeqOne and (seq1[col] in CHARS_GAP):
                continousGapCost += GAP_EXT_COST
                #print(seq1[col], seq2[col], f"EXTEND GAP s1 ({continousGapCost}) [{distance}]")
                continue

            if (not inGapExtenstionSeqTwo) and (seq2[col] in CHARS_GAP):
                continousGapCost += GAP_OPEN_COST
                inGapExtenstionSeqTwo = True
                #print(seq1[col], seq2[col], f"OPEN GAP s2 ({GAP_OPEN_COST}) [{distance}]")
                continue
            elif inGapExtenstionSeqTwo and (seq2[col] in CHARS_GAP):
                continousGapCost += GAP_EXT_COST
                #print(seq1[col], seq2[col], f"EXTEND GAP s2 ({continousGapCost}) [{distance}]")
                continue
            
            if IGNORE_AMBIG_DIFFS:
                if any(seq1[col].upper() in s and seq2[col].upper() in s for s in CODES):
                    #print(seq1[col], seq2[col], f"IGNORE AMBIG [{distance}]")
                    containsAmbig = True
                    distance += 0  # maybe add a specific cost for those aswell...
                    continue

            distance += 1
            print(seq1[col], seq2[col], col,  f"MISMATCH (+{MISMATCH}) [{distance}]")
        #else:
            #print(seq1[col], seq2[col], f"[{distance}]")

    #print(distance, containsAmbig)
    return distance, containsAmbig




# main
parser = argparse.ArgumentParser(
                    prog='SampleSummarizer',
                    description='Summarizes multiple MSA of one sample into one tsv, which will be SampleSummary.tsv in the same directory as the inputs')

parser.add_argument('files', nargs='+', help="Output tsvs from MSA")
args = parser.parse_args()

allDistance = 0
rowsToWrite = [['segment', 'distance', 'path'], ['All']] # header & all row to be replaced
rowsAmbig = [['segment', 'path']]
  
for f in args.files:
    # handle one MSA distance file at a time
    dirPath = os.path.dirname(os.path.realpath(f))
    segment = f.split('/')[-1].split('_')[-1].split('.')[0]
    distance, containsAmbig = calcDistance(f)
    #if distance >= 0:
    rowsToWrite.append([segment, distance, os.path.realpath(f)])
    allDistance += int(distance)
    if containsAmbig:
        rowsAmbig.append([segment, os.path.realpath(f)])


rowsToWrite[1] = ['All', str(allDistance), dirPath] # replace all row!   

# write summary
outfile = dirPath + '/SampleSummary.tsv'
print("Writing summary to " + outfile)
with open(outfile, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    for row in rowsToWrite:
        writer.writerow(row)
outAmbigs = dirPath + '/AmbigsSummary.tsv' # additionally added. Not part of thesis. Gives overview over how many Ambgig-mismatches are "ignored"
with open(outAmbigs, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    for row in rowsAmbig:
        writer.writerow(row)



