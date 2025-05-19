# pass as many sampleSummaries (in specific format) as you want. Tool will name the sample after the containing folder of each summary.
# this script is meant to be called as defined at the bottom of the Snakefile.
# produces three 4 output tsvs.

import argparse
import csv
import os

EXCLUDE_SAMPLES = []
"""
These samples were excluded only for VAPOR-BUGFIXED comparison. They failed due to local problems, further explained in the plotting script plotVaporVaporBugFixed_dist.py.
#EXCLUDE_SAMPLES = ['SRR27181281','SRR6117005', 'SRR28995157', 'SRR6117006', 'SRR28995158', 'SRR28995159', 'SRR6117008', 'SRR28995160', 'SRR6117009', 'SRR7598294', 'SRR28995161', 'SRR7598295', 
#                  'SRR28995162', 'SRR7598298', 'SRR6117012', 'SRR28388469', 'SRR28995164', 'SRR6117013', 'SRR28995165', 'SRR6117014', 'SRR28995166', 'SRR6117015', 'SRR28995167', 
#                  'SRR28995168', 'SRR28995169', 'SRR28995170', 'SRR28995171', 'SRR8875156', 'SRR28995172', 'SRR8875157', 'SRR28995173', 'SRR28388480', 'SRR28995175', 'SRR28995176', 
#                  'SRR28995177', 'SRR7106939', 'SRR8875164', 'SRR6364208', 'SRR6364209', 'SRR28388487', 'SRR7253124', 'SRR5428785', 'SRR5428786', 'SRR28995121', 'SRR5428787', 'SRR28995122', 
#                  'SRR28995123', 'SRR5428789', 'SRR28995124', 'SRR5930376', 'SRR28995125', 'SRR5930377', 'SRR5428791', 'SRR28995126', 'SRR30610047', 'SRR28995127', 'SRR5930379', 
#                  'SRR28995129', 'SRR28995130', 'SRR28995131', 'SRR5930383', 'SRR28995132', 'SRR28995133', 'SRR28995134', 'SRR28995135', 'SRR28995136', 'SRR5894315', 'SRR28995138', 
#                  'SRR6116991', 'SRR6116992', 'SRR6116993', 'SRR28995145', 'SRR6116995', 'SRR28995147', 'SRR28995148', 'SRR28995149', 'SRR28995155', 'SRR6117004'] 
"""

parser = argparse.ArgumentParser(
                    prog='BatchSummarizer',
                    description='Summarizes multiple sampleSummaries into one big Distant Summary TSV file')

parser.add_argument('--output', '-o', help="Specify output file name")
parser.add_argument('--threshold', '-t', help="Threshold for highDistsOut argument")
parser.add_argument('--highDistsOut', '-d', help="File that contains samples that have a higher distance than threshold t")
parser.add_argument('files', nargs='+', help="TSV output of sampleSummaries.py")
args = parser.parse_args()

rowsToWrite = [['Sample', 'Segment', 'Distance', 'Path']] # header row
rowsWithoutError = []

for f in args.files:
    with open(f, 'r') as file:
        matrix = list(csv.reader(file, delimiter="\t", quotechar='"'))

        sampleName = os.path.realpath(f).split('/')[-2]
        for segmentRow in matrix[2:]:
            if sampleName in EXCLUDE_SAMPLES:
                print("Exclude sample ", sampleName)
            elif segmentRow[1] != "0":
                rowsToWrite.append([sampleName] + segmentRow)
            else:
                rowsWithoutError.append([sampleName] + segmentRow)

print("# samples:", len(args.files))
        
# write summary for different configurations
print("Writing summary to " + args.output)
with open(args.output, 'w', newline='') as tsvfile, \
     open(args.highDistsOut, 'w') as highDistsFile: #additional output for quantitative analysis of high distant segments
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    writerHighDists = csv.writer(highDistsFile, delimiter='\t', lineterminator='\n')
    for row in rowsToWrite:
        writer.writerow(row)
        if row[2] == 'Distance': #skip header
            continue
        if int(row[2]) >= int(args.threshold): 
            writerHighDists.writerow(row)

print("and to " + "withZeroDist_" + args.output)
with open("withZeroDist_" + args.output, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    for row in rowsToWrite:
        writer.writerow(row)
    for row in rowsWithoutError:
        writer.writerow(row)
        
print("and to " + "withOutTwos_" + args.output)
with open("withOutTwos_" + args.output, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    for row in rowsToWrite:
        if row[2] == 'Distance':#skip header
            writer.writerow(row)
            continue
        if int(row[2]) != -2:
            writer.writerow(row)
    for row in rowsWithoutError:
        writer.writerow(row)
        