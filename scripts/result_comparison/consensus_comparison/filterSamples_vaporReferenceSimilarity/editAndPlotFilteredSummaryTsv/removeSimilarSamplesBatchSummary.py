import argparse

parser = argparse.ArgumentParser(description='Filter lines batch summary based on sample IDs')
parser.add_argument('batchSummary', type=str, help='The first file with the sample data')
parser.add_argument('identifierList', type=str, help='The second file with the list of samples')
#std out then produces SummaryTsv with samples of the list removed

args = parser.parse_args()

with open(args.identifierList, 'r') as f:
    sample_set = set(f.read().splitlines())

filtered_lines = []
with open(args.batchSummary, 'r') as f:
    for line in f:
        columns = line.split('\t')
        sample = columns[0]
        if sample not in sample_set:
            filtered_lines.append(line)

for line in filtered_lines:
    print(line, end='')
