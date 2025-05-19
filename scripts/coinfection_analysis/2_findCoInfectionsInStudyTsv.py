import csv
import ast

KEYWORDS = ["co-inf", "coinfection"]

sample_to_study_tsv = 'sampleRunId_to_BioProjectData.tsv' # result of "1_queryStudyDesigs.py"
outList = 'samples_studyConfirmed_coinfection.txt' # list of samples that have keyowrd appearing in title or description, ready to copy and paste into plotCoinfectionScatter.py

bioproject_metadata = {}
with open(sample_to_study_tsv, encoding='utf-8') as tsvfile, \
    open(outList, 'w') as outFile:
    reader = csv.DictReader(tsvfile, delimiter='\t')

    outFile.write("['SRR27181099', 'SRR27181748', 'SRR27181764', ...")

    for row in reader:
        for keyword in KEYWORDS:
            if keyword in row["Title"].lower() or keyword in row["Abstract"].lower():
                print(row)
                outFile.write("'" + row["Sample"] + "',")
                break

    outFile.write("]")