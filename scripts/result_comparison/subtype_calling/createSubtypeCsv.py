import os
import csv
import glob
import re
from collections import defaultdict

"""
This script takes VAPORs and IRMAs result (see very end), extracts their called subtype in different ways (first functions),
and creates a subtypeTSV result, whith these columns "sample", "irma_HA", "irma_NA", "vapor_HA", "vapor_NA"
If no subtpye could be extracted (different reasons possible listed below) a '?' is written at the corresponding position.
"""

def extract_irma_subtype(irma_folder, segment):
    irma_subtype = defaultdict(list)
    for sample_folder in glob.glob(os.path.join(irma_folder, "*")):
        sample_name = os.path.basename(sample_folder)
        fasta_file = os.path.join(sample_folder, f"{segment}.fasta")
        
        if os.path.exists(fasta_file):
            with open(fasta_file, 'r') as file:
                header = file.readline().strip()
                if segment in header:
                    parts = header.split('_')
                    if len(parts) > 2: # irmas fasta file names are of the form A_HA_H1 for example.
                        subtype = parts[2]
                        irma_subtype[sample_name] = subtype[1]
                    else:
                        irma_subtype[sample_name] = '?'  # No valid subtype found
                else:
                    irma_subtype[sample_name] = '?'  # No match for the segment in header
        else:
            irma_subtype[sample_name] = '?'  # No fasra file found in the folder

    return irma_subtype


def extract_vapor_subtype(vapor_file):
    vapor_subtype = {}

    with open(vapor_file, 'r') as file:
        for line in file.readlines():
            cols = line.split('\t')

            if cols[1] == "\n":
                vapor_subtype[cols[0]] = {"HA": '?', "NA": '?'}
                continue
            elif len(cols) == 2:
                ha = extract_type_number(cols[1], 'H')
                na = extract_type_number(cols[1], 'N')
            elif len(cols) == 3:
                if cols[1] == "":
                    ha = "?"
                    na = extract_type_number(cols[2], 'N')
                elif cols[2] == "\n":
                    na = "?"
                    ha = extract_type_number(cols[1], 'H')
                else:
                    print("Irgendwas kann nicht stimmen", len(cols))
                    continue
            else:
                print("Weird line found", cols)

            vapor_subtype[cols[0]] = {"HA": ha, "NA": na} 

    return vapor_subtype

# regex helper for VAPOR
def extract_type_number(input_string, letter): 
    pattern = rf"{letter}(\d+)"
    match = re.search(pattern, input_string)
    return int(match.group(1)) if match else None


def combine_data(irma_folder, vapor_file, output_file):
    irma_data_H = extract_irma_subtype(irma_folder, "HA")
    irma_data_N = extract_irma_subtype(irma_folder, "NA")
    vapor_data = extract_vapor_subtype(vapor_file)

    samples = sorted(set(irma_data_H.keys()).union(vapor_data.keys()))
    with open(sample_list, 'r') as file:
        allSamples = [line.strip() for line in file if line.strip()]

    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(["sample", "irma_HA", "irma_NA", "vapor_HA", "vapor_NA"])
        commonToolFails = []

        for sample in allSamples:
            if sample not in samples:
                commonToolFails.append(sample)
                continue
            irmaH = irma_data_H[sample]
            irmaN = irma_data_N[sample]
            if irmaH == []:
                irmaH = irmaN = '?'
            if sample in vapor_data.keys():
                vaporH = vapor_data[sample]["HA"]
                vaporN = vapor_data[sample]["NA"]
            else:
                vaporH = vaporN = '?'
            writer.writerow([sample, irmaH, irmaN, vaporH, vaporN])
            if not (str(irmaH) == str(vaporH) and str(irmaN) == str(vaporN)):
                print("Not equel in: ", sample, "irma:", irmaH, irmaN, "vapor:", vaporH, vaporN)
    print(commonToolFails) # makes sure every sample could be processed
    print(len(commonToolFails))

sample_list = "/home/galaxy/akolbecher/workdir/data/quality_comparison_set/Batches/all_first64/all_first64_identifiers.txt" # txt-file containing all sample accessions to be analyzed, each in a new row.
irma_folder = "../consensus_comparisom/irma/fasta_consensus" # path to directory containing IRMAs consensus results (./<sampleAccession>/<segment>.fasta)
vapor_file = "./vapor/vapor_subtype.tabular" # subtype results from AIV VAPOR workflow. Cat them together if all batches should be included in one tsv.
output_file = "./united_variant_callings.tsv" # output tsv file name 

combine_data(irma_folder, vapor_file, output_file)