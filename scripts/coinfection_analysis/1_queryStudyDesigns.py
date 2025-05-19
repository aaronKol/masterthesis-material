import csv
import xml.etree.ElementTree as ET
from Bio import Entrez

Entrez.email = "-"

sra_runTable = "/home/galaxy/akolbecher/workdir/data/quality_comparison_set/SraRunTable.csv" # csv runtable containing columns BioProject and Run with metadata of interest
accession_list = "/home/galaxy/akolbecher/workdir/data/quality_comparison_set/Batches/all_first64/all_first64_identifiers.txt" # IN: txt list (newline seperated) of used sample identifiers

output_tsv = "sampleRunId_to_BioProjectData.tsv"  # OUT: tsv outfile name

used_samples = []
with open(accession_list, 'r') as f:
    for line in f:
        used_samples.append(line.strip('\n'))

print("Searching for samples:", len(used_samples))


bioproject_metadata = {}
with open(sra_runTable, newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',')

    sample_bioproject_map = {}
    for row in reader:
        run_id = row['Run']
        if run_id in used_samples: # nur die, die wir auch prozessiert haben
            bioproject_id = row['BioProject']  
            sample_bioproject_map[run_id] = bioproject_id
            bioproject_metadata[bioproject_id] = None #query the results into this dict, then write sample-to-metadata out
    print("Done reading the sra run table. Bioprojects to query: ", len(bioproject_metadata.keys()))

for bioproject_id in bioproject_metadata.keys():
    print("Fetching: ", bioproject_id)
    try:
        handle = Entrez.efetch(db="bioproject", id=bioproject_id, rettype="xml")
        records = handle.read()
        handle.close()
        print(records)

        root = ET.fromstring(records)
        bioproject_data = {}

        title = root.find('.//Title').text
        description = root.find('.//Description').text
        bioproject_data = {"Title": title, "Abstract": description}
        print(bioproject_data)

        bioproject_metadata[bioproject_id] = bioproject_data
    except Exception as e:
        try: # print handle and Try one more time
            print(handle)
            handle = Entrez.efetch(db="bioproject", id=bioproject_id, rettype="xml")
            records = handle.read()
            handle.close()
            print(records)

            root = ET.fromstring(records)
            bioproject_data = {}

            title = root.find('.//Title').text
            description = root.find('.//Description').text
            bioproject_data = {"Title": title, "Abstract": description}
            print(bioproject_data)

            bioproject_metadata[bioproject_id] = bioproject_data
        except Exception as e:
            print("Failed to query (2times)", bioproject_id, "with eror:", e)
            bioproject_metadata[bioproject_id] = {"Title": "Error", "Abstract": "Error"}

# write out
with open(output_tsv, "w", newline='', encoding='utf-8') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(["Sample", "Title", "Abstract"])

    for run_id, bioproject_id in sample_bioproject_map.items():
        meta = bioproject_metadata.get(bioproject_id, {"Title": "N/A", "Abstract": "N/A"})
        writer.writerow([run_id, meta["Title"], meta["Abstract"]])

print("DONE")
