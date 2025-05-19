
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import re
import csv
import concurrent.futures

# script that fastly (due to multithreading) can filter potential reference candidates, based on their label to bestAlignment subtype-match. For usage, see main function below.

# extract subtype from refernce collection header convention like strings
def extract_segment_number(header, segment):
    if segment == 'N':
        match = re.search(r'\|H\d+N(\d+)\|', header)
    elif segment == 'H':
        match = re.search(r'\|H(\d+)N\d+\|', header)
    else:
        return -1 # code internal mistake (other segments not allowed)
    if match:
        return int(match.group(1))
    else: 
        print("Error detecting subtype number (VAPOR ref header)")
        return None


# computes simple global pairwise alignment and returns their scores.
def align_sequences(query_seq, ref_sequences):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    scores = []
    for ref_id, ref_seq in ref_sequences.items():
        score = aligner.score(query_seq, ref_seq)
        scores.append((ref_id, score))
    return scores


# takes a single reference candidate, and aligns it to every representative subtype sequence, returns None of best alignment with correct subtype, Sequence header and scores otherwise
def process_sequence(record, ref_sequences, segment):
    query_header = record.id
    query_seq = str(record.seq)

    extracted_subtype = extract_segment_number(query_header, segment)
    scores = align_sequences(query_seq, ref_sequences)

    best_match = max(scores, key=lambda x: x[1])
    best_ref = best_match[0]

    match = re.search(rf'{segment}(\d+)', best_ref)
    if match:
        best_number = int(match.group(1))
    else:
        print("Error detecting subtype number (representative template sequence)")
        return None

    if extracted_subtype is not None and extracted_subtype != best_number:
        print("WRONGREF! Header (extracted", segment, "):", query_header, "(", extracted_subtype, "): Score:", scores)
        return [query_header] + [s[1] for s in scores]
    else:
        print("OK: Header (extracted", segment, "):", query_header, "(", extracted_subtype, "): Score:", scores)
        return None


# takes all inputs (reference fasta file and sequence candidates) and distribute to workers
def findWrongRefs(ref_fasta, input_fasta, output_tsv, segment):
    ref_sequences = {}
    for record in SeqIO.parse(ref_fasta, "fasta"):
        ref_sequences[record.id] = str(record.seq)

    with open(output_tsv, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        header = ["Sequence_Header"] + list(ref_sequences.keys())
        writer.writerow(header)


        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as executor:
            futures = []
            for record in SeqIO.parse(input_fasta, "fasta"):
                future = executor.submit(process_sequence, record, ref_sequences, segment)
                futures.append(future)

            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                if result:
                    writer.writerow(result)


if __name__ == "__main__":
    # NA
    ref_fasta = "../irmaNARefs.fasta"   
    input_fasta = "/home/galaxy/akolbecher/workdir/data/vapor_ref_collection/Influenza_per-segment_reference_collection/6-NA.fasta"
    output_tsv = "outputWrongRefsNA.tsv"       

    # HA
    #ref_fasta = "../irmaHARefs.fasta"
    #input_fasta = "/home/galaxy/akolbecher/workdir/data/vapor_ref_collection/Influenza_per-segment_reference_collection/4-HA.fasta"
    #output_tsv = "outputWrongRefsHA.tsv"       

    # call this function with: 
    # 1. fasta of reliable subtype sequences (Header should contain unique N<int> ord H<int> for subtype detectiuon)
    # 2. fasta of the potential reference candidates, that should be checked. (only of HA/NA, corresponding to the first input)
    # 3. name of output tsv, which shows the candidates, that achieved the highest alignment score to a subtype, different to their label.
    # stdOut of this script is a well formated log of the alignment score for every input sequence
    findWrongRefs(ref_fasta, input_fasta, output_tsv, "N")
