import os
import argparse
from Bio import SeqIO, Align


def process_batch(batch_path, identity_threshold):
    consensus_path = os.path.join(batch_path, "consensus_HA.fasta")
    ref_folder_path = os.path.join(batch_path, "ref_folder")
    output_path = os.path.join(batch_path, "filtered_identifiers_" + str(identity_threshold) + ".txt")

    if not os.path.exists(consensus_path) or not os.path.isdir(ref_folder_path):
        print("Missing required files or directories in batch:", consensus_path, ref_folder_path)
        return

    identifiers_above_threshold = []

    # consensus file
    with open(consensus_path, "r") as consensus_file:
        consensus_records = list(SeqIO.parse(consensus_file, "fasta"))

    for record in consensus_records:
        identifier = record.id.split("|")[0] 

        print("Handling:", identifier)
        if len(record.seq) == 0:
            print("Consensus sequence is empty for identifier:", identifier)
            continue

        file_name = identifier + ".fasta"
        ref_file_path = os.path.join(ref_folder_path, file_name)

        if not os.path.exists(ref_file_path):
            print("Reference file not found for identifier:", identifier)
            continue

        with open(ref_file_path, "r") as ref_file:
            ref_records = list(SeqIO.parse(ref_file, "fasta"))

        if not ref_records:
            print("No sequences found in reference file for identifier:", identifier)
            continue

        aligner = Align.PairwiseAligner()
        aligner.match_score = 1
        # mismatch > open and close gap
        aligner.open_gap_score = -1
        aligner.extend_gap_score = -0.3
        aligner.mismatch_score = 0

        alignment = aligner.align(record.seq, ref_records[0].seq)[0]

        sumMatches = 0
        for i in range(len(alignment[0])):
            char0 = alignment[0][i]
            char1 = alignment[1][i]
            if char0 == char1 or char0 == 'N' or char1 == 'N':  # mismatch with N should be free
                sumMatches += 1

        identity = sumMatches / len(alignment[0])

        print("identity", identity)

        if identity >= identity_threshold:
            print("Exceeds identity threshold:", identity)
            print(alignment)
            identifiers_above_threshold.append(identifier)
        else:
            print("Check :)")

    # output file
    with open(output_path, "w") as output_file:
        output_file.write("\n".join(identifiers_above_threshold))

    print("Processed batch", batch_path, ": ", len(identifiers_above_threshold), " identifiers written to", output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process batch directories for pairwise alignment of VAPOR reference and consensus")
    parser.add_argument("--input_dir", required=True, help="Path to the folder containing batch directories") 
    # the specified folder can have multiple subfolders, corresponding to batch folders
    # the batch folders must have 2 elements: consensus_HA.fasta & ref_folder (contains <sampleId>.fasta files)
    # first file is "Per-segment-consensus" from galaxy
    # second folder "ref_folder"  has to contain the sequences chosen by vapor (with sequence!)
    # I obtained by "(step 43) seqtk>collapse>replace output" of my workflow from galaxy - possibly easier to obtain with an additional seqtk step on VAPOR output.
    # HA häufiger konstruiert als NA. Deshalb das als abgleich

    args = parser.parse_args()

    batch_directories = [os.path.join(args.input_dir, d) for d in os.listdir(args.input_dir) if os.path.isdir(os.path.join(args.input_dir, d))]

    threshs = [0.999, 0.998, 0.997, 0.994, 0.992, 0.990, 0.988, 0.986, 0.984, 0.982, 0.980, 0.995, 0.993, 0.991, 0.989, 0.987, 0.985, 0.983, 0.981]

    for batch_path in batch_directories:
        for thresh in threshs:
            process_batch(batch_path, thresh)
