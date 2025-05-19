import os
import re

# script that takes regular VAPOR scoring output files and creates a TSV of distances of score and pid of best reference to best refernce of another subtype

def extract_subtype(reference_string, segment_type):
    if segment_type == 'HA':
        match = re.search(r'H(\d+)', reference_string)
    elif segment_type == 'NA':
        match = re.search(r'N(\d+)', reference_string) 
    if match:
        return int(match.group(1)) 
    return None

def calculate_distance(file_path, segment_type):
    if not os.path.exists(file_path):
        return -1, -1, -1 
    with open(file_path, 'r') as f:
        first_line = f.readline().strip().split('\t')
        first_identity = float(first_line[0])
        first_score = float(first_line[1])
        first_reference = first_line[-1].split('|')[-2]
        first_subtype = extract_subtype(first_reference, segment_type)
        second_subtype = None
        found_second = False
        dist_identity = -1
        dist_score = -1
        lines_first = 1
        lines_second = 0

        for line in f:
            line_data = line.strip().split('\t')
            curr_identity = float(line_data[0])
            curr_score = float(line_data[1])
            curr_reference = line_data[-1].split('|')[-2]
            curr_subtype = extract_subtype(curr_reference, segment_type)

            if first_subtype == curr_subtype:
                lines_first += 1
                continue

            if found_second is False:
                # for first other hit, calculate the distance
                if curr_subtype != first_subtype:
                    #print(first_line)
                    #print(line_data)
                    #print(first_subtype, curr_subtype)
                    found_second = True
                    #dist_identity = abs(first_identity - curr_identity)
                    #dist_score = abs(first_score - curr_score) - if only positive distances are desired, use this
                    dist_identity = first_identity - curr_identity
                    dist_score = first_score - curr_score
                    if first_identity - curr_identity < -0.05: # more thatn 5% negative
                        print(file_path, dist_identity)
                        #return -3, -3, -3
                    second_subtype = curr_subtype
                    lines_second += 1
            elif curr_subtype == second_subtype:
                # for later hits, only increase counter
                lines_second += 1

        if found_second:
            return dist_identity, dist_score, lines_second / lines_first

    return -2, -2, -2  # = nothing found, but file exists = really no hits

def process_samples(input_folder, output_file):
    with open(output_file, 'w') as out_file:
        out_file.write("sample\tdist_identity_HA\tdist_identity_NA\tdist_score_HA\tdist_score_NA\trel.Freq_2nd_Subtype_HA\trel.Freq_2nd_Subtype_NA\n")

        for sample_folder in os.listdir(input_folder):
            sample_path = os.path.join(input_folder, sample_folder)

            if os.path.isdir(sample_path):
                ha_file = os.path.join(sample_path, '4-HA.fasta.tabular')
                na_file = os.path.join(sample_path, '6-NA.fasta.tabular')

                dist_ha_id, dist_ha_sc, relFreq_2ndref_ha = calculate_distance(ha_file, 'HA')
                dist_na_id, dist_na_sc, relFreq_2ndref_na = calculate_distance(na_file, 'NA')


                out_file.write(f"{sample_folder}\t{dist_ha_id}\t{dist_na_id}\t{dist_ha_sc}\t{dist_na_sc}\t{relFreq_2ndref_ha}\t{relFreq_2ndref_na}\n")
    
    print("Results written to", output_file)

if __name__ == "__main__":
    input_folder = '/home/galaxy/akolbecher/workdir/data/Download_f64_seqtk_vaporScores/VaporRefScores_filteredRef_n10000' # folder containing VAPORs scoring result files
    output_file = 'vapor2ndRefDist_filtered_n10000.tsv'
    process_samples(input_folder, output_file)
