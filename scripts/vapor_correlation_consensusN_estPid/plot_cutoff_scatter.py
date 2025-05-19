import os
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

# script that reads in VAPOR result file and VAPOR consensus sequence and plots a scatter of VAPOR pid & consensus fraction of Ns

# read VAPOR scoring result file
def parse_tabular(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip().split('\t')
        if len(first_line) < 2:
            return None, None  
        return float(first_line[0]), float(first_line[1])

# read Vapor consensus and count N 
def count_N_in_fasta(fasta_path):
    if not os.path.exists(fasta_path):
        return None, None  # no vapor consensus
    with open(fasta_path, 'r') as file:
        lines = file.readlines()
        sequence = ''.join(lines[1:]).replace('\n', '') 
        if len(sequence) == 0:
            return None, None #-1
        relative_count = sequence.count('N') 
        return relative_count, len(sequence)

def collect_and_plot(vapor_scores_path, vapor_fasta_path, batch_summary_tsv, plot3D=False):
    segment_mapping = {
        "1-PB2.fa.tabular": "PB2", "3-PA.fa.tabular": "PA", "4-HA.fa.tabular": "HA",
        "5-NP.fa.tabular": "NP", "6-NA.fa.tabular": "NA", "7-MP.fa.tabular": "MP", "8-NS.fa.tabular": "NS"
    }

    df = pd.read_csv(batch_summary_tsv, sep='\t')

    data_score = []  
    data_identity = []  
    data_all = []
    data_score_irma_identity = []
    data_identity_irma_numN = []

    for sample in os.listdir(vapor_scores_path): #start from existing vapor score files
        sample_path = os.path.join(vapor_scores_path, sample)
        if not os.path.isdir(sample_path):
            continue

        for file_name, segment in segment_mapping.items(): # look if and which consensus exist
            tabular_path = os.path.join(sample_path, file_name)
            file_name = segment + ".fasta"
            fasta_path = os.path.join(vapor_fasta_path, sample, file_name)

            if not os.path.exists(tabular_path):
                continue
            
            # now collect all desired values
            identity, score = parse_tabular(tabular_path)
            num_N, seq_length = count_N_in_fasta(fasta_path)
            
            sample_segment_row = df.loc[(df['Sample'] == sample) & (df['Segment'] == segment), 'Distance']
            if len(sample_segment_row.values > 0):
                dist_to_irma = sample_segment_row.values[0]
            else:
                dist_to_irma = None

            if score is not None and num_N is not None:
                num_N = num_N / seq_length
                score = score / seq_length

                num_not_N = 1 - num_N
                data_score.append((score, num_not_N))
                data_identity.append((identity, num_not_N))
                data_all.append((identity, score, num_not_N))
                #save sample/semgments, that have high coverage, many Ns.
                if 0 < num_not_N < 0.2 and identity > 0.8:
                    print("I", sample, segment, "NumNotN:", num_not_N, "VaporCov:", identity)
                if 0 == num_not_N and identity > 0.8:
                    print("II", sample, segment, "NumNotN:", num_not_N, "VaporCov:", identity)

                if dist_to_irma is not None:
                    data_score_irma_identity.append((score, dist_to_irma, identity))
                    data_identity_irma_numN.append((identity, num_not_N, dist_to_irma))

    # plot (3d plot didnt make sense...)
    # 1
    x_vals, y_vals, color_vals = zip(*data_all) if data_all else ([], [], [])
    plt.scatter(x_vals, y_vals, c=color_vals, cmap='viridis', alpha=0.6)
    plt.colorbar(label='Called bases in consensus (non-Ns)')
    plt.xlabel('% coverage')
    plt.ylabel('Score')
    plt.yscale('log') 
    plt.title('Vapor coverage + score leading to amount of base calls (color - relative amount of non-Ns)')
    plt.show()
    # 2
    x_vals, y_vals = zip(*data_identity) if data_identity else ([], [])
    plt.scatter(x_vals, y_vals, alpha=0.6)
    plt.xlabel('% Identity')
    plt.xlim(0.8, 1.0)
    #plt.xscale('log') für score
    plt.ylabel('Relative amount of base calls in consensus (non-Ns)')
    plt.title('Correlation: Vapor coverage (%) vs. base calls consensus')
    plt.show()

    # consider dist to irma
    # 3
    x_vals, y_vals, color_vals = zip(*data_score_irma_identity) if data_score_irma_identity else ([], [], []) #catch case empty
    plt.scatter(x_vals, y_vals, c=color_vals, cmap='viridis', alpha=0.6)
    plt.colorbar(label='Vapor coverage (%)')
    plt.xlabel('Score')
    plt.xscale('log')
    plt.ylabel('Dist To Irma')
    plt.ylim(-0.5, 20)
    #plt.yscale('log') 
    plt.title('Vapor Score (and coverage (color)) vs Dist To Irma')
    plt.show()
    # 4 
    x_vals, y_vals, irma_dists = zip(*data_identity_irma_numN) if data_identity_irma_numN else ([], [], []) # ""
    # map irma dist color values to (0-75). 75 then represent also the completely high values
    color_vals = [min(value, 15) for value in irma_dists]
    plt.scatter(x_vals, y_vals, c=color_vals, cmap='viridis', alpha=0.6)
    plt.colorbar(label='Distance to irma')
    plt.xlabel('% coverage')
    plt.xlim(0.8, 1.0)
    plt.ylabel('Relative amount of base calls in consensus (non-Ns)')
    #plt.ylim(-0.5, 20)
    #plt.yscale('log') 
    plt.title('Vapor coverage vs relative amount of non-Ns (colored by consensus distance to irma )')
    plt.show()


if __name__ == "__main__":
    vapor_scores_folder = "/home/galaxy/akolbecher/workdir/data/Download_f64_seqtk_vaporScores/VaporRefScores"  # folder of VAPOR tabular scoring output files
    vapor_consensus_folder = "/home/galaxy/akolbecher/workdir/quality_comparison/consensus_comparisom/vapor/fasta_consensus" # folder of VAPORs consensus sequences
    batch_summary_tsv = "/home/galaxy/akolbecher/workdir/quality_comparison/consensus_comparisom/64_withZeroDist_batchSummary.tsv" # SummaryTSV (distance to irma) for including this as a parameter
    collect_and_plot(vapor_scores_folder, vapor_consensus_folder, batch_summary_tsv, plot3D=True)
