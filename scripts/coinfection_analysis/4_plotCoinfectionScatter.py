import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# needs the output of 2_findCoInfectionsInStudyTsv.py to be pasted in "coinfection_samples_HA/NA. "
# and takes the path to the tsv produced by 3_calcVaporSndRefDistTsv.py as input (in bottom main method)
# this script then plots the coinfection scatter plot from thesis Fig. 4.32
# not used but kept here are code to plot histograms of pid (l 89-105) and a 3d Scatter, including the relative frequency of references of the second subtype in VAPOR output (function plot_scatter3D)

def plot_scatter(input_file, highlight_samples, logScale=False):
    df = pd.read_csv(input_file, sep='\t')
    
    ha_data = df[['sample', 'dist_identity_HA', 'dist_score_HA', 'rel.Freq_2nd_Subtype_HA']]
    na_data = df[['sample', 'dist_identity_NA', 'dist_score_NA', 'rel.Freq_2nd_Subtype_NA']]

    print("Total samples: ", len(ha_data))

    ha_data = ha_data[(ha_data['dist_identity_HA'] != -1) & (ha_data['dist_score_HA'] != -1)]
    na_data = na_data[(na_data['dist_identity_NA'] != -1) & (na_data['dist_score_NA'] != -1)]

    ha_data = ha_data[(ha_data['dist_identity_HA'] != -3) & (ha_data['dist_score_HA'] != -3)]
    na_data = na_data[(na_data['dist_identity_NA'] != -3) & (na_data['dist_score_NA'] != -3)]

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 14))


    ### HA ###
    ax1.scatter(ha_data['dist_identity_HA'], ha_data['dist_score_HA'], color='blue', label='HA')
    
    # conifection highlighted (obtained from 2_findCoinfectionsInStudyTsv.py)
    confirmedCo_HA_noSecondHit = len(ha_data[(ha_data['dist_identity_HA'] == -2) & (ha_data['dist_score_HA'] == -2) & ha_data['sample'].isin(COINFECTION_SAMPLES_HA)])
    confirmedCo_HA_WithSecondHit = len(ha_data[((ha_data['dist_identity_HA'] != -1) & (ha_data['dist_score_HA'] != -1) & ha_data['dist_identity_HA'] != -2) & (ha_data['dist_score_HA'] != -2) & ha_data['sample'].isin(coinfection_samples_HA)])
    if highlight_samples:
        highlighted_ha = ha_data[ha_data['sample'].isin(COINFECTION_SAMPLES_HA)]
        ax1.scatter(highlighted_ha['dist_identity_HA'], highlighted_ha['dist_score_HA'], color='green', label='Coinfection HA', s=100, edgecolors='black', zorder=5)
        print("HA conifection distances:")
        print(highlighted_ha['sample'], highlighted_ha['dist_identity_HA'], highlighted_ha['dist_score_HA'])
        print("Number of samples to highlight in NA:", len(highlighted_ha))

    ax1.set_title('HA segments by VAPOR coverage vs score distances (Best to next best hit with different HA call)')
    ax1.set_xlabel('Vapor query coverage difference HA')
    ax1.set_ylabel('Vapor score difference HA')
    ax1.set_xlim([-0.7, ha_data['dist_identity_HA'].max() * 1.1]) 
    ax1.set_ylim([0, ha_data['dist_score_HA'].max() * 1.1])
    if (logScale):
        ax1.set_ylim([max(ha_data['dist_score_HA'].min() * 0.1, 1), ha_data['dist_score_HA'].max() * 1.5])
        ax1.set_yscale('log')
    # display fraction of coinfections found in the plot.
    ax1.text(0.99, 0.05, f'Coinfects with no second hit: {confirmedCo_HA_noSecondHit}/{len(COINFECTION_SAMPLES_HA)} total', transform=ax1.transAxes, ha='right', va='bottom', fontsize=10)

    ax1.grid(True)
    ax1.legend()


    ### NA ### 
    ax2.scatter(na_data['dist_identity_NA'], na_data['dist_score_NA'], color='red', label='NA')
    
    # (coinfections obtained from 2_findCoinfectionsInStudyTsv.py)
    threeCoinfects = ['SRR16232028', 'SRR27181748', 'SRR27181764']
    confirmedCo_NA_noSecondHit = len(na_data[(na_data['dist_identity_NA'] == -2) & (na_data['dist_score_NA'] == -2) & na_data['sample'].isin(COINFECTION_SAMPLES_NA)])
    confirmedCo_NA_WithSecondHit = len(na_data[(na_data['dist_identity_NA'] != -2) & (na_data['dist_score_NA'] != -2) & na_data['sample'].isin(COINFECTION_SAMPLES_NA)])
    if highlight_samples:
        highlighted_na = na_data[na_data['sample'].isin(COINFECTION_SAMPLES_NA)]
        ax2.scatter(highlighted_na['dist_identity_NA'], highlighted_na['dist_score_NA'], color='green', label='Coinfection NA', s=100, edgecolors='black', zorder=5)
        print("NA conifection distances:")
        print(highlighted_na['sample'], highlighted_na['dist_identity_NA'], highlighted_na['dist_score_NA'])
        print("Number of samples to highlight in NA:", len(highlighted_na))

    ax2.set_title('NA segments by VAPOR coverage vs score distances (Best to next best hit with different NA call)')
    ax2.set_xlabel('Vapor query coverage difference NA')
    ax2.set_ylabel('Vapor score difference NA')
    ax2.set_xlim([-0.7, na_data['dist_identity_NA'].max() * 1.1])  
    ax2.set_ylim([0, na_data['dist_score_NA'].max() * 1.1])
    if (logScale):
        ax2.set_ylim([max(na_data['dist_score_NA'].min() * 0.1, 1), na_data['dist_score_NA'].max() * 1.5])
        ax2.set_yscale('log')

    ax2.text(0.99, 0.05, f'Coinfects with no second hit: {confirmedCo_NA_noSecondHit}/{len(COINFECTION_SAMPLES_NA)} total', transform=ax2.transAxes, ha='right', va='bottom', fontsize=10)

    ax2.grid(True)
    ax2.legend()


    ha_data = ha_data[(ha_data['dist_identity_HA'] != -2) & (ha_data['dist_score_HA'] != -2)]
    na_data = na_data[(na_data['dist_identity_NA'] != -2) & (na_data['dist_score_NA'] != -2)]

    # histogram HA pid distribution
    ax3.hist(ha_data['dist_identity_HA'], bins=200, color='blue', alpha=0.7, label='HA Identity') 
    ax3.set_title('Histogram for HA Identity')
    ax3.set_xlabel('dist_identity_HA')
    ax3.set_ylabel('Frequency')
    ax3.set_xlim([0, 0.3])
    ax3.grid(True)
    ax3.legend()

    # histogram NA pid distribution
    ax4.hist(na_data['dist_identity_NA'], bins=200, color='red', alpha=0.7, label='NA Identity') 
    ax4.set_title('Histogram for NA Identity')
    ax4.set_xlabel('dist_identity_NA')
    ax4.set_ylabel('Frequency')
    ax4.set_xlim([0, 0.3]) 
    ax4.grid(True)
    ax4.legend()

    ###

    plt.tight_layout()
    if highlight_samples:
        plt.savefig("Scatter_log_withCoInfectionHighlights.png")
    else:
        plt.savefig("Scatter_seconBestRef.png")
    plt.show()


def plot_scatter3D(input_file, highlight_samples, logScale=False):
    # without histograms, collapse with plot_scatter with many arguments? better than copy paste
    df = pd.read_csv(input_file, sep='\t')

    
    ha_data = df[['sample', 'dist_identity_HA', 'dist_score_HA', 'rel.Freq_2nd_Subtype_HA']]
    na_data = df[['sample', 'dist_identity_NA', 'dist_score_NA', 'rel.Freq_2nd_Subtype_NA']]

    print("Total samples: ", len(ha_data))

    ha_data = ha_data[(ha_data['dist_identity_HA'] != -1) & (ha_data['dist_score_HA'] != -1)]
    na_data = na_data[(na_data['dist_identity_NA'] != -1) & (na_data['dist_score_NA'] != -1)]
    ha_data = ha_data[(ha_data['dist_identity_HA'] != -2) & (ha_data['dist_score_HA'] != -2)]
    na_data = na_data[(na_data['dist_identity_NA'] != -2) & (na_data['dist_score_NA'] != -2)]

    fig = plt.figure(figsize=(20, 14))

    ### HA Plot ###
    ax1 = fig.add_subplot(121, projection='3d') 
    ax1.scatter(ha_data['dist_identity_HA'], ha_data['dist_score_HA'], ha_data['rel.Freq_2nd_Subtype_HA'], color='blue', label='HA')
    

    if highlight_samples:
        highlighted_ha = ha_data[ha_data['sample'].isin(COINFECTION_SAMPLES_HA)]
        ax1.scatter(highlighted_ha['dist_identity_HA'], highlighted_ha['dist_score_HA'], highlighted_ha['rel.Freq_2nd_Subtype_HA'], color='green', label='Coinfection HA', s=100, edgecolors='black', zorder=5)
        print("HA conifection distances:")
        print(highlighted_ha['sample'], highlighted_ha['dist_identity_HA'], highlighted_ha['dist_score_HA'], highlighted_ha['rel.Freq_2nd_Subtype_HA'])
        print("Number of samples to highlight in HA:", len(highlighted_ha))

    ax1.set_title('3D Scatter Plot for HA')
    ax1.set_xlabel('Vapor identity distance HA')
    ax1.set_ylabel('Vapor score distance HA')
    ax1.set_zlim([0, 20])
    ax1.set_zlabel('Relative Frequency of 2nd Subtype HA (%)')
    ax1.legend()

    ### NA Plot ###
    ax2 = fig.add_subplot(122, projection='3d') 
    ax2.scatter(na_data['dist_identity_NA'], na_data['dist_score_NA'], na_data['rel.Freq_2nd_Subtype_NA'], color='red', label='NA')
    
    if highlight_samples:
        highlighted_na = na_data[na_data['sample'].isin(COINFECTION_SAMPLES_NA)]
        ax2.scatter(highlighted_na['dist_identity_NA'], highlighted_na['dist_score_NA'], highlighted_na['rel.Freq_2nd_Subtype_NA'], color='green', label='Coinfection NA', s=100, edgecolors='black', zorder=5)
        print("NA conifection distances:")
        print(highlighted_na['sample'], highlighted_na['dist_identity_NA'], highlighted_na['dist_score_NA'], highlighted_na['rel.Freq_2nd_Subtype_NA'])
        print("Number of samples to highlight in NA:", len(highlighted_na))

    ax2.set_title('3D Scatter Plot for NA')
    ax2.set_xlabel('Vapor identity distance NA')
    ax2.set_ylabel('Vapor score distance NA')
    ax2.set_zlim([0, 20])
    ax2.set_zlabel('Relative Frequency of 2nd Subtype NA (%)')
    ax2.legend()

    plt.tight_layout()
    plt.savefig("Scatter3D_log_withCoInfectionHighlights.png")
    plt.show()

if __name__ == "__main__":
    input_file = 'vapor2ndRefDist_filtered_n10000.tsv'
    COINFECTION_SAMPLES_HA = ['SRR16232025', 'SRR16232026', 'SRR16232027', 'SRR16232028', 'SRR16232029', 'SRR16232030', 'SRR16232031', 'SRR16232032', 'SRR16232033', 'SRR16232034', 'SRR16232035', 'SRR16183372', 'SRR16183373', 'SRR16183374', 'SRR16183375', 'SRR16183376', 'SRR16183377', 'SRR16183378', 'SRR16183379', 'SRR16183381', 'SRR16183382', 'SRR16183383', 'SRR18273475', 'SRR18273476', 'SRR18273477', 'SRR18273478', 'SRR18273479', 'SRR18273480', 'SRR18273481', 'SRR18273482', 'SRR18273483', 'SRR18273484', 'SRR18273485', 'SRR18273486', 'SRR18273487', 'SRR18273488', 'SRR18273489', 'SRR18273490', 'SRR18273491', 'SRR18273492', 'SRR18273493', 'SRR18273494', 'SRR18273495', 'SRR18273496', 'SRR18273497', 'SRR18273498', 'SRR18273499', 'SRR18273500', 'SRR18273501', 'SRR18273502', 'SRR18273503', 'SRR18273504', 'SRR18273505', 'SRR18273506', 'SRR18273507', 'SRR18273508', 'SRR18273509', 'SRR18273510', 'SRR18273511', 'SRR18273512', 'SRR18273513', 'SRR18273514', 'SRR18273515', 'SRR18273516', 'SRR18273517', 'SRR18273518', 'SRR18273519', 'SRR18273520', 'SRR18273521', 'SRR18273522', 'SRR18273523', 'SRR18273524', 'SRR18273525', 'SRR18273526', 'SRR18273527', 'SRR27181099', 'SRR27181748', 'SRR27181764', 'SRR18273528', 'SRR18273529', 'SRR18273530', 'SRR18273531', 'SRR18273532', 'SRR18273533', 'SRR18273534', 'SRR18273535', 'SRR18273536', 'SRR18273537', 'SRR18273538', 'SRR18273539']
    COINFECTION_SAMPLES_NA = ['SRR16232025', 'SRR16232026', 'SRR16232027', 'SRR16232028', 'SRR16232029', 'SRR16232030', 'SRR16232031', 'SRR16232032', 'SRR16232033', 'SRR16232034', 'SRR16232035', 'SRR16183372', 'SRR16183373', 'SRR16183374', 'SRR16183375', 'SRR16183376', 'SRR16183377', 'SRR16183378', 'SRR16183379', 'SRR16183381', 'SRR16183382', 'SRR16183383', 'SRR18273475', 'SRR18273476', 'SRR18273477', 'SRR18273478', 'SRR18273479', 'SRR18273480', 'SRR18273481', 'SRR18273482', 'SRR18273483', 'SRR18273484', 'SRR18273485', 'SRR18273486', 'SRR18273487', 'SRR18273488', 'SRR18273489', 'SRR18273490', 'SRR18273491', 'SRR18273492', 'SRR18273493', 'SRR18273494', 'SRR18273495', 'SRR18273496', 'SRR18273497', 'SRR18273498', 'SRR18273499', 'SRR18273500', 'SRR18273501', 'SRR18273502', 'SRR18273503', 'SRR18273504', 'SRR18273505', 'SRR18273506', 'SRR18273507', 'SRR18273508', 'SRR18273509', 'SRR18273510', 'SRR18273511', 'SRR18273512', 'SRR18273513', 'SRR18273514', 'SRR18273515', 'SRR18273516', 'SRR18273517', 'SRR18273518', 'SRR18273519', 'SRR18273520', 'SRR18273521', 'SRR18273522', 'SRR18273523', 'SRR18273524', 'SRR18273525', 'SRR18273526', 'SRR18273527', 'SRR27181099', 'SRR27181748', 'SRR27181764', 'SRR18273528', 'SRR18273529', 'SRR18273530', 'SRR18273531', 'SRR18273532', 'SRR18273533', 'SRR18273534', 'SRR18273535', 'SRR18273536', 'SRR18273537', 'SRR18273538', 'SRR18273539']
    plot_scatter(input_file, True, True) # args: vaporScoreFile, highlight Coinfection toglle, score axis log scale toggle

