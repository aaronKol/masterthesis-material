import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
"""
example of plotting script for pie and boxplot charts.

COMMENT ONLY RELEVANT FOR VAPOR BUG FIXED COMPARISON PLOTS
906 files were not created by vaporBugFixed due to local runtime issues. they will be excluded from the analysis. Changing the order of the output file of vapor cannot lead to Not creating this file, while vapor did-
Of the remaining '-1' values, this is the breakdown of which run procudeced a result and if they are counted in the result:

galaxy issue (excluded from analysis via batchSummarizer-py):
Sample: SRR27181281, Segment: HA, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR27181281/*HA.fasta.tabular
Sample: SRR27181281, Segment: MP, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR27181281/*MP.fasta.tabular
Sample: SRR27181281, Segment: NA, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR27181281/*NA.fasta.tabular
Sample: SRR27181281, Segment: NP, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR27181281/*NP.fasta.tabular
Sample: SRR27181281, Segment: NS, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR27181281/*NS.fasta.tabular
Sample: SRR27181281, Segment: PA, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR27181281/*PA.fasta.tabular
Sample: SRR27181281, Segment: PB1, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR27181281/*PB1.fasta.tabular
Sample: SRR27181281, Segment: PB2, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR27181281/*PB2.fasta.tabular

8 one run empy, other N only - replaced by '-2' since no tool built meaningful consensus
Sample: SRR14420360, Segment: NP, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR14420360/*NP.fasta.tabular
Sample: SRR14420370, Segment: PB2, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR14420370/*PB2.fasta.tabular
Sample: SRR10972179, Segment: HA, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR10972179/*HA.fasta.tabular
Sample: SRR14420375, Segment: NA, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR14420375/*NA.fasta.tabular
Sample: SRR14420375, Segment: NP, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR14420375/*NP.fasta.tabular
Sample: SRR16359867, Segment: NP, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR16359867/*NP.fasta.tabular
Sample: SRR16359869, Segment: NP, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR16359869/*NP.fasta.tabular
Sample: SRR29478192, Segment: PB1, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/SRR29478192/*PB1.fasta.tabular

VAPOR managed to create consensus, Bugfixed didnt. keep this. real problems
Sample: ERR3474626, Segment: MP, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/ERR3474626/*MP.fasta.tabular
Sample: ERR3474634, Segment: NP, Expected pattern: /home/galaxy/akolbecher/others/VaporRefScoreBug/simultanForFailedSamples/outDir_backup/ERR3474634/*NP.fasta.tabular


-> subtract 906 from -1 set
"""

file_path = 'comment_withOutTwos_batchSummaryVAPORbugFixed.tsv'
df = pd.read_csv(file_path, sep='\t')

print(df.head()) 

""" all 3
category_counts = {
    '-1': (df['Distance'] == -1).sum(),
    '0': (df['Distance'] == 0).sum(),
    'positive': (df['Distance'] > 0).sum()
}
print(category_counts)

labels = ['Only VAPOR built segment', 'Built segments match', 'Built segments have pos. distance']
sizes = [category_counts['-1']-906, category_counts['0'], category_counts['positive']] #sizes for vaporBugRerun
colors = ['#ff9999','#66b3ff','#99ff99']  
""" 

#only 2
category_counts = {
    '0': (df['Distance'] == 0).sum(),
    'positive': (df['Distance'] > 0).sum()
}
print(category_counts)

labels = ['Built segments match', 'Built segments have pos. distance']
sizes = [category_counts['0'], category_counts['positive']] #sizes for vaporBugRerun
colors = ['#66b3ff','#99ff99']  

# pie chart
fig, ax = plt.subplots(figsize=(7, 7))
ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90, pctdistance=0.9, textprops={'fontsize': 13} )
ax.axis('equal') 
plt.title('Consensus comparison VAPOR vs VAPOR-Fixed')
plt.show()

# boxplot
positive_df = df[df['Distance'] > 0]

plt.figure(figsize=(8, 6))
sns.boxplot(x='Distance', data=positive_df)
plt.title('Consensus segmetns with Positive Distances')
positive_distances = sorted(positive_df['Distance'].unique())

plt.xticks(range(0, positive_df['Distance'].astype(int).max() + 10, 40))
plt.show()
