import matplotlib.pyplot as plt
import numpy as np


# values are obtained, by using the script "removeSimilarSamplesSummaryTsv.py" 
# and from resulting filteredSummaryTsv then counting the remaining TOTAL and NON-ZERO entries. 
# This script then creates Fig. 4.18 of the thesis (# histogram 2)
fractionZeroDistSegments = {
    # percent, nonZeros, total
    "All samples":  (92.09, 2120, 26806),
    "100.0%":       (92.24, 1974, 25439),
    "99.8%":        (92.29, 1651, 21414),
    "99.6%":        (92.63, 1331, 18069),
    "99.4%":        (93.24, 1141, 16878),
    "99.2%":        (93.36, 1094, 16471),
    "99.0%":        (93.39, 1082, 16368), 
    "98.8%":        (93.44, 1067, 16262),
    "98.6%":        (93.56, 1037, 16097),
    "98.4%":        (93.60, 1019, 15933),
    "98.2%":        (93.55, 1011, 15676),
    "98.0%":        (93.44, 998, 15204),
}


labels = list(fractionZeroDistSegments.keys())
percentages = [value[0] for value in fractionZeroDistSegments.values()]
non_zeros = [value[1] for value in fractionZeroDistSegments.values()]
totals = [value[2] for value in fractionZeroDistSegments.values()]

# histogram 1: 
plt.figure(figsize=(10, 5))
plt.bar(labels, percentages)
plt.xlabel('Vapor Reference-Consensus Similarity Kickout Threshold')
plt.ylabel('Fraction of IRMA-Vapor Consensus Matches')
plt.title('Histogram of Vapor Reference Similarity Threshold vs. Portion of Consensus Matches')
plt.ylim(80, 100)
plt.xticks(rotation=45, ha="right")
plt.show()


# histogram 2:
zeroDistPortion = np.subtract(totals, non_zeros)

fig, ax = plt.subplots(figsize=(10, 5))

for i in range(len(labels)):
    # zeroDistPortion
    ax.bar(labels[i], zeroDistPortion[i], bottom=0, color='green', label='Consensus match' if i == 0 else "")#avoid multiple legends to be shown
    # nonZeroDistPortion
    ax.bar(labels[i], totals[i] - zeroDistPortion[i], bottom=zeroDistPortion[i], color='red', label='Consensus distance > 0' if i == 0 else "")
    # perccentages
    ax.text(labels[i], zeroDistPortion[i] / 2, f'{percentages[i]:.2f}%', 
            ha='center', va='center', color='white', fontsize=10)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], title="Legend")

plt.xlabel('Vapor Reference-Consensus Similarity Kickout Threshold')
plt.ylabel('Remaining assembled segment consensi')
plt.title('Histogram of Vapor Reference Similarity Threshold vs. Remaining Segments with Portion of Consensus Matches')
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.show()