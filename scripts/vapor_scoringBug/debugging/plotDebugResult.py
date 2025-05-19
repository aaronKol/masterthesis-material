import pandas as pd
import matplotlib.pyplot as plt

# script for visualizing the scoring of wrong and correct references. 
# input is a VAPOR Search Result, which can be produced with running VAPOR in a debug mode, specifing the reference sequence identifiers of interest (The best (bugged) and second best probably)
# an example debugSearchResult for which the scoring bug occured is placed next to the script.
data = pd.read_csv('debugSearchResult.txt', sep=' ', header=None)

positions = data[0]
dataset1_bool = data[2] 
dataset1_raw = data[3]
dataset1_filled = data[4] 
dataset1_deque = data[5]  
dataset2_bool = data[9]
dataset2_raw = data[10]
dataset2_filled = data[11] 
dataset2_deque = data[12] 

dataset1_isGap = []
dataset2_isGap = []

for i in positions:
    dataset1_isGap.append('g' in data[6][i])
    dataset2_isGap.append('g' in data[13][i])


plt.figure(figsize=(10, 8))

# correct ref
plt.subplot(2, 2, 2)  
plt.plot(positions, dataset1_raw, label="Weight - raw", linestyle='-', color='b')
plt.plot(positions, dataset1_filled, label="Weight - filled", linestyle='--', color='g')
plt.plot(positions, dataset1_deque, label="Weight - deque", linestyle='-.', color='y')
plt.title("Array values of second reference hit")
plt.xlabel("Position")
plt.ylabel("Values")
plt.legend()


# bugged
plt.subplot(2, 2, 1) 
plt.plot(positions, dataset2_raw, label="Weight - raw", linestyle='-', color='b')
plt.plot(positions, dataset2_filled, label="Weight - filled", linestyle='--', color='g')
plt.plot(positions, dataset2_deque, label="Weight - deque", linestyle='-.', color='y')
plt.title("Array values of first reference hit")
plt.xlabel("Position")
plt.ylabel("Values")
plt.legend()


#  boolean value present in dbg
plt.subplot(2, 2, 4)  
plt.scatter(positions, dataset1_bool, color='b', alpha=0.7)
plt.title("Second ref. - Kmer present in dbg")
plt.xlabel("Position")
plt.ylabel("Present in dbg (True/False)")
plt.yticks([0, 1], ['False', 'True'])  
plt.legend()

plt.subplot(2, 2, 3)  
plt.scatter(positions, dataset2_bool, color='r', alpha=0.7)
plt.title("First ref. - Kmer present in dbg")
plt.xlabel("Position")
plt.ylabel("Present in dbg (True/False)")
plt.yticks([0, 1], ['False', 'True'])  
plt.legend()

plt.tight_layout()
plt.show()


# gaps plots
plt.figure(figsize=(10, 8))

plt.subplot(1, 2, 1)  
plt.scatter(positions, dataset1_isGap, color='b', alpha=0.7)
plt.title("Correct - Gaps")
plt.xlabel("Position")
plt.ylabel("Is a gap in dbg (True/False)")
plt.yticks([0, 1], ['False', 'True'])  
plt.legend()
plt.text(0.99, 0.05, f'#non-Gaps / #Positions: {dataset1_isGap.count(False)}/{len(positions)} total', ha='left', va='center_baseline', fontsize=10)


plt.subplot(1, 2, 2)  
plt.scatter(positions, dataset2_isGap, color='r', alpha=0.7)
plt.title("Bugged - Gaps")
plt.xlabel("Position")
plt.ylabel("Is a gap in dbg (True/False)")
plt.yticks([0, 1], ['False', 'True'])  
plt.legend()
plt.text(0.99, 0.05, f'#non-Gaps / #Positions: {dataset2_isGap.count(False)}/{len(positions)} total', ha='left', va='center_baseline', fontsize=10)

plt.tight_layout()
plt.show()

# verify other calculation changes scoring order in VAPOR result
dataset1_score = sum(dataset1_deque)
dataset2_score = sum(dataset2_deque)

dataset1_fracNonGaps = dataset1_isGap.count(False)/len(positions)
dataset2_fracNonGaps = dataset2_isGap.count(False)/len(positions)

print(dataset1_fracNonGaps, dataset1_score)
print(dataset2_fracNonGaps, dataset2_score)

print("Total correct:", dataset1_fracNonGaps*dataset1_score)
print("Total bugged:", dataset2_fracNonGaps*dataset2_score)

