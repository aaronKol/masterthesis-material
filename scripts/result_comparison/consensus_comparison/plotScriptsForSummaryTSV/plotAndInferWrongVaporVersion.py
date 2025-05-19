import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# script performing the majority voting analysis on IRMA, VAPOR, VAPOR-Bugfixed. (Vapor,vapor-Bugfixed agreeing is ignored)
# input the 3 distance summary tsvs. line 39-60 are plotting another plot - which I didnt found as usefull
# from line 61 on, the "incondsistent" tool is infered and plotted as in Thesis Fig. 4.31B

df_irma_vapor = pd.read_csv("../results_vapor_irma/withZeroDist_batchSummaryVAPOR_IRMA.tsv", sep="\t")  # IRMA-VAPOR distance summary tsv
df_vapor_vaporFixed = pd.read_csv("../results_vapor_vaporBugFixed/withZeroDist_batchSummaryVAPORbugFixed.tsv", sep="\t")  # VAPOR-VAPORFIXED distance summary tsv
df_vaporFixed_irma = pd.read_csv("../results_vaporBugFixed_irma/withZeroDist_batchSummaryVAPORbugFixed_irma.tsv", sep="\t")  # VAPORFIXED-IRMA distance summary tsv

# only look at samples, somehow processed by every approach
samples1 = set(df_irma_vapor['Sample'])
samples2 = set(df_vapor_vaporFixed['Sample'])
samples3 = set(df_vaporFixed_irma['Sample'])
shared_samples = samples1 & samples2 & samples3
df_irma_vapor = df_irma_vapor[df_irma_vapor['Sample'].isin(shared_samples)]
df_vapor_vaporFixed = df_vapor_vaporFixed[df_vapor_vaporFixed['Sample'].isin(shared_samples)]
df_vaporFixed_irma = df_vaporFixed_irma[df_vaporFixed_irma['Sample'].isin(shared_samples)]

# additionally remove the samples, that could not be processed locally (in Bug fixed run - Not unprocessable due to the changed VAPOR code !)
exclude_df = pd.read_csv("../results_vapor_vaporBugFixed/local_errors.txt", sep="[:,]", usecols=[0,1], names=["Sample", "Segment"], skiprows=1)
exclude_set = set(zip(exclude_df["Sample"].str.strip(), exclude_df["Segment"].str.strip()))

def filter_exclusions(df):
    mask = df.apply(lambda row: (row["Sample"], row["Segment"]) not in exclude_set, axis=1)
    return df[mask]

df_irma_vapor = filter_exclusions(df_irma_vapor)
df_vapor_vaporFixed = filter_exclusions(df_vapor_vaporFixed)
df_vaporFixed_irma = filter_exclusions(df_vaporFixed_irma)

# concatinate dataframes
df_irma_vapor['Comparison'] = 'IRMA-VAPOR'
df_vapor_vaporFixed['Comparison'] = 'VAPOR-VAPORFIXED'
df_vaporFixed_irma['Comparison'] = 'VAPORFIXED-IRMA'
df = pd.concat([df_irma_vapor, df_vapor_vaporFixed, df_vaporFixed_irma])

# map distances to categories
def categorize(distance):
    if distance == -2:
        return "Both missing (-2)"
    elif distance == -1:
        return "One missing (-1)"
    elif distance == 0:
        return "Match (0)"
    else:
        return "Mismatch (>0)"

df['Category'] = df['Distance'].apply(categorize)

# plot count of agreement categories - not used in thesis, less informative. Kept if wanted in future
sns.countplot(data=df, x='Category', hue='Comparison')
plt.title("Agreement Categories (Shared Samples Only)")
plt.ylabel("Count")
plt.xticks(rotation=15)
plt.tight_layout()
plt.show()

# merge on sample-segmen and name distance column to dist_tool1_tool2
merged = df_irma_vapor.merge(df_vapor_vaporFixed, on=["Sample", "Segment"], suffixes=('_irma_vapor', '_vapor_vaporfixed'))
merged = merged.merge(df_vaporFixed_irma, on=["Sample", "Segment"])
merged.rename(columns={"Distance": "dist_vaporfixed_irma"}) 
#some empty fields were still there
merged = merged[merged['Segment'].notna() & (merged['Segment'].str.strip() != "")] 

def infer_disagreeing_tool(row):
    d1 = row['Distance_irma_vapor']
    d2 = row['Distance_vapor_vaporfixed']
    d3 = row['dist_vaporfixed_irma']

    # nur wenn alle das segment gebaut haben
    if min(d1, d2, d3) < 0:
        return "skip"

    if d1 == d3 == 0:
        return "all agree"
    elif d1 == 0:
        return "VAPORFIXED"
    elif d3 == 0:
        return "VAPOR"
    else:
        return "unclear"


merged["Inconsistent_Tool"] = merged.apply(infer_disagreeing_tool, axis=1)

disagreements = merged[merged["Inconsistent_Tool"].isin(["VAPOR", "VAPORFIXED"])]
print("VAPOR:", disagreements['Inconsistent_Tool'].value_counts()['VAPOR'])
print("VAPOR_FIXED:", disagreements['Inconsistent_Tool'].value_counts()['VAPORFIXED'])

sns.countplot(data=disagreements, x='Inconsistent_Tool')
plt.title("Inferred Tool-Specific Errors")
plt.ylabel("Count of Likely Errors")
plt.xlabel("Tool Likely to be Incorrect")
plt.tight_layout()
plt.show()
