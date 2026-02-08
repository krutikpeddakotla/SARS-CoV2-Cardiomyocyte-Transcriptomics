import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats

RNA_seq_read = pd.read_csv("GSE151879_raw_counts_genes.Adult_human_cardiomyocytes.txt", sep='\t')
#print(RNA_seq_read.head().to_string())

#### PHASE 1: Data Preparation & Cleaning ####
## Creating a Design Matrix
design_matrix_list = RNA_seq_read.columns
split = RNA_seq_read.columns.str.split(pat="_", expand=True).to_list()
dsm = pd.DataFrame(split, index=design_matrix_list).drop(index="gene_id")
dsm.rename(columns={0:"Age", 1:"Species", 2:"Sample Type", 3:"Condition", 4:"Type_ID"}, inplace=True)
#print(dsm.to_string())

## Dropping values with 0 counts across all six samples
RNA_seq_read.replace(0,np.nan, inplace=True)
RNA_seq_read.dropna(subset=[
    "Adult_human_cardiomyocytes_Mock_1", "Adult_human_cardiomyocytes_Mock_2", "Adult_human_cardiomyocytes_Mock_3", "Adult_human_cardiomyocytes_SARS-CoV2_1", "Adult_human_cardiomyocytes_SARS-CoV2_2", "Adult_human_cardiomyocytes_SARS-CoV2_3"],
    how="all", inplace=True)
RNA_seq_read.replace(np.nan,0, inplace=True)
#print(RNA_seq_read.head().to_string())

## Normalization
gene_id = RNA_seq_read["gene_id"]
drop_col = RNA_seq_read.drop(columns="gene_id")
normalized_col = drop_col.div(drop_col.sum(axis=0), axis=1) * 1000000
RNA_seq_read_normal = pd.concat([gene_id, normalized_col], axis=1)
#print(RNA_seq_read_normal.to_string())

#### PHASE 2: Differential Expression Analysis ####
## Applying Log2 to balance skewed data
RNA_seq_read_normlog = np.log2(RNA_seq_read_normal.iloc[:,1:7] + 1)

## Calculating Means
mock_cols = [
    "Adult_human_cardiomyocytes_Mock_1",
    "Adult_human_cardiomyocytes_Mock_2",
    "Adult_human_cardiomyocytes_Mock_3"
]
sars_cols = [
    "Adult_human_cardiomyocytes_SARS-CoV2_1",
    "Adult_human_cardiomyocytes_SARS-CoV2_2",
    "Adult_human_cardiomyocytes_SARS-CoV2_3"
]
mean_mock = RNA_seq_read_normlog[mock_cols]
RNA_seq_read_normlog["Mean_Mock_Log"] = mean_mock.mean(axis=1, numeric_only=True)
mean_infected = RNA_seq_read_normlog[sars_cols]
RNA_seq_read_normlog["Mean_Infected_Log"] = mean_infected.mean(axis=1, numeric_only=True)

## Calculating Fold Change
RNA_seq_read_normlog["Fold_Change"] = (RNA_seq_read_normlog["Mean_Infected_Log"] -
                                       RNA_seq_read_normlog["Mean_Mock_Log"])
RNA_seq_read_normlog = pd.concat([gene_id, RNA_seq_read_normlog], axis=1)
#print(RNA_seq_read_normlog.sort_values(by="Fold_Change", ascending=False).head(50).to_string())

#### PHASE 3: More Data Analysis ####
RNA_seq_read_normlog = RNA_seq_read_normlog.set_index("gene_id")

## T-Test
tstat, pvalue = stats.ttest_ind(RNA_seq_read_normlog[mock_cols], RNA_seq_read_normlog[sars_cols], equal_var=False, axis=1)
RNA_seq_read_normlog["p_value"] = pvalue
RNA_seq_read_normlog["log_p_value"] = -np.log10(RNA_seq_read_normlog["p_value"])

gene_reference = pd.read_csv("Ensemble_Gene_NameDescription.txt")
gene_reference.rename(columns={"Gene stable ID":"gene_id"}, inplace=True)
gene_table = RNA_seq_read_normlog.merge(gene_reference, how="outer", on="gene_id").set_index(
    "gene_id")
print(gene_table.sort_values(by="Fold_Change", ascending=False).head(50).to_string())
#//
gene_table.to_csv("SARS-CoV2_DGE_Results")
#//

#### PHASE 4: Visualization ####
magnitude_change = gene_table["Fold_Change"]
stat_significance = gene_table["log_p_value"]

gene_table["Expression"] = "Not Significant"
gene_table.loc[(gene_table["Fold_Change"] > 1) & (gene_table["p_value"] < 0.05), "Expression"] = "Up-Regulated"
gene_table.loc[(gene_table["Fold_Change"] < -1) & (gene_table["p_value"] < 0.05), "Expression"] = "Down-Regulated"

## Volcano plot -- Differential Gene Expression in SARS-CoV-2 Infected Adult Human Cardiomyocytes
f1 = plt.figure(1)
volcano_plot = sns.scatterplot(gene_table, x=magnitude_change,y=stat_significance, hue="Expression", palette="magma")
plt.title("Differential Gene Expression in SARS-CoV-2 Infected Adult Human Cardiomyocytes")
plt.axvline(x=-1, color="black", linestyle="dashed")
plt.axvline(x=1, color="black", linestyle="dashed")
plt.axhline(y=1.301, color="black", linestyle="dashed")
plt.tight_layout()
plt.savefig("volcano_plot",dpi=300)

## Heatmap plot - Top 20 Genes by Fold Change
f2= plt.figure(2)
gene_table["Display_Name"] = gene_table["Gene name"].fillna(gene_table.index.to_series())
gene_table_1 = gene_table.sort_values(by="Fold_Change", ascending=False).head(20)
gene_table_20 = gene_table_1.set_index("Display_Name")[mock_cols + sars_cols]
heatmap_plot = sns.heatmap(gene_table_20, cmap="RdBu_r", center=0)
plt.title("Top 20 Genes by Fold Change")
#//
plt.savefig("heatmap_top20", dpi=300)
plt.show()
#//
