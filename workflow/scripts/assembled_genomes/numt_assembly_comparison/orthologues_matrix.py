import pandas as pd
import numpy as np
import seaborn as sns

# genomes=pd.read_table("config/genomes.tsv").iloc[0:-2]
genomes = pd.read_table(snakemake.input.genomes).iloc[0:-2]

# df=pd.read_csv("results/genome_comparison/all_numts_updated.csv")
df = pd.read_csv(snakemake.input.df)
df = df.replace(".", np.nan)

statuscols = [i for i in df.columns if "_status" in i]

dfs = {}

for index, row in genomes.iterrows():
    tempdf = pd.read_csv(
        f"data/genome_comparison/additional_regions_compatible_seqs/{row['Species']}/{row['Assembly Name']}/{row['Assembly Accession']}/{row['Genome Abbreviation']}/all_numts/{row['Genome Abbreviation']}.csv"
    )
    tempdf["genome"] = row["Genome Abbreviation"]
    dfs[row["Genome Abbreviation"]] = tempdf


# df[df['REF10_status'].fillna(">").str.contains("\|COMPATIBLE_SEQUENCE")][['REF10_chromosome','REF10_region_start','REF10_region_end']].to_csv("test/ref10_compatible.bed", index=False, header=None, sep='\t')

dfs = pd.concat(dfs)

dfs["region_id"] = dfs["region_id"].where(
    dfs["additional_region_id"].isna(), dfs["additional_region_id"]
)

single_genome_regions = {i: 0 for i in genomes["Genome Abbreviation"]}

for i in single_genome_regions.keys():
    tempdf = dfs[dfs["genome"] == i]
    single_genome_regions[i] = len(
        tempdf[~tempdf["region_id"].isin(df[f"{i}_region_id"])]["region_id"].unique()
    )

dfs = dfs.set_index("genome")

found_matrix = {i: {} for i in genomes["Genome Abbreviation"]}
not_found_matrix = {i: {} for i in genomes["Genome Abbreviation"]}

for a in genomes["Genome Abbreviation"]:
    found_matrix[a] = {}
    for b in genomes["Genome Abbreviation"]:
        if a != b:
            tempdf = df[[f"{a}_region_id", f"{b}_region_id"]]
            found_matrix[a][b] = len(tempdf.dropna())
            not_found_matrix[a][b] = (
                len(dfs.loc[a]["region_id"].dropna().unique()) - found_matrix[a][b]
            )
        else:
            found_matrix[a][b] = len(dfs.loc[a]["region_id"].dropna().unique())
            not_found_matrix[a][b] = 0

found_df = pd.DataFrame(found_matrix)
not_found_df = pd.DataFrame(not_found_matrix)

found_df = found_df[genomes["Genome Abbreviation"]].loc[genomes["Genome Abbreviation"]]
not_found_df = not_found_df[genomes["Genome Abbreviation"]].loc[
    genomes["Genome Abbreviation"]
]

found_matrix_lower = np.tril(found_df.astype(str))
not_found_matrix_upper = np.triu(not_found_df.transpose().astype(str) + "|") + np.triu(
    not_found_df.astype(str)
)
found_matrix = pd.DataFrame(found_matrix_lower)
not_found_matrix = pd.DataFrame(not_found_matrix_upper)
not_found_matrix = not_found_matrix.astype(str).replace("0", "")
found_matrix = found_matrix.astype(str).replace("0", "")

matrix = pd.DataFrame(found_matrix + not_found_matrix)
matrix.columns = genomes["Genome Abbreviation"].tolist()
matrix.index = genomes["Genome Abbreviation"].tolist()

np.mean(list(single_genome_regions.values())[0:-2])
np.min(list(single_genome_regions.values()))
np.max(list(single_genome_regions.values())[0:-2])


for i in single_genome_regions.keys():
    matrix.loc[i][i] = (
        str(len(dfs.loc[i]["region_id"].unique())) + "|" + str(single_genome_regions[i])
    )


matrix.to_csv(snakemake.output[0])
# matrix.to_csv("tables/genome_comparison/orthologues_matrix_updated.csv")

#####
"""
for i in matrix.columns:
    matrix[i][i]=matrix[i][i].split("|")[0]


pctg_in_11={}

for i in matrix.columns:
    pctg_in_11[i]=int(matrix.loc[i]['REF11'])/int(matrix[i][i])

pd.DataFrame([pctg_in_11]).transpose().to_csv("tables/genome_comparison/regions_in_ref11_pctg_updated.csv")


"""
